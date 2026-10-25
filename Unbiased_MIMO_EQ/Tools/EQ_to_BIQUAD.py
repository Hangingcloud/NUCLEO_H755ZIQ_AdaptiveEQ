#!/usr/bin/env python3
"""Generate CM7/Core/Inc/mic_eq_coeffs.h from the REW files in Tools/.

Usage (project root):  python3 Tools/EQ_to_BIQUAD.py [--check]    (--check also prints how well the result follows the REW files)
Needs numpy and scipy.

Inputs (Tools/):
  PCM1808_L_fit.txt, PCM1808_R_fit.txt        REW "Filter Settings" of the PCM1808 channels (PK / LS / HS), converted to biquads
  SN1 to reference.txt, SN2 to reference.txt  REW minimum phase copy of A/B = ECM8000 / SN: the correction of the SN microphone
                                              (magnitude and phase). SN_ON_PF6 / SN_ON_PA3 below select the file per input pin.

For every microphone:
  1. 16 biquads (13 PK, 2 LS, 1 HS; RBJ cookbook, shelf Q = 1/sqrt(2)) fitted to the complex response of the file:
     magnitude 20 Hz .. 20 kHz, phase 20 Hz .. 3 kHz (taper to 5 kHz). The file is averaged in 1/48 octave bins.
  2. A 128 tap FIR (delay 64 samples) for what is left above about 1.4 kHz in phase and magnitude (gain limited to +-6 dB).
The polarity is not in the files: see MIMO_MIC_POLARITY in mimo_config.h.
"""
import math, os, re, sys
import numpy as np
from scipy.optimize import least_squares

FS = 44100.0                       # sample rate the firmware runs at (44.1 kHz)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "CM7", "Core", "Inc", "mic_eq_coeffs.h")
MAX_SECTIONS = 20                  # must match MIC_EQ_MAX_SECTIONS in main.c
SHELF_Q = 1.0 / math.sqrt(2.0)

PCM = {"pcmL": "Tools/PCM1808_L_fit.txt", "pcmR": "Tools/PCM1808_R_fit.txt"}
# Which microphone is plugged into which input pin (adcL = PF6 / ADC3, adcR = PA3 / ADC2: channel names, not left / right)
SN_ON_PF6 = "SN2"                  # microphone on PF6 (ADC3, Zio A6): the right microphone
SN_ON_PA3 = "SN1"                  # microphone on PA3 (ADC2, Zio A0): the left microphone
MICS = {"adcL": "Tools/%s to reference.txt" % SN_ON_PF6,
        "adcR": "Tools/%s to reference.txt" % SN_ON_PA3}

# --- fit settings --------------------------------------------------------------------------------
TARGET_BINS = 48                   # the REW file is averaged in bins of 1/48 octave
N_PK = 13                          # peaking sections (+ 2 low shelves + 1 high shelf = 16)
F_PH_FULL, F_PH_ZERO = 3000.0, 5000.0          # phase fitted fully up to 3 kHz, tapers to 0 at 5 kHz
W_MAG, W_PH = 2.0, 0.17            # weights: 0.5 dB and 6 deg count the same
FIR_TAPS = 128                     # delay = FIR_TAPS / 2 samples, resolution about 4 * FS / FIR_TAPS (1378 Hz)
FIR_MAX_GAIN_DB = 6.0

LINE = re.compile(r"Filter\s+(\d+):\s+(ON|OFF)\s+(\w+)(.*)")

# --- REW filter files (PCM1808) ------------------------------------------------------------------
def parse(path):
    seen, out = set(), []
    for ln in open(path, encoding="utf-8", errors="replace"):
        m = LINE.match(ln.strip())
        if not m:
            continue
        idx, on, typ, rest = int(m.group(1)), m.group(2), m.group(3), m.group(4)
        if idx in seen:            # a second equaliser block starts: stop
            break
        seen.add(idx)
        if on != "ON" or typ == "None":
            continue
        fc = re.search(r"Fc\s+([\d.]+)", rest)
        g = re.search(r"Gain\s+(-?[\d.]+)", rest)
        q = re.search(r"Q\s+([\d.]+)", rest)
        if typ not in ("PK", "LS", "HS"):
            sys.exit(f"{path}: filter {idx}: type {typ} not supported (PK, LS, HS only)")
        out.append((typ, float(fc.group(1)), float(g.group(1)), float(q.group(1)) if q else SHELF_Q))
    return out

def biquad(typ, fc, gain_db, q):
    """[b0, b1, b2, a1, a2], a0 = 1"""
    A = 10.0 ** (gain_db / 40.0)
    w0 = 2.0 * math.pi * fc / FS
    c, s = math.cos(w0), math.sin(w0)
    al = s / (2.0 * q)
    if typ == "PK":
        b = [1 + al * A, -2 * c, 1 - al * A]
        a = [1 + al / A, -2 * c, 1 - al / A]
    else:
        t = 2.0 * math.sqrt(A) * al
        if typ == "LS":
            b = [A * ((A + 1) - (A - 1) * c + t), 2 * A * ((A - 1) - (A + 1) * c), A * ((A + 1) - (A - 1) * c - t)]
            a = [(A + 1) + (A - 1) * c + t, -2 * ((A - 1) + (A + 1) * c), (A + 1) + (A - 1) * c - t]
        else:
            b = [A * ((A + 1) + (A - 1) * c + t), -2 * A * ((A - 1) + (A + 1) * c), A * ((A + 1) + (A - 1) * c - t)]
            a = [(A + 1) - (A - 1) * c + t, 2 * ((A - 1) - (A + 1) * c), (A + 1) - (A - 1) * c - t]
    return [b[0] / a[0], b[1] / a[0], b[2] / a[0], a[1] / a[0], a[2] / a[0]]

def chain(secs, f):
    """complex response of biquads [b0, b1, b2, a1, a2] at the frequencies f"""
    z = np.exp(-2j * np.pi * np.asarray(f, float) / FS)
    h = np.ones_like(z)
    for b0, b1, b2, a1, a2 in secs:
        h = h * (b0 + b1 * z + b2 * z * z) / (1 + a1 * z + a2 * z * z)
    return h

def params_response(params, f):
    return chain([biquad(*p) for p in params], f)

# --- target: the REW file -------------------------------------------------------------------------
def load_target(rel):
    """complex response of the REW file, 20 Hz .. 20 kHz, averaged in 1/TARGET_BINS octave bins"""
    d = np.loadtxt(os.path.join(ROOT, rel), comments="*")
    b = (d[:, 0] >= 20.0) & (d[:, 0] <= 20000.0)
    f, m, ph = d[b, 0], d[b, 1], d[b, 2]
    H = 10 ** (m / 20) * np.exp(1j * np.radians(ph))
    k = np.floor(np.log2(f) * TARGET_BINS).astype(int)
    keys, inv = np.unique(k, return_inverse=True)
    n = np.bincount(inv)
    return 2.0 ** ((keys + 0.5) / TARGET_BINS), np.bincount(inv, H.real) / n + 1j * np.bincount(inv, H.imag) / n

# --- biquad fit ------------------------------------------------------------------------------------
# parameter vector: [LS1 fc, LS1 g, LS2 fc, LS2 g, HS fc, HS g] + [PK fc, g, Q] * N_PK
def unpack(p):
    secs = [("LS", p[0], p[1], SHELF_Q), ("LS", p[2], p[3], SHELF_Q), ("HS", p[4], p[5], SHELF_Q)]
    for k in range(N_PK):
        secs.append(("PK", p[6 + 3 * k], p[7 + 3 * k], p[8 + 3 * k]))
    return secs

def residuals(p, fg, T, wp):
    C = params_response(unpack(p), fg)
    rm = W_MAG * (20 * np.log10(np.abs(C)) - 20 * np.log10(np.abs(T)))
    rp = W_PH * wp * np.degrees(np.angle(C / T))
    pen = np.concatenate([[0.012 * p[1], 0.012 * p[3], 0.012 * p[5]], 0.12 * p[7::3]])   # gentle preference for small gains
    return np.concatenate([rm, rp, pen])

def fit_biquads(fg, T):
    wp = np.clip((F_PH_ZERO - fg) / (F_PH_ZERO - F_PH_FULL), 0.0, 1.0)
    p0 = [14.0, 0.0, 25.0, 0.0, 12000.0, 0.0]
    lo = [4.0, -40.0, 6.0, -30.0, 3000.0, -20.0]
    hi = [60.0, 20.0, 150.0, 30.0, 20000.0, 20.0]
    for fc in np.geomspace(35.0, 17000.0, N_PK):
        p0 += [fc, 0.0, 1.5]; lo += [fc / 1.6, -15.0, 0.4]; hi += [fc * 1.6, 15.0, 8.0]
    p = np.array(p0, float)
    for it in range(3):                      # the PK gains start from what is still missing
        err = 20 * np.log10(np.abs(T)) - 20 * np.log10(np.abs(params_response(unpack(p), fg)))
        for k in range(N_PK):
            p[7 + 3 * k] = float(np.clip(0.6 * np.interp(p[6 + 3 * k], fg, err), -12, 12)) + (p[7 + 3 * k] if it else 0.0)
        p = least_squares(residuals, p, bounds=(lo, hi), args=(fg, T, wp), max_nfev=400, x_scale="jac").x
    secs = sorted(unpack(p), key=lambda s: {"LS": 0, "PK": 1, "HS": 2}[s[0]] * 1e6 + s[1])
    out = []
    for typ, fc, g, q in secs:             # rounded like a REW filter file (Fc 3, Gain 4, Q 4 decimals); tiny peaks dropped
        if typ == "PK" and abs(g) < 0.02:
            continue
        out.append((typ, float(f"{fc:.3f}"), float(f"{g:.4f}"), float(f"{q:.4f}") if typ == "PK" else SHELF_Q))
    return out

# --- FIR: what the biquads leave in phase and magnitude -------------------------------------------
def design_fir(f, T, secs):
    taps, nfft = FIR_TAPS, 8192
    C = chain([biquad(*s) for s in secs], f)
    k = np.argmin(abs(f - 1000.0))
    phT = np.unwrap(np.angle(T)); phT -= np.round(phT[k] / (2 * np.pi)) * 2 * np.pi
    phC = np.unwrap(np.angle(C));  phC -= np.round(phC[k] / (2 * np.pi)) * 2 * np.pi
    resid = phC - phT                                              # G = exp(-j resid): the phase of the target minus that of the biquads
    lg = np.clip(np.log(np.abs(T)) - np.log(np.abs(C)), -np.log(10 ** (FIR_MAX_GAIN_DB / 20)), np.log(10 ** (FIR_MAX_GAIN_DB / 20)))
    fg = np.arange(nfft // 2 + 1) * FS / nfft
    fmin = 4.0 * FS / taps
    fade = np.clip((fg - 0.5 * fmin) / (0.5 * fmin), 0.0, 1.0)     # nothing below about 345 Hz, full from 689 Hz
    pg = np.interp(fg, f, resid, left=resid[0], right=resid[-1]) * fade
    mag = np.exp(np.interp(fg, f, lg, left=lg[0], right=lg[-1]) * fade)
    D = taps // 2
    G = mag * np.exp(-1j * pg) * np.exp(-2j * np.pi * fg * D / FS)
    g = np.fft.irfft(G, nfft)[:taps]
    win = np.ones(taps); r = taps // 4
    win[:r] = 0.5 * (1 - np.cos(np.pi * np.arange(r) / r)); win[-r:] = win[:r][::-1]
    return (g * win).astype(np.float32)

def fir_response(g, f):
    n = len(g)
    return (np.exp(-2j * np.pi * np.outer(f, np.arange(n)) / FS) @ g.astype(np.float64)) * np.exp(2j * np.pi * f * (n // 2) / FS)   # delay removed

# --- check against the raw files -----------------------------------------------------------------
def check(tag, rel, secs, g):
    d = np.loadtxt(os.path.join(ROOT, rel), comments="*")
    b = (d[:, 0] >= 20) & (d[:, 0] <= 20000)
    f = d[b, 0][::4]; T = (10 ** (d[b, 1] / 20) * np.exp(1j * np.radians(d[b, 2])))[::4]
    C = chain([biquad(*s) for s in secs], f) * fir_response(g, f)
    dm = 20 * np.log10(np.abs(C)) - 20 * np.log10(np.abs(T)); dp = np.degrees(np.angle(C / T))
    print(f"   {tag}: biquads x FIR minus the file {rel}   (rms / 95 % of the points)")
    for lo, hi in ((20, 50), (50, 200), (200, 1000), (1000, 3000), (3000, 10000), (10000, 20000)):
        s = (f >= lo) & (f < hi)
        print("      %5d..%-6d Hz  magnitude %5.2f / %5.2f dB   phase %5.1f / %5.1f deg" % (lo, hi, np.sqrt(np.mean(dm[s] ** 2)), np.percentile(np.abs(dm[s]), 95),
              np.sqrt(np.mean(dp[s] ** 2)), np.percentile(np.abs(dp[s]), 95)))

def main():
    do_check = "--check" in sys.argv
    lines = ["/* Generated by Tools/EQ_to_BIQUAD.py from the REW files in Tools/ - do not edit by hand. */",
             "#ifndef MIC_EQ_COEFFS_H", "#define MIC_EQ_COEFFS_H", "",
             f"#define MIC_EQ_MAX_SECTIONS {MAX_SECTIONS}U", "",
             "typedef struct { double b0, b1, b2, a1, a2; } MicEqBiquad;   /* a0 = 1 */", ""]
    freqs = [30, 50, 100, 200, 500, 1000, 2000, 3000, 5000, 8000, 10000, 15000, 20000]
    chains, firs = {}, {}
    for name, path in PCM.items():
        chains[name] = [biquad(*f) for f in parse(os.path.join(ROOT, path))]
        print(f"{name}: {path}: {len(chains[name])} sections")
    for name, path in MICS.items():
        f, T = load_target(path)
        secs = fit_biquads(f, T)
        firs[name] = design_fir(f, T, secs)
        chains[name] = [biquad(*s) for s in secs]
        resp = chain(chains[name], f)
        print(f"{name}: {path}: {len(secs)} biquads (fit to the file: magnitude rms {np.sqrt(np.mean((20*np.log10(abs(resp))-20*np.log10(abs(T)))**2)):.2f} dB), "
              f"FIR {len(firs[name])} taps (gain {20*np.log10(abs(fir_response(firs[name], f))).min():+.1f} .. {20*np.log10(abs(fir_response(firs[name], f))).max():+.1f} dB)")
        if do_check:
            check(name, path, secs, firs[name])
    for name in list(PCM) + list(MICS):
        secs = chains[name]
        if len(secs) > MAX_SECTIONS:
            sys.exit(f"{name}: {len(secs)} filters, max {MAX_SECTIONS}")
        rmax = max([max(abs(np.roots([1.0, s[3], s[4]]))) for s in secs] or [0.0])
        if rmax >= 1.0:
            sys.exit(f"{name}: unstable biquad (pole radius {rmax})")
        print(f"   {name}: " + "  ".join(f"{k}Hz {20*np.log10(abs(chain(secs, [k])[0])):+.1f}" for k in freqs) + f"   (largest pole radius {rmax:.5f})")
        lines.append(f"#define MIC_EQ_{name.upper()}_N {len(secs)}U")
        lines.append(f"static const MicEqBiquad micEq_{name}[{max(len(secs), 1)}] = {{")
        for s in (secs or [[1.0, 0.0, 0.0, 0.0, 0.0]]):
            lines.append("    { " + ", ".join(f"{v:.17g}" for v in s) + " },")
        lines.append("};")
        lines.append("")
    lines.append("/* Phase / magnitude FIR of the microphones, after the biquads (taps = N, delay = N / 2 samples), used by mic_phase.c */")
    for name, g in firs.items():
        lines.append(f"#define MIC_PH_{name.upper()}_N {len(g)}U")
        lines.append(f"static const float micPh_{name}[{len(g)}] = {{")
        for i in range(0, len(g), 4):
            lines.append("    " + ", ".join(f"{float(v):.9e}f" for v in g[i:i + 4]) + ",")
        lines.append("};")
        lines.append("")
    lines.append("#endif /* MIC_EQ_COEFFS_H */")
    open(OUT, "w").write("\n".join(lines) + "\n")
    print("wrote", OUT)

if __name__ == "__main__":
    main()
