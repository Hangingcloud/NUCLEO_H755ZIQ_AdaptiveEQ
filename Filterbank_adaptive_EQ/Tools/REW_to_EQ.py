#!/usr/bin/env python3
"""Fit PK / LS / HS biquads to the inverse of a measured REW response and write them as a REW "Filter Settings"
file (input of Tools/EQ_to_BIQUAD.py).

Usage (project root):  python3 Tools/REW_to_EQ.py MEASUREMENT.txt OUT.txt [peaking sections] [options]
  --refine OLD_FIT.txt  the measurement was taken with OLD_FIT active: new correction = old correction - measured error
  --lf-pole CAP_DB      below the measured range continue the correction as a one-pole high-pass, capped at CAP_DB
  --adjust F1-F2:dB     add dB to the correction from F1 upwards (half-octave raised-cosine ramp below F1), repeatable
Input is a REW text export (Freq, SPL, Phase), normalised to its level at 1 kHz.
"""
import math, re, sys
import numpy as np
from scipy.optimize import least_squares

FS = 44100.0

def load(path):
    f, d = [], []
    for ln in open(path, encoding="utf-8", errors="replace"):
        if ln.startswith("*") or not ln.strip():
            continue
        p = ln.split()
        f.append(float(p[0])); d.append(float(p[1]))
    f, d = np.array(f), np.array(d)
    return f, d - np.interp(1000.0, f, d)           # dB relative to 1 kHz

def sections(typ, fc, g, q):
    """RBJ biquad as (b, a) arrays, a0 = 1 (same definitions as EQ_to_BIQUAD.py)"""
    A = 10.0 ** (g / 40.0); w0 = 2 * math.pi * fc / FS
    c, s = math.cos(w0), math.sin(w0); al = s / (2 * q)
    if typ == "PK":
        b = [1 + al * A, -2 * c, 1 - al * A]; a = [1 + al / A, -2 * c, 1 - al / A]
    else:
        t = 2 * math.sqrt(A) * al
        if typ == "LS":
            b = [A * ((A + 1) - (A - 1) * c + t), 2 * A * ((A - 1) - (A + 1) * c), A * ((A + 1) - (A - 1) * c - t)]
            a = [(A + 1) + (A - 1) * c + t, -2 * ((A - 1) + (A + 1) * c), (A + 1) + (A - 1) * c - t]
        else:
            b = [A * ((A + 1) + (A - 1) * c + t), -2 * A * ((A - 1) + (A + 1) * c), A * ((A + 1) + (A - 1) * c - t)]
            a = [(A + 1) - (A - 1) * c + t, 2 * ((A - 1) - (A + 1) * c), (A + 1) - (A - 1) * c - t]
    return np.array(b) / a[0], np.array(a) / a[0]

def resp_db(freqs, filt):
    z = np.exp(-2j * np.pi * freqs / FS); h = np.ones_like(z)
    for typ, fc, g, q in filt:
        b, a = sections(typ, fc, g, q)
        h *= (b[0] + b[1] * z + b[2] * z * z) / (1 + a[1] * z + a[2] * z * z)
    return 20 * np.log10(np.abs(h))

def fit(f, target, types, start, lo, hi, n_meas=None, starts=60, seed=1):
    """Multi-start least squares, then Lawson re-weighting towards the minimax solution (smallest worst-case error).
    The first len(f)-n_meas points extend below the measured range (half weight, not counted)."""
    start, lo, hi = np.array(start, float), np.array(lo, float), np.array(hi, float)
    nx = 0 if n_meas is None else len(f) - n_meas
    base_w = np.ones(len(f)); base_w[:nx] = 0.5
    rng = np.random.default_rng(seed)

    def unpack(x):
        return [(t, x[3 * i], x[3 * i + 1], x[3 * i + 2]) for i, t in enumerate(types)]

    def run(x0, w):
        r = least_squares(lambda x: w * (resp_db(f, unpack(x)) - target), np.clip(x0, lo, hi), bounds=(lo, hi),
                          xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=400)
        return r.x

    def worst(x):
        e = resp_db(f, unpack(x)) - target
        return np.max(np.abs(e[nx:]))

    best, best_w = None, 1e9
    for k in range(starts):
        x0 = start.copy() if k == 0 else np.where(                    # jitter inside the bounds (log scale for fc, Q)
            np.tile([True, False, True], len(types)),
            np.exp(rng.uniform(np.log(np.maximum(lo, 1e-3)), np.log(hi))),
            rng.uniform(lo, hi))
        try:
            x = run(x0, base_w)
        except Exception:
            continue
        w_ = worst(x)
        if w_ < best_w: best, best_w = x, w_
    # Lawson refinement
    w = base_w.copy(); x = best
    for it in range(40):
        x = run(x, w)
        e = np.abs(resp_db(f, unpack(x)) - target)
        w = w * (0.5 + e / max(np.max(e[nx:]), 1e-9)); w = w / np.mean(w) * np.mean(base_w); w[:nx] *= 0.5
        if worst(x) < best_w: best, best_w = x, worst(x)
    return unpack(best), resp_db(f, unpack(best)) - target

def load_filters(path):
    out = []
    for ln in open(path, encoding="utf-8", errors="replace"):
        m = re.match(r"Filter\s+\d+:\s+ON\s+(PK|LS|HS)\s+Fc\s+([\d.]+)\s*Hz\s+Gain\s+(-?[\d.]+)\s*dB\s+Q\s+([\d.]+)", ln)
        if m:
            out.append((m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))))
    return out

def adjustment(f, f1, f2, db):
    """db from f1 upwards (f2 is only used for the log text), raised-cosine ramp over the half octave below f1"""
    lo = f1 / math.sqrt(2.0)
    w = np.zeros_like(f)
    w[f >= f1] = 1.0
    r = (f >= lo) & (f < f1)
    w[r] = 0.5 - 0.5 * np.cos(np.pi * (np.log(f[r]) - np.log(lo)) / (np.log(f1) - np.log(lo)))
    return db * w

def main():
    args = sys.argv[1:]
    adj = []
    while "--adjust" in args:
        i = args.index("--adjust")
        band, db = args[i + 1].split(":")
        f1, f2 = (float(x) for x in band.split("-"))
        adj.append((f1, f2, float(db)))
        del args[i:i + 2]
    old_filt = None
    if "--refine" in args:
        i = args.index("--refine"); old_filt = load_filters(args[i + 1]); del args[i:i + 2]
        if not old_filt: sys.exit("--refine: no filters found in that file")
    lf_cap = None
    if "--lf-pole" in args:
        i = args.index("--lf-pole"); lf_cap = float(args[i + 1]); del args[i:i + 2]
    src, dst = args[0], args[1]
    n_pk = int(args[2]) if len(args) > 2 else 1
    f, dev = load(src)
    target = -dev                                     # what the filters must do
    if old_filt:                                      # measured with the old filters on
        target = resp_db(f, old_filt) - dev
    for f1, f2, db in adj:
        target = target + adjustment(f, f1, f2, db)
    # below the measured range the correction is held at its 20 Hz value (a one-pole high-pass would ask for unlimited boost)
    fx = np.geomspace(5.0, f[0] * 0.98, 24)
    if old_filt:
        tx = resp_db(fx, old_filt)                     # keep the old correction
    elif lf_cap is None:
        tx = np.full_like(fx, target[0])
    else:                                              # continue as the one-pole the data follows, capped
        m = (f >= 20) & (f <= 80)
        pole = lambda fc, ff: 10 * np.log10(1 + (fc / ff) ** 2)
        fc = least_squares(lambda p: pole(p[0], f[m]) - target[m], [27.0]).x[0]
        tx = np.minimum(pole(fc, fx), lf_cap)
        print(f"one-pole fit of the droop: corner {fc:.2f} Hz")
    f_all = np.concatenate([fx, f]); t_all = np.concatenate([tx, target])
    # LS (bass), n_pk peaking sections, HS (treble)
    types = ["LS"] + ["PK"] * n_pk + ["HS"]
    start = [(30, 3, 0.7)] + [(60 * 3 ** i, 0.5, 1.0) for i in range(n_pk)] + [(9000, 0.2, 0.7)]
    lo = [(10, -2, 0.2)] + [(15, -6, 0.3)] * n_pk + [(2000, -3, 0.2)]
    hi = [(300, max(6.0, lf_cap or 0.0) + 0.5, 2.0)] + [(18000, 6, 6.0)] * n_pk + [(20000, 3, 2.0)]
    flat = lambda L: [v for t in L for v in t]
    filt, err_all = fit(f_all, t_all, types, flat(start), flat(lo), flat(hi), n_meas=len(f))
    err = err_all[len(fx):]                           # report the error on the measured range only
    note = (", refined from the measurement with the previous filters on" if old_filt else "") + "".join(f", adjust {db:+.2f} dB at {f1:.0f}-{f2:.0f} Hz" for f1, f2, db in adj)
    lines = ["Filter Settings file", "", "Generated by Tools/REW_to_EQ.py from " + src + note, "", "Equaliser: Generic"]
    for i, (t, fc, g, q) in enumerate(filt, 1):
        lines.append(f"Filter {i:2d}: ON  {t}  Fc {fc:10.4f} Hz  Gain {g:9.4f} dB  Q {q:7.4f}")
    open(dst, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines[5:]))
    print(f"max error {np.max(np.abs(err)):.3f} dB, rms {np.sqrt(np.mean(err**2)):.3f} dB (over {f[0]:.0f}-{f[-1]:.0f} Hz)")
    print("  correction below the measured range:", "  ".join(f"{x:.0f}Hz {resp_db(np.array([x]), filt)[0]:+.2f}" for x in (5, 10, 15)))
    print("  Hz     measured  target  correction  residual (dB)   [target = -measured + adjustment]"
          + ("   (refine: 'measured' = error with the old filters on)" if old_filt else ""))
    for fq in (20, 25, 30, 40, 50, 70, 100, 200, 500, 1000, 3000, 6000, 8000, 10000, 12000, 15000, 18000, 20000):
        i = int(np.argmin(np.abs(f - fq)))
        print(f"{f[i]:7.0f} {dev[i]:+8.3f} {target[i]:+8.3f} {resp_db(f[i:i+1], filt)[0]:+10.3f} {err[i]:+9.3f}")

main()
