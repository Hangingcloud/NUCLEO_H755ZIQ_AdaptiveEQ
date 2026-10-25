/**
 * @file    farina.c
 * @brief   Farina exponential sweep and impulse response identification, see farina.h.
 */
/* Always optimized, independent of the IDE's Debug / Release setting. */
#pragma GCC optimize ("O2")
#include "farina.h"
#include "mimo_config.h"
#include <math.h>
#include <string.h>

#define N           MIMO_SWEEP_N
#define PI_D        3.14159265358979323846

/* One complex buffer (interleaved re, im):
 *   re = x(n)       the sweep, generated again when the recordings are complete
 *   im = y_avg(n)   filled while the recordings arrive
 * x and y are real, so one complex FFT gives both spectra:
 *   X[k] = (Z[k] + conj(Z[N-k])) / 2        Y[k] = (Z[k] - conj(Z[N-k])) / (2j)
 * After the division H[k] = Y X* / (|X|^2 + reg) the inverse FFT is real. */
static float buf[2U * N] __attribute__((aligned(8)));
static float cosTab[N / 4U + 1U];            /* cos(2 pi t / N), t = 0 .. N/4 */

#define SWEEP_LEVEL (powf(10.0f, MIMO_SWEEP_LEVEL_DB / 20.0f))

/* ---------------------------------------------------------------------------------------------- */
void Farina_Init(void)
{
    for (uint32_t t = 0; t <= N / 4U; t++)
        cosTab[t] = (float)cos(2.0 * PI_D * (double)t / (double)N);
}

uint32_t Farina_BufferBytes(void) { return (uint32_t)sizeof(buf); }

/* ---------------------------------------------------------------------------------------------- */
/* Sweep generator                                                                                */
/* ---------------------------------------------------------------------------------------------- */
uint32_t Farina_SweepLen(void)
{
    return (uint32_t)((double)MIMO_SWEEP_T_S * (double)MIMO_FS_HZ + 0.5);
}

/* r = ln(w2/w1) / (T fs): per sample growth of the exponent; K = w1 T / ln(w2/w1); MULT = exp(r) */
static double sweepR, sweepK, sweepMult;
static uint32_t sweepLen, fadeLen;

void Farina_GenInit(FarinaGen *g, uint32_t n0)
{
    const double l = log((double)MIMO_SWEEP_F2_HZ / (double)MIMO_SWEEP_F1_HZ);
    sweepLen  = Farina_SweepLen();
    fadeLen   = (uint32_t)((double)MIMO_SWEEP_FADE_MS * 0.001 * (double)MIMO_FS_HZ);
    sweepR    = l / (double)sweepLen;
    sweepK    = 2.0 * PI_D * (double)MIMO_SWEEP_F1_HZ * (double)MIMO_SWEEP_T_S / l;
    sweepMult = exp(sweepR);
    g->n = n0;
    g->e = exp((double)n0 * sweepR);
}

float Farina_GenNext(FarinaGen *g)
{
    const uint32_t len = sweepLen;
    const uint32_t n   = g->n;
    float out = 0.0f;

    if (n < len)
    {
        double ph = sweepK * (g->e - 1.0);
        ph -= 2.0 * PI_D * floor(ph / (2.0 * PI_D));                 /* wrap: float sinf stays accurate */
        float v = sinf((float)ph);

        const uint32_t nf = fadeLen;
        if (n < nf)            v *= 0.5f * (1.0f - cosf((float)(PI_D * (double)n / (double)nf)));
        else if (n + nf >= len) v *= 0.5f * (1.0f - cosf((float)(PI_D * (double)(len - 1U - n) / (double)nf)));
        out = v * SWEEP_LEVEL;
    }

    g->n = n + 1U;
    if ((g->n & 1023U) == 0U) g->e = exp((double)g->n * sweepR);     /* re-sync, keeps the product error small */
    else                      g->e *= sweepMult;
    return out;
}

/* ---------------------------------------------------------------------------------------------- */
/* FFT                                                                                            */
/* ---------------------------------------------------------------------------------------------- */
static inline float twCos(uint32_t t)            /* cos(2 pi t / N), 0 <= t < N/2 */
{
    return (t <= N / 4U) ? cosTab[t] : -cosTab[N / 2U - t];
}
static inline float twSin(uint32_t t)            /* sin(2 pi t / N) = cos(2 pi (N/4 - t) / N) */
{
    return cosTab[(t <= N / 4U) ? (N / 4U - t) : (t - N / 4U)];
}

static void fft(float *a, int inverse)
{
    /* bit reversal permutation */
    for (uint32_t i = 1U, j = 0U; i < N; i++)
    {
        uint32_t bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j)
        {
            float tr = a[2U * i], ti = a[2U * i + 1U];
            a[2U * i] = a[2U * j];  a[2U * i + 1U] = a[2U * j + 1U];
            a[2U * j] = tr;         a[2U * j + 1U] = ti;
        }
    }

    const float sgn = inverse ? 1.0f : -1.0f;      /* forward: exp(-j ...) */
    for (uint32_t len = 2U; len <= N; len <<= 1)
    {
        const uint32_t half = len >> 1, step = N / len;
        for (uint32_t k = 0U; k < half; k++)
        {
            const float wr = twCos(k * step);
            const float wi = sgn * twSin(k * step);
            for (uint32_t i = k; i < N; i += len)
            {
                float *p = &a[2U * i];
                float *q = &a[2U * (i + half)];
                const float tr = q[0] * wr - q[1] * wi;
                const float ti = q[0] * wi + q[1] * wr;
                q[0] = p[0] - tr;  q[1] = p[1] - ti;
                p[0] += tr;        p[1] += ti;
            }
        }
    }
    if (inverse)
    {
        const float s = 1.0f / (float)N;
        for (uint32_t i = 0U; i < 2U * N; i++) a[i] *= s;
    }
}

/* ---------------------------------------------------------------------------------------------- */
void Farina_Clear(void)
{
    memset(buf, 0, sizeof(buf));
}

void Farina_Accumulate(uint32_t n0, const float *z, uint32_t cnt)
{
    const float w = 1.0f / (float)MIMO_SWEEP_AVG;
    for (uint32_t i = 0U; i < cnt; i++)
    {
        const uint32_t n = n0 + i;
        if (n >= N) break;
        buf[2U * n + 1U] += z[i] * w;
    }
}

/* Copy of the averaged recording y_avg as 16 bit integers, scale = largest |y| / 32767 (y = value * scale). Statistics of the second half of the record:
 * mean, rms, and the largest |y| of the whole record, in units of full scale. */
void Farina_StoreRecord(int16_t *dst, float *scale, float *tailMean, float *tailRms, float *peak)
{
    float mx = 0.0f;
    for (uint32_t n = 0U; n < N; n++) { const float a = fabsf(buf[2U * n + 1U]); if (a > mx) mx = a; }
    const float sc = (mx > 0.0f) ? mx / 32767.0f : 1.0f;
    for (uint32_t n = 0U; n < N; n++)
    {
        float v = buf[2U * n + 1U] / sc;
        v += (v >= 0.0f) ? 0.5f : -0.5f;
        dst[n] = (int16_t)v;
    }
    double sm = 0.0, sq = 0.0;
    for (uint32_t n = N / 2U; n < N; n++) { const double v = (double)buf[2U * n + 1U]; sm += v; sq += v * v; }
    const double cnt = (double)(N / 2U);
    *scale = sc;  *peak = mx;
    *tailMean = (float)(sm / cnt);
    *tailRms  = (float)sqrt(sq / cnt - (sm / cnt) * (sm / cnt));
}

void Farina_Deconvolve(void)
{
    FarinaGen g;
    Farina_GenInit(&g, 0U);
    for (uint32_t n = 0U; n < N; n++) buf[2U * n] = Farina_GenNext(&g);

    fft(buf, 0);

    /* regularization level: mean |X|^2 in the top octave of the sweep band (least energy per bin); outside of the band |X|^2 is about 0 */
    const double binHz = (double)MIMO_FS_HZ / (double)N;
    uint32_t kLo = (uint32_t)(0.5 * (double)MIMO_SWEEP_F2_HZ / binHz);
    uint32_t kHi = (uint32_t)(0.95 * (double)MIMO_SWEEP_F2_HZ / binHz);
    if (kHi > N / 2U - 1U) kHi = N / 2U - 1U;
    if (kLo >= kHi) kLo = kHi > 8U ? kHi - 8U : 1U;
    double acc = 0.0;
    for (uint32_t k = kLo; k <= kHi; k++)
    {
        const float *zk = &buf[2U * k], *zn = &buf[2U * (N - k)];
        const float xr = 0.5f * (zk[0] + zn[0]), xi = 0.5f * (zk[1] - zn[1]);
        acc += (double)xr * xr + (double)xi * xi;
    }
    const float reg = (float)((double)MIMO_SWEEP_REG * acc / (double)(kHi - kLo + 1U));

    for (uint32_t k = 0U; k <= N / 2U; k++)
    {
        const uint32_t k2 = (N - k) & (N - 1U);
        const float ar = buf[2U * k],  ai = buf[2U * k + 1U];
        const float br = buf[2U * k2], bi = buf[2U * k2 + 1U];
        const float xr = 0.5f * (ar + br), xi = 0.5f * (ai - bi);
        const float yr = 0.5f * (ai + bi), yi = -0.5f * (ar - br);
        const float den = xr * xr + xi * xi + reg;
        /* Outside of the sweep band the response is faded out (raised cosine below F1 / over 15 % above F2) */
        float bw = 1.0f;
        {
            const float fk = (float)((double)k * binHz);
            const float f1 = MIMO_SWEEP_F1_HZ, f2 = MIMO_SWEEP_F2_HZ;
            if (fk < MIMO_SWEEP_BW_LOW * f1)       bw = 0.0f;
            else if (fk < f1)         bw = 0.5f * (1.0f - cosf(3.14159265f * (fk - MIMO_SWEEP_BW_LOW * f1) / ((1.0f - MIMO_SWEEP_BW_LOW) * f1)));
            else if (fk > 1.15f * f2) bw = 0.0f;
            else if (fk > f2)         bw = 0.5f * (1.0f + cosf(3.14159265f * (fk - f2) / (0.15f * f2)));
        }
        const float hr = bw * (yr * xr + yi * xi) / den;      /* Y conj(X) */
        const float hi = bw * (yi * xr - yr * xi) / den;
        if (k == k2)       /* k = 0 and k = N/2: the spectrum of a real signal is real there */
        {
            buf[2U * k] = hr;  buf[2U * k + 1U] = 0.0f;
        }
        else
        {
            buf[2U * k]  = hr;  buf[2U * k + 1U]  = hi;
            buf[2U * k2] = hr;  buf[2U * k2 + 1U] = -hi;       /* Hermitian symmetry */
        }
    }

    fft(buf, 1);
}

float Farina_Ir(int32_t lag)
{
    return buf[2U * ((uint32_t)lag & (N - 1U))];                  /* circular: a negative lag is the end of the buffer */
}

void Farina_Locate(int32_t *peakLag, float *peakVal, float *snrDb)
{
    /* The response lies around lag 0 of the circular deconvolution, the direct sound can sit at a negative lag (the end of the buffer). Search -2048 .. N/2. */
    float pk = 0.0f;
    int32_t pl = 0;
    for (int32_t l = -2048; l < (int32_t)(N / 2U); l++)
    {
        const float a = fabsf(Farina_Ir(l));
        if (a > pk) { pk = a; pl = l; }
    }
    double sq = 0.0;
    for (uint32_t n = N / 4U; n < N / 2U; n++) sq += (double)buf[2U * n] * (double)buf[2U * n];
    const double noise = sqrt(sq / (double)(N / 4U)) + 1.0e-30;
    *peakLag = pl;
    *peakVal = pk;
    *snrDb   = (float)(20.0 * log10((double)pk / noise + 1.0e-30));
}

/* Lag where the energy of the response rises: the first block of 32 samples between lag -1024 and +1536 above 3 % (-15 dB) of the strongest block. */
int32_t Farina_EnergyStart(void)
{
    float e[80];
    float mx = 0.0f;
    for (int32_t b = 0; b < 80; b++)
    {
        double s2 = 0.0;
        for (int32_t k = 0; k < 32; k++) { const float v = Farina_Ir(-1024 + 32 * b + k); s2 += (double)v * (double)v; }
        e[b] = (float)s2;
        if (e[b] > mx) mx = e[b];
    }
    for (int32_t b = 0; b < 80; b++)
        if (e[b] >= 0.03f * mx) return -1024 + 32 * b;
    return 0;
}
