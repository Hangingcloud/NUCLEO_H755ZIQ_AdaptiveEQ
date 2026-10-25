/**
 * @file    room_eq.c
 * @brief   Adaptive multiple position room response equalizer, see room_eq.h.
 *
 * Plain C99, no hardware access.
 * It is built with -O2 even in the Debug configuration (the CubeIDE Debug build uses -O0).
 */
#pragma GCC optimize ("O2")

#include "room_eq.h"
#include <math.h>
#include <string.h>

#ifdef __arm__
#define RE_DTCM  __attribute__((section(".dtcm"), aligned(8)))   /* zero wait state RAM, CPU only */
#else
#define RE_DTCM
#endif

#define K        ROOMEQ_K
#define K2       (ROOMEQ_K / 2U)           /* complex FFT length, number of bins - 1 */
#define M        ROOMEQ_BANDS
#define R        ROOMEQ_MICS
#define L        ROOMEQ_FRAMES_AVG
#define NF       ROOMEQ_MAX_FILTERS
#define RING     8192U                     /* per stream, power of two */
#define STREAMS  (1U + ROOMEQ_MICS)        /* 0 = reference, 1.. = mics */
#define P        96U                       /* design grid points (log spaced, 20 Hz .. 20 kHz) */
#define PI_D     3.14159265358979323846

/* ---- band centres [Hz] and step sizes ---- */
static const double bandCenter[M] = {
#if ROOMEQ_WITH_25HZ
                                      25,
#endif
                                      47,
#if ROOMEQ_WITH_84HZ
                                      84,                                     /* closes the 1.6 octave gap between 47 and 140 Hz */
#endif
                                      140, 234, 352, 492, 633, 820, 1008, 1242, 1523,
                                      1875, 2320, 2930, 3680, 4805, 6352, 8507, 11414, 15210, 20883 };
static const float  bandStep[M]   = {
#if ROOMEQ_WITH_25HZ
                                      0.05f,                                  /* 25 Hz */
#endif
                                      0.05f,                                  /* 47 Hz */
#if ROOMEQ_WITH_84HZ
                                      0.05f,                                  /* 84 Hz */
#endif
                                      0.05f,                                  /* 140 Hz */
                                      0.05f, 0.05f,                           /* 234, 352 Hz */
                                      0.05f, 0.05f,                           /* 492, 633 Hz */
                                      0.05f, 0.05f,                           /* 820, 1008 Hz */
                                      0.05f, 0.05f,                           /* 1242, 1523 Hz */
                                      0.01f, 0.01f,                           /* 1875, 2320 Hz */
                                      0.01f, 0.01f,                           /* 2930, 3680 Hz */
                                      0.05f, 0.05f,                           /* 4805, 6352 Hz */
                                      0.05f, 0.05f,                           /* 8507, 11414 Hz */
                                      0.05f, 0.05f };                         /* 15210, 20883 Hz */

/* ---- FFT tables / work area (DTCM) ---- */
RE_DTCM static float    tabCos[K2 + 1U], tabSin[K2 + 1U];     /* cos/sin(2 pi i / K), i = 0..K/2 */
RE_DTCM static uint16_t bitrev[K2];
RE_DTCM static float    fftRe[K2], fftIm[K2];

/* ---- mic frequency response correction (linear magnitude per bin; applied in the analysis) ---- */
static float        micCorr[R][K2 + 1U];
static uint8_t      micCorrSet[R];
static volatile int micCorrOn = 1;

/* ---- input rings (AXI SRAM is fine, they are read sequentially) ---- */
static int16_t           ring[STREAMS][RING];
static volatile uint32_t total[STREAMS];                      /* samples written so far per stream */

/* ---- bands ---- */
static double   bandEdge[M + 1U];                             /* bin mapping edges */
static double   smoothEdge[M + 1U];                           /* edges used for the smoothing (extended to 0 and infinity) */
static uint32_t binLo[M], binHi[M];
static uint8_t  bandInRange[M];
static float    mu[M];

/* ---- adaptation state ---- */
static float    weightDb[R][M];
static float    errLin[R][M][L];
static uint32_t histCount, histPos;
static uint32_t frameStart;
static uint32_t framesDone;
static volatile int adaptOn = 1;
static volatile int eqActive = 1;                              /* 0 = equalizer bypassed: loop open, no adaptation */
static float    lastRefRms, lastMicRms[R];
static double   micFloor[R][M];                               /* noise floor per band (linear magnitude) */
static uint8_t  floorKnown[R];
static uint8_t  bandValid[R][M];
static float    bandSnrDb[R][M], refRelDb[M];
static float    lastErrRms[R];                                /* rms of the normalized error over the used bands, last frame (dB) */
static uint32_t lastValidCnt[R];

/* ---- design ---- */
static double   gridF[P], gridCos[P], gridCos2[P], gridW[P];
#if ROOMEQ_COMPLEX_INVERSE
/* Complex inverse of the prototype for an FFT / FIR realisation of the equalizer; the peaking filter bank only needs the magnitude.
 * Missing for the FFT / FIR route: a phase for Hp (e.g. minimum phase of the magnitude) and a linear frequency grid. */
static double   hpPhase[P];                                    /* phase of the prototype Hp on the grid (rad); 0 = no phase information */
static double   hinvRe[P], hinvIm[P];                          /* complex inverse */
#endif
RE_DTCM static double unitShape[NF][P];                        /* dB response per dB of gain of each filter */
static double   cholL[NF][NF];
static double   fFreq[NF], fCosW[NF], fSinW[NF], fQ[NF];
static double   fGainDb[NF];
static RoomEqCoefs coefBuf[2];
static volatile uint32_t coefCur, coefSeq;

/* ====================================================================== */
/*  FFT                                                                   */
/* ====================================================================== */
static void fftInit(void)
{
    for (uint32_t i = 0; i <= K2; i++)
    {
        double a = 2.0 * PI_D * (double)i / (double)K;
        tabCos[i] = (float)cos(a);
        tabSin[i] = (float)sin(a);
    }
    for (uint32_t n = 0; n < K2; n++)
    {
        uint32_t r = 0, v = n;
        for (uint32_t b = 0; b < 11U; b++) { r = (r << 1) | (v & 1U); v >>= 1; }
        bitrev[n] = (uint16_t)r;
    }
}

/* in place complex FFT of length 2048 on fftRe/fftIm (input already in bit reversed order) */
static void fft2048(void)
{
    for (uint32_t len = 2U; len <= K2; len <<= 1)
    {
        const uint32_t half = len >> 1, stride = K2 / len;
        for (uint32_t i = 0; i < K2; i += len)
        {
            for (uint32_t j = 0; j < half; j++)
            {
                const uint32_t idx = 2U * j * stride;
                const float wr = tabCos[idx], wi = -tabSin[idx];
                const uint32_t a = i + j, b = a + half;
                const float tr = wr * fftRe[b] - wi * fftIm[b];
                const float ti = wr * fftIm[b] + wi * fftRe[b];
                fftRe[b] = fftRe[a] - tr;  fftIm[b] = fftIm[a] - ti;
                fftRe[a] += tr;            fftIm[a] += ti;
            }
        }
    }
}

/* Hann window of the 4096 point frame (periodic) */
static inline float hann(uint32_t i)
{
    const float c = (i <= K2) ? tabCos[i] : tabCos[K - i];
    return 0.5f - 0.5f * c;
}

/* windowed real FFT of K samples of stream s starting at sample `start`; returns the energy sum(x^2) */
static double loadAndTransform(const int16_t *rg, uint32_t start)
{
    double energy = 0.0;
    for (uint32_t n = 0; n < K2; n++)
    {
        const float a = (float)rg[(start + 2U * n) & (RING - 1U)];
        const float b = (float)rg[(start + 2U * n + 1U) & (RING - 1U)];
        energy += (double)(a * a + b * b);
        const uint32_t r = bitrev[n];
        fftRe[r] = a * hann(2U * n);
        fftIm[r] = b * hann(2U * n + 1U);
    }
    fft2048();
    return energy;
}

/* magnitude |X(k)| of the real signal from the half length complex transform Z */
static inline float binMagnitude(uint32_t k)
{
    const uint32_t kk = (k == K2) ? 0U : k;
    const uint32_t nn = (K2 - k) & (K2 - 1U);
    const float zr = fftRe[kk], zi = fftIm[kk], nr = fftRe[nn], ni = fftIm[nn];
    const float er = 0.5f * (zr + nr), ei = 0.5f * (zi - ni);                 /* even part */
    const float tr = zr - nr, ti = zi + ni;
    const float orr = 0.5f * ti, oi = -0.5f * tr;                             /* odd part = t / (2j) */
    const float c = tabCos[k], s = tabSin[k];
    const float xr = er + (c * orr + s * oi);
    const float xi = ei + (c * oi - s * orr);
    return sqrtf(xr * xr + xi * xi);
}

/* average magnitude per band (linear); corr = optional magnitude correction per bin (mic frequency response) */
static void bandMeans(double *mean, const float *corr)
{
    for (uint32_t m = 0; m < M; m++)
    {
        double sum = 0.0;
        if (corr)
            for (uint32_t k = binLo[m]; k <= binHi[m]; k++) sum += (double)(binMagnitude(k) * corr[k]);
        else
            for (uint32_t k = binLo[m]; k <= binHi[m]; k++) sum += (double)binMagnitude(k);
        mean[m] = sum / (double)(binHi[m] - binLo[m] + 1U);
    }
}

/* ====================================================================== */
/*  Peaking filter design helpers                                         */
/* ====================================================================== */
typedef struct { double nc, n1, n2, dc, d1, d2; } PkTerms;      /* |H|^2 = (nc + n1 cos w + n2 cos 2w) / (dc + d1 cos w + d2 cos 2w) */

static void pkCoefs(double gainDb, double cosw, double sinw, double q, RoomEqBiquad *c)
{
    const double A = pow(10.0, gainDb / 40.0), al = sinw / (2.0 * q);
    const double a0 = 1.0 + al / A;
    c->b0 = (1.0 + al * A) / a0;
    c->b1 = (-2.0 * cosw) / a0;
    c->b2 = (1.0 - al * A) / a0;
    c->a1 = (-2.0 * cosw) / a0;
    c->a2 = (1.0 - al / A) / a0;
}

static void pkTerms(const RoomEqBiquad *c, PkTerms *t)
{
    t->nc = c->b0 * c->b0 + c->b1 * c->b1 + c->b2 * c->b2;
    t->n1 = 2.0 * (c->b0 * c->b1 + c->b1 * c->b2);
    t->n2 = 2.0 * c->b0 * c->b2;
    t->dc = 1.0 + c->a1 * c->a1 + c->a2 * c->a2;
    t->d1 = 2.0 * (c->a1 + c->a1 * c->a2);
    t->d2 = 2.0 * c->a2;
}

/* model response in dB on the design grid for the given filter gains */
static void modelDb(const double *gdb, double *out)
{
    PkTerms t[NF];
    for (uint32_t i = 0; i < NF; i++)
    {
        RoomEqBiquad c;
        pkCoefs(gdb[i], fCosW[i], fSinW[i], fQ[i], &c);
        pkTerms(&c, &t[i]);
    }
    for (uint32_t p = 0; p < P; p++)
    {
        double prod = 1.0;
        for (uint32_t i = 0; i < NF; i++)
        {
            const double num = t[i].nc + t[i].n1 * gridCos[p] + t[i].n2 * gridCos2[p];
            const double den = t[i].dc + t[i].d1 * gridCos[p] + t[i].d2 * gridCos2[p];
            prod *= num / den;
        }
        out[p] = 10.0 * log10(prod);
    }
}

/* one peaking filter per band centre; Q from the spacing, unit shapes and the normal matrix */
static void prepareFilters(void)
{
    const uint32_t n = NF;
    for (uint32_t i = 0; i < n; i++) fFreq[i] = bandCenter[i];

    for (uint32_t i = 0; i < n; i++)
    {
        const double fp = (i > 0U) ? fFreq[i - 1U] : fFreq[i] * fFreq[i] / fFreq[i + 1U];
        const double fn = (i + 1U < n) ? fFreq[i + 1U] : fFreq[i] * fFreq[i] / fFreq[i - 1U];
        const double bwOct = ROOMEQ_Q_BW_FACTOR * 0.5 * log2(fn / fp);
        double q = 1.0 / (2.0 * sinh(0.5 * log(2.0) * bwOct));
        if (q < 0.4) q = 0.4;
        if (q > 8.0) q = 8.0;
        fQ[i] = q;
        const double w0 = 2.0 * PI_D * fFreq[i] / ROOMEQ_FS_HZ;
        fCosW[i] = cos(w0);
        fSinW[i] = sin(w0);
        fGainDb[i] = 0.0;
    }

    /* unit shapes: response of one filter alone at +6 dB, per dB */
    for (uint32_t i = 0; i < n; i++)
    {
        RoomEqBiquad c; PkTerms t;
        pkCoefs(6.0, fCosW[i], fSinW[i], fQ[i], &c);
        pkTerms(&c, &t);
        for (uint32_t p = 0; p < P; p++)
        {
            const double num = t.nc + t.n1 * gridCos[p] + t.n2 * gridCos2[p];
            const double den = t.dc + t.d1 * gridCos[p] + t.d2 * gridCos2[p];
            unitShape[i][p] = 10.0 * log10(num / den) / 6.0;
        }
    }

    /* normal matrix A^T W A + lambda I and its Cholesky factor */
    double Nm[NF][NF], trace = 0.0;
    for (uint32_t i = 0; i < n; i++)
        for (uint32_t j = 0; j < n; j++)
        {
            double s = 0.0;
            for (uint32_t p = 0; p < P; p++) s += gridW[p] * unitShape[i][p] * unitShape[j][p];
            Nm[i][j] = s;
            if (i == j) trace += s;
        }
    const double lambda = 0.01 * trace / (double)n;
    for (uint32_t i = 0; i < n; i++) Nm[i][i] += lambda;
    for (uint32_t i = 0; i < n; i++)
        for (uint32_t j = 0; j <= i; j++)
        {
            double s = Nm[i][j];
            for (uint32_t k = 0; k < j; k++) s -= cholL[i][k] * cholL[j][k];
            if (i == j) cholL[i][i] = sqrt(s > 1e-12 ? s : 1e-12);
            else        cholL[i][j] = s / cholL[j][j];
        }
}

static void cholSolve(const double *rhs, double *x, uint32_t n)
{
    double y[NF];
    for (uint32_t i = 0; i < n; i++)
    {
        double s = rhs[i];
        for (uint32_t k = 0; k < i; k++) s -= cholL[i][k] * y[k];
        y[i] = s / cholL[i][i];
    }
    for (int32_t i = (int32_t)n - 1; i >= 0; i--)
    {
        double s = y[i];
        for (uint32_t k = (uint32_t)i + 1U; k < n; k++) s -= cholL[k][i] * x[k];
        x[i] = s / cholL[i][i];
    }
}

/* fit the filter gains to the target (dB on the grid) and publish the coefficients */
static void fitAndPublish(const double *targetDb, double *fitOut)
{
    double model[P], res[P], rhs[NF], delta[NF];
    for (int it = 0; it < 4; it++)
    {
        modelDb(fGainDb, model);
        for (uint32_t p = 0; p < P; p++) res[p] = targetDb[p] - model[p];
        for (uint32_t i = 0; i < NF; i++)
        {
            double s = 0.0;
            for (uint32_t p = 0; p < P; p++) s += gridW[p] * unitShape[i][p] * res[p];
            rhs[i] = s;
        }
        cholSolve(rhs, delta, NF);
        for (uint32_t i = 0; i < NF; i++)
        {
            double g = fGainDb[i] + delta[i];
            if (g >  ROOMEQ_FILTER_MAX_DB) g =  ROOMEQ_FILTER_MAX_DB;
            if (g < -ROOMEQ_FILTER_MAX_DB) g = -ROOMEQ_FILTER_MAX_DB;
            fGainDb[i] = g;
        }
    }
    if (fitOut) modelDb(fGainDb, fitOut);

    RoomEqCoefs *dst = &coefBuf[1U - coefCur];
    dst->n = NF;
    for (uint32_t i = 0; i < NF; i++) pkCoefs(fGainDb[i], fCosW[i], fSinW[i], fQ[i], &dst->bq[i]);
    coefCur = 1U - coefCur;
    coefSeq = coefSeq + 1U;
}

/* ====================================================================== */
/*  Target: smoothed, averaged, regularized inverse of the weights         */
/* ====================================================================== */
static void buildTarget(double *targetDb)
{
    double hLin[R][M];
    for (uint32_t r = 0; r < R; r++)
        for (uint32_t m = 0; m < M; m++) hLin[r][m] = pow(10.0, (double)weightDb[r][m] / 20.0);

    const double half = pow(2.0, 0.5 * ROOMEQ_SMOOTH_OCT);                 /* 1/6 octave = +-1/12 octave */
    for (uint32_t p = 0; p < P; p++)
    {
        const double lo = gridF[p] / half, hi = gridF[p] * half;
        double hp = 0.0;
        for (uint32_t r = 0; r < R; r++)
        {
            double acc = 0.0;                                              /* mean of the constant-per-band magnitude over [lo, hi] */
            for (uint32_t m = 0; m < M; m++)
            {
                const double a = (lo > smoothEdge[m]) ? lo : smoothEdge[m];
                const double b = (hi < smoothEdge[m + 1U]) ? hi : smoothEdge[m + 1U];
                if (b > a) acc += hLin[r][m] * (b - a);
            }
            hp += acc / (hi - lo);
        }
        hp /= (double)R;                                                   /* prototype: mean of the mic responses */
#if ROOMEQ_USE_REGULARIZATION
        const double beta = ROOMEQ_BETA;
#else
        const double beta = 0.0;                                           /* hp > 0: the weights are limited to +-WEIGHT_MAX_DB */
#endif
#if ROOMEQ_COMPLEX_INVERSE
        /* Hinv = conj(Hp) / (|Hp|^2 + beta), Hp = hp * exp(j*phase) */
        const double hcRe = hp * cos(hpPhase[p]), hcIm = hp * sin(hpPhase[p]);
        const double dn = hcRe * hcRe + hcIm * hcIm + beta;
        hinvRe[p] = hcRe / dn;
        hinvIm[p] = -hcIm / dn;
        const double hinv = sqrt(hinvRe[p] * hinvRe[p] + hinvIm[p] * hinvIm[p]);
#else
        const double hinv = hp / (hp * hp + beta);                         /* phase = 0, so conj(Hp) = Hp */
#endif
        targetDb[p] = 20.0 * log10(hinv);
    }

    /* band control: hold the value at the edges of the equalized range */
    uint32_t pLo = 0, pHi = P - 1U;
    while (pLo < P - 1U && gridF[pLo] < ROOMEQ_F_LOW_HZ) pLo++;
    while (pHi > 0U && gridF[pHi] > ROOMEQ_F_HIGH_HZ) pHi--;
    for (uint32_t p = 0; p < pLo; p++) targetDb[p] = targetDb[pLo];
    for (uint32_t p = pHi + 1U; p < P; p++) targetDb[p] = targetDb[pHi];

    /* No extra level normalization: the weights are zero-mean over the bands that were heard, bands that were
     * never heard stay at 0 dB, so the equalizer keeps the level and does not move the untouched regions. */
}

/* ====================================================================== */
/*  Setup                                                                 */
/* ====================================================================== */
static void flatCoefs(void)
{
    for (uint32_t i = 0; i < NF; i++) fGainDb[i] = 0.0;
    for (uint32_t b = 0; b < 2U; b++)
    {
        coefBuf[b].n = NF;
        for (uint32_t i = 0; i < NF; i++) pkCoefs(0.0, fCosW[i], fSinW[i], fQ[i], &coefBuf[b].bq[i]);
    }
    coefSeq = coefSeq + 1U;
}

void RoomEq_Init(void)
{
    fftInit();
    memset(ring, 0, sizeof(ring));
    for (uint32_t s = 0; s < STREAMS; s++) total[s] = 0;
    memset(micCorrSet, 0, sizeof(micCorrSet));          /* no mic correction until RoomEq_SetMicCorrection() */
    micCorrOn = 1;

    const double df = ROOMEQ_FS_HZ / (double)K;
    bandEdge[0] = bandCenter[0] * sqrt(bandCenter[0] / bandCenter[1]);
    for (uint32_t m = 1; m < M; m++) bandEdge[m] = sqrt(bandCenter[m - 1U] * bandCenter[m]);
    bandEdge[M] = ROOMEQ_FS_HZ / 2.0 + df;
    smoothEdge[0] = 0.0;
    for (uint32_t m = 1; m < M; m++) smoothEdge[m] = bandEdge[m];
    smoothEdge[M] = 1.0e9;
    for (uint32_t m = 0; m < M; m++)
    {
        uint32_t lo = (uint32_t)ceil(bandEdge[m] / df), hi = (uint32_t)ceil(bandEdge[m + 1U] / df) - 1U;
        if (lo < 1U) lo = 1U;
        if (hi > K2) hi = K2;
        if (hi < lo) hi = lo;
        binLo[m] = lo; binHi[m] = hi;
        bandInRange[m] = (bandCenter[m] >= ROOMEQ_F_LOW_HZ && bandCenter[m] <= ROOMEQ_F_HIGH_HZ) ? 1U : 0U;
        mu[m] = bandStep[m];
    }

    for (uint32_t p = 0; p < P; p++)
    {
        gridF[p] = 20.0 * pow(1000.0, (double)p / (double)(P - 1U));
        const double w = 2.0 * PI_D * gridF[p] / ROOMEQ_FS_HZ;
        gridCos[p] = cos(w);
        gridCos2[p] = cos(2.0 * w);
        gridW[p] = (gridF[p] >= ROOMEQ_F_LOW_HZ && gridF[p] <= ROOMEQ_F_HIGH_HZ) ? 1.0 : 0.25;
    }

    coefCur = 0;
    prepareFilters();
    RoomEq_Reset();
}

void RoomEq_Reset(void)
{
    memset(weightDb, 0, sizeof(weightDb));
    memset(errLin, 0, sizeof(errLin));
    histCount = 0; histPos = 0; framesDone = 0;
    frameStart = ROOMEQ_REF_DELAY;
    for (uint32_t s = 0; s < STREAMS; s++) if (total[s] > frameStart) frameStart = total[s];
    flatCoefs();
}

void RoomEq_SetMicCorrection(uint32_t mic, const double *c, uint32_t n)
{
    if (mic >= R) return;
    for (uint32_t k = 0; k <= K2; k++)
    {
        const double w = 2.0 * PI_D * (double)k / (double)K, c1 = cos(w), c2 = cos(2.0 * w);
        double prod = 1.0;
        for (uint32_t i = 0; i < n; i++)
        {
            const double b0 = c[5U * i], b1 = c[5U * i + 1U], b2 = c[5U * i + 2U], a1 = c[5U * i + 3U], a2 = c[5U * i + 4U];
            const double num = b0 * b0 + b1 * b1 + b2 * b2 + 2.0 * (b0 * b1 + b1 * b2) * c1 + 2.0 * b0 * b2 * c2;
            const double den = 1.0 + a1 * a1 + a2 * a2 + 2.0 * (a1 + a1 * a2) * c1 + 2.0 * a2 * c2;
            prod *= num / den;
        }
        micCorr[mic][k] = (float)sqrt(prod);
    }
    micCorrSet[mic] = 1U;
}

void RoomEq_SetMicCorrectionEnabled(int on) { micCorrOn = on ? 1 : 0; }

void RoomEq_ForgetNoiseFloor(void) { for (uint32_t r = 0; r < R; r++) floorKnown[r] = 0U; }
int  RoomEq_NoiseFloorKnown(void) { return floorKnown[0] && floorKnown[1]; }
int  RoomEq_BandValid(uint32_t mic, uint32_t band) { return bandValid[mic][band]; }
float RoomEq_BandSnrDb(uint32_t mic, uint32_t band) { return bandSnrDb[mic][band]; }
float RoomEq_RefRelDb(uint32_t band) { return refRelDb[band]; }
float RoomEq_ErrRmsDb(uint32_t mic) { return lastErrRms[mic]; }
uint32_t RoomEq_ValidBands(uint32_t mic) { return lastValidCnt[mic]; }

void RoomEq_SetAdapt(int on) { adaptOn = on ? 1 : 0; }
void RoomEq_SetEqActive(int on) { eqActive = on ? 1 : 0; }
int  RoomEq_GetAdapt(void)   { return adaptOn; }

/* ====================================================================== */
/*  Streams                                                               */
/* ====================================================================== */
static void push(uint32_t stream, const int16_t *s, uint32_t n)
{
    uint32_t w = total[stream];
    for (uint32_t i = 0; i < n; i++) ring[stream][(w + i) & (RING - 1U)] = s[i];
    total[stream] = w + n;
}

void RoomEq_PushRef(const int16_t *s, uint32_t n) { push(0U, s, n); }
void RoomEq_PushMic(uint32_t mic, const int16_t *s, uint32_t n) { if (mic < R) push(1U + mic, s, n); }

/* ====================================================================== */
/*  Frame analysis and adaptation                                         */
/* ====================================================================== */
int RoomEq_Poll(void)
{
    uint32_t t[STREAMS];
    for (uint32_t s = 0; s < STREAMS; s++) t[s] = total[s];

    uint32_t tmin = t[0] + ROOMEQ_REF_DELAY, tmax = tmin;
    for (uint32_t s = 1; s < STREAMS; s++)
    {
        if (t[s] < tmin) tmin = t[s];
        if (t[s] > tmax) tmax = t[s];
    }

    if (!adaptOn || !eqActive)                          /* frozen, or the loop is open (equalizer bypassed): resume with fresh data */
    {
        frameStart = tmin;
        return 0;
    }

    /* the main loop was too slow and the rings wrapped: jump forward to data that is still there */
    if (tmax > frameStart + RING - 512U) frameStart = tmax - (RING - 512U);
#if ROOMEQ_REF_DELAY > 0U
    if (frameStart < ROOMEQ_REF_DELAY) frameStart = ROOMEQ_REF_DELAY;
#endif
    if (tmin < frameStart + K) return 0;                /* not enough data yet */

    double bandMean[STREAMS][M];
    double energy[STREAMS];
    for (uint32_t s = 0; s < STREAMS; s++)
    {
        const uint32_t start = (s == 0U) ? (frameStart - ROOMEQ_REF_DELAY) : frameStart;
        energy[s] = loadAndTransform(ring[s], start);
        bandMeans(bandMean[s], (s >= 1U && micCorrOn && micCorrSet[s - 1U]) ? micCorr[s - 1U] : 0);
    }
    frameStart += ROOMEQ_HOP;

    lastRefRms = (float)sqrt(energy[0] / (double)K);
    for (uint32_t r = 0; r < R; r++) lastMicRms[r] = (float)sqrt(energy[1U + r] / (double)K);
    if (lastRefRms < ROOMEQ_REF_RMS_MIN)                /* nothing is played: the mics hear their own noise */
    {
        for (uint32_t r = 0; r < R; r++)
        {
            for (uint32_t m = 0; m < M; m++)
                micFloor[r][m] = floorKnown[r] ? 0.7 * micFloor[r][m] + 0.3 * bandMean[1U + r][m] : bandMean[1U + r][m];
            floorKnown[r] = 1U;
        }
        return 1;
    }
    for (uint32_t r = 0; r < R; r++) if (lastMicRms[r] < ROOMEQ_MIC_RMS_MIN) return 1;   /* mic dead / not connected */
    if (!floorKnown[0] || !floorKnown[1]) return 1;                                      /* wait for the noise floor */

    /* which bands carry a usable reference and are heard above the mic noise? */
    double xmax = 0.0;
    for (uint32_t m = 0; m < M; m++) if (bandMean[0][m] > xmax) xmax = bandMean[0][m];
    const double relMin = pow(10.0, -ROOMEQ_REF_BAND_REL_DB / 20.0);
    for (uint32_t m = 0; m < M; m++) refRelDb[m] = (float)(20.0 * log10((bandMean[0][m] + 1e-9) / (xmax + 1e-9)));
    for (uint32_t r = 0; r < R; r++)
        for (uint32_t m = 0; m < M; m++)
        {
            const double snr = bandMean[1U + r][m] / (micFloor[r][m] + 1e-9);
            bandSnrDb[r][m] = (float)(20.0 * log10(snr + 1e-9));
            bandValid[r][m] = (bandMean[0][m] >= relMin * xmax && snr >= ROOMEQ_SNR_MIN) ? 1U : 0U;
        }

    /* e = d - x with the common level removed. x is taken before the equalizer, so the equalizer in the signal path
     * (about -H) closes the loop and H settles at the room response. */
    const uint32_t slot = histPos;                       /* history slots 0..histCount-1 are valid */
    histPos = (histPos + 1U) % L;
    for (uint32_t r = 0; r < R; r++)
    {
        double e[M], mean = 0.0;
        uint32_t cnt = 0;
        for (uint32_t m = 0; m < M; m++)
        {
            e[m] = 20.0 * log10((bandMean[1U + r][m] + 1e-9) / (bandMean[0][m] + 1e-9));
            if (bandValid[r][m] && bandInRange[m]) { mean += e[m]; cnt++; }
        }
        {
            double ss = 0.0;
            for (uint32_t m = 0; m < M; m++)
                if (bandValid[r][m] && bandInRange[m]) ss += (e[m] - (cnt ? mean / (double)cnt : 0.0)) * (e[m] - (cnt ? mean / (double)cnt : 0.0));
            lastErrRms[r]   = cnt ? (float)sqrt(ss / (double)cnt) : 0.0f;
            lastValidCnt[r] = cnt;
        }
        const int enough = (cnt >= ROOMEQ_MIN_VALID_BANDS);
        if (cnt) mean /= (double)cnt;
        for (uint32_t m = 0; m < M; m++)                 /* invalid bands: error 0 dB = the weight is not changed */
            errLin[r][m][slot] = (enough && bandValid[r][m]) ? (float)pow(10.0, (e[m] - mean) / 20.0) : 1.0f;
    }
    if (histCount < L) histCount++;

    /* average of the last L errors, weight update, weight limit */
    for (uint32_t r = 0; r < R; r++)
        for (uint32_t m = 0; m < M; m++)
        {
            double acc = 0.0;
            for (uint32_t i = 0; i < histCount; i++) acc += (double)errLin[r][m][i];
            const double E = 20.0 * log10(acc / (double)histCount);
            double h = (double)weightDb[r][m] + (double)mu[m] * E;
            if (h >  ROOMEQ_WEIGHT_MAX_DB) h =  ROOMEQ_WEIGHT_MAX_DB;
            if (h < -ROOMEQ_WEIGHT_MAX_DB) h = -ROOMEQ_WEIGHT_MAX_DB;
            weightDb[r][m] = (float)h;
        }
    framesDone++;

    double target[P];
    buildTarget(target);
    fitAndPublish(target, 0);
    return 2;
}

/* ====================================================================== */
/*  Accessors                                                             */
/* ====================================================================== */
void RoomEq_PeakCoefs(double fc_hz, double gain_db, double q, RoomEqBiquad *out)
{
    const double w0 = 2.0 * PI_D * fc_hz / ROOMEQ_FS_HZ;
    pkCoefs(gain_db, cos(w0), sin(w0), q, out);
}

#if ROOMEQ_COMPLEX_INVERSE
uint32_t RoomEq_GetComplexInverse(double *freq, double *re, double *im)
{
    for (uint32_t p = 0; p < P; p++) { freq[p] = gridF[p]; re[p] = hinvRe[p]; im[p] = hinvIm[p]; }
    return P;
}
#endif

uint32_t RoomEq_CoefSeq(void) { return coefSeq; }
const RoomEqCoefs *RoomEq_Coefs(void) { return &coefBuf[coefCur]; }
double RoomEq_BandCenter(uint32_t m) { return bandCenter[m]; }
double RoomEq_BandStep(uint32_t m) { return mu[m]; }
float  RoomEq_Weight(uint32_t mic, uint32_t band) { return weightDb[mic][band]; }
double RoomEq_FilterFreq(uint32_t i) { return fFreq[i]; }
double RoomEq_FilterGainDb(uint32_t i) { return fGainDb[i]; }
uint32_t RoomEq_Frames(void) { return framesDone; }
float  RoomEq_RefRms(void) { return lastRefRms; }
float  RoomEq_MicRms(uint32_t mic) { return lastMicRms[mic]; }
