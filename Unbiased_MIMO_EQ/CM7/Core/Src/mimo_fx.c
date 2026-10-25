/**
 * @file    mimo_fx.c
 * @brief   Multichannel room equalizer based on IPNLMS-FX, see mimo_fx.h.
 *
 * Where the algorithm is in the code:
 *   y_j(n) = w_j^T x(n)                                      MimoFx_ProcessBlock()
 *   z_m(n) = sum_j h_mj * y_j(n)                             the real room (the microphone signal zw[m])
 *   e_m(n) = d_m(n) - z_m(n)                                 update(), "e[m]"  (d_m(n) = x(n - tau_m))
 *   x_f(m,j) = h^_mj * x                                     runBlock(), "xfw"  (h^ placed at its true lag S_m)
 *   g_j (IPNLMS gains)                                       update(), "A + Bc |w|"
 *   mu_j(n) (step size per tap)                              update(), "c"
 *   w_j(n) = w_j(n-1) + mu_j(n) sum_m e_m(n) x_f(m,j)(n)     update(), last loop
 *
 * Filters are stored reversed in memory (oldest tap first): a dot product over a window of the
 * history then runs over two ascending arrays.
 */
/* Always optimized, independent of the IDE's Debug / Release setting. */
#pragma GCC optimize ("O2")
#include "mimo_fx.h"
#include "farina.h"
#include <math.h>
#include <string.h>

#ifdef MIMO_HOST
#define DTCM
#define BARRIER() __sync_synchronize()
#else
#define DTCM __attribute__((section(".dtcm"), aligned(8)))
#define BARRIER() __asm volatile ("dmb" ::: "memory")
#endif
#define HOT  __attribute__((optimize("O2")))

#define J   MIMO_J
#define M   MIMO_M
#define LW  MIMO_LW
#define LH  MIMO_LH
#define B   MIMO_BLOCK
#define D   MIMO_ADAPT_DECIM
#define HIST MIMO_HIST
#define RMASK (MIMO_RING - 1U)

/* ---------------------------------------------------------------------------------------------- */
/* Ring buffers between the interrupts and the main loop (AXI SRAM, only accessed in blocks)       */
/* ---------------------------------------------------------------------------------------------- */
static float xRing[MIMO_RING];
static float zRing[M][MIMO_RING];
static volatile uint32_t refBlocks, micBlocks[M];

/* ---------------------------------------------------------------------------------------------- */
/* Interrupt part                                                                                  */
/* ---------------------------------------------------------------------------------------------- */
DTCM static float xd[LW - 1U + B];            /* linear history window of the input of filter 1 (the mono reference, or the left channel) */
DTCM static float wPub[2][J][LW];             /* published weights (reversed order), double buffer */
static volatile uint32_t wCur;

/* ---------------------------------------------------------------------------------------------- */
/* Main loop part                                                                                  */
/* ---------------------------------------------------------------------------------------------- */
DTCM static float xw[HIST + B];               /* x(n0 - HIST + i) */
DTCM static float zw[M][B];                   /* z_m / G_m of the block */
DTCM static float xfw[M][J][LW - 1U + B];     /* x_f(m,j): the last LW-1 samples of the previous block + this block */
DTCM static float xsq[J][LW - 1U + B];        /* sum over m of x_f(m,j)^2 */
DTCM static float hrev[M][J][LH];             /* h^_mj, newest tap last */
DTCM static float w[J][LW];                   /* working weights, reversed: w[j][LW-1-l] = w_j(l) */

static uint16_t shiftM[M], tauM[M];
static float    invG[M];
static MimoIr   irStore;
static int      haveIr;
static volatile int adaptOn = 1;
static int      needInit = 1;
static uint32_t pollBlock, blocksDone, skipped, backlog;
static float    muScale = 1.0f;               /* step size factor: halved every time a filter diverges */
static uint32_t resets;                       /* number of divergence resets */
static uint32_t checkCnt;
static int      irNew;                        /* 1 after a new h^ was installed, read once by MimoFx_IrInstalled */
static uint32_t rtBlocks;                      /* blocks processed since the adaptation started */
static float    gMu = MIMO_RT_MU, gDelta = MIMO_RT_DELTA, gEps = MIMO_RT_EPSILON;      /* the gain parameters in use (runtime adjustable) */
DTCM static float wSnap[J][LW];                   /* a recent state of the weights: where a divergence returns to */
static float    errDb[M];                     /* real microphone error e_m relative to d_m, dB, smoothed */
HOT static inline __attribute__((always_inline)) float dot(const float *a, const float *b, uint32_t n)
{
    float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f, s3 = 0.0f;
    for (uint32_t i = 0U; i < n; i += 4U)
    {
        s0 += a[i]      * b[i];
        s1 += a[i + 1U] * b[i + 1U];
        s2 += a[i + 2U] * b[i + 2U];
        s3 += a[i + 3U] * b[i + 3U];
    }
    return (s0 + s1) + (s2 + s3);
}

/* 4 outputs of two filters that share the same input window:  o0[s] = sum_k a0[k] x[k+s],  o1[s] = sum_k a1[k] x[k+s],  s = 0 .. 3.
 * Coefficient and input loads are reused (6 loads for 8 multiply-adds): the Cortex-M7 loads once per cycle. Valid because the weights / h^ are constant over the 4 outputs. */
HOT static inline __attribute__((always_inline)) void dot4x2(const float *a0, const float *a1, const float *x, uint32_t n, float *o0, float *o1)
{
    float p0 = 0.0f, p1 = 0.0f, p2 = 0.0f, p3 = 0.0f, q0 = 0.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;
    for (uint32_t k = 0U; k < n; k++)
    {
        const float x0 = x[k], x1 = x[k + 1U], x2 = x[k + 2U], x3 = x[k + 3U];
        const float u = a0[k], v = a1[k];
        p0 += u * x0;  p1 += u * x1;  p2 += u * x2;  p3 += u * x3;
        q0 += v * x0;  q1 += v * x1;  q2 += v * x2;  q3 += v * x3;
    }
    o0[0] = p0; o0[1] = p1; o0[2] = p2; o0[3] = p3;
    o1[0] = q0; o1[1] = q1; o1[2] = q2; o1[3] = q3;
}

/* 4 outputs of one filter: o[s] = sum_k a[k] x[k+s] */
HOT static inline __attribute__((always_inline)) void dot4(const float *a, const float *x, uint32_t n, float *o)
{
    float p0 = 0.0f, p1 = 0.0f, p2 = 0.0f, p3 = 0.0f;
    for (uint32_t k = 0U; k < n; k++)
    {
        const float u = a[k];
        p0 += u * x[k];  p1 += u * x[k + 1U];  p2 += u * x[k + 2U];  p3 += u * x[k + 3U];
    }
    o[0] = p0; o[1] = p1; o[2] = p2; o[3] = p3;
}

/* ---------------------------------------------------------------------------------------------- */
/* Bass of the filters stays neutral (MIMO_LF_NEUTRAL_HZ): the deviation of w from the neutral filter is kept orthogonal to cos / sin vectors of the frequencies
 * 0 .. MIMO_LF_NEUTRAL_HZ (every 20 Hz), so its frequency response is zero there (projection after every block). The vectors are orthonormalized once
 * (Gram-Schmidt, nearly dependent ones are dropped). */
#if MIMO_LF_NEUTRAL_HZ > 0
#define LF_MAXQ 8U      /* maximum number of basis vectors */
static float    lfQ[LF_MAXQ][LW];                 /* orthonormal basis, in the storage order of w (oldest tap first) */
static uint32_t lfN;
static void lfInit(void)
{
    lfN = 0U;
    for (uint32_t f = 0U; f <= (uint32_t)MIMO_LF_NEUTRAL_HZ; f += 20U)
        for (uint32_t kind = 0U; kind < 2U; kind++)
        {
            if (f == 0U && kind == 1U) continue;                          /* sin(0) = 0 */
            float v[LW];
            double n0 = 0.0;
            for (uint32_t l = 0U; l < LW; l++)
            {
                const double ph = 2.0 * 3.14159265358979 * (double)f * (double)(LW - 1U - l) / (double)MIMO_FS_HZ;     /* lag of tap l is LW-1-l */
                v[l] = (float)(kind ? sin(ph) : cos(ph));
                n0 += (double)v[l] * (double)v[l];
            }
            for (uint32_t pass = 0U; pass < 2U; pass++)                    /* twice: modified Gram-Schmidt with re-orthogonalization */
                for (uint32_t q = 0U; q < lfN; q++)
                {
                    float c = 0.0f;
                    for (uint32_t l = 0U; l < LW; l++) c += lfQ[q][l] * v[l];
                    for (uint32_t l = 0U; l < LW; l++) v[l] -= c * lfQ[q][l];
                }
            double n1 = 0.0;
            for (uint32_t l = 0U; l < LW; l++) n1 += (double)v[l] * (double)v[l];
            if (n1 < 0.03 * n0 || lfN >= LF_MAXQ) continue;               /* (nearly) dependent on the ones before: dropped */
            const float inv = 1.0f / sqrtf((float)n1);
            for (uint32_t l = 0U; l < LW; l++) lfQ[lfN][l] = v[l] * inv;
            lfN++;
        }
}
static void lfConstrain(void)
{
    uint32_t dly = (uint32_t)(((uint32_t)(tauM[0] - shiftM[0]) + (uint32_t)(tauM[1] - shiftM[1])) / 2U);
    if (dly > LW - 1U) dly = LW - 1U;
    const uint32_t ip = LW - 1U - dly;                                      /* index of the neutral pulse */
    for (uint32_t j = 0U; j < J; j++)
    {
        float dev[LW];
        for (uint32_t l = 0U; l < LW; l++) dev[l] = w[j][l] - ((l == ip) ? 1.0f / (float)J : 0.0f);
        for (uint32_t q = 0U; q < lfN; q++)
        {
            float c = 0.0f;
            for (uint32_t l = 0U; l < LW; l++) c += lfQ[q][l] * dev[l];
            for (uint32_t l = 0U; l < LW; l++) dev[l] -= c * lfQ[q][l];
        }
        for (uint32_t l = 0U; l < LW; l++) w[j][l] = dev[l] + ((l == ip) ? 1.0f / (float)J : 0.0f);
    }
}
#endif

/* ---------------------------------------------------------------------------------------------- */
#if MIMO_SMOOTH_OCT_DIV > 0 || MIMO_IR_LF_IGNORE_HZ > 0
/* Fractional octave smoothing of a measured response h (LH taps), magnitude and phase (complex smoothing): every spectrum value is replaced by the mean of the complex values within
 * +- 1/(2 MIMO_SMOOTH_OCT_DIV) octave. The direct peak is rotated to index 0 first so that its delay does not rotate the phase inside the averaging band. */
#define SM_N 512U
static float smCos[SM_N], smSin[SM_N];
static float smPreR[SM_N / 2U + 2U], smPreI[SM_N / 2U + 2U];
static float smG[SM_N], smSr[SM_N / 2U + 1U], smSi[SM_N / 2U + 1U];

static void smoothResponse(float *h)
{
    const uint32_t N = SM_N, K = SM_N / 2U;
    static uint8_t tab = 0;
    if (!tab)
    {
        for (uint32_t i = 0U; i < N; i++) { smCos[i] = cosf(6.28318531f * (float)i / (float)N); smSin[i] = sinf(6.28318531f * (float)i / (float)N); }
        tab = 1;
    }
    uint32_t pk = 0U; float pv = 0.0f;
    for (uint32_t i = 0U; i < LH; i++) if (fabsf(h[i]) > pv) { pv = fabsf(h[i]); pk = i; }
    for (uint32_t i = 0U; i < N; i++) smG[i] = 0.0f;
    for (uint32_t i = 0U; i < LH; i++) smG[(i + N - pk) % N] = h[i];          /* peak -> index 0 */

    smPreR[0] = 0.0f; smPreI[0] = 0.0f;                                       /* spectrum, running sums for the band means */
    for (uint32_t k = 0U; k <= K; k++)
    {
        float re = 0.0f, im = 0.0f;
        for (uint32_t i = 0U; i < N; i++)
        {
            const uint32_t t = (uint32_t)(((uint64_t)i * k) % N);
            re += smG[i] * smCos[t];
            im -= smG[i] * smSin[t];
        }
        smPreR[k + 1U] = smPreR[k] + re;
        smPreI[k + 1U] = smPreI[k] + im;
    }
#if MIMO_SMOOTH_OCT_DIV > 0
    const float r = powf(2.0f, 0.5f / (float)MIMO_SMOOTH_OCT_DIV);
#else
    const float r = 1.0f;
#endif
    for (uint32_t k = 0U; k <= K; k++)
    {
        uint32_t lo = (uint32_t)ceilf((float)k / r), hi = (uint32_t)floorf((float)k * r);
        if (k == 0U) { lo = 0U; hi = 0U; }
        if (hi > K) hi = K;
        if (lo > k) lo = k;
        if (hi < k) hi = k;
        const float inv = 1.0f / (float)(hi - lo + 1U);
        smSr[k] = (smPreR[hi + 1U] - smPreR[lo]) * inv;
        smSi[k] = (smPreI[hi + 1U] - smPreI[lo]) * inv;
    }
#if MIMO_IR_LF_IGNORE_HZ > 0
    {   /* below MIMO_IR_LF_IGNORE_HZ the spectrum is replaced by the mean level of 500 .. 2000 Hz with the phase of the direct sound (the peak sits at index 0),
         * faded back to the measurement over one octave above it */
        const float binHz = (float)MIMO_FS_HZ / (float)N;
        float acc = 0.0f; uint32_t cnt = 0U;
        for (uint32_t k = 1U; k <= K; k++)
        {
            const float f = (float)k * binHz;
            if (f >= 500.0f && f <= 2000.0f) { acc += sqrtf(smSr[k] * smSr[k] + smSi[k] * smSi[k]); cnt++; }
        }
        const float amid = (cnt > 0U) ? acc / (float)cnt : 0.0f;
        for (uint32_t k = 0U; k <= K; k++)
        {
            const float f = (float)k * binHz;
            float c = 0.0f;
            if (f <= (float)MIMO_IR_LF_IGNORE_HZ) c = 1.0f;
            else if (f < 2.0f * (float)MIMO_IR_LF_IGNORE_HZ) c = 0.5f * (1.0f + cosf(3.14159265f * (f / (float)MIMO_IR_LF_IGNORE_HZ - 1.0f)));
            if (c > 0.0f) { smSr[k] = c * amid + (1.0f - c) * smSr[k]; smSi[k] = (1.0f - c) * smSi[k]; }
        }
    }
#endif
    for (uint32_t i = 0U; i < LH; i++)                                         /* inverse DFT, only the LH wanted samples, rotated back */
    {
        const uint32_t n = (i + N - pk) % N;
        float acc = smSr[0] + smSr[K] * ((n & 1U) ? -1.0f : 1.0f);
        for (uint32_t k = 1U; k < K; k++)
        {
            const uint32_t t = (uint32_t)(((uint64_t)n * k) % N);
            acc += 2.0f * (smSr[k] * smCos[t] - smSi[k] * smSin[t]);
        }
        h[i] = acc / (float)N;
    }
}
#endif

int MimoIr_Build(MimoIr *ir, const float *const win[M][J])
{
    memset(ir, 0, sizeof(*ir));
    for (uint32_t m = 0U; m < M; m++)
    {
        /* the window of every path starts at the same lag (MIMO_IR_START) so the delay between the loudspeakers stays in h^; the desired signal is delayed by MIMO_MODEL_DELAY after it */
        ir->shift[m] = (uint16_t)MIMO_IR_START;
        ir->tau[m]   = (uint16_t)(MIMO_IR_START + MIMO_MODEL_DELAY);

        double energy = 0.0;
        for (uint32_t j = 0U; j < J; j++)
        {
            float pk = 0.0f; int8_t sg = 1;                       /* polarity of the strongest sample of the window */
            for (uint32_t k = 0U; k < LH; k++)
            {
                float v = win[m][j][k];
                if (fabsf(v) > pk) { pk = fabsf(v); sg = (v < 0.0f) ? -1 : 1; }
                if (k + MIMO_IR_FADE >= LH)                       /* fade out the truncated end */
                {
                    const float u = (float)(LH - 1U - k) / (float)MIMO_IR_FADE;
                    v *= 0.5f * (1.0f - cosf(3.14159265f * u));
                }
                ir->h[m][j][k] = v;
                energy += (double)v * (double)v;
            }
            ir->polarity[m][j] = sg;
        }
#if MIMO_SMOOTH_OCT_DIV > 0 || MIMO_IR_LF_IGNORE_HZ > 0
        energy = 0.0;
        for (uint32_t j = 0U; j < J; j++)
        {
            smoothResponse(ir->h[m][j]);
            for (uint32_t k = 0U; k < LH; k++) energy += (double)ir->h[m][j][k] * (double)ir->h[m][j][k];
        }
#endif
        if (energy < 1.0e-12) return -1;
        ir->gain[m] = MIMO_IR_NORMALIZE ? (float)sqrt(energy / (double)J) : 1.0f;
        for (uint32_t j = 0U; j < J; j++)
            for (uint32_t k = 0U; k < LH; k++) ir->h[m][j][k] /= ir->gain[m];
    }
    return 0;
}

/* ---------------------------------------------------------------------------------------------- */
static float jumpPeak;                      /* largest relative change of the played weights between two publishes since it was read */
static void publish(void)
{
    const uint32_t idle = 1U - wCur;
    memcpy(wPub[idle], w, sizeof(w));
    {
        float dn = 0.0f, wn = 1.0e-12f;
        for (uint32_t j = 0U; j < J; j++)
            for (uint32_t l = 0U; l < LW; l++) { const float d = wPub[idle][j][l] - wPub[wCur][j][l]; dn += d * d; wn += wPub[wCur][j][l] * wPub[wCur][j][l]; }
        const float r = sqrtf(dn / wn);
        if (r > jumpPeak) jumpPeak = r;
    }
    BARRIER();
    wCur = idle;
}

static void resetWeights(void)
{
    memset(w, 0, sizeof(w));
    for (uint32_t j = 0U; j < J; j++) w[j][LW - 1U] = 1.0f / (float)J;     /* unit pulse, 1/J for each loudspeaker */
    needInit = 1;
}

void MimoFx_Init(void)
{
    memset(xRing, 0, sizeof(xRing));
    memset(zRing, 0, sizeof(zRing));
    memset(xd, 0, sizeof(xd));
    memset(wPub, 0, sizeof(wPub));
    memset(xw, 0, sizeof(xw));
    memset(zw, 0, sizeof(zw));
    memset(xfw, 0, sizeof(xfw));
    memset(xsq, 0, sizeof(xsq));
    memset(hrev, 0, sizeof(hrev));
    refBlocks = 0U;
    for (uint32_t m = 0U; m < M; m++) { micBlocks[m] = 0U; invG[m] = 1.0f; errDb[m] = 0.0f; }
#if MIMO_LF_NEUTRAL_HZ > 0
    lfInit();
#endif
    haveIr = 0;
    pollBlock = blocksDone = skipped = backlog = 0U;
    wCur = 0U;
    resetWeights();
    publish();
    publish();                                   /* both buffers hold the start weights */
}

void MimoFx_Reset(void)
{
    resetWeights();
    publish();
    publish();
}

int MimoFx_SetIr(const MimoIr *ir)
{
    muScale = 1.0f;
    irStore = *ir;
    for (uint32_t m = 0U; m < M; m++)
    {
        shiftM[m] = ir->shift[m];
        tauM[m]   = ir->tau[m];
        invG[m]   = 1.0f / ir->gain[m];
        for (uint32_t j = 0U; j < J; j++)
            for (uint32_t k = 0U; k < LH; k++) hrev[m][j][LH - 1U - k] = ir->h[m][j][k];
    }
    haveIr = 1;
    resetWeights();
    {   /* start from the neutral filter: unit pulse with the delay of the desired signal */
        uint32_t dly = (uint32_t)(((uint32_t)(tauM[0] - shiftM[0]) + (uint32_t)(tauM[1] - shiftM[1])) / 2U);
        if (dly > LW - 1U) dly = LW - 1U;
        memset(w, 0, sizeof(w));
        for (uint32_t j = 0U; j < J; j++) w[j][LW - 1U - dly] = 1.0f / (float)J;
        memcpy(wSnap, w, sizeof(w));
        publish();
        publish();
        gMu = MIMO_RT_MU; gDelta = MIMO_RT_DELTA; gEps = MIMO_RT_EPSILON;
        rtBlocks = 0U;
        irNew = 1;
        needInit = 1;
        pollBlock = refBlocks;
    }
    return 0;
}

int  MimoFx_HaveIr(void)                { return haveIr; }
const MimoIr *MimoFx_GetIr(void)        { return &irStore; }
void MimoFx_SetAdapt(int on)            { adaptOn = on; needInit = 1; }
int  MimoFx_GetAdapt(void)              { return adaptOn; }
float MimoFx_Weight(uint32_t j, uint32_t l) { return w[j][LW - 1U - l]; }
float MimoFx_ErrDb(uint32_t m)          { return errDb[m]; }
uint32_t MimoFx_BlocksDone(void)        { return blocksDone; }
uint32_t MimoFx_Skipped(void)           { return skipped; }
uint32_t MimoFx_Backlog(void)           { return backlog; }
uint32_t MimoFx_RefBlocks(void)         { return refBlocks; }
uint32_t MimoFx_MicBlocks(uint32_t m)   { return micBlocks[m]; }

/* largest and smallest gain (dB) of the working filter between 20 Hz and 20 kHz */
void MimoFx_GainRangeDb(uint32_t j, float *maxDb, float *minDb)
{
    float hi = 0.0f, lo = 1.0e30f;
    for (uint32_t p = 0U; p < 64U; p++)
    {
        const float f  = 20.0f * powf(900.0f, (float)p / 63.0f);               /* 20 Hz .. 18 kHz */
        const float wr = 2.0f * 3.14159265f * f / (float)MIMO_FS_HZ;
        const float cr = cosf(wr), ci = -sinf(wr);
        float zr = 1.0f, zi = 0.0f, sr = 0.0f, si = 0.0f;
        for (uint32_t l = 0U; l < LW; l++)
        {
            const float v = w[j][LW - 1U - l];
            sr += v * zr;  si += v * zi;
            const float t = zr * cr - zi * ci;  zi = zr * ci + zi * cr;  zr = t;
        }
        const float m = sr * sr + si * si;
        if (m > hi) hi = m;
        if (m < lo) lo = m;
    }
    if (maxDb) *maxDb = 10.0f * log10f(hi + 1.0e-20f);
    if (minDb) *minDb = 10.0f * log10f(lo + 1.0e-20f);
}

float MimoFx_PeakGainDb(uint32_t j, float *freqHz)
{
    float best = 0.0f, bf = 0.0f;
    for (uint32_t p = 0U; p < 96U; p++)                       /* 96 log spaced points, 20 Hz .. 20 kHz */
    {
        const float f  = 20.0f * powf(1000.0f, (float)p / 95.0f);
        const float wr = 2.0f * 3.14159265f * f / (float)MIMO_FS_HZ;
        const float cr = cosf(wr), ci = -sinf(wr);              /* e^{-j w} */
        float zr = 1.0f, zi = 0.0f, sr = 0.0f, si = 0.0f;
        for (uint32_t l = 0U; l < LW; l++)                      /* sum_l w_j(l) e^{-j w l}; w[] is stored reversed */
        {
            const float v = w[j][LW - 1U - l];
            sr += v * zr;  si += v * zi;
            const float t = zr * cr - zi * ci;  zi = zr * ci + zi * cr;  zr = t;
        }
        const float m = sr * sr + si * si;
        if (m > best) { best = m; bf = f; }
    }
    if (freqHz) *freqHz = bf;
    return 10.0f * log10f(best + 1.0e-20f);
}

void MimoFx_FilterResponse(uint32_t j, float freqHz, float *magDb, float *phaseDeg)
{
    const float wr = 2.0f * 3.14159265f * freqHz / (float)MIMO_FS_HZ;
    float sr = 0.0f, si = 0.0f;
    for (uint32_t l = 0U; l < LW; l++)                       /* sum_l w_j(l) e^{-j w l}; w[] is stored reversed */
    {
        const float v = wPub[wCur][j][LW - 1U - l];          /* the published filter (what the loudspeaker plays) */
        sr += v * cosf(wr * (float)l);
        si -= v * sinf(wr * (float)l);
    }
    if (magDb)    *magDb    = 10.0f * log10f(sr * sr + si * si + 1.0e-20f);
    if (phaseDeg) *phaseDeg = atan2f(si, sr) * (180.0f / 3.14159265f);
}

float MimoFx_WeightRms(uint32_t j)
{
    double s = 0.0;
    for (uint32_t l = 0U; l < LW; l++) s += (double)w[j][l] * (double)w[j][l];
    return (float)sqrt(s);
}

/* ---------------------------------------------------------------------------------------------- */
/* Audio interrupt                                                                                 */
/* ---------------------------------------------------------------------------------------------- */
HOT void MimoFx_ProcessBlock(const float *x0, float out[J][B])
{
    const uint32_t cur = wCur;

    memmove(xd, xd + B, (LW - 1U) * sizeof(float));
    memcpy(xd + (LW - 1U), x0, B * sizeof(float));

    const uint32_t base = (refBlocks * B) & RMASK;            /* the ring for the adaptation gets the reference */
    memcpy(&xRing[base], x0, B * sizeof(float));

    for (uint32_t i = 0U; i < B; i += 4U)                      /* one input, both filters, weights are constant over the block */
    {
        float o0[4], o1[4];
        dot4x2(wPub[cur][0], wPub[cur][1], &xd[i], LW, o0, o1);
        for (uint32_t t = 0U; t < 4U; t++) { out[0][i + t] = o0[t]; out[1][i + t] = o1[t]; }
    }
    BARRIER();
    refBlocks = refBlocks + 1U;
}

void MimoFx_PushMic(uint32_t m, const float *z)
{
    static uint8_t micSynced[M];
    if (!micSynced[m]) { micBlocks[m] = refBlocks; micSynced[m] = 1U; }    /* the first microphone block is paired with the reference block that is running now */
    memcpy(&zRing[m][(micBlocks[m] * B) & RMASK], z, B * sizeof(float));
    BARRIER();
    micBlocks[m] = micBlocks[m] + 1U;
}

/* ---------------------------------------------------------------------------------------------- */
/* Main loop                                                                                       */
/* ---------------------------------------------------------------------------------------------- */
/* One update step of the adaptation at sample i of the block.
 *
 * e_m(n) = d_m(n) - z_m(n) is the error of the real microphone: d_m(n) = x(n - tau_m), z_m(n) the microphone signal divided by G_m. The room, the audio chain and the microphone
 * filters delay the loudspeaker signal by the bulk delay S_m (MIMO_IR_START), so the filtered reference x_f(m,j)(n) = sum_k h^_mj[k] x(n - S_m - k) is computed with h^ at its true
 * lag (see runBlock): e_m(n) and x_f(n - l) belong together in time. This delay is inside the adaptation loop (S_m / D updates) and limits the step size mu. */
HOT static void update(uint32_t i, float *sumE2, float *sumD2)
{
    float e[M];
    for (uint32_t m = 0U; m < M; m++)
    {
        const float d = xw[HIST + i - ((tauM[m] < HIST) ? tauM[m] : HIST)];     /* d_m(n) = x(n - tau_m) */
        e[m] = d - zw[m][i];
        *sumE2 += e[m] * e[m];
        *sumD2 += d * d;
    }

    for (uint32_t j = 0U; j < J; j++)
    {
        float *wj = w[j];
        const float *xf0 = &xfw[0][j][i];                    /* x_f(m,j)(n - l), oldest first */
        const float *xf1 = &xfw[1][j][i];
        const float *xs  = &xsq[j][i];

        float sumAbs = 0.0f, s2 = 0.0f, wxs = 0.0f;
        for (uint32_t l = 0U; l < LW; l++)
        {
            const float aw = fabsf(wj[l]);
            sumAbs += aw;
            s2     += xs[l];
            wxs    += aw * xs[l];
        }
        const float A  = (1.0f - MIMO_KAPPA) / (2.0f * (float)LW);                    /* IPNLMS gain: constant part */
        const float Bc = (1.0f + MIMO_KAPPA) / (gEps + 2.0f * sumAbs);
        float c  = muScale * gMu / (gDelta + A * s2 + Bc * wxs);   /* step size normalization */
        const float cA = c * A, cB = c * Bc;
        const float e0 = e[0], e1 = e[1];

        for (uint32_t l = 0U; l < LW; l++)
        {
            const float wl = wj[l];
            const float t  = e0 * xf0[l] + e1 * xf1[l];                                 /* sum_m e_m x_f(m,j)(n-l) */
            wj[l] = wl + (cA + cB * fabsf(wl)) * t;
        }
    }
}

/* x_f(m,j) of the LW-1 samples in front of the block (after a reset, a new h^ or a gap) */
HOT static void initTails(void)
{
    for (uint32_t m = 0U; m < M; m++)
    {
        for (uint32_t j = 0U; j < J; j++)
        {
            for (uint32_t t = 0U; t < LW - 1U; t++)
            {
                /* sample index relative to the block start: i = t - (LW-1) */
                const int32_t i = (int32_t)t - (int32_t)(LW - 1U);
                xfw[m][j][t] = dot(hrev[m][j], &xw[(int32_t)HIST + i - (int32_t)shiftM[m] - (int32_t)LH + 1], LH);
            }
        }
    }
    for (uint32_t j = 0U; j < J; j++)
        for (uint32_t t = 0U; t < LW - 1U; t++)
        {
            float s = 0.0f;
            for (uint32_t m = 0U; m < M; m++) s += xfw[m][j][t] * xfw[m][j][t];
            xsq[j][t] = s;
        }
}

/* A filter diverged: back to a recent state of the filters, with a half as large step size. */
static void diverged(void)
{
    resets++;
    if (muScale > 1.0f / 64.0f) muScale *= 0.5f;
    memcpy(w, wSnap, sizeof(w));
    needInit = 1;
    publish();
    publish();
}

/* One block of the adaptation. xw[] (history + this block) and zw[] must be filled. */
HOT static void runBlock(void)
{
    if (needInit) { initTails(); needInit = 0; }

    float e2 = 0.0f, d2 = 0.0f;

    for (uint32_t g0 = 0U; g0 < B; g0 += D)
    {
        /* the weights change only in update(), once per group of D samples; x_f of the group does not depend on them */
        for (uint32_t sdx = 0U; sdx < D; sdx += 4U)
        {
            const uint32_t i = g0 + sdx;
            float o0[4], o1[4];

            /* x_f(m,j)(n) = sum_k h^_mj[k] x(n - S_m - k)   (h^ at its true lag) */
            for (uint32_t m = 0U; m < M; m++)
            {
                dot4x2(hrev[m][0], hrev[m][1], &xw[HIST + i - shiftM[m] - LH + 1U], LH, o0, o1);
                for (uint32_t t = 0U; t < 4U; t++) { xfw[m][0][LW - 1U + i + t] = o0[t]; xfw[m][1][LW - 1U + i + t] = o1[t]; }
            }
            for (uint32_t t = 0U; t < 4U; t++)
                for (uint32_t j = 0U; j < J; j++)
                {
                    float sq = 0.0f;
                    for (uint32_t m = 0U; m < M; m++) sq += xfw[m][j][LW - 1U + i + t] * xfw[m][j][LW - 1U + i + t];
                    xsq[j][LW - 1U + i + t] = sq;
                }
        }

        update(g0 + D - 1U, &e2, &d2);
    }

    /* error figures (both mics together, smoothed in dB) */
    {
        const float r = 10.0f * log10f((e2 + 1.0e-20f) / (d2 + 1.0e-20f));
        for (uint32_t m = 0U; m < M; m++) errDb[m] = 0.9f * errDb[m] + 0.1f * r;
    }

    for (uint32_t m = 0U; m < M; m++)
        for (uint32_t j = 0U; j < J; j++)
            memmove(xfw[m][j], &xfw[m][j][B], (LW - 1U) * sizeof(float));
    for (uint32_t j = 0U; j < J; j++) memmove(xsq[j], &xsq[j][B], (LW - 1U) * sizeof(float));

#if MIMO_LF_NEUTRAL_HZ > 0
    lfConstrain();                                           /* the bass of the filters goes back to the neutral filter */
#endif

    /* NaN / runaway protection */
    for (uint32_t j = 0U; j < J; j++)
        for (uint32_t l = 0U; l < LW; l++)
            if (!(w[j][l] > -100.0f && w[j][l] < 100.0f)) { diverged(); return; }

    /* divergence: the largest gain of a filter between 20 Hz and 20 kHz (checked every block) */
    if (MIMO_GAIN_LIMIT_DB > 0.0f && ((++checkCnt) & 3U) == 0U)
        if (MimoFx_PeakGainDb(0U, NULL) > MIMO_GAIN_LIMIT_DB || MimoFx_PeakGainDb(1U, NULL) > MIMO_GAIN_LIMIT_DB) { diverged(); return; }

#if MIMO_GAIN_MAX_DB > 0
    /* limit of the filter gain, +- MIMO_GAIN_MAX_DB: while a filter exceeds it, its weights are pulled towards the neutral filter (delayed unit pulse), harder the further it is outside */
    {
        uint32_t dly = (uint32_t)(((uint32_t)(tauM[0] - shiftM[0]) + (uint32_t)(tauM[1] - shiftM[1])) / 2U);
        if (dly > LW - 1U) dly = LW - 1U;
        for (uint32_t j = 0U; j < J; j++)
        {
            float hi, lo;
            MimoFx_GainRangeDb(j, &hi, &lo);
            const float over = fmaxf(hi - (float)MIMO_GAIN_MAX_DB, (-(float)MIMO_GAIN_MAX_DB) - lo);
            if (over > 0.0f)
            {
                const float a = fminf(0.1f, 0.01f * over + 0.005f);
                for (uint32_t l = 0U; l < LW; l++) w[j][l] *= (1.0f - a);
                w[j][LW - 1U - dly] += a * (1.0f / (float)J);
            }
        }
    }
#endif

    publish();
}

/* Fill the work windows from the rings (filled by the audio interrupts) and run one block */
HOT static void processBlock(uint32_t blk)
{
    const uint32_t n0 = blk * B;

    for (uint32_t i = 0U; i < HIST + B; i++)
    {
        const uint32_t idx = (n0 - HIST + i) & RMASK;
        xw[i] = xRing[idx];
    }
    for (uint32_t m = 0U; m < M; m++)
        for (uint32_t i = 0U; i < B; i++) zw[m][i] = zRing[m][(n0 + i) & RMASK] * invG[m];
    runBlock();
}

int MimoFx_IrInstalled(void)               /* 1 once after a new h^ was installed */
{
    const int d = irNew;
    irNew = 0;
    return d;
}
float MimoFx_JumpDb(void) { const float r = jumpPeak; jumpPeak = 0.0f; return 20.0f * log10f(r + 1.0e-9f); }
uint32_t MimoFx_RtBlocks(void) { return rtBlocks; }
void  MimoFx_SetGains(float mu, float delta, float eps) { gMu = mu; gDelta = delta; gEps = eps; }
void  MimoFx_GetGains(float *mu, float *delta, float *eps) { *mu = gMu; *delta = gDelta; *eps = gEps; }
uint32_t MimoFx_Resets(void)                { return resets; }
float    MimoFx_MuScale(void)               { return muScale; }

uint32_t MimoFx_Poll(uint32_t maxBlocks)
{
    uint32_t done = 0U;
    while (done < maxBlocks)
    {
        uint32_t ready = refBlocks;
        for (uint32_t m = 0U; m < M; m++) if (micBlocks[m] < ready) ready = micBlocks[m];

        if (!haveIr || !adaptOn)                 /* nothing to adapt: do not build up a backlog */
        {
            pollBlock = ready;
            needInit  = 1;
            backlog   = 0U;
            break;
        }
        if (ready - pollBlock > (MIMO_RING - HIST) / B - 4U)         /* history would be overwritten */
        {
            pollBlock = ready - 1U;
            skipped++;
            needInit = 1;
        }
        backlog = ready - pollBlock;
        if (pollBlock >= ready) break;

        rtBlocks++;
        if ((rtBlocks & 255U) == 0U)
        {
            memcpy(wSnap, w, sizeof(w));                     /* where a divergence returns to */
        }
        processBlock(pollBlock);
        pollBlock++;
        blocksDone++;
        done++;
    }
    return done;
}
