/**
 * @file    mic_phase.c
 * @brief   Phase correction FIR of the microphones, see mic_phase.h.
 */
/* Real time code: always compiled with -O2, independent of the IDE's optimization setting. */
#pragma GCC optimize ("O2")
#include "mic_phase.h"
#include <string.h>

#ifdef MIMO_HOST
#define DTCM
#else
#define DTCM __attribute__((section(".dtcm"), aligned(8)))
#endif

#define MAXN 1024U

static struct { const float *c; uint32_t n; } table[MIC_PH_CHANNELS];

DTCM static float coef[MIC_PH_CHANNELS][MIC_PH_MAX_TAPS];                 /* newest tap last */
DTCM static float hist[MIC_PH_CHANNELS][MIC_PH_MAX_TAPS - 1U + MAXN];    /* previous samples + block */

void MicPhase_Init(const float *c0, uint32_t n0, const float *c1, uint32_t n1)
{
    table[0].c = c0;  table[0].n = (n0 > MIC_PH_MAX_TAPS) ? MIC_PH_MAX_TAPS : n0;
    table[1].c = c1;  table[1].n = (n1 > MIC_PH_MAX_TAPS) ? MIC_PH_MAX_TAPS : n1;
    memset(hist, 0, sizeof(hist));
    memset(coef, 0, sizeof(coef));
    for (uint32_t ch = 0U; ch < MIC_PH_CHANNELS; ch++)
        for (uint32_t k = 0U; k < table[ch].n; k++)
            coef[ch][table[ch].n - 1U - k] = table[ch].c[k];     /* reversed: oldest sample first */
}

uint32_t MicPhase_Taps(uint32_t ch) { return table[ch].n; }

void MicPhase_Process(uint32_t ch, float *buf, uint32_t n)
{
    const uint32_t taps = table[ch].n;
    if (taps == 0U || n > MAXN) return;                      /* pass-through */

    float *h = hist[ch];
    const uint32_t keep = taps - 1U;
    const float *c = coef[ch];
    memcpy(&h[keep], buf, n * sizeof(float));
    uint32_t i = 0U;
    for (; i + 8U <= n; i += 8U)                              /* 8 outputs per pass: every coefficient is loaded once for 8 multiply-adds */
    {
        const float *x = &h[i];
        float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f, s3 = 0.0f, s4 = 0.0f, s5 = 0.0f, s6 = 0.0f, s7 = 0.0f;
        for (uint32_t k = 0U; k < taps; k++)
        {
            const float ck = c[k];
            s0 += ck * x[k];       s1 += ck * x[k + 1U];  s2 += ck * x[k + 2U];  s3 += ck * x[k + 3U];
            s4 += ck * x[k + 4U];  s5 += ck * x[k + 5U];  s6 += ck * x[k + 6U];  s7 += ck * x[k + 7U];
        }
        buf[i] = s0;  buf[i + 1U] = s1;  buf[i + 2U] = s2;  buf[i + 3U] = s3;
        buf[i + 4U] = s4;  buf[i + 5U] = s5;  buf[i + 6U] = s6;  buf[i + 7U] = s7;
    }
    for (; i < n; i++)
    {
        const float *x = &h[i];
        float s0 = 0.0f;
        for (uint32_t k = 0U; k < taps; k++) s0 += c[k] * x[k];
        buf[i] = s0;
    }
    memmove(h, &h[n], keep * sizeof(float));
}
