/**
 * @file    mimo_fx.h
 * @brief   Multichannel room equalizer based on IPNLMS-FX, 2 loudspeakers x 2 microphones, adapted with the real microphone error.
 *
 * Signals (all float, +-1 = full scale):
 *   x(n)        mono reference (PCM1808)
 *   y_j(n)      = w_j^T x      output of adaptive filter j = the signal of loudspeaker j
 *   z_m(n)      corrected microphone signal m, divided by the gain G_m of the measured responses
 *   d_m(n)      = x(n - tau_m)                        desired signal
 *   e_m(n)      = d_m(n) - z_m(n)                     error of the real microphone
 *   x_f(m,j)    = h^_mj * x                           filtered reference (h^ at its true lag S_m)
 *
 * Two parts keep the interrupt short:
 *   - MimoFx_ProcessBlock()  audio interrupt (SAI block): 2 FIR filters with the published weights.
 *   - MimoFx_Poll()          main loop: x_f stream, adaptation, publishes the new weights.
 * Plain C without hardware access.
 */
#ifndef MIMO_FX_H
#define MIMO_FX_H

#include <stdint.h>
#include "mimo_config.h"

/* Measured room responses h^_mj = response from loudspeaker j to microphone m, cut out of the sweep measurement */
typedef struct
{
    float    h[MIMO_M][MIMO_J][MIMO_LH];    /* window of the response, divided by gain[m] */
    uint16_t shift[MIMO_M];                 /* S_m: bulk delay in front of the window (samples) */
    uint16_t tau[MIMO_M];                   /* tau_m: delay of the desired signal (samples) */
    float    gain[MIMO_M];                  /* G_m */
    int8_t   polarity[MIMO_M][MIMO_J];      /* sign (+1 / -1) of the strongest peak of the direct sound of every path */
} MimoIr;

/* Build MimoIr from the 4 measured responses: win[m][j] = MIMO_LH samples of the response starting at the lag MIMO_IR_START.
 * Returns 0 on success, -1 response too weak. */
int      MimoIr_Build(MimoIr *ir, const float *const win[MIMO_M][MIMO_J]);

void     MimoFx_Init(void);
int      MimoFx_SetIr(const MimoIr *ir);    /* install new h^, restart the adaptation. Main loop only. */
int      MimoFx_HaveIr(void);
void     MimoFx_FilterResponse(uint32_t j, float freqHz, float *magDb, float *phaseDeg);   /* response of w_j (main loop only) */
void     MimoFx_Reset(void);                /* weights back to the start value. Main loop only. */

/* Audio interrupt: one block of the mono reference x in, the loudspeaker signals out. */
void     MimoFx_ProcessBlock(const float *x, float out[MIMO_J][MIMO_BLOCK]);
/* Audio interrupt: one block of corrected microphone m (m = 0 .. M-1, not yet divided by G_m). */
void     MimoFx_PushMic(uint32_t m, const float *z);
/* Next block counter of the reference / microphone. */
uint32_t MimoFx_RefBlocks(void);
uint32_t MimoFx_MicBlocks(uint32_t m);

/* Main loop: process as many finished blocks as are available (at most maxBlocks). Returns the number processed. */
uint32_t MimoFx_Poll(uint32_t maxBlocks);

/* Control / status */
void     MimoFx_SetAdapt(int on);
int      MimoFx_GetAdapt(void);
float    MimoFx_Weight(uint32_t j, uint32_t l);          /* w_j[l], l = 0 .. MIMO_LW-1 */
float    MimoFx_WeightRms(uint32_t j);
float    MimoFx_PeakGainDb(uint32_t j, float *freqHz);   /* largest |W_j(f)| between 20 Hz and 20 kHz in dB and where (main loop only) */
void     MimoFx_GainRangeDb(uint32_t j, float *maxDb, float *minDb);   /* largest / smallest gain of the played filter, 20 Hz .. 20 kHz (main loop only) */
float    MimoFx_ErrDb(uint32_t m);          /* real mic error e_m relative to d_m, dB (smoothed) */
uint32_t MimoFx_BlocksDone(void);
int      MimoFx_IrInstalled(void);           /* 1 once after a new h^ was installed */
uint32_t MimoFx_Resets(void);                /* divergence resets so far (each halves the step size) */
void     MimoFx_SetGains(float mu, float delta, float eps);   /* gain parameters of the IPNLMS, effective immediately */
void     MimoFx_GetGains(float *mu, float *delta, float *eps);
float    MimoFx_JumpDb(void);                /* largest change of the played weights between two publishes since the last call, 20 log10(|dw| / |w|) dB */
uint32_t MimoFx_RtBlocks(void);              /* blocks adapted so far */
float    MimoFx_MuScale(void);               /* step size factor (reduced after a divergence reset) */
uint32_t MimoFx_Skipped(void);              /* blocks dropped because the main loop was too slow */
uint32_t MimoFx_Backlog(void);              /* blocks waiting right now */
const MimoIr *MimoFx_GetIr(void);

#endif /* MIMO_FX_H */
