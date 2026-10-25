/**
 * @file    mic_phase.h
 * @brief   Phase (and residual magnitude) correction of the analog microphones: FIR after the magnitude biquads, delay = taps / 2 samples.
 *          Coefficients: micPh_adcL / micPh_adcR in mic_eq_coeffs.h, handed over by MicPhase_Init(). With 0 taps the stage is a pass-through.
 *          Plain C without hardware access.
 */
#ifndef MIC_PHASE_H
#define MIC_PHASE_H

#include <stdint.h>

#define MIC_PH_CHANNELS 2U          /* 0 = PF6 / ADC3 ("adcL"), 1 = PA3 / ADC2 ("adcR") */

#define MIC_PH_MAX_TAPS 256U

/* c0 / n0: FIR of channel 0 (PF6 / ADC3), c1 / n1: channel 1 (PA3 / ADC2); newest tap last is NOT required: the first coefficient is the first tap. Copies them. */
void MicPhase_Init(const float *c0, uint32_t n0, const float *c1, uint32_t n1);
/* Filter one block in place (n samples, n <= 1024) of channel ch. */
void MicPhase_Process(uint32_t ch, float *buf, uint32_t n);
/* Number of taps and the delay (samples) of the stage, for the printout */
uint32_t MicPhase_Taps(uint32_t ch);

#endif /* MIC_PHASE_H */
