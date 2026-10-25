/**
 * @file    farina.h
 * @brief   Exponential sweep (Farina) generator and impulse response identification.
 *
 *   x(t) = sin[ w1 T / ln(w2/w1) * ( exp( t/T * ln(w2/w1) ) - 1 ) ]
 *   h(t) = IFFT[ FFT(y_avg(t)) / FFT(x(t)) ]
 *
 * y_avg is the average of several recordings of the same sweep, formed in the time domain while the samples
 * arrive (one buffer). The division is regularized outside of the sweep band. Plain C without hardware access.
 */
#ifndef FARINA_H
#define FARINA_H

#include <stdint.h>

typedef struct
{
    double   e;         /* exp(n * r), advanced by multiplication */
    uint32_t n;         /* sample index of the next call */
} FarinaGen;

/* Build the cosine table of the FFT (once, at start-up). */
void     Farina_Init(void);

/* Sweep length in samples (T * fs, rounded) */
uint32_t Farina_SweepLen(void);

/* Sequential generator: returns x(n) * level, with fade in / out; 0 after the end of the sweep. */
void     Farina_GenInit(FarinaGen *g, uint32_t n0);
float    Farina_GenNext(FarinaGen *g);

/* Averaging of the recordings: Clear -> Accumulate (blocks, any order) x MIMO_SWEEP_AVG -> Deconvolve. */
void     Farina_Clear(void);
void     Farina_Accumulate(uint32_t n0, const float *z, uint32_t cnt);   /* y_avg[n0 + i] += z[i] / MIMO_SWEEP_AVG */

/* Before the deconvolution: copy of y_avg as int16 (y = value * scale) and statistics of the second half of the record (mean, rms) and the peak, in full scale units. */
void     Farina_StoreRecord(int16_t *dst, float *scale, float *tailMean, float *tailRms, float *peak);

/* x -> FFT, divide, IFFT. */
void     Farina_Deconvolve(void);

/* Result of the deconvolution: h(n), n = 0 .. MIMO_SWEEP_N - 1 (circular: the response follows the sweep latency). */
float    Farina_Ir(int32_t lag);

/* Locate the strongest response in the record: largest |h| at lags -2048 .. N/2 -> lag and value, and the signal to noise ratio =
 * peak / rms of the quarter 1/4 .. 1/2 of the record. */
void     Farina_Locate(int32_t *peakLag, float *peakVal, float *snrDb);

/* Lag of the first 32 sample block whose energy is above 3 % of the strongest block (lags -1024 .. +1536). */
int32_t  Farina_EnergyStart(void);

/* Approximate size (bytes) of the static buffer, for printouts. */
uint32_t Farina_BufferBytes(void);

#endif /* FARINA_H */
