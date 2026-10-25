/**
 * @file    room_eq.h
 * @brief   Adaptive multiple position room response equalizer.
 *
 * Signal flow
 *   reference x(n)  = input signal before the room equalizer (PCM1808 signal after the corrections, mono sum, with the output gain)
 *   microphones d_r = R = 2 measurement microphones
 *
 *   1. STFT (frame K = 4096, hop 3072, Hann window); the magnitude spectra of x and d_r are averaged in M bands.
 *   2. Per mic r and band m the weights H_r(m) [dB] are adapted frame by frame:
 *          e(m)   = d(m) - x(m)                       (dB, common level removed; the room equalizer in the signal path,
 *                                                      about -H, closes the loop)
 *          E{e}   = 20 log10( 1/L * sum_{i=n-L+1..n} 10^(e_i/20) )      (L = ROOMEQ_FRAMES_AVG)
 *          H(m)  += mu(m) * E{e}
 *      Band control: weights outside the equalized range are held. Gain control: |H| <= ROOMEQ_WEIGHT_MAX_DB.
 *   3. The weights give a constant magnitude per band; it is smoothed in 1/6 octave and averaged over the mics (linear
 *      magnitude). Regularized inverse: H_inv = H_p / (H_p^2 + beta).
 *   4. One peaking biquad per band centre realises the inverse magnitude; the gains come from a least-squares fit in dB
 *      (Gauss-Newton).
 *
 * Threading: RoomEq_Push*() are called from the audio interrupts, RoomEq_Poll() and the setters from the main loop;
 * the audio interrupt only reads the coefficients with RoomEq_CoefSeq() / RoomEq_Coefs().
 */
#ifndef ROOM_EQ_H
#define ROOM_EQ_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- settings ---- */
#ifndef ROOMEQ_FS_HZ
#define ROOMEQ_FS_HZ            44100.0
#endif
#define ROOMEQ_K                4096U         /* STFT frame length */
#define ROOMEQ_HOP              3072U         /* 25 % overlap */
#ifndef ROOMEQ_WITH_25HZ
#define ROOMEQ_WITH_25HZ        0             /* 1 = extra band / filter at 25 Hz (below 47 Hz) */
#endif
#ifndef ROOMEQ_WITH_84HZ
#define ROOMEQ_WITH_84HZ        1             /* 1 = extra band / filter at 84 Hz (between 47 and 140 Hz) */
#endif
#define ROOMEQ_BANDS            (20U + ROOMEQ_WITH_25HZ + ROOMEQ_WITH_84HZ)   /* M: 20 base centres + optional 25 Hz and 84 Hz */
#define ROOMEQ_MICS             2U            /* R */
#ifndef ROOMEQ_FRAMES_AVG
#define ROOMEQ_FRAMES_AVG       10U           /* N: frames in the error average */
#endif
#ifndef ROOMEQ_USE_REGULARIZATION
#define ROOMEQ_USE_REGULARIZATION 1          /* 1 = regularized inverse H/(H^2+beta), 0 = plain 1/H (beta is not used) */
#endif
#ifndef ROOMEQ_COMPLEX_INVERSE
#define ROOMEQ_COMPLEX_INVERSE  0             /* Foundation for a future FFT / FIR realisation of the equalizer instead of the peaking filter
                                               * bank, which only needs the magnitude. 1 = complex inverse conj(Hp)/(|Hp|^2+beta) with a phase
                                               * for Hp (hpPhase[], zero until a phase source exists); the magnitude is identical. */
#endif
#ifndef ROOMEQ_BETA
#define ROOMEQ_BETA             0.00001       /* regularization factor */
#endif
#define ROOMEQ_SMOOTH_OCT       (1.0 / 6.0)   /* magnitude smoothing width */
#define ROOMEQ_MAX_FILTERS      ROOMEQ_BANDS  /* one peaking filter per band centre */

#ifndef ROOMEQ_F_LOW_HZ
#define ROOMEQ_F_LOW_HZ         20.0          /* band control: equalized range (weights outside are held) */
#endif
#ifndef ROOMEQ_F_HIGH_HZ
#define ROOMEQ_F_HIGH_HZ        21000.0
#endif
#ifndef ROOMEQ_WEIGHT_MAX_DB
#define ROOMEQ_WEIGHT_MAX_DB    12.0          /* gain control: |H(m)| limit */
#endif
#ifndef ROOMEQ_FILTER_MAX_DB
#define ROOMEQ_FILTER_MAX_DB    15.0          /* limit of a single peaking filter */
#endif
#ifndef ROOMEQ_REF_RMS_MIN
#define ROOMEQ_REF_RMS_MIN      30.0          /* no adaptation below this reference level (16 bit LSB rms) */
#endif
#ifndef ROOMEQ_MIC_RMS_MIN
#define ROOMEQ_MIC_RMS_MIN      4.0           /* ... or if the mic hears only its own noise */
#endif
#ifndef ROOMEQ_REF_DELAY
#define ROOMEQ_REF_DELAY        0U            /* samples the reference is shifted into the past (acoustic + converter latency) */
#endif
#ifndef ROOMEQ_REF_BAND_REL_DB
#define ROOMEQ_REF_BAND_REL_DB  55.0          /* a band is only used if the reference has at least this level (dB) re its strongest band */
#endif
#ifndef ROOMEQ_SNR_MIN
#define ROOMEQ_SNR_MIN          3.0           /* ... and the mic hears it at least this much (linear) above its measured noise floor */
#endif
#ifndef ROOMEQ_MIN_VALID_BANDS
#define ROOMEQ_MIN_VALID_BANDS  4U            /* no adaptation of a mic in a frame with fewer valid bands */
#endif
#ifndef ROOMEQ_Q_BW_FACTOR
#define ROOMEQ_Q_BW_FACTOR      1.5           /* peaking filter bandwidth = factor * band spacing */
#endif

typedef struct { double b0, b1, b2, a1, a2; } RoomEqBiquad;       /* a0 = 1, transposed direct form II */
typedef struct { uint32_t n; RoomEqBiquad bq[ROOMEQ_MAX_FILTERS]; } RoomEqCoefs;

void     RoomEq_Init(void);                        /* once, at start-up (builds tables, flat equalizer) */
void     RoomEq_Reset(void);                       /* forget the room estimate, equalizer becomes flat */
void     RoomEq_SetAdapt(int on);
void     RoomEq_SetEqActive(int on);               /* 0 while the equalizer is bypassed: the loop is open, no adaptation (the weights would run away) */
int      RoomEq_GetAdapt(void);

/* Frequency response correction of a measurement mic, applied to the magnitude spectrum in the analysis
 * (equivalent to filtering the mic signal, but it costs nothing in the audio interrupt).
 * coef5 = n biquads as { b0, b1, b2, a1, a2 } (a0 = 1), e.g. the generated micEq_* arrays. */
void     RoomEq_SetMicCorrection(uint32_t mic, const double *coef5, uint32_t n);
void     RoomEq_SetMicCorrectionEnabled(int on);

/* Mic noise floor per band: learned from frames in which the reference is silent (output muted for a moment by
 * the caller, or a pause in the music). Without a known floor of both mics there is no adaptation. Bands where
 * a mic does not hear the reference clearly above its floor are not adapted (their weight stays as it is,
 * 0 dB = no correction at the start) - otherwise noise in weak bands would be "equalized" (treble loss). */
void     RoomEq_ForgetNoiseFloor(void);
int      RoomEq_NoiseFloorKnown(void);                 /* 1 when both mics have a floor */
int      RoomEq_BandValid(uint32_t mic, uint32_t band);
float    RoomEq_BandSnrDb(uint32_t mic, uint32_t band);
float    RoomEq_RefRelDb(uint32_t band);               /* reference level of the band re its strongest band */
float    RoomEq_ErrRmsDb(uint32_t mic);                /* how far the room estimate is off in the last frame: rms of the error over the used bands (dB);
                                                        * jumps when the room changes, then falls with the adaptation time */
uint32_t RoomEq_ValidBands(uint32_t mic);              /* bands used in the last frame; fewer than ROOMEQ_MIN_VALID_BANDS = no adaptation */

void     RoomEq_PushRef(const int16_t *s, uint32_t n);
void     RoomEq_PushMic(uint32_t mic, const int16_t *s, uint32_t n);
int      RoomEq_Poll(void);                        /* 0 = nothing to do, 1 = frame analysed but too quiet to learn, 2 = adapted (new coefficients) */

/* RBJ peaking filter (the same design as the equalizer filters), e.g. for manual notches */
void     RoomEq_PeakCoefs(double fc_hz, double gain_db, double q, RoomEqBiquad *out);

uint32_t            RoomEq_CoefSeq(void);          /* changes whenever new coefficients are available */
const RoomEqCoefs  *RoomEq_Coefs(void);

/* diagnostics */
double   RoomEq_BandCenter(uint32_t m);
double   RoomEq_BandStep(uint32_t m);
float    RoomEq_Weight(uint32_t mic, uint32_t band);   /* dB */
double   RoomEq_FilterFreq(uint32_t i);
double   RoomEq_FilterGainDb(uint32_t i);
uint32_t RoomEq_Frames(void);                      /* frames adapted so far */
float    RoomEq_RefRms(void);
float    RoomEq_MicRms(uint32_t mic);

#if ROOMEQ_COMPLEX_INVERSE
/* complex inverse on the design grid (log spaced), as calculated in the last adaptation step; returns the number of points */
uint32_t RoomEq_GetComplexInverse(double *freq, double *re, double *im);
#endif

#ifdef __cplusplus
}
#endif

#endif /* ROOM_EQ_H */
