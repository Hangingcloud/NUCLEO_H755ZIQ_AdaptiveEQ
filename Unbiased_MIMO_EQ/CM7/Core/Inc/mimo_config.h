/**
 * @file    mimo_config.h
 * @brief   Settings of the multichannel room equalizer (2 loudspeakers, 2 microphones).
 *
 * The impulse responses h^ are not adapted, they are measured with an exponential (Farina) sweep.
 * Plain C, no hardware access.
 */
#ifndef MIMO_CONFIG_H
#define MIMO_CONFIG_H

#ifndef MIMO_FS_HZ
#define MIMO_FS_HZ              44100U
#endif

/* ------------------------------------------------------------------------------------------------
 * System size
 * ---------------------------------------------------------------------------------------------- */
#define MIMO_J                  2U      /* loudspeakers = adaptive filters w_j (I2S left = 1, I2S right = 2) */
#define MIMO_M                  2U      /* microphones m = 1, 2 */

/* Block length of the audio callbacks (frames per DMA half buffer) */
#define MIMO_BLOCK              256U

/* ------------------------------------------------------------------------------------------------
 * Filter lengths, both must be multiples of 4.
 * ---------------------------------------------------------------------------------------------- */
#ifndef MIMO_LW
#define MIMO_LW                 256U    /* adaptive filters w_j (2 x LW multiply-adds per sample in the audio interrupt) */
#endif
#ifndef MIMO_LH
#define MIMO_LH                 256U    /* estimated impulse responses h^_mj (window of the measured response) */
#endif

/* ------------------------------------------------------------------------------------------------
 * Adaptation
 * ---------------------------------------------------------------------------------------------- */
#ifndef MIMO_KAPPA
#define MIMO_KAPPA              (-0.5f) /* -1 = NLMS ... +1 = fully proportionate */
#endif

/* The weights are updated every MIMO_ADAPT_DECIM samples; the x_f stream and the filtering run at every sample.
 * Larger values save CPU and keep the loop stable for strongly low-passed input. mu is not scaled by the decimation. */
#ifndef MIMO_ADAPT_DECIM
#define MIMO_ADAPT_DECIM        8U      /* must divide MIMO_BLOCK */
#endif

/* Button B1: hold it for MIMO_CAL_HOLD_MS to start the calibration; it starts MIMO_CAL_COUNTDOWN_S seconds after the button is released
 * (LD2 blinks slowly), then the sweeps run (LD2 blinks fast), then LD2 stays lit while the adaptation runs.
 * A rejected measurement is shown by LD3 blinking. A short press pauses / continues the adaptation. */
#ifndef MIMO_CAL_HOLD_MS
#define MIMO_CAL_HOLD_MS        3000U
#endif
#ifndef MIMO_CAL_COUNTDOWN_S
#define MIMO_CAL_COUNTDOWN_S    3U
#endif

/* Limit of the filter gain (20 Hz .. 20 kHz): kept between -MIMO_GAIN_MAX_DB and +MIMO_GAIN_MAX_DB (absolute, the neutral filter sits at -6 dB).
 * A filter outside the limit is pulled smoothly towards the neutral filter. 0 = no limit. Must be smaller than MIMO_GAIN_LIMIT_DB. */
#ifndef MIMO_GAIN_MAX_DB
#define MIMO_GAIN_MAX_DB        0
#endif

/* Safety: if the largest gain of a filter w_j (20 Hz .. 20 kHz, checked once per second) exceeds this, the weights are reset to the start value
 * and a message is printed. 0 = off. */
#ifndef MIMO_GAIN_LIMIT_DB
#define MIMO_GAIN_LIMIT_DB      20.0f
#endif

/* ------------------------------------------------------------------------------------------------
 * Signals
 * ---------------------------------------------------------------------------------------------- */
/* Reference x(n) from the PCM1808: 0 = (L+R)/2, 1 = left only, 2 = right only */
#define MIMO_REF_MODE           0

/* Which ADC input is microphone 1 / 2. Microphone 1 pairs with loudspeaker 1 (left), microphone 2 with loudspeaker 2 (right).
 * ADC index 0 = PF6 (ADC3, Zio pin A6), index 1 = PA3 (ADC2, Zio pin A0). */
#define MIMO_MIC1_ADC           1U
#define MIMO_MIC2_ADC           0U

/* Polarity of the microphone signals: -1 inverts both microphones. Set -1 if the first peak of the measured responses is negative. */
#ifndef MIMO_MIC_POLARITY
#define MIMO_MIC_POLARITY       (1.0f)
#endif

/* Desired signal d_m(n) = x(n - tau_m), tau_m = MIMO_IR_START + MIMO_MODEL_DELAY (the same for both microphones). */
#ifndef MIMO_MODEL_DELAY
#define MIMO_MODEL_DELAY        128U    /* delay of the desired signal after the start of the window of h^ (samples) */
#endif

/* ------------------------------------------------------------------------------------------------
 * Farina sweep measurement of h^
 * x(t) = sin[ w1 T / ln(w2/w1) * (exp(t/T ln(w2/w1)) - 1) ],  h(t) = IFFT[ FFT(y_avg(t)) / FFT(x(t)) ],
 * y_avg = average of MIMO_SWEEP_AVG recordings. Every loudspeaker/microphone combination is measured on its own.
 * One recording is 2^MIMO_SWEEP_LOG2N samples long (sweep + decay of the room). The FFT buffer takes
 * 8 x 2^LOG2N bytes: 256 KB at 15, 512 KB at 16 does not fit into the AXI SRAM.
 * ---------------------------------------------------------------------------------------------- */
#define MIMO_SWEEP_F1_HZ        20.0f   /* start of the sweep */
#define MIMO_SWEEP_F2_HZ        20000.0f
#define MIMO_SWEEP_T_S          0.5f    /* sweep duration T in seconds (shorter than the record) */
#define MIMO_SWEEP_LOG2N        15U     /* record / FFT length N = 2^15 = 32768 samples (0.74 s) */
#ifndef MIMO_SWEEP_AVG
#define MIMO_SWEEP_AVG          4U      /* recordings averaged per combination */
#endif
#define MIMO_SWEEP_LEVEL_DB     (-20.0f) /* sweep level re full scale at the reference (before the output gain) */
#ifndef MIMO_SWEEP_BW_LOW
#define MIMO_SWEEP_BW_LOW       0.5f    /* below F1 the response is faded out between MIMO_SWEEP_BW_LOW * F1 (zero) and F1 (full) */
#endif
#define MIMO_SWEEP_FADE_MS      10.0f   /* raised cosine fade in / out of the sweep */
#define MIMO_SWEEP_GAP_BLOCKS   172U    /* pause before every recording in blocks of MIMO_BLOCK (172 blocks = 1.0 s) */
#define MIMO_SWEEP_REG          1.0e-2f /* regularization of the division: |X|^2 + REG * (|X|^2 at the top of the band) */

/* Window of the measured response used as h^: MIMO_LH samples starting at the lag MIMO_IR_START of the (circular) deconvolution, the same for
 * all four paths so that the delay between the loudspeakers stays in h^. MIMO_IR_START must contain the audio chain delay and the sound travel
 * time, with a small margin in front. After a calibration the printout "response energy starts at lag ..." shows where the response begins:
 * set MIMO_IR_START about 32 samples before that, as a multiple of MIMO_ADAPT_DECIM. */
#ifndef MIMO_IR_START
#define MIMO_IR_START           864U
#endif
/* A measured response is only installed if its peak is at least this far above the noise of the record (every path). */
#ifndef MIMO_IR_MIN_SNR_DB
#define MIMO_IR_MIN_SNR_DB      20.0f
#endif
/* Below this frequency (Hz) the measured h^ is replaced by the mid band level (500 .. 2000 Hz) with the phase of the direct sound,
 * faded back over one octave above it. 0 = off (the measured bass is used). */
#ifndef MIMO_IR_LF_IGNORE_HZ
#define MIMO_IR_LF_IGNORE_HZ    150
#endif

/* Smoothing of the measured responses h^ (complex smoothing) in 1/N octave: 6 = 1/6 octave, 3 = 1/3, 1 = one octave, 0 = off. */
#ifndef MIMO_SMOOTH_OCT_DIV
#define MIMO_SMOOTH_OCT_DIV     6
#endif

/* Frequencies below this (Hz) are not corrected (0 = off). After every adapted block the deviation of the weights from the neutral filter
 * (the delayed pulse the adaptation starts from) is projected so that its frequency response is zero at 0, 20, 40 ... up to this frequency
 * (linearly constrained LMS). */
#ifndef MIMO_LF_NEUTRAL_HZ
#define MIMO_LF_NEUTRAL_HZ      150
#endif

/* Gain parameters of the adaptation, start values (changeable on the serial port while it runs). */
#ifndef MIMO_RT_MU
#define MIMO_RT_MU              0.029f
#endif
#ifndef MIMO_RT_DELTA
#define MIMO_RT_DELTA           2.58e-12f
#endif
#ifndef MIMO_RT_EPSILON
#define MIMO_RT_EPSILON         1.0e-3f
#endif

#define MIMO_IR_FADE            32U     /* the last taps of h^ are faded out (truncation) */
#ifndef MIMO_MAX_SHIFT
#define MIMO_MAX_SHIFT          1536U   /* largest bulk delay S_m that the history can hold */
#endif
#ifndef MIMO_MAX_TAU
#define MIMO_MAX_TAU            2048U   /* largest tau_m */
#endif

/* Gain normalization: h^ is divided by G_m (RMS of the measured responses of microphone m) and the mic signal is divided by the same G_m,
 * so that d and the mic signal are comparable regardless of speaker volume and mic gain. 0 = off. */
#define MIMO_IR_NORMALIZE       1

/* ------------------------------------------------------------------------------------------------
 * Derived sizes, do not edit
 * ---------------------------------------------------------------------------------------------- */
#define MIMO_SWEEP_N            (1UL << MIMO_SWEEP_LOG2N)
#define MIMO_HIST               3072U   /* samples of history in the work windows of the adaptation (>= max(S + LH, tau)) */
#define MIMO_RING               8192U   /* ring buffers between interrupts and main loop (power of two) */

#if (MIMO_LW % 4U) || (MIMO_LH % 4U)
#error "MIMO_LW and MIMO_LH must be multiples of 4"
#endif
#if (MIMO_ADAPT_DECIM % 4U) || (MIMO_J != 2U) || (MIMO_M != 2U)
#error "the kernels need MIMO_ADAPT_DECIM a multiple of 4 and J = M = 2"
#endif
#if (MIMO_BLOCK % MIMO_ADAPT_DECIM)
#error "MIMO_ADAPT_DECIM must divide MIMO_BLOCK"
#endif
#if (MIMO_MAX_SHIFT + MIMO_LH > MIMO_HIST) || (MIMO_MAX_TAU > MIMO_HIST)
#error "MIMO_HIST too small for MIMO_MAX_SHIFT / MIMO_LH / MIMO_MAX_TAU"
#endif
#if (MIMO_IR_START % MIMO_ADAPT_DECIM) || (MIMO_IR_START > MIMO_MAX_SHIFT)
#error "MIMO_IR_START must be a multiple of MIMO_ADAPT_DECIM and not larger than MIMO_MAX_SHIFT"
#endif

#endif /* MIMO_CONFIG_H */
