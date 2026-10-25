/* USER CODE BEGIN Header */
/*
 * @file           : main.c
 * @brief          : MIMO room equalizer based on IPNLMS-FX
 *
 * Signal path:
 *   PCM1808 (SAI1 A) -> DC block -> PCM1808 correction (biquads) -> mono reference x(n)
 *       -> adaptive filters w_1, w_2 -> I2S output: left = loudspeaker 1, right = loudspeaker 2
 *   mic 1 (PF6 / ADC3), mic 2 (PA3 / ADC2) -> DC block -> mic correction (biquads, magnitude) -> phase FIR
 *       -> z_m(n)
 *   The impulse responses h^_mj are measured with an exponential sweep (button B1 or UART key 's'), every loudspeaker / microphone
 *   combination on its own. The adaptation (main loop) is in mimo_fx.c, the sweep in farina.c.
 */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "adc.h"
#include "bdma.h"
#include "dma.h"
#include "i2s.h"
#include "sai.h"
#include "usart.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mic_eq_coeffs.h"   /* generated filter coefficients */
#include "mic_phase.h"       /* phase / magnitude FIR of the microphones (coefficients: mic_eq_coeffs.h) */
#include "mimo_config.h"
#include "mimo_fx.h"
#include "farina.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */
/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
/* I2S/SAI capture buffer sizes: 16-bit words, circular DMA fills it continuously.
 * Half-complete callback fires when the first half is full, full-complete
 * when the second half is full -- classic double-buffer / ping-pong. */
#define I2S_HALF_FRAMES         256U                         /* stereo frames per half-buffer */
#define I2S_HALF_SAMPLES        (I2S_HALF_FRAMES * 2U)       /* L+R samples per half-buffer */
#define I2S_BUFFER_SAMPLES      (I2S_HALF_SAMPLES * 2U)      /* both halves */
#if (I2S_HALF_FRAMES != MIMO_BLOCK)
#error "MIMO_BLOCK must be equal to I2S_HALF_FRAMES"
#endif

/* Gain of the PCM1808 input in dB (applied to the reference x), any value (e.g. 17.5f, negative = quieter). */
#define MIC_GAIN_DB             2.0f        /* PCM1808, SAI1 A */
/* Roles: the PCM1808 (SAI1 A: PE2 SCKI, PE4 LRCK, PE5 BCK, PE6 DOUT) is the audio SOURCE (mono reference x), the two
 * analog mics (PF6 = ADC3_INP8, PA3 / A0 = ADC2_INP15) are the measurement microphones of the model.
 * B1 (PC13): short press = the adaptation pauses / continues (LD2 lit = running),
 *            hold 3 s and let go = the sweep calibration without a PC: countdown, sweeps (LD2 shows the state, see Cal_Leds()).
 * Serial terminal (USART3, 115200): see Gain_HandleUart(). */
#define MIC_SRC_I2S             0U       /* PCM1808: source, played */
#define MIC_SRC_ADC             1U       /* analog mics: only measured */
#define MIC_SRC_COUNT           2U
#define ADC_HALF_SAMPLES        I2S_HALF_FRAMES              /* per ADC channel, same block length as SAI */
/* Master output volume in dB (after the filters, it is part of the measured loop).
 * DAC full scale = 2.1 Vrms per channel. */
#define OUTPUT_GAIN_DB          -5.0f
/* DAC output: I2S1 (PB3 BCK, PA15 LRCK, PB5 DIN), plain stereo: left channel = loudspeaker 1, right channel = loudspeaker 2. */
#define hi2sOut                 hi2s1
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */
/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
/* Buffer filled by the PCM1808 (SAI RX DMA), stereo interleaved L,R */
int16_t micRxBuffer[I2S_BUFFER_SAMPLES];

/* Buffer read by the Speaker/DAC (I2S TX DMA), stereo interleaved L,R */
int16_t i2sTxBuffer[I2S_BUFFER_SAMPLES];

volatile uint32_t rxCount = 0, txCount = 0;
volatile int16_t  micPeak = 0;              /* peak |sample| of last block, for UART debug */
volatile int16_t  rawMin[2], rawMax[2];     /* raw SAI data per slot (L,R) of last block */
volatile uint32_t saiSr = 0;                /* SAI1 Block A status register, sticky errors */

/* Analog mics: channel 0 = ADC3 / PF6, channel 1 = ADC2 / PA3.
 * ADC3 is served by the BDMA, which can only reach SRAM4 -> .sram4 section (see linker script).
 * ADC2 is served by DMA1 (stream 3), which reaches the normal AXI SRAM. */
__attribute__((section(".sram4"), aligned(32))) uint16_t adcRxBuffer[ADC_HALF_SAMPLES * 2U];
__attribute__((aligned(32))) uint16_t adc2RxBuffer[ADC_HALF_SAMPLES * 2U];
/* The two ADCs deliver their halves in separate callbacks: collect them here, process when both are in */
static int16_t  adcPcm[2][2][ADC_HALF_SAMPLES];     /* [half][channel][sample], signed */
static uint8_t  adcReady[2];                        /* per half: bit 0 = ADC3 in, bit 1 = ADC2 in */
volatile uint32_t adcCount = 0;
volatile uint16_t adcMin[2], adcMax[2];             /* raw ADC range of last block per channel, for UART debug */

static uint32_t lastPrint = 0;
/* Per source (index = MIC_SRC_xxx) */
static int32_t        srcDc[MIC_SRC_COUNT][2];           /* DC blocker states per channel, Q8 */
static float          pcmGain;                           /* linear gain of the PCM1808 reference, from MIC_GAIN_DB */
static float          outGain = 1.0f;                    /* master volume, from OUTPUT_GAIN_DB */
static volatile float lvlSumSq[MIC_SRC_COUNT][2];        /* sum of squares BEFORE gain per channel, for the level meter */
static volatile uint32_t lvlCount[MIC_SRC_COUNT];

/* Frequency-response correction: one biquad cascade per source and channel (coefficients in mic_eq_coeffs.h).
 * Sections with poles close to the unit circle run in double precision. Both the PCM1808 and the analog mics are corrected in the time domain. */
static const struct { const MicEqBiquad *bq; uint32_t n; } eqTable[MIC_SRC_COUNT][2] = {
    { { micEq_pcmL, MIC_EQ_PCML_N }, { micEq_pcmR, MIC_EQ_PCMR_N } },    /* MIC_SRC_I2S */
    { { micEq_adcL, MIC_EQ_ADCL_N }, { micEq_adcR, MIC_EQ_ADCR_N } },    /* MIC_SRC_ADC */
};
/* The D-cache is off: coefficients and states live in DTCM (zero wait states, .dtcm is NOLOAD, the coefficients are copied from flash in Eq_Init()).
 * Sections whose poles are very close to the unit circle run in double precision, all others in single precision (a double biquad step costs about 40 cycles on the Cortex-M7 FPU). */
#define EQ_DOUBLE_RADIUS        0.99        /* pole radius above which a section stays in double precision */
__attribute__((section(".dtcm"), aligned(8))) static MicEqBiquad eqCoef[MIC_SRC_COUNT][2][MIC_EQ_MAX_SECTIONS];                    /* double sections */
__attribute__((section(".dtcm"), aligned(8))) static double      eqState[MIC_SRC_COUNT][2][MIC_EQ_MAX_SECTIONS][2];   /* transposed DF2: z1, z2 */
__attribute__((section(".dtcm"), aligned(8))) static float       eqCoefF[MIC_SRC_COUNT][2][MIC_EQ_MAX_SECTIONS][5];   /* single precision sections: b0 b1 b2 a1 a2 */
__attribute__((section(".dtcm"), aligned(8))) static float       eqStateF[MIC_SRC_COUNT][2][MIC_EQ_MAX_SECTIONS][2];
static uint32_t eqNd[MIC_SRC_COUNT][2], eqNf[MIC_SRC_COUNT][2];                                                                   /* number of double / float sections */
volatile uint32_t procCyclesMax = 0;    /* worst time of the PCM1808 block (incl. the output filters) */
volatile uint32_t procMicCyclesMax = 0; /* worst time of the mic block */
volatile uint32_t pollCyclesMax = 0;    /* worst time of one adaptation block (main loop) */

static volatile uint8_t errSai = 0, errAdc = 0, errI2s = 0;   /* set by the error callbacks, printed in the main loop */
static volatile uint8_t bypass = 0;     /* 1 = equalizer off: the corrected PCM1808 stereo signal is played (B1 hold / 'q' / overload) */
static uint8_t  calPending = 0U;                             /* the calibration (sweep + training) starts when calStartTick is reached */
static uint32_t calStartTick = 0U;
static uint32_t failBlinkUntil = 0U;                         /* LD3 blinks until then: the last measurement was rejected */
static uint32_t ackUntil = 0U;                               /* LD3 lit until then: B1 was held long enough */
static volatile uint8_t shed = 0;       /* set when a block took too long: the equalizer was switched off */
static volatile uint32_t shedCycles = 0;
static volatile uint32_t clipCount = 0;                     /* output samples that were clipped */
static volatile int16_t  outPeak[2];                        /* peak |sample| to the loudspeakers since the last printout */

/* Sweep measurement (shared between the audio interrupts and the main loop). One combination = loudspeaker j,
 * microphone m, MIMO_SWEEP_AVG recordings. Both interrupts count blocks the same way (MimoFx_RefBlocks / MicBlocks):
 * the sweep starts at the first sample of reference block startBlock and the recording at the first sample of microphone block startBlock. */
typedef struct
{
    volatile uint8_t  armed;        /* a sweep is waiting / playing / recording */
    volatile uint8_t  playing;      /* SAI interrupt: sweep still running */
    volatile uint8_t  recDone;      /* mic interrupt: record complete */
    volatile uint32_t startBlock;
    volatile uint32_t recPos;       /* samples recorded */
    uint8_t           j, m;         /* loudspeaker, microphone (logical, 0 .. 1) */
    FarinaGen         gen;
} MeasCtl;
static MeasCtl meas;

static struct
{
    uint8_t  active;
    uint8_t  step;                  /* 0 arm, 1 wait for the recording, 2 deconvolve */
    uint8_t  combo, rep;
    uint32_t t0;
} ms;
/* Averaged sweep recordings of the 4 combinations (int16, 64 KB each) in the D2 SRAM (see the linker script), for the 'y' dump */
__attribute__((section(".ram_d2"), aligned(8))) static int16_t recStore[MIMO_J * MIMO_M][MIMO_SWEEP_N];
static float recScale[MIMO_J * MIMO_M];
static uint8_t recValid[MIMO_J * MIMO_M];
static uint8_t recNext;
static float fullDb[MIMO_M][MIMO_J][16];                  /* spectrum of the long measured response (0.37 s from the window start MIMO_IR_START), dB at the 16 print frequencies */
static const float fullFreq[16] = { 20.f, 31.f, 50.f, 80.f, 125.f, 200.f, 315.f, 500.f, 800.f, 1250.f, 2000.f, 3150.f, 5000.f, 8000.f, 12500.f, 19000.f };
static float irWin[MIMO_M][MIMO_J][MIMO_LH];                    /* measured responses of the 4 combinations: the window of MIMO_LH samples at the lag MIMO_IR_START */
static int32_t irPeakLag[MIMO_M][MIMO_J];                        /* information: lag of the strongest sample, lag where the energy starts */
static int32_t irEnergyStart[MIMO_M][MIMO_J];
static float    irPeakVal[MIMO_M][MIMO_J], irSnr[MIMO_M][MIMO_J];
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
void PeriphCommonClock_Config(void);
/* USER CODE BEGIN PFP */
void UART_Print(char *str);
void Error_Handler(void);
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/**
 * @brief  Return the half of i2sTxBuffer that the I2S DMA is NOT playing right now.
 *
 * SAI1, I2S and TIM6 (ADC trigger) are all derived from PLL1, so every source
 * delivers exactly one block per TX half-buffer. The sources are started
 * while TX is in the middle of a half (Audio_WaitTxMidHalf), so this decision
 * always has half a block of margin.
 */
static int16_t *Audio_TxHalfToFill(void)
{
    uint32_t txRemaining = __HAL_DMA_GET_COUNTER(hi2sOut.hdmatx);        /* counts down */
    return (txRemaining > I2S_HALF_SAMPLES) ? &i2sTxBuffer[I2S_HALF_SAMPLES]  /* TX plays half 0 */
                                            : &i2sTxBuffer[0];                /* TX plays half 1 */
}

static inline __attribute__((always_inline)) float Audio_Clip(float v)
{
    if (v >  32767.0f) return  32767.0f;
    if (v < -32768.0f) return -32768.0f;
    return v;
}

/**
 * @brief  DC blocker, level measurement and frequency-response correction of one block of a stereo source.
 * @param  source  MIC_SRC_xxx
 * @param  srcL    first left sample (PCM: left slot, ADC: PF6 / ADC3)
 * @param  srcR    first right sample (PCM: right slot, ADC: PA3 / ADC2)
 * @param  stride  distance between samples of one channel (2 = interleaved stereo, 1 = separate buffers)
 * @param  outL    corrected samples in int16 units (float), not clipped
 */
__attribute__((optimize("O2")))
static void Eq_Block(uint32_t source, const int16_t *srcL, const int16_t *srcR, uint32_t stride,
                     float *outL, float *outR)
{
    int32_t dcL = srcDc[source][0];
    int32_t dcR = srcDc[source][1];
    float   sumL = 0.0f, sumR = 0.0f;
    const MicEqBiquad *bqL = eqCoef[source][0];
    const MicEqBiquad *bqR = eqCoef[source][1];
    const uint32_t     nL  = eqNd[source][0];
    const uint32_t     nR  = eqNd[source][1];
    double (*zL)[2]        = eqState[source][0];
    double (*zR)[2]        = eqState[source][1];
    const float (*cFL)[5]  = eqCoefF[source][0];
    const float (*cFR)[5]  = eqCoefF[source][1];
    const uint32_t     nFL = eqNf[source][0];
    const uint32_t     nFR = eqNf[source][1];
    float (*zFL)[2]        = eqStateF[source][0];
    float (*zFR)[2]        = eqStateF[source][1];

    for (uint32_t i = 0; i < I2S_HALF_FRAMES; i++)
    {
        int32_t xl = srcL[i * stride];
        int32_t xr = srcR[i * stride];

        /* one-pole DC blocker: dc += (x - dc) / 256 */
        dcL += ((xl << 8) - dcL) >> 8;
        dcR += ((xr << 8) - dcR) >> 8;
        float dl = (float)(xl - (dcL >> 8));
        float dr = (float)(xr - (dcR >> 8));
        sumL += dl * dl;
        sumR += dr * dr;

        double ql = (double)dl;
        double qr = (double)dr;
        for (uint32_t k = 0; k < nL; k++)
        {
            double y = bqL[k].b0 * ql + zL[k][0];
            zL[k][0] = bqL[k].b1 * ql - bqL[k].a1 * y + zL[k][1];
            zL[k][1] = bqL[k].b2 * ql - bqL[k].a2 * y;
            ql = y;
        }
        for (uint32_t k = 0; k < nR; k++)
        {
            double y = bqR[k].b0 * qr + zR[k][0];
            zR[k][0] = bqR[k].b1 * qr - bqR[k].a1 * y + zR[k][1];
            zR[k][1] = bqR[k].b2 * qr - bqR[k].a2 * y;
            qr = y;
        }
        float fl = (float)ql;
        float fr = (float)qr;
        for (uint32_t k = 0; k < nFL; k++)                       /* single precision sections */
        {
            const float y = cFL[k][0] * fl + zFL[k][0];
            zFL[k][0] = cFL[k][1] * fl - cFL[k][3] * y + zFL[k][1];
            zFL[k][1] = cFL[k][2] * fl - cFL[k][4] * y;
            fl = y;
        }
        for (uint32_t k = 0; k < nFR; k++)
        {
            const float y = cFR[k][0] * fr + zFR[k][0];
            zFR[k][0] = cFR[k][1] * fr - cFR[k][3] * y + zFR[k][1];
            zFR[k][1] = cFR[k][2] * fr - cFR[k][4] * y;
            fr = y;
        }
        outL[i] = fl;
        outR[i] = fr;
    }
    srcDc[source][0] = dcL;
    srcDc[source][1] = dcR;
    lvlSumSq[source][0] += sumL;
    lvlSumSq[source][1] += sumR;
    lvlCount[source]    += I2S_HALF_FRAMES;
}

/* One finished half of micRxBuffer (PCM1808, stereo slots: slot 0 = left, slot 1 = right):
 * correction -> mono reference x -> adaptive filters -> loudspeakers. */
__attribute__((optimize("O2")))
static void Audio_SaiBlock(const int16_t *src)
{
    const uint32_t t0  = DWT->CYCCNT;
    int16_t       *dst = Audio_TxHalfToFill();

    /* raw data, both slots, before DC blocker / gain */
    for (uint32_t ch = 0; ch < 2U; ch++)
    {
        int16_t mn = INT16_MAX, mx = INT16_MIN;
        for (uint32_t i = 0; i < I2S_HALF_FRAMES; i++)
        {
            int16_t v = src[2U * i + ch];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        rawMin[ch] = mn;
        rawMax[ch] = mx;
    }
    saiSr |= SAI1_Block_A->SR;
    rxCount++;

    float pl[I2S_HALF_FRAMES], pr[I2S_HALF_FRAMES], x[I2S_HALF_FRAMES];
    float out[MIMO_J][MIMO_BLOCK];
    Eq_Block(MIC_SRC_I2S, &src[0], &src[1], 2U, pl, pr);

    const float g = pcmGain * (1.0f / 32768.0f);
    for (uint32_t i = 0; i < I2S_HALF_FRAMES; i++)
    {
        pl[i] *= g;
        pr[i] *= g;
#if (MIMO_REF_MODE == 1)
        x[i] = pl[i];
#elif (MIMO_REF_MODE == 2)
        x[i] = pr[i];
#else
        x[i] = 0.5f * (pl[i] + pr[i]);
#endif
    }

    /* sweep measurement: the reference is replaced by the sweep, which goes to ONE loudspeaker only */
    const uint32_t blk = MimoFx_RefBlocks();
    float sw[I2S_HALF_FRAMES];
    memset(sw, 0, sizeof(sw));
    if (meas.armed && meas.playing && blk >= meas.startBlock)
    {
        for (uint32_t i = 0; i < I2S_HALF_FRAMES; i++) sw[i] = Farina_GenNext(&meas.gen);
        if (meas.gen.n >= MIMO_SWEEP_N) meas.playing = 0U;
    }
    const uint8_t measuring = ms.active;

    MimoFx_ProcessBlock(measuring ? sw : x, out);                    /* the sweep (one loudspeaker) during a measurement, otherwise the reference -> w1 and w2 */
    const uint8_t useEq = (uint8_t)(MimoFx_HaveIr() && !bypass);  /* until a calibration has finished the plain stereo signal is played */

    int16_t peak = 0;
    const float og = outGain * 32768.0f;
    for (uint32_t i = 0; i < I2S_HALF_FRAMES; i++)
    {
        float vl, vr;
        if (measuring)                                              /* only loudspeaker j plays the sweep */
        {
            vl = (meas.j == 0U) ? sw[i] : 0.0f;
            vr = (meas.j == 1U) ? sw[i] : 0.0f;
        }
        else if (!useEq)
        {
            vl = pl[i];
            vr = pr[i];
        }
        else
        {
            vl = out[0][i];
            vr = out[1][i];
        }
        const int16_t yl = (int16_t)Audio_Clip(vl * og);
        const int16_t yr = (int16_t)Audio_Clip(vr * og);
        dst[2U * i]      = yl;         /* loudspeaker 1 */
        dst[2U * i + 1U] = yr;         /* loudspeaker 2 */
        const int16_t al = (int16_t)((yl < 0) ? -yl : yl);
        const int16_t ar = (int16_t)((yr < 0) ? -yr : yr);
        if (al > peak) peak = al;
        if (ar > peak) peak = ar;
        if (al > outPeak[0]) outPeak[0] = al;
        if (ar > outPeak[1]) outPeak[1] = ar;
        if ((vl * og > 32767.0f) || (vl * og < -32768.0f) || (vr * og > 32767.0f) || (vr * og < -32768.0f)) clipCount++;
    }
    micPeak = peak;

    const uint32_t dt = DWT->CYCCNT - t0;
    if (dt > procCyclesMax) procCyclesMax = dt;
    /* load protection: more than 75 % of the block time for this call alone would starve the main loop
     * (adaptation, serial port, button) and make the audio stutter -> bypass and report it */
    if ((dt > (SystemCoreClock / MIMO_FS_HZ) * I2S_HALF_FRAMES * 3U / 4U) && !bypass && !measuring)
    {
        bypass     = 1U;
        shedCycles = dt;
        shed       = 1U;
    }
}

/**
 * @brief  One finished half of the ADC buffers, both microphones: DC blocker, magnitude correction (biquads), phase
 *         correction (FIR), then they go to the adaptation (model microphone 1 / 2 as set by MIMO_MIC1_ADC / MIMO_MIC2_ADC)
 *         and, while a sweep is recorded, into the averaging buffer of the measurement.
 */
__attribute__((optimize("O2")))
static void Audio_MicBlock(const int16_t *adc0, const int16_t *adc1)
{
    const uint32_t t0 = DWT->CYCCNT;
    float z0[I2S_HALF_FRAMES], z1[I2S_HALF_FRAMES];

    Eq_Block(MIC_SRC_ADC, adc0, adc1, 1U, z0, z1);
    for (uint32_t i = 0; i < I2S_HALF_FRAMES; i++)
    {
        z0[i] *= (MIMO_MIC_POLARITY / 32768.0f);       /* MIMO_MIC_POLARITY = -1 inverts the microphones */
        z1[i] *= (MIMO_MIC_POLARITY / 32768.0f);
    }
    MicPhase_Process(0U, z0, I2S_HALF_FRAMES);
    MicPhase_Process(1U, z1, I2S_HALF_FRAMES);

    float *zc[2] = { z0, z1 };
    float *m1 = zc[MIMO_MIC1_ADC];
    float *m2 = zc[MIMO_MIC2_ADC];

    const uint32_t blk = MimoFx_MicBlocks(0U);
    if (meas.armed && !meas.recDone && blk >= meas.startBlock)
    {
        Farina_Accumulate(meas.recPos, (meas.m == 0U) ? m1 : m2, I2S_HALF_FRAMES);
        meas.recPos += I2S_HALF_FRAMES;
        if (meas.recPos >= MIMO_SWEEP_N) meas.recDone = 1U;
    }
    MimoFx_PushMic(0U, m1);
    MimoFx_PushMic(1U, m2);

    const uint32_t dt = DWT->CYCCNT - t0;
    if (dt > procMicCyclesMax) procMicCyclesMax = dt;
}

/* The two ADCs deliver their halves in separate callbacks: the block is processed as soon as both are in */
static void Audio_AdcBlock(uint32_t ch, uint32_t half, const uint16_t *src)
{
    int16_t  *pcm = adcPcm[half][ch];
    uint16_t  mn = UINT16_MAX, mx = 0;

    for (uint32_t i = 0; i < ADC_HALF_SAMPLES; i++)
    {
        uint16_t v = src[i];
        if (v < mn) mn = v;
        if (v > mx) mx = v;
        pcm[i] = (int16_t)((int32_t)v - 32768);   /* 0..65535 -> signed, mic bias ~ VDD/2 ~ 0 */
    }
    adcMin[ch] = mn;
    adcMax[ch] = mx;

    adcReady[half] |= (uint8_t)(1U << ch);
    if (adcReady[half] != 3U) return;             /* wait for the other mic */
    adcReady[half] = 0U;
    adcCount++;

    Audio_MicBlock(adcPcm[half][0], adcPcm[half][1]);
}

/* Fired by the HAL SAI IRQ handler when the first half of micRxBuffer is full */
void HAL_SAI_RxHalfCpltCallback(SAI_HandleTypeDef *hsai)
{
    if (hsai->Instance == hsai_BlockA1.Instance)
    {
        Audio_SaiBlock(&micRxBuffer[0]);
    }
}

/* Fired when the second half of micRxBuffer is full */
void HAL_SAI_RxCpltCallback(SAI_HandleTypeDef *hsai)
{
    if (hsai->Instance == hsai_BlockA1.Instance)
    {
        Audio_SaiBlock(&micRxBuffer[I2S_HALF_SAMPLES]);

        /* Toggle LD1 to show audio is flowing */
        HAL_GPIO_TogglePin(GPIOB, GPIO_PIN_0);
    }
}

/* Fired by the BDMA / DMA1 IRQ handler when the first half of an ADC buffer is full */
void HAL_ADC_ConvHalfCpltCallback(ADC_HandleTypeDef *hadc)
{
    if (hadc->Instance == ADC3)
    {
        Audio_AdcBlock(0U, 0U, &adcRxBuffer[0]);
    }
    else if (hadc->Instance == ADC2)
    {
        Audio_AdcBlock(1U, 0U, &adc2RxBuffer[0]);
    }
}

/* Fired when the second half of an ADC buffer is full */
void HAL_ADC_ConvCpltCallback(ADC_HandleTypeDef *hadc)
{
    if (hadc->Instance == ADC3)
    {
        Audio_AdcBlock(0U, 1U, &adcRxBuffer[ADC_HALF_SAMPLES]);
    }
    else if (hadc->Instance == ADC2)
    {
        Audio_AdcBlock(1U, 1U, &adc2RxBuffer[ADC_HALF_SAMPLES]);
    }
}

void HAL_ADC_ErrorCallback(ADC_HandleTypeDef *hadc)
{
    errAdc = 1U;                        /* never print (blocking) from an interrupt */
}

/**
 * @brief  Busy-wait until the I2S TX DMA is in the middle of a half-buffer.
 *         A source started right after this delivers its blocks half a block
 *         away from the TX half switch, and stays there (common PLL1 clock).
 */
static void Audio_WaitTxMidHalf(void)
{
    uint32_t pos;
    do
    {
        pos = __HAL_DMA_GET_COUNTER(hi2sOut.hdmatx) % I2S_HALF_SAMPLES;
    } while ((pos < I2S_HALF_SAMPLES / 2U - 16U) || (pos > I2S_HALF_SAMPLES / 2U + 16U));
}

/**
 * @brief  Reconfigure an ADC (CubeMX: software start, one-shot DMA) for audio:
 *         TIM6-triggered at fs, circular DMA, 16 bit, 8x oversampling.
 *         Rev.V silicon halves the ADC kernel clock internally: 80 MHz / 4 / 2 = 10 MHz.
 *         One sample: 8 x (16.5 + 8.5) cycles / 10 MHz = 20 us < 22.7 us (1/fs), little margin.
 *         If this ever exceeds 1/fs, triggers get skipped -> half rate -> "robot voice".
 *         Used for ADC3 (PF6, channel 8) and ADC2 (PA3, channel 15).
 */
static void Audio_ADC_Config(ADC_HandleTypeDef *hadc, uint32_t channel)
{
    ADC_ChannelConfTypeDef ch = {0};

    hadc->Init.ClockPrescaler           = ADC_CLOCK_ASYNC_DIV4;           /* PLL2P 80 MHz / 4 (/2 rev.V) */
    hadc->Init.Resolution               = ADC_RESOLUTION_16B;
    hadc->Init.ContinuousConvMode       = DISABLE;
    hadc->Init.ExternalTrigConv         = ADC_EXTERNALTRIG_T6_TRGO;
    hadc->Init.ExternalTrigConvEdge     = ADC_EXTERNALTRIGCONVEDGE_RISING;
    hadc->Init.ConversionDataManagement = ADC_CONVERSIONDATA_DMA_CIRCULAR;
    hadc->Init.Overrun                  = ADC_OVR_DATA_OVERWRITTEN;
    hadc->Init.OversamplingMode         = ENABLE;
    hadc->Init.Oversampling.Ratio                 = 8;
    hadc->Init.Oversampling.RightBitShift         = ADC_RIGHTBITSHIFT_3;
    hadc->Init.Oversampling.TriggeredMode         = ADC_TRIGGEREDMODE_SINGLE_TRIGGER;
    hadc->Init.Oversampling.OversamplingStopReset = ADC_REGOVERSAMPLING_CONTINUED_MODE;
    if (HAL_ADC_Init(hadc) != HAL_OK) Error_Handler();

    ch.Channel      = channel;
    ch.Rank         = ADC_REGULAR_RANK_1;
    ch.SamplingTime = ADC_SAMPLETIME_16CYCLES_5;
    ch.SingleDiff   = ADC_SINGLE_ENDED;
    ch.OffsetNumber = ADC_OFFSET_NONE;
    if (HAL_ADC_ConfigChannel(hadc, &ch) != HAL_OK) Error_Handler();

    if (HAL_ADCEx_Calibration_Start(hadc, ADC_CALIB_OFFSET_LINEARITY, ADC_SINGLE_ENDED) != HAL_OK)
        Error_Handler();
}

/**
 * @brief  Configure SAI1 Block A for the PCM1808 (slave, I2S, MD0/MD1/FMT = GND):
 *         master RX, 16 bit data in 32 bit slots (top 16 of the 24 PCM1808 bits) and the
 *         256 x fs master clock on PE2 = PCM1808 SCKI. Set here as well as in CubeMX so a
 *         regenerated sai.c cannot silently turn the block back into a transmitter
 *         (CubeMX resets the audio mode when the block is re-created).
 */
static void Audio_SAI_Config(void)
{
    hsai_BlockA1.Init.AudioMode = SAI_MODEMASTER_RX;
    hsai_BlockA1.Init.MckOutput = SAI_MCK_OUTPUT_ENABLE;
    if (HAL_SAI_InitProtocol(&hsai_BlockA1, SAI_I2S_STANDARD, SAI_PROTOCOL_DATASIZE_16BITEXTENDED, 2) != HAL_OK)
    {
        Error_Handler();
    }

    /* The generated pin setup uses the LOW speed class, which is specified for ~12 MHz at best:
     * the 11.29 MHz MCLK (PE2) would come out as a rounded, weak clock. Output pins to the PCM1808
     * (MCLK PE2, LRCK PE4, BCK PE5) get the fastest edges; PE6 (DOUT) is an input. */
    GPIO_InitTypeDef io = {0};
    io.Pin       = GPIO_PIN_2 | GPIO_PIN_4 | GPIO_PIN_5;
    io.Mode      = GPIO_MODE_AF_PP;
    io.Pull      = GPIO_NOPULL;
    io.Speed     = GPIO_SPEED_FREQ_VERY_HIGH;
    io.Alternate = GPIO_AF6_SAI1;
    HAL_GPIO_Init(GPIOE, &io);
}

/**
 * @brief  TIM6 as ADC trigger at exactly the I2S/SAI sample rate.
 *         fs = PLL1Q / 256 and the timer clock also comes from PLL1, so
 *         ARR = timclk * 256 / PLL1Q (= 1600 here) gives zero drift.
 */
static void Audio_TIM6_Config(void)
{
    PLL1_ClocksTypeDef pll;
    HAL_RCCEx_GetPLL1ClockFreq(&pll);
    uint32_t timClk = 2U * HAL_RCC_GetPCLK1Freq();                        /* APB1 prescaler 4 -> x2 */
    uint32_t div    = (uint32_t)(((uint64_t)timClk * 256U + pll.PLL1_Q_Frequency / 2U)
                                 / pll.PLL1_Q_Frequency);

    __HAL_RCC_TIM6_CLK_ENABLE();
    TIM6->CR1 = 0;
    TIM6->PSC = 0;
    TIM6->ARR = div - 1U;
    TIM6->EGR = TIM_EGR_UG;             /* load PSC/ARR */
    TIM6->CR2 = TIM_CR2_MMS_1;          /* TRGO on update event */
}

void HAL_I2S_TxCpltCallback(I2S_HandleTypeDef *hi2s)
{
    if (hi2s->Instance == hi2sOut.Instance) { txCount++; }
}

void HAL_SAI_ErrorCallback(SAI_HandleTypeDef *hsai)
{
    errSai = 1U;
}

void HAL_I2S_ErrorCallback(I2S_HandleTypeDef *hi2s)
{
    errI2s = 1U;
}

static float Audio_Db(float lin)
{
    return (lin > 1e-9f) ? 20.0f * log10f(lin) : -200.0f;
}

/* ---------------------------------------------------------------------------------------------------------------
 * Control
 * ------------------------------------------------------------------------------------------------------------- */
/* LEDs without a PC.  LD2 (yellow): off = not calibrated, blinking 2 Hz = countdown before the sweeps (step away, quiet room), very fast = sweeps running,
 * lit = adaptation running, off = paused.  LD3 (red): lit = equalizer off (bypass, key 'q'), blinking = the last measurement was rejected, short flash = B1 held long enough. */
static void Cal_Leds(void)
{
    const uint32_t t = HAL_GetTick();
    uint8_t ld2, ld3;
    if (calPending)                                   ld2 = (uint8_t)((t / 250U) & 1U);
    else if (ms.active)                               ld2 = (uint8_t)((t / 50U) & 1U);
    else                                              ld2 = (MimoFx_HaveIr() && MimoFx_GetAdapt()) ? 1U : 0U;
    if ((int32_t)(failBlinkUntil - t) > 0)       ld3 = (uint8_t)((t / 100U) & 1U);
    else if ((int32_t)(ackUntil - t) > 0)             ld3 = 1U;
    else                                              ld3 = bypass;
    HAL_GPIO_WritePin(GPIOE, GPIO_PIN_1, ld2 ? GPIO_PIN_SET : GPIO_PIN_RESET);
    HAL_GPIO_WritePin(GPIOB, GPIO_PIN_14, ld3 ? GPIO_PIN_SET : GPIO_PIN_RESET);
}

static void Mimo_ShowAdapt(void)
{
    Cal_Leds();
}

static void Mimo_ToggleAdapt(void)
{
    if (ms.active) { UART_Print("Measurement running, 'x' aborts it\r\n"); return; }
    if (!MimoFx_HaveIr()) { UART_Print("No impulse response yet: measure first (key 's')\r\n"); return; }
    MimoFx_SetAdapt(!MimoFx_GetAdapt());
    Mimo_ShowAdapt();
    UART_Print(MimoFx_GetAdapt() ? "Adaptation: RUNNING (LD2 lit)\r\n" : "Adaptation: FROZEN (filters kept)\r\n");
}

/* equalizer on / off for an A/B comparison (the adaptation keeps running) */
static void Mimo_ToggleBypass(void)
{
    bypass = (uint8_t)!bypass;
    Cal_Leds();
    UART_Print(bypass ? "Equalizer: OFF - corrected PCM1808 stereo signal (LD3 lit)\r\n" : "Equalizer: ON\r\n");
}

/* ---------------------------------------------------------------------------------------------------------------
 * Sweep measurement of h^ (main loop part, the interrupts do the playing / recording)
 * ------------------------------------------------------------------------------------------------------------- */
static void Meas_Start(void)
{
    if (ms.active) { UART_Print("Measurement already running\r\n"); return; }
    MimoFx_SetAdapt(0);
    Mimo_ShowAdapt();
    memset(irWin, 0, sizeof(irWin));
    Farina_Clear();
    meas.armed = 0U;
    ms.combo = 0U;  ms.rep = 0U;  ms.step = 0U;
    ms.active = 1U;
    char m[160];
    snprintf(m, sizeof(m), "Sweep measurement: %.0f Hz .. %.0f Hz, T = %.2f s, record %lu samples, %u averages, 4 combinations "
             "(about %.0f s). Keep the room quiet.\r\n", (double)MIMO_SWEEP_F1_HZ, (double)MIMO_SWEEP_F2_HZ,
             (double)MIMO_SWEEP_T_S, (unsigned long)MIMO_SWEEP_N, (unsigned)MIMO_SWEEP_AVG,
             (double)(4U * MIMO_SWEEP_AVG) * ((double)MIMO_SWEEP_N / MIMO_FS_HZ + 0.1));
    UART_Print(m);
}

static void Meas_Abort(const char *why)
{
    meas.armed = 0U;
    meas.playing = 0U;
    ms.active = 0U;
    UART_Print((char *)why);
}

static void Meas_Finish(void)
{
    char m[340];
    float minSnr = 1.0e9f;
    ms.active = 0U;
    for (uint32_t mi = 0; mi < MIMO_M; mi++)
        for (uint32_t j = 0; j < MIMO_J; j++)
        {
            snprintf(m, sizeof(m), "  speaker %u -> mic %u: strongest peak at lag %ld (%.1f ms), response energy starts at lag %ld, %.1f dB above the noise of the record\r\n", (unsigned)(j + 1U), (unsigned)(mi + 1U),
                     (long)irPeakLag[mi][j], (double)irPeakLag[mi][j] * 1000.0 / (double)MIMO_FS_HZ, (long)irEnergyStart[mi][j], (double)irSnr[mi][j]);
            UART_Print(m);
            if (irSnr[mi][j] < minSnr) minSnr = irSnr[mi][j];
        }
    if (minSnr < MIMO_IR_MIN_SNR_DB)
    {
        snprintf(m, sizeof(m), "Measurement REJECTED: the weakest response is only %.1f dB above the noise (minimum %.0f dB): the sweep does not reach the microphone well enough. "
                 "Raise the amplifier volume / MIMO_SWEEP_LEVEL_DB or move the microphones closer; nothing was changed.\r\n", (double)minSnr, (double)MIMO_IR_MIN_SNR_DB);
        UART_Print(m);
        failBlinkUntil = HAL_GetTick() + 2500U;
        return;
    }
    const float *full[MIMO_M][MIMO_J];
    for (uint32_t mi = 0; mi < MIMO_M; mi++) for (uint32_t j = 0; j < MIMO_J; j++) full[mi][j] = irWin[mi][j];
    static MimoIr ir;
    const int rc = MimoIr_Build(&ir, full);
    if (rc != 0)
    {
        UART_Print("Measurement FAILED: a response is empty (is the amplifier on, the mics connected?)\r\n");
        return;
    }
    for (uint32_t mi = 0; mi < MIMO_M; mi++)
    {
        snprintf(m, sizeof(m), "mic %u: window of h^ starts at lag S = %u, desired delay tau = %u, G = %.5f\r\n", (unsigned)(mi + 1U), (unsigned)ir.shift[mi], (unsigned)ir.tau[mi], (double)ir.gain[mi]);
        UART_Print(m);
    }
    snprintf(m, sizeof(m), "polarity of the strongest sample of the window (MIMO_MIC_POLARITY = %+d): spk1->mic1 %c  spk1->mic2 %c  spk2->mic1 %c  spk2->mic2 %c\r\n",
             (int)MIMO_MIC_POLARITY, ir.polarity[0][0] > 0 ? '+' : '-', ir.polarity[1][0] > 0 ? '+' : '-',
             ir.polarity[0][1] > 0 ? '+' : '-', ir.polarity[1][1] > 0 ? '+' : '-');
    UART_Print(m);
    MimoFx_SetIr(&ir);              /* restarts the adaptation from the start weights */
    MimoFx_SetAdapt(1);
    Mimo_ShowAdapt();
    UART_Print("Impulse responses installed, adaptation RUNNING\r\n");
}

static void Meas_Poll(void)
{
    if (!ms.active) return;
    const uint8_t j = (uint8_t)(ms.combo / MIMO_M);
    const uint8_t mi = (uint8_t)(ms.combo % MIMO_M);
    char m[100];

    switch (ms.step)
    {
    case 0:                                                  /* arm: both interrupts start 4 blocks from now */
        meas.j = j;  meas.m = mi;
        meas.recPos = 0U;
        meas.recDone = 0U;
        Farina_GenInit(&meas.gen, 0U);
        meas.startBlock = MimoFx_RefBlocks() + 4U + MIMO_SWEEP_GAP_BLOCKS;
        meas.playing = 1U;
        meas.armed = 1U;
        ms.t0 = HAL_GetTick();
        ms.step = 1U;
        break;
    case 1:
        if (meas.recDone && !meas.playing)
        {
            meas.armed = 0U;
            ms.rep++;
            snprintf(m, sizeof(m), "  loudspeaker %u -> mic %u: recording %u / %u done\r\n", (unsigned)(j + 1U), (unsigned)(mi + 1U),
                     (unsigned)ms.rep, (unsigned)MIMO_SWEEP_AVG);
            UART_Print(m);
            ms.step = (ms.rep >= MIMO_SWEEP_AVG) ? 2U : 0U;
        }
        else if (HAL_GetTick() - ms.t0 > 5000U + (uint32_t)(1000.0 * (double)MIMO_SWEEP_N / MIMO_FS_HZ))
        {
            Meas_Abort("Measurement ABORTED: timeout (audio blocks not running?)\r\n");
        }
        break;
    case 2:                                                  /* all recordings of this combination are in: x -> FFT, Y/X, IFFT */
    {
        float tm, tr, pk;
        char rm[160];
        const uint32_t ri = (uint32_t)j * MIMO_M + mi;
        Farina_StoreRecord(recStore[ri], &recScale[ri], &tm, &tr, &pk);
        recValid[ri] = 1U;
        snprintf(rm, sizeof(rm), "  record loudspeaker %u -> mic %u: peak %.4f, second half: mean %+.2e, rms %.2e (full scale 1.0)\r\n", (unsigned)(j + 1U), (unsigned)(mi + 1U), (double)pk, (double)tm, (double)tr);
        UART_Print(rm);
    }
        Farina_Deconvolve();
        Farina_Locate(&irPeakLag[mi][j], &irPeakVal[mi][j], &irSnr[mi][j]);          /* information and the signal to noise ratio */
        irEnergyStart[mi][j] = Farina_EnergyStart();                                 /* information: where the response begins */
        for (uint32_t k = 0; k < MIMO_LH; k++) irWin[mi][j][k] = Farina_Ir((int32_t)MIMO_IR_START + (int32_t)k);   /* fixed window */
        {   /* spectrum of the long response, frequency resolution about 2.7 Hz */
            const uint32_t L = 16384U, FADE = 2048U;
            for (uint32_t p = 0U; p < 16U; p++)
            {
                const double w = 2.0 * 3.14159265358979 * (double)fullFreq[p] / (double)MIMO_FS_HZ;
                const double cr = cos(w), ci = -sin(w);
                double zr = 1.0, zi = 0.0, re = 0.0, im = 0.0;
                for (uint32_t n = 0U; n < L; n++)
                {
                    double v = (double)Farina_Ir((int32_t)MIMO_IR_START + (int32_t)n);
                    if (n >= L - FADE) v *= 0.5 * (1.0 + cos(3.14159265358979 * (double)(n - (L - FADE)) / (double)FADE));
                    re += v * zr;  im += v * zi;
                    const double t = zr * cr - zi * ci;  zi = zr * ci + zi * cr;  zr = t;
                }
                fullDb[mi][j][p] = (float)(10.0 * log10(re * re + im * im + 1.0e-30));
            }
        }
        Farina_Clear();
        ms.combo++;
        ms.rep = 0U;
        ms.step = 0U;
        if (ms.combo >= MIMO_J * MIMO_M) Meas_Finish();
        break;
    default:
        break;
    }
}

/* ---------------------------------------------------------------------------------------------------------------
 * Printouts
 * ------------------------------------------------------------------------------------------------------------- */
/* 'h': spectrum of the measured responses h^ (for the four paths) and the numbers of the measurement; 'H': the taps themselves */
static void Mimo_PrintIr(int taps)
{
    if (!MimoFx_HaveIr()) { UART_Print("No impulse response yet\r\n"); return; }
    const MimoIr *ir = MimoFx_GetIr();
    char m[220];
    for (uint32_t mi = 0; mi < MIMO_M; mi++)
    {
        snprintf(m, sizeof(m), "mic %u: window start S = %u, tau = %u, G = %.6f, polarity %c %c\r\n", (unsigned)(mi + 1U), (unsigned)ir->shift[mi], (unsigned)ir->tau[mi], (double)ir->gain[mi],
                 ir->polarity[mi][0] > 0 ? '+' : '-', ir->polarity[mi][1] > 0 ? '+' : '-');
        UART_Print(m);
    }
    if (!taps)
    {
        static const float fr[16] = { 20.f, 31.f, 50.f, 80.f, 125.f, 200.f, 315.f, 500.f, 800.f, 1250.f, 2000.f, 3150.f, 5000.f, 8000.f, 12500.f, 19000.f };
        UART_Print("|h^(f)| in dB at 20 31 50 80 125 200 315 500 800 1250 2000 3150 5000 8000 12500 19000 Hz (normalized: all four paths together have an energy of 2 per mic):\r\n");
        for (uint32_t mi = 0; mi < MIMO_M; mi++)
            for (uint32_t j = 0; j < MIMO_J; j++)
            {
                int n = snprintf(m, sizeof(m), "m%u<-spk%u:", (unsigned)(mi + 1U), (unsigned)(j + 1U));
                for (uint32_t p = 0; p < 16U; p++)
                {
                    const float wr = 2.0f * 3.14159265f * fr[p] / (float)MIMO_FS_HZ;
                    float re = 0.0f, im = 0.0f;
                    for (uint32_t k = 0; k < MIMO_LH; k++) { re += ir->h[mi][j][k] * cosf(wr * (float)k); im -= ir->h[mi][j][k] * sinf(wr * (float)k); }
                    n += snprintf(&m[n], sizeof(m) - (size_t)n, " %+5.1f", (double)(10.0f * log10f(re * re + im * im + 1.0e-12f)));
                }
                UART_Print(m); UART_Print("\r\n");
            }
        /* the same in dB relative to the mean of 800 .. 3150 Hz, for the 256 sample window and for the long response (0.37 s) */
        UART_Print("relative to the mean of 800 .. 3150 Hz, 'win' = the 256 sample window used by the filters, 'full' = the long measured response:\r\n");
        for (uint32_t mi = 0; mi < MIMO_M; mi++)
            for (uint32_t j = 0; j < MIMO_J; j++)
            {
                float win[16];
                for (uint32_t p = 0; p < 16U; p++)
                {
                    const float wr = 2.0f * 3.14159265f * fullFreq[p] / (float)MIMO_FS_HZ;
                    float re = 0.0f, im = 0.0f;
                    for (uint32_t k = 0; k < MIMO_LH; k++) { re += ir->h[mi][j][k] * cosf(wr * (float)k); im -= ir->h[mi][j][k] * sinf(wr * (float)k); }
                    win[p] = 10.0f * log10f(re * re + im * im + 1.0e-12f);
                }
                const float rw = 0.25f * (win[8] + win[9] + win[10] + win[11]);
                const float rf = 0.25f * (fullDb[mi][j][8] + fullDb[mi][j][9] + fullDb[mi][j][10] + fullDb[mi][j][11]);
                int n = snprintf(m, sizeof(m), "m%u<-spk%u win :", (unsigned)(mi + 1U), (unsigned)(j + 1U));
                for (uint32_t p = 0; p < 16U; p++) n += snprintf(&m[n], sizeof(m) - (size_t)n, " %+5.1f", (double)(win[p] - rw));
                UART_Print(m); UART_Print("\r\n");
                n = snprintf(m, sizeof(m), "m%u<-spk%u full:", (unsigned)(mi + 1U), (unsigned)(j + 1U));
                for (uint32_t p = 0; p < 16U; p++) n += snprintf(&m[n], sizeof(m) - (size_t)n, " %+5.1f", (double)(fullDb[mi][j][p] - rf));
                UART_Print(m); UART_Print("\r\n");
            }
        /* what the model expects at each microphone when both loudspeakers play the same signal: the sum of the paths through the filters, dB relative to the mean of
         * 800 .. 3150 Hz. 'plain' = unit filters (no correction), 'corrected' = with the published filters w_j. */
        UART_Print("model prediction at the microphones, both loudspeakers playing the same signal, relative to 800 .. 3150 Hz:\r\n");
        for (uint32_t mi = 0; mi < MIMO_M; mi++)
        {
            float plain[16], corr[16];
            for (uint32_t p = 0; p < 16U; p++)
            {
                const float wr = 2.0f * 3.14159265f * fullFreq[p] / (float)MIMO_FS_HZ;
                float pr = 0.0f, pi = 0.0f, cr = 0.0f, ci = 0.0f;
                for (uint32_t j = 0; j < MIMO_J; j++)
                {
                    float hr = 0.0f, hi = 0.0f;
                    for (uint32_t k = 0; k < MIMO_LH; k++) { hr += ir->h[mi][j][k] * cosf(wr * (float)k); hi -= ir->h[mi][j][k] * sinf(wr * (float)k); }
                    float md, pd;
                    MimoFx_FilterResponse(j, fullFreq[p], &md, &pd);
                    const float wm = powf(10.0f, md * 0.05f), wa = pd * (3.14159265f / 180.0f);
                    const float wre = wm * cosf(wa), wim = wm * sinf(wa);
                    pr += hr;                       pi += hi;                       /* unit filters (the common delay does not matter for the level) */
                    cr += hr * wre - hi * wim;      ci += hr * wim + hi * wre;      /* with the filters */
                }
                plain[p] = 10.0f * log10f(pr * pr + pi * pi + 1.0e-12f);
                corr[p]  = 10.0f * log10f(cr * cr + ci * ci + 1.0e-12f);
            }
            const float r0 = 0.25f * (plain[8] + plain[9] + plain[10] + plain[11]);
            const float r1 = 0.25f * (corr[8] + corr[9] + corr[10] + corr[11]);
            int n = snprintf(m, sizeof(m), "m%u plain    :", (unsigned)(mi + 1U));
            for (uint32_t p = 0; p < 16U; p++) n += snprintf(&m[n], sizeof(m) - (size_t)n, " %+5.1f", (double)(plain[p] - r0));
            UART_Print(m); UART_Print("\r\n");
            n = snprintf(m, sizeof(m), "m%u corrected:", (unsigned)(mi + 1U));
            for (uint32_t p = 0; p < 16U; p++) n += snprintf(&m[n], sizeof(m) - (size_t)n, " %+5.1f", (double)(corr[p] - r1));
            UART_Print(m); UART_Print("\r\n");
        }
        /* conditioning of the 2 x 2 system [h^_mj(f)]: ratio of the largest to the smallest singular value in dB. Where it is large (> 20 dB) the two microphones see the two loudspeakers almost alike. */
        {
            int n = snprintf(m, sizeof(m), "conditioning of the 2x2 system, largest / smallest singular value in dB:");
            UART_Print(m); UART_Print("\r\n");
            n = snprintf(m, sizeof(m), "cond dB    :");
            for (uint32_t p = 0; p < 16U; p++)
            {
                const float wr = 2.0f * 3.14159265f * fullFreq[p] / (float)MIMO_FS_HZ;
                float hr[2][2], hi[2][2];
                for (uint32_t mi = 0; mi < 2U; mi++)
                    for (uint32_t j = 0; j < 2U; j++)
                    {
                        float a = 0.0f, b = 0.0f;
                        for (uint32_t k = 0; k < MIMO_LH; k++) { a += ir->h[mi][j][k] * cosf(wr * (float)k); b -= ir->h[mi][j][k] * sinf(wr * (float)k); }
                        hr[mi][j] = a; hi[mi][j] = b;
                    }
                const float T  = hr[0][0] * hr[0][0] + hi[0][0] * hi[0][0] + hr[0][1] * hr[0][1] + hi[0][1] * hi[0][1]
                               + hr[1][0] * hr[1][0] + hi[1][0] * hi[1][0] + hr[1][1] * hr[1][1] + hi[1][1] * hi[1][1];
                const float dr = hr[0][0] * hr[1][1] - hi[0][0] * hi[1][1] - (hr[0][1] * hr[1][0] - hi[0][1] * hi[1][0]);
                const float di = hr[0][0] * hi[1][1] + hi[0][0] * hr[1][1] - (hr[0][1] * hi[1][0] + hi[0][1] * hr[1][0]);
                const float D2 = dr * dr + di * di;                                   /* |det|^2 = (s1 s2)^2 */
                float disc = T * T - 4.0f * D2;
                if (disc < 0.0f) disc = 0.0f;
                const float s2max = 0.5f * (T + sqrtf(disc)), s2min = 0.5f * (T - sqrtf(disc));
                n += snprintf(&m[n], sizeof(m) - (size_t)n, " %+5.1f", (double)(10.0f * log10f((s2max + 1.0e-12f) / (s2min + 1.0e-12f))));
            }
            UART_Print(m); UART_Print("\r\n");
        }
        return;
    }
    for (uint32_t mi = 0; mi < MIMO_M; mi++)
        for (uint32_t j = 0; j < MIMO_J; j++)
        {
            snprintf(m, sizeof(m), "H m%u spk%u\r\n", (unsigned)(mi + 1U), (unsigned)(j + 1U));
            UART_Print(m);
            for (uint32_t k = 0; k < MIMO_LH; k += 16U)
            {
                int n = 0;
                for (uint32_t t = 0; t < 16U && k + t < MIMO_LH; t++) n += snprintf(&m[n], sizeof(m) - (size_t)n, "%.6f ", (double)ir->h[mi][j][k + t]);
                UART_Print(m); UART_Print("\r\n");
            }
        }
    UART_Print("END H\r\n");
}

/* w1 and w2: level of each and the phase of w2 relative to w1. Around 180 deg = one loudspeaker is inverted against the other at that frequency (the phantom centre disappears). */
static void Mimo_PrintFilterRelation(void)
{
    static const float fr[10] = { 50.f, 100.f, 200.f, 500.f, 1000.f, 2000.f, 4000.f, 8000.f, 12000.f, 16000.f };
    char m[200];
    UART_Print("  Hz        : "); for (uint32_t p = 0; p < 10U; p++) { snprintf(m, sizeof(m), "%7.0f", (double)fr[p]); UART_Print(m); } UART_Print("\r\n");
    float mag[2][10], ph[2][10];
    for (uint32_t j = 0; j < 2U; j++) for (uint32_t p = 0; p < 10U; p++) MimoFx_FilterResponse(j, fr[p], &mag[j][p], &ph[j][p]);
    for (uint32_t j = 0; j < 2U; j++) { snprintf(m, sizeof(m), "  |w%u| dB   : ", (unsigned)(j + 1U)); UART_Print(m); for (uint32_t p = 0; p < 10U; p++) { snprintf(m, sizeof(m), "%+7.1f", (double)mag[j][p]); UART_Print(m); } UART_Print("\r\n"); }
    UART_Print("  w2 vs w1 deg: ");
    for (uint32_t p = 0; p < 10U; p++)
    {
        float d = ph[1][p] - ph[0][p];
        while (d > 180.0f) d -= 360.0f;
        while (d < -180.0f) d += 360.0f;
        snprintf(m, sizeof(m), "%+7.0f", (double)d);
        UART_Print(m);
    }
    UART_Print("\r\n");
}

static void Mimo_PrintState(void)
{
    char m[200];
    snprintf(m, sizeof(m), "MIMO (unbiased): ir=%s adapt=%s eq=%s | blocks=%lu backlog=%lu skipped=%lu\r\n",
             MimoFx_HaveIr() ? "yes" : "NO", MimoFx_GetAdapt() ? "run" : "frozen", bypass ? "OFF" : "on", (unsigned long)MimoFx_BlocksDone(),
             (unsigned long)MimoFx_Backlog(), (unsigned long)MimoFx_Skipped());
    UART_Print(m);
    for (uint32_t j = 0; j < MIMO_J; j++)
    {
        uint32_t pk = 0; float pv = 0.0f;
        for (uint32_t l = 0; l < MIMO_LW; l++)
        {
            const float v = fabsf(MimoFx_Weight(j, l));
            if (v > pv) { pv = v; pk = l; }
        }
        snprintf(m, sizeof(m), "  w%u: rms %.3f, peak %+.3f at tap %lu\r\n", (unsigned)(j + 1U),
                 (double)MimoFx_WeightRms(j), (double)MimoFx_Weight(j, pk), (unsigned long)pk);
        UART_Print(m);
    }
    snprintf(m, sizeof(m), "  mic error e %.1f / %.1f dB (re d, smoothed)\r\n",
             (double)MimoFx_ErrDb(0), (double)MimoFx_ErrDb(1));
    UART_Print(m);
    Mimo_PrintFilterRelation();
    snprintf(m, sizeof(m), "  mic phase FIR taps: %lu / %lu, CPU %lu MHz\r\n", (unsigned long)MicPhase_Taps(0U),
             (unsigned long)MicPhase_Taps(1U), (unsigned long)(SystemCoreClock / 1000000U));
    UART_Print(m);
}

/**
 * @brief  Control from the ST-LINK serial terminal (USART3, 115200), non-blocking, one key per call.
 *           's'        measure the impulse responses with the Farina sweep (4 combinations), then start the adaptation
 *           'x'        abort a running measurement
 *           'a'        adaptation run / freeze (same as a short B1 press, LD2 lit = running)
 *           'q'        equalizer on / off (LD3 red lit = off)
 *           'u' 'j'    IPNLMS step size mu up / down (x1.5)      'i' 'k'  delta up / down (x3)      'o' 'l'  epsilon up / down (x3)   (upper key = increase, key below = decrease)
 *           'y'        dump an averaged sweep recording (next combination per press, 20 s each)
 *           'z'        weights back to the start value (the measured responses are kept)
 *           'w'        print the state          'h' / 'H'  print the measured responses h^ (spectrum / the taps)
 *           '+' / '-'  PCM1808 gain +1 dB / -1 dB      ']' / '['  +0.1 dB / -0.1 dB
 *           'g'        print the gain                   'r'        reset the gain to MIC_GAIN_DB
 *         The gain is kept until reset/power-off; copy it into MIC_GAIN_DB to make it permanent.
 */
static void Gain_HandleUart(void)
{
    if (__HAL_UART_GET_FLAG(&huart3, UART_FLAG_ORE)) __HAL_UART_CLEAR_OREFLAG(&huart3);
    if (!__HAL_UART_GET_FLAG(&huart3, UART_FLAG_RXNE)) return;

    char  c    = (char)(huart3.Instance->RDR & 0xFFU);
    float step = 0.0f;

    switch (c)
    {
        case '+': step =  1.0f; break;
        case '-': step = -1.0f; break;
        case ']': step =  0.1f; break;
        case '[': step = -0.1f; break;
        case 'r': pcmGain = powf(10.0f, MIC_GAIN_DB / 20.0f); break;
        case 'g': break;
        case 's': Meas_Start(); return;
        case 'u': case 'j': case 'i': case 'k': case 'o': case 'l':   /* gain parameters of the IPNLMS, live */
        {
            float mu, de, ep;
            MimoFx_GetGains(&mu, &de, &ep);
            if (c == 'u') { mu *= 1.5f; }
            if (c == 'j') { mu /= 1.5f; }
            if (c == 'o') { ep *= 3.0f; }
            if (c == 'l') { ep /= 3.0f; }
            if (c == 'i') { de *= 3.0f; }
            if (c == 'k') { de /= 3.0f; }
            MimoFx_SetGains(mu, de, ep);
            char m[260];
            snprintf(m, sizeof(m), "-----------------------------------------------------------------------------\r\ngains: mu %.4f  delta %.2e  epsilon %.2e\r\n(up / down:  u/j mu,  i/k delta,  o/l epsilon)\r\n-----------------------------------------------------------------------------\r\n", (double)mu, (double)de, (double)ep);
            UART_Print(m);
            return;
        }
        case 'x': if (ms.active) Meas_Abort("Measurement aborted\r\n"); return;
        case 'a': Mimo_ToggleAdapt(); return;
        case 'q': Mimo_ToggleBypass(); return;
        case 'z':
            MimoFx_Reset();
            UART_Print("Weights reset to the start value\r\n");
            return;
        case 'w': Mimo_PrintState(); return;
        case 'y':   /* dump of an averaged sweep recording (next combination at every press) */
        {
            uint8_t ri = recNext;
            for (uint8_t t = 0U; t < MIMO_J * MIMO_M && !recValid[ri]; t++) ri = (uint8_t)((ri + 1U) % (MIMO_J * MIMO_M));
            if (!recValid[ri]) { UART_Print("No recording stored: calibrate first\r\n"); return; }
            recNext = (uint8_t)((ri + 1U) % (MIMO_J * MIMO_M));
            char mm[200];
            snprintf(mm, sizeof(mm), "RECORD loudspeaker %u -> mic %u, N %lu, fs %u Hz, y = value * %.6e (full scale 1.0), sweep %.0f .. %.0f Hz, T %.2f s, level %.1f dB, %u averages\r\n",
                     (unsigned)(ri / MIMO_M + 1U), (unsigned)(ri % MIMO_M + 1U), (unsigned long)MIMO_SWEEP_N, (unsigned)MIMO_FS_HZ, (double)recScale[ri],
                     (double)MIMO_SWEEP_F1_HZ, (double)MIMO_SWEEP_F2_HZ, (double)MIMO_SWEEP_T_S, (double)MIMO_SWEEP_LEVEL_DB, (unsigned)MIMO_SWEEP_AVG);
            UART_Print(mm);
            for (uint32_t n0 = 0U; n0 < MIMO_SWEEP_N; n0 += 64U)
            {
                char line[64 * 8 + 4];
                uint32_t len = 0U;
                for (uint32_t n = n0; n < n0 + 64U && n < MIMO_SWEEP_N; n++)
                    len += (uint32_t)snprintf(&line[len], sizeof(line) - len, "%d\r\n", (int)recStore[ri][n]);
                HAL_UART_Transmit(&huart3, (uint8_t *)line, len, HAL_MAX_DELAY);
            }
            UART_Print("END RECORD\r\n");
            return;
        }
        case 'h': Mimo_PrintIr(0); return;
        case 'H': Mimo_PrintIr(1); return;
        default:  return;
    }

    if (step != 0.0f)
    {
        float db = Audio_Db(pcmGain) + step;
        if (db < -40.0f) db = -40.0f;
        if (db >  60.0f) db =  60.0f;
        pcmGain = powf(10.0f, db / 20.0f);   /* single float store: safe vs. the ISRs */
    }

    char m[60];
    snprintf(m, sizeof(m), "Gain PCM1808 = %.1f dB\r\n", (double)Audio_Db(pcmGain));
    UART_Print(m);
}

/**
 * @brief  Copy the correction coefficients from flash into DTCM and clear the filter states.
 */
static void Eq_Init(void)
{
    memset(eqState, 0, sizeof(eqState));
    memset(eqStateF, 0, sizeof(eqStateF));
    for (uint32_t k = 0; k < MIC_SRC_COUNT; k++)
    {
        for (uint32_t c = 0; c < 2U; c++)
        {
            eqNd[k][c] = 0U;
            eqNf[k][c] = 0U;
            for (uint32_t n = 0; n < eqTable[k][c].n; n++)
            {
                const MicEqBiquad *b = &eqTable[k][c].bq[n];
                /* pole radius of z^2 + a1 z + a2 */
                const double disc = b->a1 * b->a1 - 4.0 * b->a2;
                const double r = (disc < 0.0) ? sqrt(b->a2) : (fabs(b->a1) + sqrt(disc)) / 2.0;
                if (r > EQ_DOUBLE_RADIUS)
                {
                    eqCoef[k][c][eqNd[k][c]++] = *b;
                }
                else
                {
                    float *f = eqCoefF[k][c][eqNf[k][c]++];
                    f[0] = (float)b->b0;  f[1] = (float)b->b1;  f[2] = (float)b->b2;  f[3] = (float)b->a1;  f[4] = (float)b->a2;
                }
            }
        }
    }
}

/**
 * @brief  Blocking UART print helper over USART3
 */
void UART_Print(char *str)
{
    HAL_UART_Transmit(&huart3, (uint8_t *)str, strlen(str), HAL_MAX_DELAY);
}
/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */
  SCB_EnableICache();     /* code runs from flash: instruction cache is safe (the D-cache stays off, DMA buffers) */
  /* USER CODE END 1 */
/* USER CODE BEGIN Boot_Mode_Sequence_0 */
#if defined(DUAL_CORE_BOOT_SYNC_SEQUENCE)
  int32_t timeout;
#endif /* DUAL_CORE_BOOT_SYNC_SEQUENCE */
/* USER CODE END Boot_Mode_Sequence_0 */

/* USER CODE BEGIN Boot_Mode_Sequence_1 */
#if defined(DUAL_CORE_BOOT_SYNC_SEQUENCE)
  /* Wait until CPU2 boots and enters in stop mode or timeout*/
  timeout = 0xFFFF;
  while((__HAL_RCC_GET_FLAG(RCC_FLAG_D2CKRDY) != RESET) && (timeout-- > 0));
  if ( timeout < 0 )
  {
  Error_Handler();
  }
#endif /* DUAL_CORE_BOOT_SYNC_SEQUENCE */
/* USER CODE END Boot_Mode_Sequence_1 */
  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* Configure the peripherals common clocks */
  PeriphCommonClock_Config();
/* USER CODE BEGIN Boot_Mode_Sequence_2 */
#if defined(DUAL_CORE_BOOT_SYNC_SEQUENCE)
/* When system initialization is finished, Cortex-M7 will release Cortex-M4 by means of
HSEM notification */
/*HW semaphore Clock enable*/
__HAL_RCC_HSEM_CLK_ENABLE();
/*Take HSEM */
HAL_HSEM_FastTake(HSEM_ID_0);
/*Release HSEM in order to notify the CPU2(CM4)*/
HAL_HSEM_Release(HSEM_ID_0,0);
/* wait until CPU2 wakes up from stop mode */
timeout = 0xFFFF;
while((__HAL_RCC_GET_FLAG(RCC_FLAG_D2CKRDY) == RESET) && (timeout-- > 0));
if ( timeout < 0 )
{
Error_Handler();
}
#endif /* DUAL_CORE_BOOT_SYNC_SEQUENCE */
/* USER CODE END Boot_Mode_Sequence_2 */

  /* USER CODE BEGIN SysInit */
  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_BDMA_Init();
  MX_SAI1_Init();
  MX_USART3_UART_Init();
  MX_ADC3_Init();
  MX_I2S1_Init();
  MX_ADC2_Init();
  /* USER CODE BEGIN 2 */

    /* D2 SRAM (SRAM1..3, 288 KB): the clocks are off after reset; the sweep recordings live there */
    SET_BIT(RCC->AHB2ENR, RCC_AHB2ENR_D2SRAM1EN | RCC_AHB2ENR_D2SRAM2EN | RCC_AHB2ENR_D2SRAM3EN);
    (void)READ_BIT(RCC->AHB2ENR, RCC_AHB2ENR_D2SRAM1EN);

    /* Initialize Debug LED on PB0 */
    __HAL_RCC_GPIOB_CLK_ENABLE();
    GPIO_InitTypeDef GPIO_InitStruct = {0};
    GPIO_InitStruct.Pin = GPIO_PIN_0;
    GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
    GPIO_InitStruct.Pull = GPIO_NOPULL;
    GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
    HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);
    HAL_GPIO_WritePin(GPIOB, GPIO_PIN_0, GPIO_PIN_RESET); // Turn OFF LD1 initially

    /* start gains (dB -> linear) */
    pcmGain = powf(10.0f, MIC_GAIN_DB / 20.0f);
    outGain = powf(10.0f, OUTPUT_GAIN_DB / 20.0f);

    /* DWT cycle counter: measures the processing time per block */
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->LAR = 0xC5ACCE55;                      /* unlock (Cortex-M7) */
    DWT->CYCCNT = 0;
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;

    Eq_Init();
    MicPhase_Init(micPh_adcL, MIC_PH_ADCL_N, micPh_adcR, MIC_PH_ADCR_N);   /* FIR arrays are in mic_eq_coeffs.h */
    Farina_Init();
    MimoFx_Init();                      /* .dtcm is NOLOAD: everything in it is cleared here (a NaN in an IIR state would never go away) */
    MimoFx_SetAdapt(0);                 /* nothing to adapt before the first sweep measurement */

    Audio_SAI_Config();                 /* PCM1808: SAI1 A master RX + MCLK on PE2 */

    /* 1. Start I2S TX DMA in CIRCULAR mode */
    /* This will continuously play whatever is inside i2sTxBuffer */

    memset(i2sTxBuffer, 0, sizeof(i2sTxBuffer));   /* silence until the first block */

    if (HAL_I2S_Transmit_DMA(&hi2sOut, (uint16_t *)i2sTxBuffer, I2S_BUFFER_SAMPLES) != HAL_OK)
    {
        UART_Print("I2S DMA TX failed!\r\n");
        Error_Handler();
    }

    /* 2. Start SAI RX DMA in CIRCULAR mode: PCM1808 (slave, I2S 24 bit, MD0/MD1/FMT = GND).
     * Block A is master: it supplies LRCK (PE4), BCK (PE5) and the 256 x fs system clock SCKI (PE2),
     * and fills micRxBuffer with the L/R samples from DOUT (PE6). */
    Audio_WaitTxMidHalf();
    if (HAL_SAI_Receive_DMA(&hsai_BlockA1, (uint8_t *)micRxBuffer, I2S_BUFFER_SAMPLES) != HAL_OK)
    {
        UART_Print("SAI DMA RX failed!\r\n");
        Error_Handler();
    }

    /* 3. Analog mics: ADC3 (PF6) and ADC2 (PA3) in circular DMA,
     *    conversions of both triggered by TIM6 at fs */
    Audio_ADC_Config(&hadc3, ADC_CHANNEL_8);
    Audio_ADC_Config(&hadc2, ADC_CHANNEL_15);
    Audio_TIM6_Config();
    if (HAL_ADC_Start_DMA(&hadc3, (uint32_t *)adcRxBuffer, ADC_HALF_SAMPLES * 2U) != HAL_OK)
    {
        UART_Print("ADC3 DMA failed!\r\n");
        Error_Handler();
    }
    if (HAL_ADC_Start_DMA(&hadc2, (uint32_t *)adc2RxBuffer, ADC_HALF_SAMPLES * 2U) != HAL_OK)
    {
        UART_Print("ADC2 DMA failed!\r\n");
        Error_Handler();
    }
    Audio_WaitTxMidHalf();
    TIM6->CR1 |= TIM_CR1_CEN;           /* first ADC trigger -> analog blocks start now */

    /* 4. User button B1 (PC13, high when pressed, external pull-down), LD2 (PE1), LD3 (PB14) */
    __HAL_RCC_GPIOC_CLK_ENABLE();
    __HAL_RCC_GPIOE_CLK_ENABLE();
    GPIO_InitStruct.Pin  = GPIO_PIN_13;
    GPIO_InitStruct.Mode = GPIO_MODE_INPUT;
    GPIO_InitStruct.Pull = GPIO_NOPULL;
    HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);
    GPIO_InitStruct.Pin   = GPIO_PIN_1;
    GPIO_InitStruct.Mode  = GPIO_MODE_OUTPUT_PP;
    GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
    HAL_GPIO_Init(GPIOE, &GPIO_InitStruct);
    GPIO_InitStruct.Pin   = GPIO_PIN_14;                 /* LD3 red: lit = equalizer bypassed */
    HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);
    HAL_GPIO_WritePin(GPIOB, GPIO_PIN_14, bypass ? GPIO_PIN_SET : GPIO_PIN_RESET);
    Mimo_ShowAdapt();

    UART_Print("MIMO room equalizer (unbiased MIPNLMS-FX, 2 loudspeakers x 2 mics) started. PCM1808 = mono reference, I2S L/R = loudspeaker 1/2.\r\n");
    UART_Print("Keys: s = sweep measurement, x abort, a adapt run/freeze, q equalizer on/off, l bias on output, z reset weights, w state, +/-/]/[ gain\r\n");
    UART_Print("Until the first calibration has finished the plain stereo signal is played. Calibrate: hold B1 for 3 s and let go (or key 's'); a short press of B1 pauses / continues the adaptation.\r\n");
  /* USER CODE END 2 */

  /* Infinite loop */
    /* USER CODE BEGIN WHILE */
    while (1)
    {
    	/* B1 (works without a PC): short press = the adaptation pauses / continues (LD2 lit = running); hold for MIMO_CAL_HOLD_MS (LD3 flashes when it is long enough) and let go = the sweep calibration starts
    	 * after MIMO_CAL_COUNTDOWN_S seconds (a short press during the countdown cancels it). The progress is shown by LD2 (see Cal_Leds). */
    	{
    	    static GPIO_PinState btnPrev = GPIO_PIN_RESET;
    	    static uint32_t btnTime = 0, pressStart = 0;
    	    static uint8_t  longDone = 0;
    	    GPIO_PinState btn = HAL_GPIO_ReadPin(GPIOC, GPIO_PIN_13);

    	    if ((btn != btnPrev) && (HAL_GetTick() - btnTime >= 50U))
    	    {
    	        btnTime = HAL_GetTick();
    	        btnPrev = btn;
    	        if (btn == GPIO_PIN_SET)
    	        {
    	            pressStart = HAL_GetTick();
    	            longDone   = 0;
    	        }
    	        else if (longDone)                                     /* released after a long hold: start the calibration */
    	        {
    	            if (!ms.active)
    	            {
    	                calPending   = 1U;
    	                calStartTick = HAL_GetTick() + MIMO_CAL_COUNTDOWN_S * 1000U;
    	                UART_Print("Calibration starts after the countdown: keep the room quiet, step away (LD2 blinks). A short press of B1 cancels.\r\n");
    	            }
    	        }
    	        else if (calPending)                                   /* short press during the countdown */
    	        {
    	            calPending = 0U;
    	            UART_Print("Calibration cancelled\r\n");
    	        }
    	        else if (HAL_GetTick() - pressStart < 1000U)           /* short press: the adaptation pauses / continues (the filters are kept) */
    	        {
    	            Mimo_ToggleAdapt();
    	        }
    	    }
    	    if ((btnPrev == GPIO_PIN_SET) && !longDone && (HAL_GetTick() - pressStart >= MIMO_CAL_HOLD_MS))
    	    {
    	        longDone = 1;
    	        ackUntil = HAL_GetTick() + 400U;                       /* LD3 flashes: let go now */
    	    }
    	}
    	if (calPending && ((int32_t)(HAL_GetTick() - calStartTick) >= 0))
    	{
    	    calPending = 0U;
    	    Meas_Start();
    	}
    	Cal_Leds();

    	Gain_HandleUart();                      /* control keys from the serial terminal */
    	Meas_Poll();                            /* sweep measurement state machine */

    	if (shed)
    	{
    	    char m[100];
    	    snprintf(m, sizeof(m), "OVERLOAD: block took %lu us, equalizer switched off ('q' to try again)\r\n",
    	             (unsigned long)(shedCycles / (SystemCoreClock / 1000000U)));
    	    shed = 0U;
    	    HAL_GPIO_WritePin(GPIOB, GPIO_PIN_14, GPIO_PIN_SET);        /* LD3: equalizer is off */
    	    UART_Print(m);
    	}
    	if (errSai) { errSai = 0U; UART_Print("SAI error\r\n"); }
    	if (errAdc) { errAdc = 0U; UART_Print("ADC error\r\n"); }
    	if (errI2s) { errI2s = 0U; UART_Print("I2S error\r\n"); }

    	/* adaptation: one block (256 samples) of x_f stream + IPNLMS-FX update per call */
    	{
    	    const uint32_t t0 = DWT->CYCCNT;
    	    if (MimoFx_Poll(1U) != 0U)
    	    {
    	        const uint32_t dt = DWT->CYCCNT - t0;
    	        if (dt > pollCyclesMax) pollCyclesMax = dt;
    	    }
    	}

    	if (MimoFx_IrInstalled())
    	{
    	    float f1 = 0.0f, f2 = 0.0f;
    	    const float g1 = MimoFx_PeakGainDb(0U, &f1), g2 = MimoFx_PeakGainDb(1U, &f2);
    	    char m[300];
    	    Mimo_ShowAdapt();
    	    snprintf(m, sizeof(m), "Calibration finished: the filters start neutral (filter gain max w1 %+.1f dB @ %.0f Hz, w2 %+.1f dB @ %.0f Hz), the adaptation runs in real time with the music ('a' freezes / resumes)\r\n",
    	             (double)g1, (double)f1, (double)g2, (double)f2);
    	    UART_Print(m);
    	    Mimo_PrintFilterRelation();
    	    UART_Print("The equalizer is now ON and adapting (LD2 lit). Short press of B1 = adaptation pause / continue, hold 3 s = new calibration.\r\n");
    	}

    	/* Audio is handled entirely in the DMA callbacks; this is debug only.
    	 * Expected once per 2 seconds: rx and adc ~ +344, tx ~ +172 (44.1 kHz / 256 frames). */
    	if (!ms.active && (HAL_GetTick() - lastPrint >= 2000))
    	{
    	    lastPrint = HAL_GetTick();

    	    char m[280];
    	    const uint32_t us = SystemCoreClock / 1000000U;
    	    snprintf(m, sizeof(m), "PCM L[%d..%d] R[%d..%d] mics[%u..%u][%u..%u] rx=%lu adc=%lu tx=%lu | pcm=%luus mic=%luus adapt=%luus (block 5805us)\r\n",
    	             (int)rawMin[0], (int)rawMax[0], (int)rawMin[1], (int)rawMax[1],
    	             (unsigned)adcMin[0], (unsigned)adcMax[0], (unsigned)adcMin[1], (unsigned)adcMax[1],
    	             (unsigned long)rxCount, (unsigned long)adcCount, (unsigned long)txCount,
    	             (unsigned long)(procCyclesMax / us), (unsigned long)(procMicCyclesMax / us), (unsigned long)(pollCyclesMax / us));
    	    procCyclesMax = 0;
    	    procMicCyclesMax = 0;
    	    pollCyclesMax = 0;
    	    UART_Print(m);

    	    /* level meter: RMS of the last 2 s per channel (raw, before the gains) */
    	    float out[MIC_SRC_COUNT][2];
    	    __disable_irq();
    	    for (uint32_t k = 0; k < MIC_SRC_COUNT; k++)
    	    {
    	        for (uint32_t c = 0; c < 2U; c++)
    	        {
    	            float rms = (lvlCount[k] > 0U) ? sqrtf(lvlSumSq[k][c] / (float)lvlCount[k]) : 0.0f;
    	            out[k][c] = Audio_Db(rms / 32768.0f);
    	            lvlSumSq[k][c] = 0.0f;
    	        }
    	        lvlCount[k] = 0;
    	    }
    	    __enable_irq();
    	    {   /* loudspeaker signals and the filters: where does the level come from? */
    	        float f1 = 0.0f, f2 = 0.0f;
    	        const float g1 = MimoFx_PeakGainDb(0U, &f1), g2 = MimoFx_PeakGainDb(1U, &f2);
    	        __disable_irq();
    	        const float pl = (float)outPeak[0], pr = (float)outPeak[1];
    	        outPeak[0] = 0;  outPeak[1] = 0;
    	        __enable_irq();
    	        UART_Print("-----------------------------------------------------------------------------\r\n");
    	        snprintf(m, sizeof(m), "speakers: peak L %.1f R %.1f dBFS, clipped samples %lu | filter gain max w1 %+.1f dB @ %.0f Hz, w2 %+.1f dB @ %.0f Hz | adaptation backlog %lu skipped %lu | step x%.3f resets %lu\r\n",
    	                 (double)Audio_Db(pl / 32768.0f), (double)Audio_Db(pr / 32768.0f), (unsigned long)clipCount,
    	                 (double)g1, (double)f1, (double)g2, (double)f2, (unsigned long)MimoFx_Backlog(), (unsigned long)MimoFx_Skipped(),
    	                 (double)MimoFx_MuScale(), (unsigned long)MimoFx_Resets());
    	        UART_Print(m);
    	        {
    	            static uint32_t lastResets = 0U;
    	            if (MimoFx_Resets() != lastResets)
    	            {
    	                lastResets = MimoFx_Resets();
    	                snprintf(m, sizeof(m), "DIVERGED: a filter exceeded %.0f dB gain, weights reset, training restarted with the step size x%.3f (reset no. %lu). Send 'h' to look at h^\r\n",
    	                         (double)MIMO_GAIN_LIMIT_DB, (double)MimoFx_MuScale(), (unsigned long)lastResets);
    	                UART_Print(m);
    	            }
    	        }
    	    }
    	    UART_Print("- variables ------------------------------------------------------------\r\n");
    	    if (MimoFx_HaveIr())
    	    {
    	        float mu, de, ep;
    	        MimoFx_GetGains(&mu, &de, &ep);
    	        snprintf(m, sizeof(m), "real time adaptation (%s): %lu blocks adapted | mu %.4g delta %.2e eps %.2e | largest weight jump %.1f dB | mic error e %.1f dB\r\n",
    	                 MimoFx_GetAdapt() ? "running" : "frozen", (unsigned long)MimoFx_RtBlocks(), (double)mu, (double)de, (double)ep, (double)MimoFx_JumpDb(), (double)MimoFx_ErrDb(0));
    	        UART_Print(m);
    	    }
    	    UART_Print("- levels ---------------------------------------------------------------\r\n");
    	    snprintf(m, sizeof(m), "level dBFS: PCM L %6.1f R %6.1f  mics PF6 %6.1f PA3 %6.1f | PCM gain %.1f dB | adapt=%s eq=%s ir=%s e %.1f dB backlog %lu\r\n",
    	             (double)out[0][0], (double)out[0][1], (double)out[1][0], (double)out[1][1], (double)Audio_Db(pcmGain),
    	             MimoFx_GetAdapt() ? "run" : "frozen", bypass ? "OFF" : "on", MimoFx_HaveIr() ? "yes" : "NO",
    	             (double)MimoFx_ErrDb(0),
    	             (unsigned long)MimoFx_Backlog());
    	    UART_Print(m);
    	}


    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
    }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Supply configuration update enable
  */
  HAL_PWREx_ConfigSupply(PWR_DIRECT_SMPS_SUPPLY);

  /** Configure the main internal regulator output voltage
  */
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE2);

  while(!__HAL_PWR_GET_FLAG(PWR_FLAG_VOSRDY)) {}

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_DIV1;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 4;
  /* The CPU runs at 282 MHz (VOS2 allows 300 MHz) while HCLK / AXI / all buses and the audio clock PLL1Q
   * (11.29 MHz = 256 x 44.1 kHz) stay at 141 MHz: VCO 282.25 -> 564.5 MHz, PLLP 2 -> CPU 282 MHz with AHB /2 = 141 MHz, PLLQ 25 -> 50.
   * Set MIMO_CPU_BOOST to 0 for the 141 MHz clock tree. */
#define MIMO_CPU_BOOST 1
#if MIMO_CPU_BOOST
  RCC_OscInitStruct.PLL.PLLN = 35;
  RCC_OscInitStruct.PLL.PLLP = 2;
  RCC_OscInitStruct.PLL.PLLQ = 50;
  RCC_OscInitStruct.PLL.PLLR = 2;
  RCC_OscInitStruct.PLL.PLLRGE = RCC_PLL1VCIRANGE_3;
  RCC_OscInitStruct.PLL.PLLVCOSEL = RCC_PLL1VCOWIDE;
  RCC_OscInitStruct.PLL.PLLFRACN = 2304;
#else
  RCC_OscInitStruct.PLL.PLLN = 17;
  RCC_OscInitStruct.PLL.PLLP = 2;
  RCC_OscInitStruct.PLL.PLLQ = 25;
  RCC_OscInitStruct.PLL.PLLR = 2;
  RCC_OscInitStruct.PLL.PLLRGE = RCC_PLL1VCIRANGE_3;
  RCC_OscInitStruct.PLL.PLLVCOSEL = RCC_PLL1VCOWIDE;
  RCC_OscInitStruct.PLL.PLLFRACN = 5248;
#endif
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2
                              |RCC_CLOCKTYPE_D3PCLK1|RCC_CLOCKTYPE_D1PCLK1;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.SYSCLKDivider = RCC_SYSCLK_DIV1;
#if MIMO_CPU_BOOST
  RCC_ClkInitStruct.AHBCLKDivider = RCC_HCLK_DIV2;
#else
  RCC_ClkInitStruct.AHBCLKDivider = RCC_HCLK_DIV1;
#endif
  RCC_ClkInitStruct.APB3CLKDivider = RCC_APB3_DIV2;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_APB1_DIV4;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_APB2_DIV2;
  RCC_ClkInitStruct.APB4CLKDivider = RCC_APB4_DIV2;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief Peripherals Common Clock Configuration
  * @retval None
  */
void PeriphCommonClock_Config(void)
{
  RCC_PeriphCLKInitTypeDef PeriphClkInitStruct = {0};

  /** Initializes the peripherals clock
  */
  PeriphClkInitStruct.PeriphClockSelection = RCC_PERIPHCLK_ADC;
  PeriphClkInitStruct.PLL2.PLL2M = 4;
  PeriphClkInitStruct.PLL2.PLL2N = 10;
  PeriphClkInitStruct.PLL2.PLL2P = 2;
  PeriphClkInitStruct.PLL2.PLL2Q = 2;
  PeriphClkInitStruct.PLL2.PLL2R = 2;
  PeriphClkInitStruct.PLL2.PLL2RGE = RCC_PLL2VCIRANGE_3;
  PeriphClkInitStruct.PLL2.PLL2VCOSEL = RCC_PLL2VCOMEDIUM;
  PeriphClkInitStruct.PLL2.PLL2FRACN = 0;
  PeriphClkInitStruct.AdcClockSelection = RCC_ADCCLKSOURCE_PLL2;
  if (HAL_RCCEx_PeriphCLKConfig(&PeriphClkInitStruct) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
    __disable_irq();
    while (1)
    {
    }
  /* USER CODE END Error_Handler_Debug */
}
#ifdef USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
