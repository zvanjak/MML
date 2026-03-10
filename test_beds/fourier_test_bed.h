#if !defined __MML_FOURIER_TEST_BED_H
#define __MML_FOURIER_TEST_BED_H

#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                       FOURIER TEST DATA STRUCTURES                                     //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Represents a test signal with known spectral properties
     */
    struct TestSignal
    {
        std::string name;                               ///< Descriptive name
        std::vector<Real> data;                         ///< Time-domain signal
        Real sampleRate = 1.0;                          ///< Samples per second
        
        // Known spectral properties
        std::vector<Real> expectedPeakFrequencies;      ///< Frequencies where peaks should appear
        std::vector<Real> expectedPeakAmplitudes;       ///< Amplitudes of frequency peaks
        Real fundamentalFrequency = 0.0;                ///< Fundamental frequency (if periodic)
        
        // Signal characteristics
        std::string category;                           ///< pure-tone, multi-tone, chirp, noise, etc.
        std::string description;                        ///< Detailed description
        bool isRealValued = true;                       ///< Whether signal is real or complex
        bool isPeriodic = false;                        ///< Whether signal is periodic
        Real signalEnergy = 0.0;                        ///< Total energy (for Parseval test)
    };

    /**
     * @brief Test case for convolution operations
     */
    struct ConvolutionTest
    {
        std::string name;
        std::vector<Real> signal;                       ///< Input signal
        std::vector<Real> kernel;                       ///< Convolution kernel
        std::vector<Real> expectedLinear;               ///< Expected linear convolution result
        std::vector<Real> expectedCircular;             ///< Expected circular convolution result
        std::string description;
    };

    /**
     * @brief Test case for correlation operations
     */
    struct CorrelationTest
    {
        std::string name;
        std::vector<Real> signal1;
        std::vector<Real> signal2;
        std::vector<Real> expectedCross;                ///< Expected cross-correlation
        std::vector<Real> expectedAuto;                 ///< Expected auto-correlation (signal1)
        int expectedPeakLag = 0;                        ///< Lag where peak occurs
        std::string description;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 1: PURE TONE SIGNALS                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Single sinusoid at known frequency
     */
    inline TestSignal getPureToneSignal()
    {
        const int N = 128;
        const Real freq = 5.0;      // 5 Hz
        const Real sampleRate = 64.0;
        const Real amplitude = 2.0;
        
        TestSignal test;
        test.name = "Pure Tone (5 Hz)";
        test.sampleRate = sampleRate;
        test.data.resize(N);
        
        Real energy = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real t = i / sampleRate;
            test.data[i] = amplitude * std::sin(2.0 * M_PI * freq * t);
            energy += test.data[i] * test.data[i];
        }
        
        test.expectedPeakFrequencies = {freq};
        test.expectedPeakAmplitudes = {amplitude};
        test.fundamentalFrequency = freq;
        test.category = "pure-tone";
        test.description = "Single sinusoid at 5 Hz. Should have single spectral peak.";
        test.isPeriodic = true;
        test.signalEnergy = energy;
        
        return test;
    }

    /**
     * @brief Cosine wave (tests real FFT conjugate symmetry)
     */
    inline TestSignal getCosineToneSignal()
    {
        const int N = 256;
        const Real freq = 8.0;
        const Real sampleRate = 128.0;
        const Real amplitude = 1.5;
        
        TestSignal test;
        test.name = "Cosine Tone (8 Hz)";
        test.sampleRate = sampleRate;
        test.data.resize(N);
        
        Real energy = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real t = i / sampleRate;
            test.data[i] = amplitude * std::cos(2.0 * M_PI * freq * t);
            energy += test.data[i] * test.data[i];
        }
        
        test.expectedPeakFrequencies = {freq};
        test.expectedPeakAmplitudes = {amplitude};
        test.fundamentalFrequency = freq;
        test.category = "pure-tone";
        test.description = "Cosine at 8 Hz. Tests real FFT conjugate symmetry.";
        test.isPeriodic = true;
        test.signalEnergy = energy;
        
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 2: MULTI-TONE SIGNALS                                      //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Sum of two sinusoids
     */
    inline TestSignal getDualToneSignal()
    {
        const int N = 512;
        const Real freq1 = 10.0;
        const Real freq2 = 25.0;
        const Real sampleRate = 256.0;
        const Real amp1 = 1.0;
        const Real amp2 = 0.5;
        
        TestSignal test;
        test.name = "Dual Tone (10 Hz + 25 Hz)";
        test.sampleRate = sampleRate;
        test.data.resize(N);
        
        Real energy = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real t = i / sampleRate;
            test.data[i] = amp1 * std::sin(2.0 * M_PI * freq1 * t) +
                          amp2 * std::sin(2.0 * M_PI * freq2 * t);
            energy += test.data[i] * test.data[i];
        }
        
        test.expectedPeakFrequencies = {freq1, freq2};
        test.expectedPeakAmplitudes = {amp1, amp2};
        test.fundamentalFrequency = 5.0;  // GCD of 10 and 25
        test.category = "multi-tone";
        test.description = "Sum of two sinusoids. Tests frequency resolution.";
        test.isPeriodic = true;
        test.signalEnergy = energy;
        
        return test;
    }

    /**
     * @brief Harmonic series (fundamental + overtones)
     */
    inline TestSignal getHarmonicSignal()
    {
        const int N = 512;
        const Real fundamental = 5.0;
        const Real sampleRate = 256.0;
        const int numHarmonics = 5;
        
        TestSignal test;
        test.name = "Harmonic Series (5 Hz fundamental)";
        test.sampleRate = sampleRate;
        test.data.resize(N);
        
        std::vector<Real> amplitudes = {1.0, 0.5, 0.33, 0.25, 0.2};
        
        Real energy = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real t = i / sampleRate;
            test.data[i] = 0.0;
            for (int h = 1; h <= numHarmonics; h++)
            {
                test.data[i] += amplitudes[h-1] * std::sin(2.0 * M_PI * h * fundamental * t);
            }
            energy += test.data[i] * test.data[i];
        }
        
        for (int h = 1; h <= numHarmonics; h++)
        {
            test.expectedPeakFrequencies.push_back(h * fundamental);
            test.expectedPeakAmplitudes.push_back(amplitudes[h-1]);
        }
        
        test.fundamentalFrequency = fundamental;
        test.category = "multi-tone";
        test.description = "Harmonic series with 1/n amplitude decay. Tests overtone analysis.";
        test.isPeriodic = true;
        test.signalEnergy = energy;
        
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 3: CHIRP SIGNALS                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Linear frequency sweep (chirp)
     */
    inline TestSignal getLinearChirpSignal()
    {
        const int N = 1024;
        const Real f0 = 5.0;        // Start frequency
        const Real f1 = 50.0;       // End frequency
        const Real sampleRate = 512.0;
        const Real duration = N / sampleRate;
        const Real k = (f1 - f0) / duration;  // Chirp rate
        
        TestSignal test;
        test.name = "Linear Chirp (5-50 Hz)";
        test.sampleRate = sampleRate;
        test.data.resize(N);
        
        Real energy = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real t = i / sampleRate;
            Real phase = 2.0 * M_PI * (f0 * t + 0.5 * k * t * t);
            test.data[i] = std::sin(phase);
            energy += test.data[i] * test.data[i];
        }
        
        test.category = "chirp";
        test.description = "Linear frequency sweep from 5 to 50 Hz. Tests STFT/spectrogram.";
        test.isPeriodic = false;
        test.signalEnergy = energy;
        
        return test;
    }

    /**
     * @brief Exponential frequency sweep
     */
    inline TestSignal getExponentialChirpSignal()
    {
        const int N = 1024;
        const Real f0 = 2.0;        // Start frequency
        const Real f1 = 64.0;       // End frequency
        const Real sampleRate = 512.0;
        const Real duration = N / sampleRate;
        const Real k = std::pow(f1/f0, 1.0/duration);
        
        TestSignal test;
        test.name = "Exponential Chirp (2-64 Hz)";
        test.sampleRate = sampleRate;
        test.data.resize(N);
        
        Real energy = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real t = i / sampleRate;
            Real freq = f0 * std::pow(k, t);
            Real phase = 2.0 * M_PI * f0 * (std::pow(k, t) - 1.0) / std::log(k);
            test.data[i] = std::sin(phase);
            energy += test.data[i] * test.data[i];
        }
        
        test.category = "chirp";
        test.description = "Exponential frequency sweep. Tests time-frequency analysis.";
        test.isPeriodic = false;
        test.signalEnergy = energy;
        
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 4: IMPULSE AND STEP SIGNALS                                //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Unit impulse (Kronecker delta)
     */
    inline TestSignal getImpulseSignal()
    {
        const int N = 128;
        const int impulsePos = N / 2;
        
        TestSignal test;
        test.name = "Unit Impulse";
        test.sampleRate = 64.0;
        test.data.resize(N, 0.0);
        test.data[impulsePos] = 1.0;
        
        test.category = "impulse";
        test.description = "Single impulse at center. Spectrum should be flat (white).";
        test.isPeriodic = false;
        test.signalEnergy = 1.0;
        
        return test;
    }

    /**
     * @brief Unit step function
     */
    inline TestSignal getStepSignal()
    {
        const int N = 256;
        const int stepPos = N / 2;
        
        TestSignal test;
        test.name = "Unit Step";
        test.sampleRate = 128.0;
        test.data.resize(N);
        
        for (int i = 0; i < N; i++)
        {
            test.data[i] = (i >= stepPos) ? 1.0 : 0.0;
        }
        
        test.category = "step";
        test.description = "Step function. Tests handling of discontinuities.";
        test.isPeriodic = false;
        test.signalEnergy = N - stepPos;
        
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 5: NOISE SIGNALS                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief White noise (synthetic - deterministic for testing)
     */
    inline TestSignal getWhiteNoiseSignal()
    {
        const int N = 512;
        
        TestSignal test;
        test.name = "Pseudo White Noise";
        test.sampleRate = 256.0;
        test.data.resize(N);
        
        // Deterministic "noise" using sine functions at many frequencies
        Real energy = 0.0;
        for (int i = 0; i < N; i++)
        {
            test.data[i] = 0.0;
            for (int k = 1; k <= 20; k++)
            {
                test.data[i] += std::sin(2.0 * M_PI * k * i / N + k) / (20.0);
            }
            energy += test.data[i] * test.data[i];
        }
        
        test.category = "noise";
        test.description = "Synthetic white noise. Spectrum should be relatively flat.";
        test.isPeriodic = false;
        test.signalEnergy = energy;
        
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 6: WINDOWED SIGNALS                                        //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Rectangular pulse
     */
    inline TestSignal getRectangularPulseSignal()
    {
        const int N = 256;
        const int pulseStart = 64;
        const int pulseEnd = 192;
        
        TestSignal test;
        test.name = "Rectangular Pulse";
        test.sampleRate = 128.0;
        test.data.resize(N, 0.0);
        
        for (int i = pulseStart; i < pulseEnd; i++)
        {
            test.data[i] = 1.0;
        }
        
        test.category = "windowed";
        test.description = "Rectangular pulse. Spectrum should be sinc-like.";
        test.isPeriodic = false;
        test.signalEnergy = pulseEnd - pulseStart;
        
        return test;
    }

    /**
     * @brief Gaussian pulse
     */
    inline TestSignal getGaussianPulseSignal()
    {
        const int N = 256;
        const Real sigma = 20.0;
        const Real center = N / 2.0;
        
        TestSignal test;
        test.name = "Gaussian Pulse";
        test.sampleRate = 128.0;
        test.data.resize(N);
        
        Real energy = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real x = (i - center) / sigma;
            test.data[i] = std::exp(-0.5 * x * x);
            energy += test.data[i] * test.data[i];
        }
        
        test.category = "windowed";
        test.description = "Gaussian pulse. Spectrum should also be Gaussian.";
        test.isPeriodic = false;
        test.signalEnergy = energy;
        
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CONVOLUTION TEST CASES                                              //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Simple convolution with known result
     */
    inline ConvolutionTest getSimpleConvolutionTest()
    {
        ConvolutionTest test;
        test.name = "Simple 3-point convolution";
        test.signal = {1.0, 2.0, 3.0};
        test.kernel = {0.5, 1.0, 0.5};
        
        // Linear: [0.5, 2.0, 4.0, 4.0, 1.5]
        // Calculation: signal[i]*kernel[0] + signal[i+1]*kernel[1] + signal[i+2]*kernel[2]
        // i=-2: 0*0 + 0*0 + 1*0.5 = 0.5
        // i=-1: 0*0 + 1*0.5 + 2*1.0 = 2.0
        // i=0:  1*0.5 + 2*1.0 + 3*0.5 = 4.0
        // i=1:  2*0.5 + 3*1.0 + 0*0.5 = 4.0
        // i=2:  3*0.5 + 0*0 + 0*0 = 1.5
        test.expectedLinear = {0.5, 2.0, 4.0, 4.0, 1.5};
        
        // Circular (same length as signal)
        test.expectedCircular = {2.0, 4.0, 4.5};
        
        test.description = "Simple 3-point signals. Tests basic convolution correctness.";
        
        return test;
    }

    /**
     * @brief Impulse response test (convolution with impulse = identity)
     */
    inline ConvolutionTest getImpulseResponseTest()
    {
        ConvolutionTest test;
        test.name = "Impulse response";
        test.signal = {0.0, 0.0, 1.0, 0.0, 0.0};
        test.kernel = {1.0, 2.0, 3.0, 2.0, 1.0};
        
        // Result should be kernel shifted to impulse position
        test.expectedLinear = {0.0, 0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 0.0};
        
        test.description = "Convolution with impulse. Result should be kernel shifted.";
        
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CORRELATION TEST CASES                                              //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Auto-correlation of periodic signal
     */
    inline CorrelationTest getPeriodicAutoCorrelationTest()
    {
        const int N = 16;
        const Real freq = 1.0;
        const Real sampleRate = 8.0;
        
        CorrelationTest test;
        test.name = "Auto-correlation of sine";
        test.signal1.resize(N);
        
        for (int i = 0; i < N; i++)
        {
            Real t = i / sampleRate;
            test.signal1[i] = std::sin(2.0 * M_PI * freq * t);
        }
        
        test.signal2 = test.signal1;
        // For full correlation (length 2N-1), zero lag is at index N-1
        test.expectedPeakLag = N - 1;  // Index 15 for N=16
        test.description = "Auto-correlation of sinusoid. Peak at zero lag (center of full correlation).";
        
        return test;
    }

    /**
     * @brief Cross-correlation with time delay
     */
    inline CorrelationTest getDelayedCrossCorrelationTest()
    {
        const int N = 32;
        const int delay = 5;
        
        CorrelationTest test;
        test.name = "Cross-correlation with delay";
        test.signal1.resize(N);
        test.signal2.resize(N, 0.0);
        
        // Signal 1: Gaussian pulse at center
        Real center = N / 2.0;
        Real sigma = 3.0;
        for (int i = 0; i < N; i++)
        {
            Real x = (i - center) / sigma;
            test.signal1[i] = std::exp(-0.5 * x * x);
        }
        
        // Signal 2: Same pulse delayed
        for (int i = delay; i < N; i++)
        {
            test.signal2[i] = test.signal1[i - delay];
        }
        
        // For full correlation of length 2N-1, center is at N-1
        // A delay of 'delay' samples appears at index (N-1) - delay
        test.expectedPeakLag = N - 1 - delay;  // Index 26 for N=32, delay=5
        test.description = "Cross-correlation of delayed signals. Peak offset from center indicates delay.";
        
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                          TEST COLLECTION GENERATORS                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Get all test signals
     */
    inline std::vector<TestSignal> getAllSignals()
    {
        return {
            getPureToneSignal(),
            getCosineToneSignal(),
            getDualToneSignal(),
            getHarmonicSignal(),
            getLinearChirpSignal(),
            getExponentialChirpSignal(),
            getImpulseSignal(),
            getStepSignal(),
            getWhiteNoiseSignal(),
            getRectangularPulseSignal(),
            getGaussianPulseSignal()
        };
    }

    /**
     * @brief Get signals by category
     */
    inline std::vector<TestSignal> getSignalsByCategory(const std::string& category)
    {
        std::vector<TestSignal> result;
        auto all = getAllSignals();
        for (const auto& sig : all)
        {
            if (sig.category == category)
                result.push_back(sig);
        }
        return result;
    }

    /**
     * @brief Get all convolution tests
     */
    inline std::vector<ConvolutionTest> getAllConvolutionTests()
    {
        return {
            getSimpleConvolutionTest(),
            getImpulseResponseTest()
        };
    }

    /**
     * @brief Get all correlation tests
     */
    inline std::vector<CorrelationTest> getAllCorrelationTests()
    {
        return {
            getPeriodicAutoCorrelationTest(),
            getDelayedCrossCorrelationTest()
        };
    }

} // namespace MML::TestBeds

#endif // __MML_FOURIER_TEST_BED_H
