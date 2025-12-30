///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Spectrum.h                                                          ///
///  Description: Power spectrum estimation and spectral analysis                     ///
///               Periodogram, Welch's method, and spectral density estimation        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_SPECTRUM_H
#define MML_SPECTRUM_H

#include "base/Vector.h"
#include "base/Matrix.h"

#include "FFT.h"
#include "RealFFT.h"
#include "Windowing.h"

#include <cmath>
#include <stdexcept>

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // PowerSpectrum - Power spectral density estimation
    //
    // PURPOSE:
    //   Compute power spectrum (squared magnitude of FFT) for signal analysis.
    //   Estimates how signal power distributes across frequency components.
    //
    // METHODS:
    //   1. Periodogram: Direct FFT -> |X[k]|² / N
    //   2. Welch: Averaged periodogram with overlapping segments (reduces variance)
    //
    // APPLICATIONS:
    //   - Identify dominant frequencies in signals
    //   - Noise analysis and characterization
    //   - System identification
    //   - Vibration analysis
    //
    // REFERENCE: Numerical Recipes spectrum.h, Welch (1967)
    ///////////////////////////////////////////////////////////////////////////////////////////
    class PowerSpectrum
    {
    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Compute - Simple periodogram (direct FFT method)
        //
        // INPUT:
        //   - data: Time-domain signal (real-valued)
        //   - window: Window function (optional, default = Hann)
        //            If empty, Hann window of data.size() is used
        //
        // OUTPUT:
        //   - Power spectrum: P[k] = |X[k]|² / N
        //   - Length: N/2 + 1 (DC to Nyquist frequency)
        //
        // NOTES:
        //   - Returns only positive frequencies (real signal -> symmetric spectrum)
        //   - DC component at index 0, Nyquist at index N/2
        //   - For two-sided spectrum, use FFT directly
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> Compute(const Vector<Real>& data, const Vector<Real>& window = Vector<Real>())
        {
            int N = data.size();
            if (N == 0) {
                throw std::invalid_argument("PowerSpectrum::Compute - empty input vector");
            }

            // Apply window (default to Hann if not provided)
            Vector<Real> windowed = data;
            Vector<Real> win = window;
            if (win.size() == 0) {
                win = Windows::Hann(N);
            }
            if (win.size() != N) {
                throw std::invalid_argument("PowerSpectrum::Compute - window size mismatch");
            }

            // Window normalization factor (compensate for window energy loss)
            Real windowNorm = 0.0;
            for (int i = 0; i < N; i++) {
                windowed[i] *= win[i];
                windowNorm += win[i] * win[i];
            }
            windowNorm /= N;

            // Compute FFT
            Vector<Complex> spectrum = RealFFT::Forward(windowed);

            // Compute power spectrum: |X[k]|² / (N * windowNorm)
            int halfN = N / 2 + 1;
            Vector<Real> power(halfN);
            
            Real scale = 1.0 / (N * windowNorm);
            for (int i = 0; i < halfN; i++) {
                Real mag = std::abs(spectrum[i]);
                power[i] = mag * mag * scale;
            }

            return power;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Welch - Averaged periodogram method (reduced variance)
        //
        // INPUT:
        //   - data: Time-domain signal (real-valued)
        //   - segmentSize: Length of each segment (should be power of 2)
        //   - overlap: Number of overlapping samples between segments
        //   - window: Window function (optional, default = Hann)
        //
        // OUTPUT:
        //   - Averaged power spectrum with reduced variance
        //   - Length: segmentSize/2 + 1
        //
        // ALGORITHM:
        //   1. Divide signal into overlapping segments
        //   2. Apply window to each segment
        //   3. Compute periodogram for each segment
        //   4. Average all periodograms
        //
        // BENEFITS:
        //   - Reduced variance compared to simple periodogram
        //   - Trade-off: Reduced frequency resolution (shorter segments)
        //   - Typical overlap: 50% (segmentSize/2)
        //
        // EXAMPLE:
        //   auto spectrum = PowerSpectrum::Welch(signal, 1024, 512);
        //   // 1024-point segments with 50% overlap
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> Welch(const Vector<Real>& data, int segmentSize, int overlap, 
                                  const Vector<Real>& window = Vector<Real>())
        {
            int N = data.size();
            if (N == 0) {
                throw std::invalid_argument("PowerSpectrum::Welch - empty input vector");
            }
            if (segmentSize <= 0 || segmentSize > N) {
                throw std::invalid_argument("PowerSpectrum::Welch - invalid segment size");
            }
            if (overlap < 0 || overlap >= segmentSize) {
                throw std::invalid_argument("PowerSpectrum::Welch - invalid overlap");
            }

            // Generate window (default to Hann)
            Vector<Real> win = window;
            if (win.size() == 0) {
                win = Windows::Hann(segmentSize);
            }
            if (win.size() != segmentSize) {
                throw std::invalid_argument("PowerSpectrum::Welch - window size must match segment size");
            }

            // Calculate number of segments
            int hop = segmentSize - overlap;
            int numSegments = 1 + (N - segmentSize) / hop;
            
            if (numSegments < 1) {
                throw std::invalid_argument("PowerSpectrum::Welch - signal too short for given parameters");
            }

            // Initialize averaged spectrum
            int halfSize = segmentSize / 2 + 1;
            Vector<Real> avgSpectrum(halfSize, 0.0);

            // Process each segment
            for (int seg = 0; seg < numSegments; seg++) {
                int start = seg * hop;
                
                // Extract segment
                Vector<Real> segment(segmentSize);
                for (int i = 0; i < segmentSize; i++) {
                    segment[i] = data[start + i];
                }

                // Compute periodogram for this segment
                Vector<Real> segmentSpectrum = Compute(segment, win);

                // Accumulate
                for (int i = 0; i < halfSize; i++) {
                    avgSpectrum[i] += segmentSpectrum[i];
                }
            }

            // Average
            for (int i = 0; i < halfSize; i++) {
                avgSpectrum[i] /= numSegments;
            }

            return avgSpectrum;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // FrequencyAxis - Generate frequency axis for power spectrum
        //
        // INPUT:
        //   - n: FFT size (number of time-domain samples)
        //   - sampleRate: Sampling rate in Hz
        //
        // OUTPUT:
        //   - Frequency values in Hz for each bin
        //   - Length: n/2 + 1 (DC to Nyquist)
        //
        // NOTES:
        //   - Frequency resolution: Δf = sampleRate / n
        //   - Nyquist frequency: fₙ = sampleRate / 2
        //   - Frequencies: [0, Δf, 2Δf, ..., fₙ]
        //
        // EXAMPLE:
        //   auto freqs = PowerSpectrum::FrequencyAxis(1024, 44100);
        //   // freqs[0] = 0 Hz, freqs[512] = 22050 Hz (Nyquist)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> FrequencyAxis(int n, Real sampleRate)
        {
            if (n <= 0) {
                throw std::invalid_argument("PowerSpectrum::FrequencyAxis - n must be positive");
            }
            if (sampleRate <= 0) {
                throw std::invalid_argument("PowerSpectrum::FrequencyAxis - sampleRate must be positive");
            }

            int halfN = n / 2 + 1;
            Vector<Real> freqs(halfN);
            Real df = sampleRate / n;

            for (int i = 0; i < halfN; i++) {
                freqs[i] = i * df;
            }

            return freqs;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // ToDB - Convert power to decibels
        //
        // INPUT:
        //   - power: Power spectrum (linear scale)
        //   - reference: Reference power level (default = 1.0)
        //
        // OUTPUT:
        //   - Power in dB: 10 * log₁₀(power / reference)
        //
        // NOTES:
        //   - Common references:
        //     * 1.0 for normalized power
        //     * Maximum power value for relative dB
        //   - Zero/negative values clamped to -∞ dB (in practice, very negative value)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> ToDB(const Vector<Real>& power, Real reference = 1.0)
        {
            if (reference <= 0) {
                throw std::invalid_argument("PowerSpectrum::ToDB - reference must be positive");
            }

            int n = power.size();
            Vector<Real> dB(n);
            const Real minDB = -200.0;  // Practical floor for dB

            for (int i = 0; i < n; i++) {
                if (power[i] > 0) {
                    dB[i] = 10.0 * std::log10(power[i] / reference);
                } else {
                    dB[i] = minDB;
                }
            }

            return dB;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    // Spectrogram - Short-Time Fourier Transform (STFT)
    //
    // PURPOSE:
    //   Compute time-frequency representation of signal.
    //   Shows how frequency content changes over time.
    //
    // ALGORITHM:
    //   1. Slide window across signal
    //   2. Compute FFT at each window position
    //   3. Stack results to create 2D time-frequency matrix
    //
    // APPLICATIONS:
    //   - Audio analysis (speech, music)
    //   - Transient detection
    //   - Time-varying system analysis
    //   - Chirp and modulated signal analysis
    //
    // REFERENCE: Numerical Recipes, Spectral analysis texts
    ///////////////////////////////////////////////////////////////////////////////////////////
    class Spectrogram
    {
    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STFT - Short-Time Fourier Transform
        //
        // INPUT:
        //   - data: Time-domain signal (real-valued)
        //   - windowSize: FFT window size (should be power of 2)
        //   - hopSize: Number of samples to advance between windows
        //   - window: Window function (optional, default = Hann)
        //
        // OUTPUT:
        //   - Matrix of complex spectra: rows = time frames, cols = frequency bins
        //   - Rows: Number of time windows
        //   - Cols: windowSize/2 + 1 (positive frequencies only)
        //
        // NOTES:
        //   - Typical overlap: 50-75% (hopSize = windowSize/4 to windowSize/2)
        //   - Time resolution: hopSize / sampleRate
        //   - Frequency resolution: sampleRate / windowSize
        //   - Trade-off: Larger window = better frequency res, worse time res
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Matrix<Complex> STFT(const Vector<Real>& data, int windowSize, int hopSize,
                                     const Vector<Real>& window = Vector<Real>())
        {
            int N = data.size();
            if (N == 0) {
                throw std::invalid_argument("Spectrogram::STFT - empty input vector");
            }
            if (windowSize <= 0 || windowSize > N) {
                throw std::invalid_argument("Spectrogram::STFT - invalid window size");
            }
            if (hopSize <= 0 || hopSize > windowSize) {
                throw std::invalid_argument("Spectrogram::STFT - invalid hop size");
            }

            // Generate window (default to Hann)
            Vector<Real> win = window;
            if (win.size() == 0) {
                win = Windows::Hann(windowSize);
            }
            if (win.size() != windowSize) {
                throw std::invalid_argument("Spectrogram::STFT - window size mismatch");
            }

            // Calculate number of frames
            int numFrames = 1 + (N - windowSize) / hopSize;
            int numFreqs = windowSize / 2 + 1;

            // Initialize result matrix
            Matrix<Complex> stft(numFrames, numFreqs);

            // Process each frame
            for (int frame = 0; frame < numFrames; frame++) {
                int start = frame * hopSize;

                // Extract and window segment
                Vector<Real> segment(windowSize);
                for (int i = 0; i < windowSize; i++) {
                    segment[i] = data[start + i] * win[i];
                }

                // Compute FFT
                Vector<Complex> spectrum = RealFFT::Forward(segment);

                // Store in matrix
                for (int k = 0; k < numFreqs; k++) {
                    stft(frame, k) = spectrum[k];
                }
            }

            return stft;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // PowerSpectrogram - Magnitude-squared STFT (power over time-frequency)
        //
        // INPUT:
        //   - data: Time-domain signal
        //   - windowSize, hopSize: Same as STFT
        //
        // OUTPUT:
        //   - Matrix of power values: |STFT|²
        //   - Suitable for visualization (log scale recommended)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Matrix<Real> PowerSpectrogram(const Vector<Real>& data, int windowSize, int hopSize)
        {
            Matrix<Complex> stft = STFT(data, windowSize, hopSize);
            
            int numFrames = stft.RowNum();
            int numFreqs = stft.ColNum();
            Matrix<Real> power(numFrames, numFreqs);

            for (int i = 0; i < numFrames; i++) {
                for (int j = 0; j < numFreqs; j++) {
                    Real mag = std::abs(stft(i, j));
                    power(i, j) = mag * mag;
                }
            }

            return power;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // TimeAxis - Generate time axis for spectrogram
        //
        // INPUT:
        //   - numFrames: Number of time frames in spectrogram
        //   - hopSize: Hop size used in STFT
        //   - sampleRate: Sampling rate in Hz
        //
        // OUTPUT:
        //   - Time values in seconds for each frame center
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> TimeAxis(int numFrames, int hopSize, Real sampleRate)
        {
            if (numFrames <= 0) {
                throw std::invalid_argument("Spectrogram::TimeAxis - numFrames must be positive");
            }
            if (hopSize <= 0) {
                throw std::invalid_argument("Spectrogram::TimeAxis - hopSize must be positive");
            }
            if (sampleRate <= 0) {
                throw std::invalid_argument("Spectrogram::TimeAxis - sampleRate must be positive");
            }

            Vector<Real> times(numFrames);
            for (int i = 0; i < numFrames; i++) {
                times[i] = (i * hopSize) / sampleRate;
            }

            return times;
        }
    };

} // namespace MML

#endif // MML_SPECTRUM_H
