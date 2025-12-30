///////////////////////////////////////////////////////////////////////////////////////////
// fourier_spectrum_tests.cpp - Tests for Power Spectrum and Spectrogram
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../mml/algorithms/Fourier/Spectrum.h"
#include "../../mml/algorithms/Fourier/FFT.h"
#include <cmath>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
// Test helpers
///////////////////////////////////////////////////////////////////////////////////////////

// Generate sinusoidal signal
static Vector<Real> GenerateSine(int n, Real frequency, Real sampleRate, Real amplitude = 1.0)
{
    Vector<Real> signal(n);
    Real omega = 2.0 * Constants::PI * frequency;
    for (int i = 0; i < n; i++) {
        Real t = i / sampleRate;
        signal[i] = amplitude * std::sin(omega * t);
    }
    return signal;
}

// Generate white noise
static Vector<Real> GenerateNoise(int n, Real amplitude = 1.0)
{
    Vector<Real> noise(n);
    for (int i = 0; i < n; i++) {
        noise[i] = amplitude * (2.0 * rand() / RAND_MAX - 1.0);
    }
    return noise;
}

// Find peak frequency in spectrum
static int FindPeakIndex(const Vector<Real>& spectrum)
{
    int peakIdx = 0;
    Real peakVal = spectrum[0];
    for (int i = 1; i < spectrum.size(); i++) {
        if (spectrum[i] > peakVal) {
            peakVal = spectrum[i];
            peakIdx = i;
        }
    }
    return peakIdx;
}

///////////////////////////////////////////////////////////////////////////////////////////
// PowerSpectrum::Compute tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("PowerSpectrum::Compute - Single frequency detection", "[fourier][spectrum]")
{
    Real sampleRate = 1000.0;  // 1000 Hz
    Real frequency = 50.0;     // 50 Hz signal
    int n = 1024;
    
    Vector<Real> signal = GenerateSine(n, frequency, sampleRate);
    Vector<Real> spectrum = PowerSpectrum::Compute(signal);
    
    REQUIRE(spectrum.size() == n / 2 + 1);
    
    // Find peak in spectrum
    Vector<Real> freqs = PowerSpectrum::FrequencyAxis(n, sampleRate);
    int peakIdx = FindPeakIndex(spectrum);
    Real peakFreq = freqs[peakIdx];
    
    // Peak should be at 50 Hz (±1 bin tolerance)
    REQUIRE(std::abs(peakFreq - frequency) < (sampleRate / n));
    
    // Peak should be significantly stronger than noise floor
    Real peakPower = spectrum[peakIdx];
    Real avgPower = 0.0;
    for (int i = 0; i < spectrum.size(); i++) {
        if (i != peakIdx && i != peakIdx - 1 && i != peakIdx + 1) {
            avgPower += spectrum[i];
        }
    }
    avgPower /= (spectrum.size() - 3);
    
    REQUIRE(peakPower > 100.0 * avgPower);  // SNR > 20 dB
}

TEST_CASE("PowerSpectrum::Compute - Multiple frequencies", "[fourier][spectrum]")
{
    Real sampleRate = 1000.0;
    int n = 2048;
    
    // Create signal with two frequencies
    Vector<Real> signal1 = GenerateSine(n, 100.0, sampleRate, 1.0);
    Vector<Real> signal2 = GenerateSine(n, 300.0, sampleRate, 0.5);
    Vector<Real> signal(n);
    for (int i = 0; i < n; i++) {
        signal[i] = signal1[i] + signal2[i];
    }
    
    Vector<Real> spectrum = PowerSpectrum::Compute(signal);
    Vector<Real> freqs = PowerSpectrum::FrequencyAxis(n, sampleRate);
    
    // Find two largest peaks
    int peak1Idx = FindPeakIndex(spectrum);
    Real peak1Power = spectrum[peak1Idx];
    
    // Zero out first peak region to find second
    Vector<Real> spectrum2 = spectrum;
    for (int i = std::max(0, peak1Idx - 2); i <= std::min((int)spectrum.size() - 1, peak1Idx + 2); i++) {
        spectrum2[i] = 0;
    }
    int peak2Idx = FindPeakIndex(spectrum2);
    
    Real freq1 = freqs[peak1Idx];
    Real freq2 = freqs[peak2Idx];
    
    // Check that peaks are at correct frequencies
    bool has100Hz = std::abs(freq1 - 100.0) < 2.0 || std::abs(freq2 - 100.0) < 2.0;
    bool has300Hz = std::abs(freq1 - 300.0) < 2.0 || std::abs(freq2 - 300.0) < 2.0;
    
    REQUIRE(has100Hz);
    REQUIRE(has300Hz);
}

TEST_CASE("PowerSpectrum::Compute - Custom window", "[fourier][spectrum]")
{
    Real sampleRate = 1000.0;
    Real frequency = 100.0;
    int n = 512;
    
    Vector<Real> signal = GenerateSine(n, frequency, sampleRate);
    
    // Compare different windows
    Vector<Real> spectrumHann = PowerSpectrum::Compute(signal, Windows::Hann(n));
    Vector<Real> spectrumHamming = PowerSpectrum::Compute(signal, Windows::Hamming(n));
    Vector<Real> spectrumRect = PowerSpectrum::Compute(signal, Windows::Rectangular(n));
    
    // All should detect the same peak
    Vector<Real> freqs = PowerSpectrum::FrequencyAxis(n, sampleRate);
    int peakHann = FindPeakIndex(spectrumHann);
    int peakHamming = FindPeakIndex(spectrumHamming);
    int peakRect = FindPeakIndex(spectrumRect);
    
    REQUIRE(std::abs(freqs[peakHann] - frequency) < 5.0);
    REQUIRE(std::abs(freqs[peakHamming] - frequency) < 5.0);
    REQUIRE(std::abs(freqs[peakRect] - frequency) < 5.0);
}

TEST_CASE("PowerSpectrum::Compute - Edge cases", "[fourier][spectrum]")
{
    SECTION("Empty input")
    {
        Vector<Real> empty;
        REQUIRE_THROWS_AS(PowerSpectrum::Compute(empty), std::invalid_argument);
    }
    
    SECTION("Window size mismatch")
    {
        Vector<Real> signal(1024);
        Vector<Real> wrongWindow(512);
        REQUIRE_THROWS_AS(PowerSpectrum::Compute(signal, wrongWindow), std::invalid_argument);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// PowerSpectrum::Welch tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("PowerSpectrum::Welch - Variance reduction", "[fourier][spectrum]")
{
    Real sampleRate = 1000.0;
    Real frequency = 50.0;
    int n = 8192;  // Long signal
    
    // Create noisy signal
    Vector<Real> signal = GenerateSine(n, frequency, sampleRate);
    Vector<Real> noise = GenerateNoise(n, 0.5);
    for (int i = 0; i < n; i++) {
        signal[i] += noise[i];
    }
    
    // Welch with 50% overlap
    int segmentSize = 1024;
    int overlap = 512;
    Vector<Real> welchSpectrum = PowerSpectrum::Welch(signal, segmentSize, overlap);
    
    REQUIRE(welchSpectrum.size() == segmentSize / 2 + 1);
    
    // Find peak
    Vector<Real> freqs = PowerSpectrum::FrequencyAxis(segmentSize, sampleRate);
    int peakIdx = FindPeakIndex(welchSpectrum);
    Real peakFreq = freqs[peakIdx];
    
    REQUIRE(std::abs(peakFreq - frequency) < (sampleRate / segmentSize + 1.0));
}

TEST_CASE("PowerSpectrum::Welch - Different overlaps", "[fourier][spectrum]")
{
    Real sampleRate = 1000.0;
    Real frequency = 100.0;
    int n = 4096;
    int segmentSize = 512;
    
    Vector<Real> signal = GenerateSine(n, frequency, sampleRate);
    
    // Test different overlap amounts
    Vector<Real> welch0 = PowerSpectrum::Welch(signal, segmentSize, 0);     // No overlap
    Vector<Real> welch50 = PowerSpectrum::Welch(signal, segmentSize, 256);  // 50%
    Vector<Real> welch75 = PowerSpectrum::Welch(signal, segmentSize, 384);  // 75%
    
    Vector<Real> freqs = PowerSpectrum::FrequencyAxis(segmentSize, sampleRate);
    
    // All should find the peak
    int peak0 = FindPeakIndex(welch0);
    int peak50 = FindPeakIndex(welch50);
    int peak75 = FindPeakIndex(welch75);
    
    REQUIRE(std::abs(freqs[peak0] - frequency) < 5.0);
    REQUIRE(std::abs(freqs[peak50] - frequency) < 5.0);
    REQUIRE(std::abs(freqs[peak75] - frequency) < 5.0);
}

TEST_CASE("PowerSpectrum::Welch - Edge cases", "[fourier][spectrum]")
{
    Vector<Real> signal(2048);
    
    SECTION("Invalid segment size")
    {
        REQUIRE_THROWS_AS(PowerSpectrum::Welch(signal, 0, 0), std::invalid_argument);
        REQUIRE_THROWS_AS(PowerSpectrum::Welch(signal, 3000, 0), std::invalid_argument);
    }
    
    SECTION("Invalid overlap")
    {
        REQUIRE_THROWS_AS(PowerSpectrum::Welch(signal, 512, -1), std::invalid_argument);
        REQUIRE_THROWS_AS(PowerSpectrum::Welch(signal, 512, 512), std::invalid_argument);
    }
    
    SECTION("Signal too short")
    {
        Vector<Real> shortSignal(100);
        REQUIRE_THROWS_AS(PowerSpectrum::Welch(shortSignal, 512, 256), std::invalid_argument);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// PowerSpectrum utility functions
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("PowerSpectrum::FrequencyAxis", "[fourier][spectrum]")
{
    Real sampleRate = 44100.0;  // CD quality
    int n = 1024;
    
    Vector<Real> freqs = PowerSpectrum::FrequencyAxis(n, sampleRate);
    
    REQUIRE(freqs.size() == n / 2 + 1);
    REQUIRE(freqs[0] == 0.0);  // DC
    REQUIRE(std::abs(freqs[freqs.size() - 1] - (sampleRate / 2.0)) < 0.01 * (sampleRate / 2.0));  // Nyquist
    
    // Check uniform spacing
    Real df = sampleRate / n;
    for (int i = 1; i < freqs.size(); i++) {
        REQUIRE(std::abs((freqs[i] - freqs[i - 1]) - df) < 0.01 * df);
    }
}

TEST_CASE("PowerSpectrum::ToDB", "[fourier][spectrum]")
{
    Vector<Real> power(5);
    power[0] = 1.0;
    power[1] = 10.0;
    power[2] = 100.0;
    power[3] = 0.1;
    power[4] = 0.01;
    
    Vector<Real> dB = PowerSpectrum::ToDB(power);
    
    REQUIRE(std::abs(dB[0] - 0.0) < 0.001);      // 10*log10(1) = 0 dB
    REQUIRE(std::abs(dB[1] - 10.0) < 0.001);     // 10*log10(10) = 10 dB
    REQUIRE(std::abs(dB[2] - 20.0) < 0.001);     // 10*log10(100) = 20 dB
    REQUIRE(std::abs(dB[3] - (-10.0)) < 0.001);    // 10*log10(0.1) = -10 dB
    REQUIRE(std::abs(dB[4] - (-20.0)) < 0.001);    // 10*log10(0.01) = -20 dB
}

///////////////////////////////////////////////////////////////////////////////////////////
// Spectrogram::STFT tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Spectrogram::STFT - Chirp signal", "[fourier][spectrum][spectrogram]")
{
    Real sampleRate = 1000.0;
    int n = 4096;
    int windowSize = 256;
    int hopSize = 128;
    
    // Generate chirp signal (frequency increasing linearly)
    Vector<Real> chirp(n);
    Real f0 = 50.0;   // Start frequency
    Real f1 = 200.0;  // End frequency
    Real duration = n / sampleRate;
    
    for (int i = 0; i < n; i++) {
        Real t = i / sampleRate;
        Real freq = f0 + (f1 - f0) * t / duration;
        Real phase = 2.0 * Constants::PI * freq * t;
        chirp[i] = std::sin(phase);
    }
    
    // Compute STFT
    Matrix<Complex> stft = Spectrogram::STFT(chirp, windowSize, hopSize);
    
    int numFrames = stft.RowNum();
    int numFreqs = stft.ColNum();
    
    REQUIRE(numFreqs == windowSize / 2 + 1);
    REQUIRE(numFrames > 1);
    
    // Check that peak frequency increases over time
    Vector<Real> freqs = PowerSpectrum::FrequencyAxis(windowSize, sampleRate);
    Vector<Real> peakFreqs(numFrames);
    
    for (int frame = 0; frame < numFrames; frame++) {
        // Find peak in this frame
        int peakIdx = 0;
        Real peakMag = std::abs(stft(frame, 0));
        for (int k = 1; k < numFreqs; k++) {
            Real mag = std::abs(stft(frame, k));
            if (mag > peakMag) {
                peakMag = mag;
                peakIdx = k;
            }
        }
        peakFreqs[frame] = freqs[peakIdx];
    }
    
    // Verify frequency increases (allowing some noise)
    int increasingCount = 0;
    for (int i = 1; i < numFrames; i++) {
        if (peakFreqs[i] >= peakFreqs[i - 1] - 5.0) {  // Allow small decreases
            increasingCount++;
        }
    }
    
    REQUIRE(increasingCount > numFrames * 0.8);  // At least 80% increasing
}

TEST_CASE("Spectrogram::STFT - Stationary signal", "[fourier][spectrum][spectrogram]")
{
    Real sampleRate = 1000.0;
    Real frequency = 100.0;
    int n = 2048;
    int windowSize = 256;
    int hopSize = 128;
    
    Vector<Real> signal = GenerateSine(n, frequency, sampleRate);
    Matrix<Complex> stft = Spectrogram::STFT(signal, windowSize, hopSize);
    
    int numFrames = stft.RowNum();
    int numFreqs = stft.ColNum();
    
    // For stationary signal, peak frequency should be constant
    Vector<Real> freqs = PowerSpectrum::FrequencyAxis(windowSize, sampleRate);
    Vector<int> peakIndices(numFrames);
    
    for (int frame = 0; frame < numFrames; frame++) {
        int peakIdx = 0;
        Real peakMag = std::abs(stft(frame, 0));
        for (int k = 1; k < numFreqs; k++) {
            Real mag = std::abs(stft(frame, k));
            if (mag > peakMag) {
                peakMag = mag;
                peakIdx = k;
            }
        }
        peakIndices[frame] = peakIdx;
    }
    
    // All peaks should be at approximately the same frequency bin
    int mostCommonPeak = peakIndices[0];
    int matchCount = 0;
    for (int i = 0; i < numFrames; i++) {
        if (std::abs(peakIndices[i] - mostCommonPeak) <= 1) {  // ±1 bin tolerance
            matchCount++;
        }
    }
    
    REQUIRE(matchCount > numFrames * 0.9);  // 90% should match
}

TEST_CASE("Spectrogram::PowerSpectrogram", "[fourier][spectrum][spectrogram]")
{
    Real sampleRate = 1000.0;
    Real frequency = 100.0;
    int n = 2048;
    int windowSize = 256;
    int hopSize = 128;
    
    Vector<Real> signal = GenerateSine(n, frequency, sampleRate);
    Matrix<Real> powerSpec = Spectrogram::PowerSpectrogram(signal, windowSize, hopSize);
    
    REQUIRE(powerSpec.RowNum() > 1);
    REQUIRE(powerSpec.ColNum() == windowSize / 2 + 1);
    
    // All power values should be non-negative
    for (int i = 0; i < powerSpec.RowNum(); i++) {
        for (int j = 0; j < powerSpec.ColNum(); j++) {
            REQUIRE(powerSpec(i, j) >= 0.0);
        }
    }
}

TEST_CASE("Spectrogram::TimeAxis", "[fourier][spectrum][spectrogram]")
{
    int numFrames = 10;
    int hopSize = 128;
    Real sampleRate = 1000.0;
    
    Vector<Real> times = Spectrogram::TimeAxis(numFrames, hopSize, sampleRate);
    
    REQUIRE(times.size() == numFrames);
    REQUIRE(times[0] == 0.0);
    
    // Check uniform spacing
    Real dt = hopSize / sampleRate;
    for (int i = 1; i < times.size(); i++) {
        REQUIRE(std::abs((times[i] - times[i - 1]) - dt) < 0.001 * dt);
    }
}

TEST_CASE("Spectrogram edge cases", "[fourier][spectrum][spectrogram]")
{
    Vector<Real> signal(2048);
    
    SECTION("Empty input")
    {
        Vector<Real> empty;
        REQUIRE_THROWS_AS(Spectrogram::STFT(empty, 256, 128), std::invalid_argument);
    }
    
    SECTION("Invalid window size")
    {
        REQUIRE_THROWS_AS(Spectrogram::STFT(signal, 0, 128), std::invalid_argument);
        REQUIRE_THROWS_AS(Spectrogram::STFT(signal, 3000, 128), std::invalid_argument);
    }
    
    SECTION("Invalid hop size")
    {
        REQUIRE_THROWS_AS(Spectrogram::STFT(signal, 256, 0), std::invalid_argument);
        REQUIRE_THROWS_AS(Spectrogram::STFT(signal, 256, 300), std::invalid_argument);
    }
}
