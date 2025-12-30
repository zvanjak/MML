#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../mml/algorithms/Fourier/RealFFT.h"
#include "../../mml/algorithms/Fourier/DFT.h"
#include "../../mml/algorithms/Fourier/FFT.h"

using namespace MML;

// Convenience constant
constexpr double PI = Constants::PI;

///////////////////////////////////////////////////////////////////////////////////////////
// Test 1: Basic Forward/Inverse Round Trip
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT - Forward/Inverse Round Trip", "[realfft]")
{
    SECTION("Size 8 - impulse signal")
    {
        Vector<Real> data{1, 0, 0, 0, 0, 0, 0, 0};
        Vector<Complex> spectrum = RealFFT::Forward(data);
        Vector<Real> recovered = RealFFT::Inverse(spectrum);
        
        REQUIRE(recovered.size() == data.size());
        for (int i = 0; i < data.size(); i++)
            CHECK(std::abs(recovered[i] - data[i]) < 1e-10);
    }
    
    SECTION("Size 16 - random signal")
    {
        Vector<Real> data{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                         7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5};
        Vector<Complex> spectrum = RealFFT::Forward(data);
        Vector<Real> recovered = RealFFT::Inverse(spectrum);
        
        REQUIRE(recovered.size() == data.size());
        for (int i = 0; i < data.size(); i++)
            CHECK(std::abs(recovered[i] - data[i]) < 1e-10);
    }
    
    SECTION("Size 32 - sinusoid")
    {
        int n = 32;
        Vector<Real> data(n);
        for (int i = 0; i < n; i++)
            data[i] = std::sin(2 * PI * 3 * i / n);  // 3 cycles
        
        Vector<Complex> spectrum = RealFFT::Forward(data);
        Vector<Real> recovered = RealFFT::Inverse(spectrum);
        
        for (int i = 0; i < n; i++)
            CHECK(std::abs(recovered[i] - data[i]) < 1e-9);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 2: Comparison with Complex FFT on Real Data
//
// CRITICAL: RealFFT should match complex FFT for real inputs
// RealFFT outputs only positive frequencies (exploiting conjugate symmetry)
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT vs Complex FFT - Matching Results", "[realfft]")
{
    auto test_size = [](int n) {
        INFO("Testing size " << n);
        
        // Create real signal
        Vector<Real> real_data(n);
        for (int i = 0; i < n; i++)
            real_data[i] = std::sin(2 * PI * 2 * i / n) + 0.5 * std::cos(2 * PI * 5 * i / n);
        
        // RealFFT transform
        Vector<Complex> real_spectrum = RealFFT::Forward(real_data);
        
        // Complex FFT transform (convert real to complex)
        Vector<Complex> complex_data(n);
        for (int i = 0; i < n; i++)
            complex_data[i] = Complex(real_data[i], 0);
        Vector<Complex> complex_spectrum = FFT::Forward(complex_data);
        
        // Compare: RealFFT gives n/2+1 components, matching first half of complex FFT
        REQUIRE(real_spectrum.size() == n/2 + 1);
        for (int k = 0; k <= n/2; k++)
        {
            INFO("Frequency bin k=" << k);
            CHECK(std::abs(real_spectrum[k].real() - complex_spectrum[k].real()) < 1e-9);
            CHECK(std::abs(real_spectrum[k].imag() - complex_spectrum[k].imag()) < 1e-9);
        }
        
        // Verify conjugate symmetry in complex FFT
        for (int k = 1; k < n/2; k++)
        {
            Complex expected = std::conj(complex_spectrum[n - k]);
            CHECK(std::abs(complex_spectrum[k].real() - expected.real()) < 1e-9);
            CHECK(std::abs(complex_spectrum[k].imag() - expected.imag()) < 1e-9);
        }
    };
    
    SECTION("Various sizes") {
        test_size(8);
        test_size(16);
        test_size(32);
        test_size(64);
        test_size(128);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 3: Spectrum Properties of Real Signals
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT - Conjugate Symmetry and Real Properties", "[realfft]")
{
    SECTION("DC component is real")
    {
        Vector<Real> data{1, 2, 3, 4, 5, 6, 7, 8};
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        // DC (k=0) must be real
        CHECK(std::abs(spectrum[0].imag()) < 1e-10);
    }
    
    SECTION("Nyquist component is real")
    {
        Vector<Real> data{1, 2, 3, 4, 5, 6, 7, 8};
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        // Nyquist (k=N/2) must be real
        int nyquist_idx = data.size() / 2;
        CHECK(std::abs(spectrum[nyquist_idx].imag()) < 1e-10);
    }
    
    SECTION("Pure cosine has zero imaginary part")
    {
        int n = 16;
        Vector<Real> data(n);
        for (int i = 0; i < n; i++)
            data[i] = std::cos(2 * PI * 2 * i / n);  // 2 cycles of cosine
        
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        // Cosine should have zero imaginary part at k=2
        CHECK(std::abs(spectrum[2].imag()) < 1e-9);
    }
    
    SECTION("Pure sine has zero real part")
    {
        int n = 16;
        Vector<Real> data(n);
        for (int i = 0; i < n; i++)
            data[i] = std::sin(2 * PI * 2 * i / n);  // 2 cycles of sine
        
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        // Sine should have zero real part at k=2
        CHECK(std::abs(spectrum[2].real()) < 1e-9);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 4: Parseval's Theorem (Energy Conservation)
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT - Parseval's Theorem", "[realfft]")
{
    auto test_parseval = [](int n) {
        INFO("Testing size " << n);
        
        Vector<Real> data(n);
        for (int i = 0; i < n; i++)
            data[i] = std::sin(2 * PI * 3 * i / n) + std::cos(2 * PI * 5 * i / n);
        
        // Time domain energy
        Real time_energy = 0;
        for (int i = 0; i < n; i++)
            time_energy += data[i] * data[i];
        
        // Frequency domain energy
        Vector<Complex> spectrum = RealFFT::Forward(data);
        Real freq_energy = 0;
        
        // DC term
        freq_energy += std::norm(spectrum[0]);
        
        // Middle frequencies (count twice for conjugate pairs)
        for (int k = 1; k < n/2; k++)
            freq_energy += 2 * std::norm(spectrum[k]);
        
        // Nyquist term
        freq_energy += std::norm(spectrum[n/2]);
        
        freq_energy /= n;  // Normalize
        
        CHECK(std::abs(time_energy - freq_energy) < 1e-8);
    };
    
    SECTION("Various sizes") {
        test_parseval(8);
        test_parseval(16);
        test_parseval(32);
        test_parseval(64);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 5: Sinusoid Peak Detection
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT - Sinusoid Peak Detection", "[realfft]")
{
    SECTION("Single frequency")
    {
        int n = 32;
        int freq = 5;  // 5 cycles in n samples
        
        Vector<Real> data(n);
        for (int i = 0; i < n; i++)
            data[i] = std::sin(2 * PI * freq * i / n);
        
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        // Find peak
        int peak_k = 0;
        Real peak_magnitude = 0;
        for (int k = 0; k <= n/2; k++)
        {
            Real mag = std::abs(spectrum[k]);
            if (mag > peak_magnitude)
            {
                peak_magnitude = mag;
                peak_k = k;
            }
        }
        
        CHECK(peak_k == freq);
        CHECK(peak_magnitude > n/2 * 0.9);  // Should be close to n/2
    }
    
    SECTION("Multiple frequencies")
    {
        int n = 64;
        Vector<Real> data(n);
        for (int i = 0; i < n; i++)
            data[i] = std::sin(2 * PI * 3 * i / n) + 0.5 * std::sin(2 * PI * 7 * i / n);
        
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        // Check that k=3 and k=7 have significant magnitude
        Real mag3 = std::abs(spectrum[3]);
        Real mag7 = std::abs(spectrum[7]);
        
        CHECK(mag3 > 10);  // Should be prominent
        CHECK(mag7 > 5);   // Half amplitude
        CHECK(mag7 < mag3);  // Proper ratio
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 6: In-Place Transform (Packed Format)
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT - In-Place Packed Transform", "[realfft]")
{
    SECTION("Forward in-place matches out-of-place")
    {
        Vector<Real> data{1, 2, 3, 4, 5, 6, 7, 8};
        
        // Out-of-place
        Vector<Complex> spectrum_oop = RealFFT::Forward(data);
        
        // In-place
        Vector<Real> data_copy = data;
        RealFFT::Transform(data_copy, 1);
        
        // Compare packed format
        CHECK(std::abs(data_copy[0] - spectrum_oop[0].real()) < 1e-10);  // DC
        CHECK(std::abs(data_copy[1] - spectrum_oop[4].real()) < 1e-10);  // Nyquist
        
        for (int k = 1; k < 4; k++)
        {
            CHECK(std::abs(data_copy[2*k] - spectrum_oop[k].real()) < 1e-10);
            CHECK(std::abs(data_copy[2*k + 1] - spectrum_oop[k].imag()) < 1e-10);
        }
    }
    
    SECTION("Round-trip in-place")
    {
        Vector<Real> original{1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0};
        Vector<Real> data = original;
        
        RealFFT::Transform(data, 1);   // Forward
        RealFFT::Transform(data, -1);  // Inverse
        
        for (int i = 0; i < original.size(); i++)
            CHECK(std::abs(data[i] - original[i]) < 1e-9);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 7: Edge Cases
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT - Edge Cases", "[realfft]")
{
    SECTION("Size 2")
    {
        Vector<Real> data{1, 2};
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        REQUIRE(spectrum.size() == 2);  // DC + Nyquist
        CHECK(std::abs(spectrum[0].real() - 3.0) < 1e-10);  // DC = sum
        CHECK(std::abs(spectrum[1].real() - (-1.0)) < 1e-10);  // Nyquist = diff
        
        Vector<Real> recovered = RealFFT::Inverse(spectrum);
        CHECK(std::abs(recovered[0] - 1.0) < 1e-10);
        CHECK(std::abs(recovered[1] - 2.0) < 1e-10);
    }
    
    SECTION("Size 4")
    {
        Vector<Real> data{1, 2, 3, 4};
        Vector<Complex> spectrum = RealFFT::Forward(data);
        Vector<Real> recovered = RealFFT::Inverse(spectrum);
        
        for (int i = 0; i < 4; i++)
            CHECK(std::abs(recovered[i] - data[i]) < 1e-10);
    }
    
    SECTION("Constant signal")
    {
        int n = 16;
        Vector<Real> data(n, 5.0);
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        // DC should be n*5, all others zero
        CHECK(std::abs(spectrum[0].real() - n * 5.0) < 1e-9);
        for (int k = 1; k <= n/2; k++)
            CHECK(std::abs(spectrum[k]) < 1e-9);
    }
    
    SECTION("Alternating signal")
    {
        Vector<Real> data{1, -1, 1, -1, 1, -1, 1, -1};
        Vector<Complex> spectrum = RealFFT::Forward(data);
        
        // Should have energy only at Nyquist frequency (k=4)
        CHECK(std::abs(spectrum[0]) < 1e-10);  // No DC
        for (int k = 1; k < 4; k++)
            CHECK(std::abs(spectrum[k]) < 1e-10);
        CHECK(std::abs(spectrum[4]) > 5);  // Strong Nyquist
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 8: Performance Comparison
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT - Performance Check", "[realfft][.performance]")
{
    // Note: This is a sanity check, not a rigorous benchmark
    SECTION("Large transform")
    {
        int n = 1024;
        Vector<Real> data(n);
        for (int i = 0; i < n; i++)
            data[i] = std::sin(2 * PI * 10 * i / n) + std::cos(2 * PI * 23 * i / n);
        
        // Should complete quickly (RealFFT is O(n log n))
        Vector<Complex> spectrum = RealFFT::Forward(data);
        Vector<Real> recovered = RealFFT::Inverse(spectrum);
        
        // Verify correctness
        Real max_error = 0;
        for (int i = 0; i < n; i++)
            max_error = std::max(max_error, std::abs(recovered[i] - data[i]));
        
        CHECK(max_error < 1e-8);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 9: Linearity Property
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("RealFFT - Linearity", "[realfft]")
{
    int n = 32;
    
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; i++)
    {
        x[i] = std::sin(2 * PI * 2 * i / n);
        y[i] = std::cos(2 * PI * 3 * i / n);
    }
    
    Real alpha = 2.0, beta = 3.0;
    Vector<Real> combo(n);
    for (int i = 0; i < n; i++)
        combo[i] = alpha * x[i] + beta * y[i];
    
    Vector<Complex> X = RealFFT::Forward(x);
    Vector<Complex> Y = RealFFT::Forward(y);
    Vector<Complex> Z = RealFFT::Forward(combo);
    
    // FFT(ax + by) = a*FFT(x) + b*FFT(y)
    for (int k = 0; k <= n/2; k++)
    {
        Complex expected = alpha * X[k] + beta * Y[k];
        CHECK(std::abs(Z[k].real() - expected.real()) < 1e-8);
        CHECK(std::abs(Z[k].imag() - expected.imag()) < 1e-8);
    }
}
