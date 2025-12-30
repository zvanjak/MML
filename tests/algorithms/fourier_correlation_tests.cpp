///////////////////////////////////////////////////////////////////////////////////////////
// Fourier Correlation Tests - Comprehensive test suite for Correlation class
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "../../mml/algorithms/Fourier/Correlation.h"
#include "../../mml/base/Vector.h"
#include <cmath>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
// Helper: Naive cross-correlation for validation
///////////////////////////////////////////////////////////////////////////////////////////
static Vector<Real> NaiveCorrelation(const Vector<Real>& x, const Vector<Real>& y)
{
    int n = x.size();
    int m = y.size();
    int result_len = n + m - 1;
    Vector<Real> result(result_len, 0.0);
    
    // Correlation formula: corr[k] = Σ x[i] * y[i - lag]
    // where lag = k - (m-1), ranging from -(m-1) to (n-1)
    // This is equivalent to convolving x with time-reversed y
    for (int k = 0; k < result_len; k++)
    {
        int lag = k - (m - 1);  // Convert index to lag
        for (int i = 0; i < n; i++)
        {
            int j = i - lag;  // Index into y
            if (j >= 0 && j < m)
                result[k] += x[i] * y[j];
        }
    }
    
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 1: Basic Correlation
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Basic Cross-Correlation", "[correlation]")
{
    SECTION("Impulse correlation")
    {
        Vector<Real> x{1, 0, 0, 0};
        Vector<Real> y{0, 0, 1, 0};
        
        Vector<Real> result = Correlation::Cross(x, y);
        
        REQUIRE(result.size() == 7);  // 4 + 4 - 1
        
        // Peak should be at lag where x[0] aligns with y[2]
        // Zero lag is at index n-1 = 3
        // Peak at lag -2: index 3-2 = 1
        CHECK(std::abs(result[1] - 1.0) < 1e-9);
        
        // All other values should be zero
        for (int i = 0; i < result.size(); i++)
        {
            if (i != 1)
                CHECK(std::abs(result[i]) < 1e-9);
        }
    }
    
    SECTION("Simple sequences")
    {
        Vector<Real> x{1, 2, 3};
        Vector<Real> y{1, 1, 1};
        
        Vector<Real> result = Correlation::Cross(x, y);
        
        REQUIRE(result.size() == 5);  // 3 + 3 - 1
        
        // Check specific values (validated against numpy.correlate(x,y,'full'))
        CHECK(std::abs(result[0] - 1.0) < 1e-9);  // lag=-2
        CHECK(std::abs(result[1] - 3.0) < 1e-9);  // lag=-1
        CHECK(std::abs(result[2] - 6.0) < 1e-9);  // lag=0 (peak)
        CHECK(std::abs(result[3] - 5.0) < 1e-9);  // lag=+1
        CHECK(std::abs(result[4] - 3.0) < 1e-9);  // lag=+2
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 2: vs Naive Implementation
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - FFT vs Naive", "[correlation]")
{
    SECTION("Various sizes")
    {
        for (int n : {4, 7, 10, 16})
        {
            for (int m : {4, 7, 10, 16})
            {
                // Generate random signals
                Vector<Real> x(n);
                Vector<Real> y(m);
                for (int i = 0; i < n; i++)
                    x[i] = std::sin(i * 0.5) + std::cos(i * 0.3);
                for (int i = 0; i < m; i++)
                    y[i] = std::sin(i * 0.7) + 0.5 * std::cos(i * 0.4);
                
                Vector<Real> fft_corr = Correlation::Cross(x, y);
                Vector<Real> naive_corr = NaiveCorrelation(x, y);
                
                REQUIRE(fft_corr.size() == naive_corr.size());
                
                for (int i = 0; i < fft_corr.size(); i++)
                {
                    INFO("n=" << n << ", m=" << m << ", i=" << i);
                    CHECK(std::abs(fft_corr[i] - naive_corr[i]) < 1e-8);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 3: Auto-Correlation
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Auto-Correlation", "[correlation]")
{
    SECTION("Auto-correlation is symmetric")
    {
        Vector<Real> x{1, 2, 3, 2, 1};
        
        Vector<Real> result = Correlation::Auto(x);
        
        REQUIRE(result.size() == 9);  // 5 + 5 - 1
        
        // Should be symmetric around center (zero lag at index 4)
        int center = 4;
        for (int i = 0; i < center; i++)
        {
            CHECK(std::abs(result[center + i] - result[center - i]) < 1e-9);
        }
    }
    
    SECTION("Auto-correlation peak at zero lag")
    {
        Vector<Real> x{1, 3, -2, 4, 1};
        
        Vector<Real> result = Correlation::Auto(x);
        
        // Zero lag is at index n-1 = 4
        int zero_lag_idx = 4;
        
        // Peak should be at zero lag
        for (int i = 0; i < result.size(); i++)
        {
            if (i != zero_lag_idx)
                CHECK(result[i] <= result[zero_lag_idx] + 1e-9);
        }
        
        // Value at zero lag equals signal energy
        Real energy = 0.0;
        for (int i = 0; i < x.size(); i++)
            energy += x[i] * x[i];
        
        CHECK(std::abs(result[zero_lag_idx] - energy) < 1e-8);
    }
    
    SECTION("Periodic signal auto-correlation")
    {
        // Create periodic signal (sine wave)
        int n = 32;
        Vector<Real> x(n);
        Real freq = 2.0 * M_PI / 8.0;  // Period = 8
        for (int i = 0; i < n; i++)
            x[i] = std::sin(i * freq);
        
        Vector<Real> result = Correlation::Auto(x);
        
        // Should have peaks at multiples of period
        int zero_lag = n - 1;
        int period = 8;
        
        // Check peak at zero lag
        CHECK(result[zero_lag] > result[zero_lag + 1]);
        CHECK(result[zero_lag] > result[zero_lag - 1]);
        
        // Check peak at lag = ±period
        CHECK(result[zero_lag + period] > result[zero_lag + period + 1]);
        CHECK(result[zero_lag - period] > result[zero_lag - period - 1]);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 4: Normalized Correlation
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Normalized Cross-Correlation", "[correlation]")
{
    SECTION("Range is [-1, 1]")
    {
        Vector<Real> x{1, 2, 3, 4, 5};
        Vector<Real> y{2, 4, 6, 8, 10};  // y = 2*x (perfect positive correlation)
        
        Vector<Real> result = Correlation::CrossNormalized(x, y);
        
        // All values should be in [-1, 1]
        for (int i = 0; i < result.size(); i++)
        {
            CHECK(result[i] >= -1.0 - 1e-9);
            CHECK(result[i] <= 1.0 + 1e-9);
        }
    }
    
    SECTION("Identical signals → correlation = 1")
    {
        Vector<Real> x{1, 3, 2, 4, 3};
        
        Vector<Real> result = Correlation::CrossNormalized(x, x);
        
        // At zero lag (index n-1), correlation should be 1
        int zero_lag = x.size() - 1;
        CHECK(std::abs(result[zero_lag] - 1.0) < 1e-8);
    }
    
    SECTION("Anti-correlated signals → correlation = -1")
    {
        Vector<Real> x{1, 2, 3, 4, 5};
        Vector<Real> y{-1, -2, -3, -4, -5};  // y = -x
        
        Vector<Real> result = Correlation::CrossNormalized(x, y);
        
        // At zero lag, should be -1
        int zero_lag = x.size() - 1;
        CHECK(std::abs(result[zero_lag] + 1.0) < 1e-8);
    }
    
    SECTION("Orthogonal signals → correlation ≈ 0")
    {
        // Sine and cosine are orthogonal
        int n = 64;
        Vector<Real> x(n);
        Vector<Real> y(n);
        Real freq = 2.0 * M_PI / 16.0;
        for (int i = 0; i < n; i++)
        {
            x[i] = std::sin(i * freq);
            y[i] = std::cos(i * freq);
        }
        
        Vector<Real> result = Correlation::CrossNormalized(x, y);
        
        // At zero lag, should be close to 0
        int zero_lag = n - 1;
        CHECK(std::abs(result[zero_lag]) < 0.1);  // Allow some numerical error
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 5: Template Matching
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Template Matching", "[correlation]")
{
    SECTION("Find template in signal")
    {
        // Signal with template embedded at known position
        Vector<Real> template_signal{1, 2, 3, 2, 1};
        
        // Embed template in longer signal
        Vector<Real> signal(20, 0.0);
        int template_start = 7;  // Embed at position 7
        for (int i = 0; i < template_signal.size(); i++)
            signal[template_start + i] = template_signal[i];
        
        // Correlate template with signal
        Vector<Real> result = Correlation::Cross(signal, template_signal);
        
        // Find peak
        int peak_idx = Correlation::FindPeakLag(result);
        
        // Convert peak index to lag
        int lag = peak_idx - (signal.size() - 1);
        
        // Expected lag: template starts at position 7
        // For cross-correlation, lag relates to where y aligns best with x
        // Peak should occur at lag = -template_start
        CHECK(std::abs(lag + template_start) <= 1);  // Allow ±1 index tolerance
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 6: Complex Signals
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Complex Signals", "[correlation]")
{
    SECTION("Complex cross-correlation")
    {
        Vector<Complex> x{Complex(1,0), Complex(0,1), Complex(-1,0), Complex(0,-1)};
        Vector<Complex> y{Complex(1,1), Complex(-1,1)};
        
        Vector<Complex> result = Correlation::Cross(x, y);
        
        REQUIRE(result.size() == 5);  // 4 + 2 - 1
        
        // Just verify it computes without error and returns correct size
        // (specific values depend on complex arithmetic)
    }
    
    SECTION("Modulated signals")
    {
        // Create complex exponential (carrier wave)
        int n = 16;
        Vector<Complex> x(n);
        Vector<Complex> y(n);
        
        Real freq1 = 2.0 * M_PI / 8.0;
        Real freq2 = 2.0 * M_PI / 4.0;
        
        for (int i = 0; i < n; i++)
        {
            x[i] = std::exp(Complex(0, i * freq1));
            y[i] = std::exp(Complex(0, i * freq2));
        }
        
        Vector<Complex> result = Correlation::Cross(x, y);
        
        REQUIRE(result.size() == 2*n - 1);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 7: Lag Utilities
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Lag Conversion", "[correlation]")
{
    SECTION("ConvertLagToIndex")
    {
        int x_size = 10;
        
        // Zero lag
        CHECK(Correlation::ConvertLagToIndex(0, x_size) == 9);  // n-1
        
        // Positive lag
        CHECK(Correlation::ConvertLagToIndex(1, x_size) == 10);
        CHECK(Correlation::ConvertLagToIndex(5, x_size) == 14);
        
        // Negative lag
        CHECK(Correlation::ConvertLagToIndex(-1, x_size) == 8);
        CHECK(Correlation::ConvertLagToIndex(-5, x_size) == 4);
    }
    
    SECTION("FindPeakLag")
    {
        Vector<Real> corr{1, 2, 5, 3, 4, 2, 1};  // Peak at index 2
        
        int peak = Correlation::FindPeakLag(corr);
        CHECK(peak == 2);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 8: Signal Delay Detection
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Signal Delay Detection", "[correlation]")
{
    SECTION("Detect time delay between signals")
    {
        // Original signal
        int n = 32;
        Vector<Real> original(n);
        for (int i = 0; i < n; i++)
            original[i] = std::sin(i * 0.3) + 0.5 * std::cos(i * 0.7);
        
        // Delayed version (shift by 5 samples)
        int delay = 5;
        Vector<Real> delayed(n, 0.0);
        for (int i = delay; i < n; i++)
            delayed[i] = original[i - delay];
        
        // Correlate
        Vector<Real> result = Correlation::Cross(original, delayed);
        
        // Find peak
        int peak_idx = Correlation::FindPeakLag(result);
        int detected_lag = peak_idx - (original.size() - 1);
        
        // Should detect the delay
        CHECK(detected_lag == -delay);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 9: Edge Cases
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Edge Cases", "[correlation]")
{
    SECTION("Single element signals")
    {
        Vector<Real> x{5.0};
        Vector<Real> y{3.0};
        
        Vector<Real> result = Correlation::Cross(x, y);
        
        REQUIRE(result.size() == 1);
        CHECK(std::abs(result[0] - 15.0) < 1e-9);
    }
    
    SECTION("Zero signals")
    {
        Vector<Real> x(10, 0.0);
        Vector<Real> y(10, 0.0);
        
        Vector<Real> result = Correlation::Cross(x, y);
        
        for (int i = 0; i < result.size(); i++)
            CHECK(std::abs(result[i]) < 1e-9);
    }
    
    SECTION("Different sizes")
    {
        Vector<Real> x{1, 2, 3, 4, 5, 6, 7, 8};  // n=8
        Vector<Real> y{1, 1, 1};                  // m=3
        
        Vector<Real> result = Correlation::Cross(x, y);
        
        REQUIRE(result.size() == 10);  // 8 + 3 - 1
    }
    
    SECTION("Normalized correlation with zero energy")
    {
        Vector<Real> x(5, 0.0);
        Vector<Real> y{1, 2, 3};
        
        Vector<Real> result = Correlation::CrossNormalized(x, y);
        
        // Should not crash, all values should be zero
        for (int i = 0; i < result.size(); i++)
            CHECK(std::abs(result[i]) < 1e-9);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Test 10: Performance (optional benchmark)
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Correlation - Performance", "[correlation][performance]")
{
    SECTION("Large signals")
    {
        int n = 1024;
        Vector<Real> x(n);
        Vector<Real> y(n);
        
        for (int i = 0; i < n; i++)
        {
            x[i] = std::sin(i * 0.1) + std::cos(i * 0.05);
            y[i] = std::cos(i * 0.15) + std::sin(i * 0.08);
        }
        
        // Should complete quickly using FFT
        Vector<Real> result = Correlation::Cross(x, y);
        
        REQUIRE(result.size() == 2*n - 1);
        
        // Verify it's not all zeros (basic sanity check)
        bool has_nonzero = false;
        for (int i = 0; i < result.size(); i++)
        {
            if (std::abs(result[i]) > 1e-6)
            {
                has_nonzero = true;
                break;
            }
        }
        CHECK(has_nonzero);
    }
}
