///////////////////////////////////////////////////////////////////////////////////////////
// Fourier Windowing Tests - Comprehensive test suite for window functions
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../mml/algorithms/Fourier/Windowing.h"
#include "../../mml/base/Vector.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace MML;
using Catch::Matchers::WithinAbs;

TEST_CASE("Windowing - Rectangular Window", "[windowing][rectangular]")
{
    SECTION("Size and values")
    {
        auto window = Windows::Rectangular(8);
        
        REQUIRE(window.size() == 8);
        for (int i = 0; i < 8; i++)
            CHECK(window[i] == 1.0);
    }
    
    SECTION("Zero input throws")
    {
        REQUIRE_THROWS_AS(Windows::Rectangular(0), std::invalid_argument);
        REQUIRE_THROWS_AS(Windows::Rectangular(-5), std::invalid_argument);
    }
}

TEST_CASE("Windowing - Hann Window", "[windowing][hann]")
{
    SECTION("Endpoints are zero")
    {
        auto window = Windows::Hann(10);
        
        CHECK(std::abs(window[0]) < 1e-9);
        CHECK(std::abs(window[9]) < 1e-9);
    }
    
    SECTION("Symmetry")
    {
        auto window = Windows::Hann(16);
        
        for (int i = 0; i < 8; i++)
            CHECK(std::abs(window[i] - window[15 - i]) < 1e-9);
    }
    
    SECTION("Peak at center")
    {
        auto window = Windows::Hann(17);  // Odd size for exact center
        
        int center = 8;
        CHECK(std::abs(window[center] - 1.0) < 1e-9);
        
        // Values decrease from center
        for (int i = 0; i < center; i++)
            CHECK(window[i] < window[center]);
    }
    
    SECTION("Known values")
    {
        auto window = Windows::Hann(4);
        
        // Hann(4): w[n] = 0.5*(1 - cos(2πn/3))
        CHECK(std::abs(window[0] - 0.0) < 1e-9);
        CHECK(std::abs(window[1] - 0.75) < 1e-9);
        CHECK(std::abs(window[2] - 0.75) < 1e-9);
        CHECK(std::abs(window[3] - 0.0) < 1e-9);
    }
}

TEST_CASE("Windowing - Hamming Window", "[windowing][hamming]")
{
    SECTION("Endpoints are non-zero (unlike Hann)")
    {
        auto window = Windows::Hamming(10);
        
        CHECK(window[0] > 0.07);  // ≈ 0.08
        CHECK(window[9] > 0.07);
    }
    
    SECTION("Symmetry")
    {
        auto window = Windows::Hamming(16);
        
        for (int i = 0; i < 8; i++)
            CHECK(std::abs(window[i] - window[15 - i]) < 1e-9);
    }
    
    SECTION("Peak at center")
    {
        auto window = Windows::Hamming(17);
        
        int center = 8;
        CHECK(std::abs(window[center] - 1.0) < 1e-9);
    }
    
    SECTION("Known values")
    {
        auto window = Windows::Hamming(4);
        
        // Hamming(4): w[n] = 0.54 - 0.46*cos(2πn/3)
        CHECK(std::abs(window[0] - 0.08) < 1e-9);
        CHECK(std::abs(window[1] - 0.77) < 1e-9);
        CHECK(std::abs(window[2] - 0.77) < 1e-9);
        CHECK(std::abs(window[3] - 0.08) < 1e-9);
    }
}

TEST_CASE("Windowing - Blackman Window", "[windowing][blackman]")
{
    SECTION("Endpoints near zero")
    {
        auto window = Windows::Blackman(10);
        
        CHECK(std::abs(window[0]) < 0.001);  // Very small but not exactly zero
        CHECK(std::abs(window[9]) < 0.001);
    }
    
    SECTION("Symmetry")
    {
        auto window = Windows::Blackman(16);
        
        for (int i = 0; i < 8; i++)
            CHECK(std::abs(window[i] - window[15 - i]) < 1e-9);
    }
    
    SECTION("Peak at center")
    {
        auto window = Windows::Blackman(17);
        
        int center = 8;
        CHECK(std::abs(window[center] - 1.0) < 1e-9);
    }
    
    SECTION("All positive values")
    {
        auto window = Windows::Blackman(32);
        
        for (int i = 0; i < 32; i++)
            CHECK(window[i] > -1e-12);  // Allow tiny numerical errors
    }
}

TEST_CASE("Windowing - Bartlett Window", "[windowing][bartlett]")
{
    SECTION("Endpoints are zero")
    {
        auto window = Windows::Bartlett(10);
        
        CHECK(std::abs(window[0]) < 1e-9);
        CHECK(std::abs(window[9]) < 1e-9);
    }
    
    SECTION("Symmetry")
    {
        auto window = Windows::Bartlett(16);
        
        for (int i = 0; i < 8; i++)
            CHECK(std::abs(window[i] - window[15 - i]) < 1e-9);
    }
    
    SECTION("Linear ramp")
    {
        auto window = Windows::Bartlett(9);
        
        // Should be triangular: 0, 0.25, 0.5, 0.75, 1.0, 0.75, 0.5, 0.25, 0
        CHECK(std::abs(window[0] - 0.0) < 1e-9);
        CHECK(std::abs(window[1] - 0.25) < 1e-9);
        CHECK(std::abs(window[2] - 0.5) < 1e-9);
        CHECK(std::abs(window[3] - 0.75) < 1e-9);
        CHECK(std::abs(window[4] - 1.0) < 1e-9);
        CHECK(std::abs(window[5] - 0.75) < 1e-9);
        CHECK(std::abs(window[6] - 0.5) < 1e-9);
        CHECK(std::abs(window[7] - 0.25) < 1e-9);
        CHECK(std::abs(window[8] - 0.0) < 1e-9);
    }
}

TEST_CASE("Windowing - Welch Window", "[windowing][welch]")
{
    SECTION("Endpoints are zero")
    {
        auto window = Windows::Welch(10);
        
        CHECK(std::abs(window[0]) < 1e-9);
        CHECK(std::abs(window[9]) < 1e-9);
    }
    
    SECTION("Symmetry")
    {
        auto window = Windows::Welch(16);
        
        for (int i = 0; i < 8; i++)
            CHECK(std::abs(window[i] - window[15 - i]) < 1e-9);
    }
    
    SECTION("Parabolic shape (smoother than Bartlett)")
    {
        auto window = Windows::Welch(5);
        
        // Welch: w[n] = 1 - (2n/(N-1) - 1)²
        // n=0: 1 - (0-1)² = 0
        // n=1: 1 - (0.5-1)² = 0.75
        // n=2: 1 - (1-1)² = 1
        // n=3: 1 - (1.5-1)² = 0.75
        // n=4: 1 - (2-1)² = 0
        CHECK(std::abs(window[0] - 0.0) < 1e-9);
        CHECK(std::abs(window[1] - 0.75) < 1e-9);
        CHECK(std::abs(window[2] - 1.0) < 1e-9);
        CHECK(std::abs(window[3] - 0.75) < 1e-9);
        CHECK(std::abs(window[4] - 0.0) < 1e-9);
    }
}

TEST_CASE("Windowing - Kaiser Window", "[windowing][kaiser]")
{
    SECTION("Beta = 0 gives rectangular")
    {
        auto kaiser = Windows::Kaiser(8, 0.0);
        auto rect = Windows::Rectangular(8);
        
        for (int i = 0; i < 8; i++)
            CHECK(std::abs(kaiser[i] - rect[i]) < 1e-6);
    }
    
    SECTION("Symmetry")
    {
        auto window = Windows::Kaiser(16, 5.0);
        
        for (int i = 0; i < 8; i++)
            CHECK(std::abs(window[i] - window[15 - i]) < 1e-9);
    }
    
    SECTION("Peak at center")
    {
        auto window = Windows::Kaiser(17, 8.6);
        
        int center = 8;
        CHECK(std::abs(window[center] - 1.0) < 1e-9);
    }
    
    SECTION("Higher beta narrows main lobe")
    {
        auto kaiser_low = Windows::Kaiser(32, 2.0);
        auto kaiser_high = Windows::Kaiser(32, 10.0);
        
        // Higher beta → endpoints closer to zero
        CHECK(kaiser_high[0] < kaiser_low[0]);
        CHECK(kaiser_high[1] < kaiser_low[1]);
    }
    
    SECTION("Negative beta throws")
    {
        REQUIRE_THROWS_AS(Windows::Kaiser(8, -1.0), std::invalid_argument);
    }
}

TEST_CASE("Windowing - Gaussian Window", "[windowing][gaussian]")
{
    SECTION("Symmetry")
    {
        auto window = Windows::Gaussian(16, 0.4);
        
        for (int i = 0; i < 8; i++)
            CHECK(std::abs(window[i] - window[15 - i]) < 1e-9);
    }
    
    SECTION("Peak at center")
    {
        auto window = Windows::Gaussian(17, 0.5);
        
        int center = 8;
        CHECK(std::abs(window[center] - 1.0) < 1e-9);
    }
    
    SECTION("All positive values")
    {
        auto window = Windows::Gaussian(32, 0.4);
        
        for (int i = 0; i < 32; i++)
            CHECK(window[i] > 0.0);
    }
    
    SECTION("Larger sigma widens window")
    {
        auto gauss_narrow = Windows::Gaussian(32, 0.2);  // Narrower spread
        auto gauss_wide = Windows::Gaussian(32, 0.6);    // Wider spread
        
        // Smaller sigma → steeper falloff → edges have lower values
        CHECK(gauss_wide[0] > gauss_narrow[0]);
        CHECK(gauss_wide[1] > gauss_narrow[1]);
    }
    
    SECTION("Zero or negative sigma throws")
    {
        REQUIRE_THROWS_AS(Windows::Gaussian(8, 0.0), std::invalid_argument);
        REQUIRE_THROWS_AS(Windows::Gaussian(8, -0.5), std::invalid_argument);
    }
}

TEST_CASE("Windowing - Window Comparison", "[windowing][comparison]")
{
    SECTION("All windows have same size")
    {
        int n = 16;
        
        auto rect = Windows::Rectangular(n);
        auto hann = Windows::Hann(n);
        auto hamming = Windows::Hamming(n);
        auto blackman = Windows::Blackman(n);
        auto bartlett = Windows::Bartlett(n);
        auto welch = Windows::Welch(n);
        auto kaiser = Windows::Kaiser(n, 5.0);
        auto gaussian = Windows::Gaussian(n, 0.4);
        
        CHECK(rect.size() == n);
        CHECK(hann.size() == n);
        CHECK(hamming.size() == n);
        CHECK(blackman.size() == n);
        CHECK(bartlett.size() == n);
        CHECK(welch.size() == n);
        CHECK(kaiser.size() == n);
        CHECK(gaussian.size() == n);
    }
    
    SECTION("Energy reduction (except Rectangular)")
    {
        int n = 128;
        
        auto rect = Windows::Rectangular(n);
        auto hann = Windows::Hann(n);
        auto hamming = Windows::Hamming(n);
        auto blackman = Windows::Blackman(n);
        
        // Compute sum (energy)
        Real sum_rect = 0, sum_hann = 0, sum_hamming = 0, sum_blackman = 0;
        for (int i = 0; i < n; i++)
        {
            sum_rect += rect[i];
            sum_hann += hann[i];
            sum_hamming += hamming[i];
            sum_blackman += blackman[i];
        }
        
        // Windows with tapering have less energy than rectangular
        CHECK(sum_hann < sum_rect);
        CHECK(sum_hamming < sum_rect);
        CHECK(sum_blackman < sum_rect);
        
        // Blackman has widest main lobe → least energy
        CHECK(sum_blackman < sum_hann);
        CHECK(sum_blackman < sum_hamming);
    }
}

TEST_CASE("Windowing - Edge Cases", "[windowing][edge-cases]")
{
    SECTION("Size = 1")
    {
        auto rect = Windows::Rectangular(1);
        auto hann = Windows::Hann(1);
        auto hamming = Windows::Hamming(1);
        
        REQUIRE(rect.size() == 1);
        REQUIRE(hann.size() == 1);
        REQUIRE(hamming.size() == 1);
        
        // All should have value 1.0 at single point
        CHECK(std::abs(rect[0] - 1.0) < 1e-9);
        CHECK(std::abs(hann[0] - 1.0) < 1e-9);
        CHECK(std::abs(hamming[0] - 1.0) < 1e-9);
    }
    
    SECTION("Size = 2")
    {
        auto hann = Windows::Hann(2);
        auto bartlett = Windows::Bartlett(2);
        
        // Hann(2): endpoints at zero
        CHECK(std::abs(hann[0] - 0.0) < 1e-9);
        CHECK(std::abs(hann[1] - 0.0) < 1e-9);
        
        // Bartlett(2): linear
        CHECK(std::abs(bartlett[0] - 0.0) < 1e-9);
        CHECK(std::abs(bartlett[1] - 0.0) < 1e-9);
    }
}

TEST_CASE("Windowing - Performance", "[windowing][performance]")
{
    SECTION("Large window creation")
    {
        int n = 8192;
        
        auto rect = Windows::Rectangular(n);
        auto hann = Windows::Hann(n);
        auto kaiser = Windows::Kaiser(n, 8.6);
        auto gaussian = Windows::Gaussian(n, 0.4);
        
        REQUIRE(rect.size() == n);
        REQUIRE(hann.size() == n);
        REQUIRE(kaiser.size() == n);
        REQUIRE(gaussian.size() == n);
    }
}
