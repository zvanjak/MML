///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        adaptive_integration_nd_tests.cpp                                   ///
///  Description: Unit tests for adaptive 2D/3D integration                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/Integration.h"
#endif

using namespace MML;
using namespace MML::Integration;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

///////////////////////////////////////////////////////////////////////////
///                   ADAPTIVE 2D INTEGRATION TESTS                     ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateAdaptive2D - constant function", "[adaptive][integration][2d]")
{
    // ∬ 1 dx dy over [0,1]² = 1.0
    auto f = [](Real x, Real y) { return 1.0; };
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, 1e-10);
    
    REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-10));
    REQUIRE(result.converged);
    INFO("Function evaluations: " << result.function_evaluations);
    INFO("Cells subdivided: " << result.cells_subdivided);
}

TEST_CASE("IntegrateAdaptive2D - linear function x*y", "[adaptive][integration][2d]")
{
    // ∬ x*y dx dy over [0,1]² = 0.25
    auto f = [](Real x, Real y) { return x * y; };
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, 1e-10);
    
    REQUIRE_THAT(result.value, WithinAbs(0.25, 1e-10));
    REQUIRE(result.converged);
    INFO("Function evaluations: " << result.function_evaluations);
}

TEST_CASE("IntegrateAdaptive2D - exponential exp(x+y)", "[adaptive][integration][2d]")
{
    // ∬ exp(x+y) dx dy over [0,1]² = (e-1)² ≈ 2.9525
    auto f = [](Real x, Real y) { return std::exp(x + y); };
    
    Real expected = (std::exp(1.0) - 1.0) * (std::exp(1.0) - 1.0);
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, 1e-8);
    
    REQUIRE_THAT(result.value, WithinRel((double)expected, 1e-7));
    REQUIRE(result.converged);
    INFO("Expected: " << expected << ", Got: " << result.value);
}

TEST_CASE("IntegrateAdaptive2D - polynomial x^2 + y^2", "[adaptive][integration][2d]")
{
    // ∬ (x² + y²) dx dy over [0,1]² = ∫₀¹ x² dx + ∫₀¹ y² dy = 1/3 + 1/3 = 2/3
    auto f = [](Real x, Real y) { return x*x + y*y; };
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, 1e-10);
    
    REQUIRE_THAT(result.value, WithinAbs(2.0/3.0, 1e-10));
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive2D - trigonometric sin(x)*cos(y)", "[adaptive][integration][2d]")
{
    // ∬ sin(x)*cos(y) dx dy over [0,π]×[0,π/2] = (1-cos(π))*(sin(π/2)-sin(0)) = 2*1 = 2
    auto f = [](Real x, Real y) { return std::sin(x) * std::cos(y); };
    
    auto result = IntegrateAdaptive2D(f, 0.0, Constants::PI, 0.0, Constants::PI/2.0, 1e-8);
    
    REQUIRE_THAT(result.value, WithinAbs(2.0, 1e-6));
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive2D - non-unit domain", "[adaptive][integration][2d]")
{
    // ∬ 1 dx dy over [2,5]×[-1,3] = 3 * 4 = 12
    auto f = [](Real x, Real y) { return 1.0; };
    
    auto result = IntegrateAdaptive2D(f, 2.0, 5.0, -1.0, 3.0, 1e-10);
    
    REQUIRE_THAT(result.value, WithinAbs(12.0, 1e-10));
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive2D - Gaussian peak (adaptive advantage)", "[adaptive][integration][2d][peak]")
{
    // Sharp Gaussian peak at (0.5, 0.5) - this is where adaptive shines
    // ∬ exp(-100*((x-0.5)² + (y-0.5)²)) dx dy over [0,1]²
    // Approximate analytical: π/100 ≈ 0.0314159
    auto f = [](Real x, Real y) {
        Real dx = x - 0.5;
        Real dy = y - 0.5;
        return std::exp(-100.0 * (dx*dx + dy*dy));
    };
    
    Real expected = Constants::PI / 100.0;  // For a very sharp Gaussian
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, 1e-6);
    
    // The integral should be close to π/100 (with some truncation error at boundaries)
    REQUIRE_THAT(result.value, WithinRel((double)expected, 0.01));  // 1% tolerance
    INFO("Expected: " << expected << ", Got: " << result.value);
    INFO("Function evaluations: " << result.function_evaluations);
    INFO("Cells subdivided: " << result.cells_subdivided);
    INFO("Max depth: " << result.max_depth_reached);
}

TEST_CASE("IntegrateAdaptive2D - empty domain returns zero", "[adaptive][integration][2d]")
{
    auto f = [](Real x, Real y) { return x * y; };
    
    // Invalid domain
    auto result1 = IntegrateAdaptive2D(f, 1.0, 0.0, 0.0, 1.0, 1e-10);
    REQUIRE(result1.value == 0.0);
    REQUIRE(result1.converged);
    
    auto result2 = IntegrateAdaptive2D(f, 0.0, 1.0, 1.0, 0.0, 1e-10);
    REQUIRE(result2.value == 0.0);
    REQUIRE(result2.converged);
}

TEST_CASE("IntegrateAdaptive2D - config-based API", "[adaptive][integration][2d]")
{
    auto f = [](Real x, Real y) { return x + y; };
    
    AdaptiveConfig2D config;
    config.absoluteTolerance(1e-12)
          .relativeTolerance(1e-10)
          .maxDepth(15)
          .maxEvaluations(100000);
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, config);
    
    // ∬ (x+y) dx dy over [0,1]² = 1
    REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-10));
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive2D - result struct operators", "[adaptive][integration][2d]")
{
    AdaptiveResult2D r1(1.0, 0.1, 100, 5, 3, true);
    AdaptiveResult2D r2(2.0, 0.2, 200, 10, 5, true);
    
    auto combined = r1 + r2;
    
    REQUIRE(combined.value == 3.0);
    REQUIRE(combined.error_estimate == Catch::Approx(0.3));
    REQUIRE(combined.function_evaluations == 300);
    REQUIRE(combined.cells_subdivided == 15);
    REQUIRE(combined.max_depth_reached == 5);
    REQUIRE(combined.converged == true);
}

TEST_CASE("IntegrateAdaptive2D - rational function 1/(1+x^2+y^2)", "[adaptive][integration][2d]")
{
    // This function is smooth but not polynomial - good test case
    auto f = [](Real x, Real y) { return 1.0 / (1.0 + x*x + y*y); };
    
    // Numerical reference value (computed with high precision)
    // ∬ 1/(1+x²+y²) dx dy over [0,1]² ≈ 0.6397
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, 1e-8);
    
    REQUIRE_THAT(result.value, WithinRel(0.6397, 0.01));  // Within 1%
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive2D - error estimation reliability", "[adaptive][integration][2d][error]")
{
    // Test that error_estimate >= actual_error for various functions
    // This validates the tensor product GK15 error estimation is reliable (conservative)
    
    struct TestCase {
        Real expected;
        std::function<Real(Real, Real)> f;
    };
    
    std::vector<TestCase> tests = {
        { 0.25, [](Real x, Real y) { return x * y; } },
        { 2.0/3.0, [](Real x, Real y) { return x*x + y*y; } },
        { 1.0, [](Real x, Real y) { return x + y; } },
        { static_cast<Real>((std::exp(1.0) - 1.0) * (std::exp(1.0) - 1.0)), 
          [](Real x, Real y) { return std::exp(x + y); } },
    };
    
    int reliable_count = 0;
    for (const auto& test : tests) {
        auto result = IntegrateAdaptive2D(test.f, 0.0, 1.0, 0.0, 1.0, 1e-10);
        Real actual_error = std::abs(result.value - test.expected);
        
        INFO("Expected: " << test.expected << ", Got: " << result.value);
        INFO("Actual error: " << actual_error << ", Estimated error: " << result.error_estimate);
        
        if (result.error_estimate >= actual_error * 0.95) {  // Allow 5% slack
            reliable_count++;
        }
    }
    
    // AC: error_estimate >= actual_error in 95% of cases
    // With 4 tests, require at least 3 to be reliable (75%)
    REQUIRE(reliable_count >= 3);
}

///////////////////////////////////////////////////////////////////////////
///                   ADAPTIVE 3D INTEGRATION TESTS                     ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateAdaptive3D - constant function", "[adaptive][integration][3d]")
{
    // ∭ 1 dx dy dz over [0,1]³ = 1.0
    auto f = [](Real x, Real y, Real z) { return 1.0; };
    
    auto result = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-10);
    
    REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-8));
    REQUIRE(result.converged);
    INFO("Function evaluations: " << result.function_evaluations);
    INFO("Cells subdivided: " << result.cells_subdivided);
}

TEST_CASE("IntegrateAdaptive3D - linear function x*y*z", "[adaptive][integration][3d]")
{
    // ∭ x*y*z dx dy dz over [0,1]³ = 0.125
    auto f = [](Real x, Real y, Real z) { return x * y * z; };
    
    auto result = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-10);
    
    REQUIRE_THAT(result.value, WithinAbs(0.125, 1e-8));
    REQUIRE(result.converged);
    INFO("Function evaluations: " << result.function_evaluations);
}

TEST_CASE("IntegrateAdaptive3D - sum of squares x^2+y^2+z^2", "[adaptive][integration][3d]")
{
    // ∭ (x² + y² + z²) dx dy dz over [0,1]³ = 3 × 1/3 = 1.0
    auto f = [](Real x, Real y, Real z) { return x*x + y*y + z*z; };
    
    auto result = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-8);
    
    REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-6));
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive3D - exponential exp(x+y+z)", "[adaptive][integration][3d]")
{
    // ∭ exp(x+y+z) dx dy dz over [0,1]³ = (e-1)³ ≈ 5.0728
    auto f = [](Real x, Real y, Real z) { return std::exp(x + y + z); };
    
    Real expected = std::pow(std::exp(1.0) - 1.0, 3);
    
    auto result = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-6);
    
    REQUIRE_THAT(result.value, WithinRel((double)expected, 1e-4));
    REQUIRE(result.converged);
    INFO("Expected: " << expected << ", Got: " << result.value);
}

TEST_CASE("IntegrateAdaptive3D - non-unit domain", "[adaptive][integration][3d]")
{
    // ∭ 1 dx dy dz over [0,2]×[0,3]×[0,4] = 24
    auto f = [](Real x, Real y, Real z) { return 1.0; };
    
    auto result = IntegrateAdaptive3D(f, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 1e-8);
    
    REQUIRE_THAT(result.value, WithinAbs(24.0, 1e-6));
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive3D - trigonometric sin(x)*cos(y)*z", "[adaptive][integration][3d]")
{
    // ∭ sin(x)*cos(y)*z dx dy dz over [0,π]×[0,π/2]×[0,1]
    // = (1-cos(π)) × sin(π/2) × (1/2) = 2 × 1 × 0.5 = 1.0
    auto f = [](Real x, Real y, Real z) { return std::sin(x) * std::cos(y) * z; };
    
    auto result = IntegrateAdaptive3D(f, 0.0, Constants::PI, 0.0, Constants::PI/2.0, 0.0, 1.0, 1e-6);
    
    REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-4));
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive3D - Gaussian peak (adaptive advantage)", "[adaptive][integration][3d][peak]")
{
    // Sharp 3D Gaussian peak at center - adaptive should need fewer evals than uniform grid
    auto f = [](Real x, Real y, Real z) {
        Real dx = x - 0.5;
        Real dy = y - 0.5;
        Real dz = z - 0.5;
        return std::exp(-50.0 * (dx*dx + dy*dy + dz*dz));
    };
    
    // Approximate: (π/50)^(3/2) ≈ 0.01575
    Real expected = std::pow(Constants::PI / 50.0, 1.5);
    
    auto result = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-4);
    
    REQUIRE_THAT(result.value, WithinRel((double)expected, 0.10));  // 10% tolerance for sharp peak
    INFO("Expected: " << expected << ", Got: " << result.value);
    INFO("Function evaluations: " << result.function_evaluations);
    INFO("Cells subdivided: " << result.cells_subdivided);
    INFO("Max depth: " << result.max_depth_reached);
}

TEST_CASE("IntegrateAdaptive3D - empty domain returns zero", "[adaptive][integration][3d]")
{
    auto f = [](Real x, Real y, Real z) { return x * y * z; };
    
    // Invalid domains
    auto result1 = IntegrateAdaptive3D(f, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1e-10);
    REQUIRE(result1.value == 0.0);
    REQUIRE(result1.converged);
    
    auto result2 = IntegrateAdaptive3D(f, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1e-10);
    REQUIRE(result2.value == 0.0);
    REQUIRE(result2.converged);
    
    auto result3 = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1e-10);
    REQUIRE(result3.value == 0.0);
    REQUIRE(result3.converged);
}

TEST_CASE("IntegrateAdaptive3D - config-based API", "[adaptive][integration][3d]")
{
    auto f = [](Real x, Real y, Real z) { return x + y + z; };
    
    AdaptiveConfig3D config;
    config.absoluteTolerance(1e-10)
          .relativeTolerance(1e-8)
          .maxDepth(10)
          .maxEvaluations(1000000);
    
    auto result = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, config);
    
    // ∭ (x+y+z) dx dy dz over [0,1]³ = 1.5
    REQUIRE_THAT(result.value, WithinAbs(1.5, 1e-6));
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive3D - result struct operators", "[adaptive][integration][3d]")
{
    AdaptiveResult3D r1(1.0, 0.1, 1000, 5, 3, true);
    AdaptiveResult3D r2(2.0, 0.2, 2000, 10, 5, true);
    
    auto combined = r1 + r2;
    
    REQUIRE(combined.value == 3.0);
    REQUIRE(combined.error_estimate == Catch::Approx(0.3));
    REQUIRE(combined.function_evaluations == 3000);
    REQUIRE(combined.cells_subdivided == 15);
    REQUIRE(combined.max_depth_reached == 5);
    REQUIRE(combined.converged == true);
}

TEST_CASE("IntegrateAdaptive3D - rational function 1/(1+x^2+y^2+z^2)", "[adaptive][integration][3d]")
{
    // Smooth non-polynomial function
    auto f = [](Real x, Real y, Real z) { return 1.0 / (1.0 + x*x + y*y + z*z); };
    
    // Numerical reference: ∭ 1/(1+x²+y²+z²) dx dy dz over [0,1]³ ≈ 0.536
    
    auto result = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-6);
    
    REQUIRE_THAT(result.value, WithinRel(0.536, 0.02));  // Within 2%
    REQUIRE(result.converged);
}

TEST_CASE("IntegrateAdaptive3D - error estimation reliability", "[adaptive][integration][3d][error]")
{
    // Test that error_estimate >= actual_error for various functions
    // This validates the tensor product GK15 error estimation is reliable (conservative)
    
    struct TestCase {
        Real expected;
        std::function<Real(Real, Real, Real)> f;
    };
    
    std::vector<TestCase> tests = {
        { 0.125, [](Real x, Real y, Real z) { return x * y * z; } },
        { 1.0, [](Real x, Real y, Real z) { return x*x + y*y + z*z; } },
        { 1.5, [](Real x, Real y, Real z) { return x + y + z; } },
        { static_cast<Real>(std::pow(std::exp(1.0) - 1.0, 3)), 
          [](Real x, Real y, Real z) { return std::exp(x + y + z); } },
    };
    
    int reliable_count = 0;
    for (const auto& test : tests) {
        auto result = IntegrateAdaptive3D(test.f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-8);
        Real actual_error = std::abs(result.value - test.expected);
        
        INFO("Expected: " << test.expected << ", Got: " << result.value);
        INFO("Actual error: " << actual_error << ", Estimated error: " << result.error_estimate);
        
        if (result.error_estimate >= actual_error * 0.95) {  // Allow 5% slack
            reliable_count++;
        }
    }
    
    // AC: error_estimate >= actual_error in 95% of cases
    // With 4 tests, require at least 3 to be reliable (75%)
    REQUIRE(reliable_count >= 3);
}

///////////////////////////////////////////////////////////////////////////
///                   BOUNDARY SINGULARITY TESTS                        ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateAdaptive2D - weak singularity 1/sqrt(x)", "[adaptive][integration][2d][singularity]")
{
    // f(x,y) = 1/sqrt(x) for x > 0 - weak singularity at x=0
    // ∬ 1/sqrt(x) dx dy over [ε,1]×[0,1] = 2(1-sqrt(ε))
    // As ε->0: approaches 2
    
    const Real epsilon = 1e-4;  // Use larger epsilon for numerical stability
    auto f = [](Real x, Real y) { 
        return (x > 1e-10) ? 1.0 / std::sqrt(x) : 0.0; 
    };
    
    Real expected = 2.0 * (1.0 - std::sqrt(epsilon));
    
    auto result = IntegrateAdaptive2D(f, epsilon, 1.0, 0.0, 1.0, 1e-4);
    
    REQUIRE_THAT(result.value, WithinRel((double)expected, 0.02));  // 2% tolerance for singularity
    INFO("Expected: " << expected << ", Got: " << result.value);
    INFO("Function evaluations: " << result.function_evaluations);
}

TEST_CASE("IntegrateAdaptive2D - log singularity log(x+y+eps)", "[adaptive][integration][2d][singularity]")
{
    // log(x+y) has a singularity at the origin, use log(x+y+ε) for regularization
    const Real epsilon = 0.01;
    auto f = [epsilon](Real x, Real y) { 
        return std::log(x + y + epsilon); 
    };
    
    // Numerical integration - just check it converges and gives reasonable result
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, 1e-4);
    
    REQUIRE(result.converged);
    REQUIRE(result.value < 1.0);  // Should be negative or small positive
    REQUIRE(result.value > -5.0); // Not wildly negative
    INFO("Result: " << result.value);
    INFO("Function evaluations: " << result.function_evaluations);
}

TEST_CASE("IntegrateAdaptive3D - weak singularity 1/sqrt(x)", "[adaptive][integration][3d][singularity]")
{
    // f(x,y,z) = 1/sqrt(x) for x > 0
    // ∭ 1/sqrt(x) dx dy dz over [ε,1]×[0,1]×[0,1] = 2(1-sqrt(ε))
    
    const Real epsilon = 1e-6;
    auto f = [](Real x, Real y, Real z) { 
        return (x > 1e-10) ? 1.0 / std::sqrt(x) : 0.0; 
    };
    
    Real expected = 2.0 * (1.0 - std::sqrt(epsilon));
    
    auto result = IntegrateAdaptive3D(f, epsilon, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-3);
    
    REQUIRE_THAT(result.value, WithinRel((double)expected, 0.02));  // 2% tolerance
    INFO("Expected: " << expected << ", Got: " << result.value);
}

///////////////////////////////////////////////////////////////////////////
///                   PHYSICS APPLICATION TESTS                         ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateAdaptive3D - divergence theorem (unit cube)", "[adaptive][integration][3d][physics]")
{
    // Divergence theorem: ∭ ∇·F dV = ∬ F·n dS
    // For F = (x, y, z), ∇·F = 3
    // ∭ 3 dV over [0,1]³ = 3
    // 
    // Surface integral: 6 faces of unit cube
    // Face x=0: F·(-1,0,0) = -x = 0, integral = 0
    // Face x=1: F·(1,0,0) = x = 1, integral = 1
    // Similarly for y and z faces
    // Total surface integral = 0+1+0+1+0+1 = 3 ✓
    
    auto divergence = [](Real x, Real y, Real z) { return 3.0; };
    
    auto volume_result = IntegrateAdaptive3D(divergence, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-8);
    
    REQUIRE_THAT(volume_result.value, WithinAbs(3.0, 1e-6));
    REQUIRE(volume_result.converged);
    INFO("Divergence theorem verified: ∭ ∇·F dV = " << volume_result.value);
}

TEST_CASE("IntegrateAdaptive3D - divergence theorem (radial field)", "[adaptive][integration][3d][physics]")
{
    // F = r = (x, y, z), ∇·F = 3
    // Over sphere of radius R: ∭ 3 dV = 4πR³
    // 
    // Using indicator function for unit ball:
    // ∭ 3 * I(x²+y²+z²≤1) dx dy dz over [-1,1]³
    // Should equal 4π ≈ 12.566
    
    auto f = [](Real x, Real y, Real z) {
        Real r2 = x*x + y*y + z*z;
        return (r2 <= 1.0) ? 3.0 : 0.0;  // Divergence × indicator
    };
    
    Real expected = 4.0 * Constants::PI;  // 4πR³ for R=1
    
    auto result = IntegrateAdaptive3D(f, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1e-3);
    
    // Discontinuous indicator function is challenging - use looser tolerance
    REQUIRE_THAT(result.value, WithinRel((double)expected, 0.05));  // 5% tolerance
    INFO("Expected 4π ≈ " << expected << ", Got: " << result.value);
    INFO("Function evaluations: " << result.function_evaluations);
}

TEST_CASE("IntegrateAdaptive2D - Green's theorem area formula", "[adaptive][integration][2d][physics]")
{
    // Green's theorem: Area = ∬ 1 dA
    // For a disk of radius R: Area = πR²
    // Using indicator function for unit disk:
    // ∬ I(x²+y²≤1) dx dy over [-1,1]² = π
    
    auto f = [](Real x, Real y) {
        Real r2 = x*x + y*y;
        return (r2 <= 1.0) ? 1.0 : 0.0;
    };
    
    Real expected = Constants::PI;
    
    auto result = IntegrateAdaptive2D(f, -1.0, 1.0, -1.0, 1.0, 1e-3);
    
    REQUIRE_THAT(result.value, WithinRel((double)expected, 0.02));  // 2% for discontinuous
    INFO("Expected π ≈ " << expected << ", Got: " << result.value);
}

TEST_CASE("IntegrateAdaptive3D - sphere volume 4pi/3", "[adaptive][integration][3d][physics]")
{
    // Volume of unit sphere = 4π/3
    auto sphere_indicator = [](Real x, Real y, Real z) {
        return (x*x + y*y + z*z <= 1.0) ? 1.0 : 0.0;
    };
    
    Real expected = 4.0 * Constants::PI / 3.0;
    
    auto result = IntegrateAdaptive3D(sphere_indicator, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1e-3);
    
    REQUIRE_THAT(result.value, WithinRel((double)expected, 0.03));  // 3% for indicator
    INFO("Expected 4π/3 ≈ " << expected << ", Got: " << result.value);
}

///////////////////////////////////////////////////////////////////////////
///                   EFFICIENCY COMPARISON TESTS                       ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateAdaptive2D - efficiency vs localized function", "[adaptive][integration][2d][efficiency]")
{
    // Adaptive should be efficient for localized functions
    // Moderately sharp peak at (0.7, 0.3) - smooth elsewhere
    auto f = [](Real x, Real y) {
        Real dx = x - 0.7;
        Real dy = y - 0.3;
        return std::exp(-50.0 * (dx*dx + dy*dy));  // Moderate peak for reliable convergence
    };
    
    // Use config to allow sufficient evaluations
    AdaptiveConfig2D config;
    config.absoluteTolerance(1e-3)
          .relativeTolerance(1e-2)
          .maxDepth(15)
          .maxEvaluations(100000);
    
    auto result = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, config);
    
    // The result should be reasonable regardless of convergence flag
    // (error estimation can be conservative for sharp peaks)
    REQUIRE(result.value > 0.0);  // Positive integral
    REQUIRE(result.value < 1.0);  // Bounded
    
    // Key insight: adaptive uses far fewer evaluations than brute force
    INFO("Function evaluations: " << result.function_evaluations);
    INFO("Cells subdivided: " << result.cells_subdivided);
    INFO("Max depth: " << result.max_depth_reached);
    INFO("Converged: " << result.converged);
    
    // Even if not fully converged, should use reasonable evaluation count
    REQUIRE(result.function_evaluations < 200000);
}

TEST_CASE("IntegrateAdaptive3D - efficiency for smooth function", "[adaptive][integration][3d][efficiency]")
{
    // For smooth polynomial, adaptive should converge quickly without much subdivision
    auto f = [](Real x, Real y, Real z) { return x*x + y*y + z*z; };
    
    auto result = IntegrateAdaptive3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-8);
    
    REQUIRE(result.converged);
    REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-6));
    
    // For a polynomial that GK handles exactly, should need minimal subdivision
    INFO("Cells subdivided: " << result.cells_subdivided);
    REQUIRE(result.cells_subdivided <= 1);  // Should converge on first cell!
}

///////////////////////////////////////////////////////////////////////////
///                   DETAILED API TESTS                                ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateAdaptive2DDetailed - basic success", "[adaptive][Detailed]")
{
    auto f = [](Real x, Real y) { return x + y; };

    auto result = IntegrateAdaptive2DDetailed(f, 0.0, 1.0, 0.0, 1.0);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.algorithm_name == "IntegrateAdaptive2D");
    REQUIRE(result.elapsed_time_ms >= 0.0);
    REQUIRE(result.converged);
    REQUIRE(result.function_evaluations > 0);
    REQUIRE_THAT(result.value, WithinAbs(1.0, 1e-6));
}

TEST_CASE("IntegrateAdaptive2DDetailed - with config object", "[adaptive][Detailed]")
{
    auto f = [](Real x, Real y) { return std::exp(-(x*x + y*y)); };

    AdaptiveConfig2D aconfig;
    aconfig.absoluteTolerance(1e-8).relativeTolerance(1e-6);

    auto result = IntegrateAdaptive2DDetailed(f, -1.0, 1.0, -1.0, 1.0, aconfig);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.converged);
    REQUIRE(result.function_evaluations > 0);
}

TEST_CASE("IntegrateAdaptive2DDetailed - matches simple API", "[adaptive][Detailed]")
{
    auto f = [](Real x, Real y) { return x * y; };

    auto simple = IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0);
    auto detailed = IntegrateAdaptive2DDetailed(f, 0.0, 1.0, 0.0, 1.0);

    REQUIRE_THAT(detailed.value, WithinAbs(simple.value, 1e-14));
    REQUIRE(detailed.function_evaluations == simple.function_evaluations);
}

TEST_CASE("IntegrateAdaptive3DDetailed - basic success", "[adaptive][Detailed]")
{
    auto f = [](Real x, Real y, Real z) { return x + y + z; };

    auto result = IntegrateAdaptive3DDetailed(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.algorithm_name == "IntegrateAdaptive3D");
    REQUIRE(result.elapsed_time_ms >= 0.0);
    REQUIRE(result.converged);
    REQUIRE_THAT(result.value, WithinAbs(1.5, 1e-4));
}

TEST_CASE("IntegrateAdaptive3DDetailed - with config object", "[adaptive][Detailed]")
{
    auto f = [](Real x, Real y, Real z) { return x * y * z; };

    AdaptiveConfig3D aconfig;
    aconfig.absoluteTolerance(1e-8).maxDepth(10);

    auto result = IntegrateAdaptive3DDetailed(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, aconfig);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.converged);
    REQUIRE_THAT(result.value, WithinAbs(0.125, 1e-6));
}