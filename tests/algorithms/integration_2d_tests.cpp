/**
 * @file integration_2d_tests.cpp
 * @brief Comprehensive tests for 2D (double) integration routines
 * 
 * Tests Integrate2D function from Integration2D.h with:
 * - Rectangular domains with constant bounds
 * - Variable bounds (triangular, circular regions)
 * - All integration methods (TRAP, SIMPSON, GAUSS10)
 * - Polynomial, trigonometric, and transcendental integrands
 * - Known analytical results for verification
 */

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/Integration.h"
#include "interfaces/IFunction.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::Integration2DTests
{
	/////////////////////////////////////////////////////////////////////////////
	///                   HELPER BOUND FUNCTIONS                              ///
	/////////////////////////////////////////////////////////////////////////////

	// Constant bounds for rectangular domain [0,1] x [0,1]
	inline Real y_zero(Real) { return REAL(0.0); }
	inline Real y_one(Real) { return REAL(1.0); }
	inline Real y_two(Real) { return REAL(2.0); }
	inline Real y_three(Real) { return REAL(3.0); }
	inline Real y_pi(Real) { return Constants::PI; }
	inline Real y_half_pi(Real) { return Constants::PI / REAL(2.0); }

	// Variable bounds for triangular domain: x ∈ [0,1], y ∈ [0, 1-x]
	inline Real y_1_minus_x(Real x) { return REAL(1.0) - x; }
	
	// Variable bounds for triangular domain: x ∈ [0,1], y ∈ [0, x]
	inline Real y_equals_x(Real x) { return x; }
	
	// Variable bounds for semicircular domain: x ∈ [0,1], y ∈ [0, √(1-x²)]
	inline Real y_sqrt_1_minus_x2(Real x) { return std::sqrt(REAL(1.0) - x * x); }
	
	// Variable bounds for quarter circle: x ∈ [-1,1], y ∈ [0, √(1-x²)]
	inline Real y_neg_sqrt_1_minus_x2(Real x) { return -std::sqrt(REAL(1.0) - x * x); }

	/////////////////////////////////////////////////////////////////////////////
	///                   HELPER SCALAR FUNCTIONS                             ///
	/////////////////////////////////////////////////////////////////////////////

	// f(x,y) = 1 - for area/volume calculations
	class ConstantFunction2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>&) const override { return REAL(1.0); }
	};

	// f(x,y) = x
	class LinearX2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override { return v[0]; }
	};

	// f(x,y) = y
	class LinearY2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override { return v[1]; }
	};

	// f(x,y) = x + y
	class LinearSum2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override { return v[0] + v[1]; }
	};

	// f(x,y) = x * y
	class Product2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override { return v[0] * v[1]; }
	};

	// f(x,y) = x² + y²
	class SumOfSquares2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override 
		{ 
			return v[0] * v[0] + v[1] * v[1]; 
		}
	};

	// f(x,y) = x² * y²
	class ProductOfSquares2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override 
		{ 
			return v[0] * v[0] * v[1] * v[1]; 
		}
	};

	// f(x,y) = e^(x+y)
	class ExponentialSum2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override 
		{ 
			return std::exp(v[0] + v[1]); 
		}
	};

	// f(x,y) = sin(x) * sin(y)
	class SinProduct2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override 
		{ 
			return std::sin(v[0]) * std::sin(v[1]); 
		}
	};

	// f(x,y) = e^(-(x²+y²)) - Gaussian
	class Gaussian2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override 
		{ 
			return std::exp(-(v[0] * v[0] + v[1] * v[1])); 
		}
	};

	// f(x,y) = x³y² - cubic product
	class CubicProduct2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override 
		{ 
			return v[0] * v[0] * v[0] * v[1] * v[1]; 
		}
	};

	// f(x,y) = cos(x)cos(y) - cosine product
	class CosProduct2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override 
		{ 
			return std::cos(v[0]) * std::cos(v[1]); 
		}
	};

	// f(x,y) = e^(-x-y) - exponential negative sum
	class ExpNegSum2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& v) const override 
		{ 
			return std::exp(-v[0] - v[1]); 
		}
	};

	/////////////////////////////////////////////////////////////////////////////
	///           RECTANGULAR DOMAIN TESTS - CONSTANT BOUNDS                  ///
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integrate2D_Constant_UnitSquare_Trap", "[integration][2d][trap]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] 1 dA = 1 (unit square area)
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, TRAP, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_Constant_UnitSquare_Simpson", "[integration][2d][simpson]")
	{
		TEST_PRECISION_INFO();
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-8)));
	}

	TEST_CASE("Integrate2D_Constant_UnitSquare_Gauss10", "[integration][2d][gauss]")
	{
		TEST_PRECISION_INFO();
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, GAUSS10, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Integrate2D_Constant_Rectangle", "[integration][2d]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,2]x[0,3] 1 dA = 6 (rectangle area)
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(2.0), y_zero, y_three);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(6.0), REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_LinearX_UnitSquare", "[integration][2d][polynomial]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] x dA = ∫[0,1] x dx * ∫[0,1] dy = 1/2 * 1 = 0.5
		LinearX2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.5), REAL(1e-8)));
	}

	TEST_CASE("Integrate2D_LinearY_UnitSquare", "[integration][2d][polynomial]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] y dA = ∫[0,1] dx * ∫[0,1] y dy = 1 * 1/2 = 0.5
		LinearY2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.5), REAL(1e-8)));
	}

	TEST_CASE("Integrate2D_LinearSum_UnitSquare", "[integration][2d][polynomial]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] (x+y) dA = 0.5 + 0.5 = 1.0
		LinearSum2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-8)));
	}

	TEST_CASE("Integrate2D_Product_UnitSquare", "[integration][2d][polynomial]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] xy dA = ∫[0,1] x dx * ∫[0,1] y dy = 1/2 * 1/2 = 0.25
		Product2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.25), REAL(1e-8)));
	}

	TEST_CASE("Integrate2D_SumOfSquares_UnitSquare", "[integration][2d][polynomial]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] (x²+y²) dA = ∫[0,1] x² dx + ∫[0,1] y² dy = 1/3 + 1/3 = 2/3
		SumOfSquares2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(2.0)/REAL(3.0), REAL(1e-8)));
	}

	TEST_CASE("Integrate2D_ProductOfSquares_UnitSquare", "[integration][2d][polynomial]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] x²y² dA = ∫[0,1] x² dx * ∫[0,1] y² dy = 1/3 * 1/3 = 1/9
		ProductOfSquares2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(9.0), REAL(1e-8)));
	}

	/////////////////////////////////////////////////////////////////////////////
	///            TRANSCENDENTAL FUNCTIONS ON RECTANGULAR DOMAIN             ///
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integrate2D_Exponential_UnitSquare", "[integration][2d][transcendental]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] e^(x+y) dA = ∫[0,1] e^x dx * ∫[0,1] e^y dy = (e-1)² ≈ 2.9525
		ExponentialSum2D func;
		Real expected = (std::exp(REAL(1.0)) - REAL(1.0)) * (std::exp(REAL(1.0)) - REAL(1.0));
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_SinProduct_UnitSquare", "[integration][2d][transcendental]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] sin(x)sin(y) dA = ∫[0,1] sin(x) dx * ∫[0,1] sin(y) dy 
		// = (1-cos(1))² ≈ 0.2123
		SinProduct2D func;
		Real sinIntegral = REAL(1.0) - std::cos(REAL(1.0));
		Real expected = sinIntegral * sinIntegral;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_SinProduct_PiSquare", "[integration][2d][transcendental]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,π]x[0,π] sin(x)sin(y) dA = 4 (both integrals = 2)
		SinProduct2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), Constants::PI, y_zero, y_pi);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(4.0), REAL(1e-5)));
	}

	TEST_CASE("Integrate2D_Gaussian_UnitSquare", "[integration][2d][transcendental]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1] e^(-(x²+y²)) dA 
		// = ∫[0,1] e^(-x²) dx * ∫[0,1] e^(-y²) dy = (erf(1) * √π/2)² ≈ 0.5577
		Gaussian2D func;
		Real erfIntegral = std::sqrt(Constants::PI) / REAL(2.0) * std::erf(REAL(1.0));
		Real expected = erfIntegral * erfIntegral;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(1e-6)));
	}

	/////////////////////////////////////////////////////////////////////////////
	///              VARIABLE BOUNDS - TRIANGULAR DOMAINS                     ///
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integrate2D_Constant_Triangle", "[integration][2d][triangle]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1-x] 1 dA = area of triangle = 1/2
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_1_minus_x);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.5), REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_LinearX_Triangle", "[integration][2d][triangle]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1-x] x dA
		// Inner: ∫[0,1-x] x dy = x(1-x)
		// Outer: ∫[0,1] x(1-x) dx = ∫[0,1] (x - x²) dx = 1/2 - 1/3 = 1/6
		LinearX2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_1_minus_x);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(6.0), REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_LinearY_Triangle", "[integration][2d][triangle]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1-x] y dA
		// Inner: ∫[0,1-x] y dy = (1-x)²/2
		// Outer: ∫[0,1] (1-x)²/2 dx = [-(1-x)³/6]_0^1 = 1/6
		LinearY2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_1_minus_x);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(6.0), REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_Product_Triangle", "[integration][2d][triangle]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1-x] xy dA
		// Inner: ∫[0,1-x] xy dy = x(1-x)²/2
		// Outer: ∫[0,1] x(1-x)²/2 dx = 1/2 * ∫[0,1] x(1-2x+x²) dx 
		//      = 1/2 * [x²/2 - 2x³/3 + x⁴/4]_0^1 = 1/2 * (1/2 - 2/3 + 1/4) = 1/24
		Product2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_1_minus_x);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(24.0), REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_Constant_LowerTriangle", "[integration][2d][triangle]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,x] 1 dA = area of lower triangle = 1/2
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_equals_x);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.5), REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_SumOfSquares_Triangle", "[integration][2d][triangle]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]x[0,1-x] (x²+y²) dA
		// = ∫[0,1] [x²(1-x) + (1-x)³/3] dx
		// = ∫[0,1] x²(1-x) dx + 1/3 * ∫[0,1] (1-x)³ dx
		// = [x³/3 - x⁴/4]_0^1 + 1/3 * [-(1-x)⁴/4]_0^1
		// = 1/3 - 1/4 + 1/3 * 1/4 = 1/12 + 1/12 = 1/6
		SumOfSquares2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_1_minus_x);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(6.0), REAL(1e-5)));
	}

	/////////////////////////////////////////////////////////////////////////////
	///              VARIABLE BOUNDS - CIRCULAR DOMAINS                       ///
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integrate2D_Constant_QuarterDisk", "[integration][2d][circular]")
	{
		TEST_PRECISION_INFO();
		// ∫∫ over quarter disk (first quadrant): x ∈ [0,1], y ∈ [0, √(1-x²)]
		// Area = π/4 ≈ 0.7854
		// Note: Circular bounds are harder to integrate - relaxed tolerance
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, GAUSS10, REAL(0.0), REAL(1.0), y_zero, y_sqrt_1_minus_x2);
		
		Real expected = Constants::PI / REAL(4.0);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(2e-3)));  // ~0.25% error for circular bounds
	}

	TEST_CASE("Integrate2D_Constant_Semicircle", "[integration][2d][circular]")
	{
		TEST_PRECISION_INFO();
		// ∫∫ over semicircle (upper half): x ∈ [-1,1], y ∈ [0, √(1-x²)]
		// Area = π/2
		// Note: Circular bounds are harder to integrate - relaxed tolerance
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, GAUSS10, -REAL(1.0), REAL(1.0), y_zero, y_sqrt_1_minus_x2);
		
		Real expected = Constants::PI / REAL(2.0);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(2e-3)));  // ~0.13% error for circular bounds
	}

	TEST_CASE("Integrate2D_Constant_FullDisk", "[integration][2d][circular]")
	{
		TEST_PRECISION_INFO();
		// ∫∫ over full unit disk: x ∈ [-1,1], y ∈ [-√(1-x²), √(1-x²)]
		// Area = π
		// Note: Circular bounds are harder to integrate - relaxed tolerance
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, GAUSS10, -REAL(1.0), REAL(1.0), y_neg_sqrt_1_minus_x2, y_sqrt_1_minus_x2);
		
		Real expected = Constants::PI;
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(4e-3)));  // ~0.13% error for full circular bounds
	}

	TEST_CASE("Integrate2D_SumOfSquares_QuarterDisk", "[integration][2d][circular]")
	{
		TEST_PRECISION_INFO();
		// ∫∫ (x²+y²) dA over quarter disk = ∫∫ r² dA
		// In polar: ∫[0,π/2] ∫[0,1] r² * r dr dθ = π/2 * 1/4 = π/8
		SumOfSquares2D func;
		
		auto result = Integrate2D(func, GAUSS10, REAL(0.0), REAL(1.0), y_zero, y_sqrt_1_minus_x2);
		
		Real expected = Constants::PI / REAL(8.0);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(1e-3)));
	}

	TEST_CASE("Integrate2D_Gaussian_QuarterDisk", "[integration][2d][circular]")
	{
		TEST_PRECISION_INFO();
		// ∫∫ e^(-(x²+y²)) dA over quarter disk (radius 1)
		// In polar: ∫[0,π/2] ∫[0,1] e^(-r²) r dr dθ = π/2 * [-e^(-r²)/2]_0^1
		// = π/2 * (1 - e^(-1))/2 = π(1-1/e)/4
		Gaussian2D func;
		
		auto result = Integrate2D(func, GAUSS10, REAL(0.0), REAL(1.0), y_zero, y_sqrt_1_minus_x2);
		
		Real expected = Constants::PI * (REAL(1.0) - std::exp(-REAL(1.0))) / REAL(4.0);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(1e-3)));
	}

	/////////////////////////////////////////////////////////////////////////////
	///              METHOD COMPARISON TESTS                                  ///
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integrate2D_MethodComparison_Polynomial", "[integration][2d][comparison]")
	{
		TEST_PRECISION_INFO();
		// Compare TRAP, SIMPSON, GAUSS10 on polynomial integral
		Product2D func;
		Real expected = REAL(0.25);  // ∫∫[0,1]² xy dA = 1/4
		
		auto trap = Integrate2D(func, TRAP, REAL(0.0), REAL(1.0), y_zero, y_one);
		auto simp = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		auto gauss = Integrate2D(func, GAUSS10, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		// All should converge
		REQUIRE(trap.converged == true);
		REQUIRE(simp.converged == true);
		
		// Gauss10 should be most accurate for polynomial
		REQUIRE_THAT(gauss.value, WithinAbs(expected, REAL(1e-12)));
		REQUIRE_THAT(simp.value, WithinAbs(expected, REAL(1e-8)));
		REQUIRE_THAT(trap.value, WithinAbs(expected, REAL(1e-4)));
	}

	TEST_CASE("Integrate2D_MethodComparison_Transcendental", "[integration][2d][comparison]")
	{
		TEST_PRECISION_INFO();
		// Compare methods on e^(x+y)
		ExponentialSum2D func;
		Real expected = (std::exp(REAL(1.0)) - REAL(1.0)) * (std::exp(REAL(1.0)) - REAL(1.0));
		
		auto trap = Integrate2D(func, TRAP, REAL(0.0), REAL(1.0), y_zero, y_one);
		auto simp = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		auto gauss = Integrate2D(func, GAUSS10, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		// All should be reasonably close
		REQUIRE_THAT(trap.value, WithinAbs(expected, REAL(1e-3)));
		REQUIRE_THAT(simp.value, WithinAbs(expected, REAL(1e-6)));
		REQUIRE_THAT(gauss.value, WithinAbs(expected, REAL(1e-8)));
	}

	/////////////////////////////////////////////////////////////////////////////
	///              EDGE CASES                                               ///
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integrate2D_EmptyDomain", "[integration][2d][edge]")
	{
		TEST_PRECISION_INFO();
		// When x1 = x2, result should be 0
		ConstantFunction2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(1.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.0), REAL(1e-12)));
	}

	TEST_CASE("Integrate2D_ReversedBounds", "[integration][2d][edge]")
	{
		TEST_PRECISION_INFO();
		// Reversed x bounds should give negative of forward integral
		Product2D func;
		
		auto forward = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		auto reverse = Integrate2D(func, SIMPSON, REAL(1.0), REAL(0.0), y_zero, y_one);
		
		REQUIRE_THAT(reverse.value, WithinAbs(-forward.value, REAL(1e-10)));
	}

	TEST_CASE("Integrate2D_NarrowStrip", "[integration][2d][edge]")
	{
		TEST_PRECISION_INFO();
		// Very narrow strip [0,1] x [0, 0.001]
		// ∫∫ 1 dA ≈ 0.001
		ConstantFunction2D func;
		auto y_small = [](Real) { return REAL(0.001); };
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_small);
		
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.001), REAL(1e-8)));
	}

	/////////////////////////////////////////////////////////////////////////////
	///              VERIFICATION WITH KNOWN INTEGRALS                        ///
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integrate2D_KnownIntegral_CubicProduct", "[integration][2d][verification]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,1]² x³y² dA = ∫[0,1] x³ dx * ∫[0,1] y² dy = 1/4 * 1/3 = 1/12
		CubicProduct2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_one);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(12.0), REAL(1e-8)));
	}

	TEST_CASE("Integrate2D_KnownIntegral_CosProduct", "[integration][2d][verification]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,π/2]² cos(x)cos(y) dA = ∫[0,π/2] cos(x) dx * ∫[0,π/2] cos(y) dy = 1 * 1 = 1
		CosProduct2D func;
		
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), Constants::PI/REAL(2.0), y_zero, y_half_pi);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-6)));
	}

	TEST_CASE("Integrate2D_KnownIntegral_ExpNegSum", "[integration][2d][verification]")
	{
		TEST_PRECISION_INFO();
		// ∫∫[0,∞]² e^(-x-y) dA = ∫[0,∞] e^(-x) dx * ∫[0,∞] e^(-y) dy = 1 * 1 = 1
		// We approximate with [0,10] x [0,10] since e^(-10) ≈ 4.5e-5
		ExpNegSum2D func;
		auto y_ten = [](Real) { return REAL(10.0); };
		
		auto result = Integrate2D(func, GAUSS10, REAL(0.0), REAL(10.0), y_zero, y_ten);
		
		// Should be very close to 1 (within e^(-10) error)
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-4)));
	}

	TEST_CASE("Integrate2D_CenterOfMass_Triangle", "[integration][2d][physics]")
	{
		TEST_PRECISION_INFO();
		// Center of mass of uniform triangle with vertices (0,0), (1,0), (0,1)
		// x_cm = ∫∫ x dA / Area = (1/6) / (1/2) = 1/3
		// y_cm = ∫∫ y dA / Area = (1/6) / (1/2) = 1/3
		
		ConstantFunction2D unit;
		LinearX2D xFunc;
		LinearY2D yFunc;
		
		auto area = Integrate2D(unit, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_1_minus_x);
		auto xMoment = Integrate2D(xFunc, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_1_minus_x);
		auto yMoment = Integrate2D(yFunc, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_1_minus_x);
		
		Real x_cm = xMoment.value / area.value;
		Real y_cm = yMoment.value / area.value;
		
		REQUIRE_THAT(area.value, WithinAbs(REAL(0.5), REAL(1e-6)));
		REQUIRE_THAT(x_cm, WithinAbs(REAL(1.0)/REAL(3.0), REAL(1e-5)));
		REQUIRE_THAT(y_cm, WithinAbs(REAL(1.0)/REAL(3.0), REAL(1e-5)));
	}

	TEST_CASE("Integrate2D_MomentOfInertia_Rectangle", "[integration][2d][physics]")
	{
		TEST_PRECISION_INFO();
		// Moment of inertia of rectangle [0,a]x[0,b] about origin
		// I = ∫∫ (x²+y²) dA = ∫[0,a] x² dx * b + a * ∫[0,b] y² dy
		//   = a³b/3 + ab³/3 = ab(a² + b²)/3
		// For a=1, b=2: I = 1*2*(1+4)/3 = 10/3
		
		SumOfSquares2D func;
		auto result = Integrate2D(func, SIMPSON, REAL(0.0), REAL(1.0), y_zero, y_two);
		
		Real expected = REAL(10.0) / REAL(3.0);
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(1e-6)));
	}
}
