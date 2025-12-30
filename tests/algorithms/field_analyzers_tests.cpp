///////////////////////////////////////////////////////////////////////////////////////////
// FieldAnalyzers Tests - Tests for ScalarFieldAnalyzer and VectorFieldAnalyzer
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "MMLBase.h"
#include "algorithms/FieldAnalyzers.h"
#include "core/Fields.h"

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace FieldAnalyzersTests
{
	/////////////////////////////////////////////////////////////////////////////////////
	///                       TEST SCALAR FIELDS                                      ///
	/////////////////////////////////////////////////////////////////////////////////////

	// Harmonic function: f(x,y,z) = x² - y² (Laplacian = 2 - 2 = 0)
	class HarmonicField2D : public IScalarFunction<3>
	{
	public:
		Real operator()(const VectorN<Real, 3>& pos) const override
		{
			return pos[0] * pos[0] - pos[1] * pos[1];
		}
	};

	// Non-harmonic function: f(x,y,z) = x² + y² + z² (Laplacian = 6)
	class NonHarmonicField : public IScalarFunction<3>
	{
	public:
		Real operator()(const VectorN<Real, 3>& pos) const override
		{
			return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
		}
	};

	// Field with critical point at origin: f = x² + y² + z²
	class ParaboloidField : public IScalarFunction<3>
	{
	public:
		Real operator()(const VectorN<Real, 3>& pos) const override
		{
			return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
		}
	};

	/////////////////////////////////////////////////////////////////////////////////////
	///                       TEST VECTOR FIELDS                                      ///
	/////////////////////////////////////////////////////////////////////////////////////

	// Solenoidal field: F = (-y, x, 0) - rotation around z-axis
	// div(F) = ∂(-y)/∂x + ∂x/∂y + ∂0/∂z = 0 + 0 + 0 = 0
	class RotationField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>({-pos[1], pos[0], REAL(0.0)});
		}
	};

	// Non-solenoidal field: F = (x, y, z) - radial outward
	// div(F) = 1 + 1 + 1 = 3
	class RadialOutwardField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>({pos[0], pos[1], pos[2]});
		}
	};

	// Irrotational (conservative) field: F = (x, y, z) = ∇(r²/2)
	// curl(F) = 0
	class GradientField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>({pos[0], pos[1], pos[2]});
		}
	};

	// Non-irrotational field: F = (-y, x, 0)
	// curl(F) = (0, 0, 2)
	class VortexField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>({-pos[1], pos[0], REAL(0.0)});
		}
	};

	// Constant field: F = (1, 0, 0)
	// div = 0, curl = 0 (both solenoidal and irrotational)
	class ConstantField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& /*pos*/) const override
		{
			return VectorN<Real, 3>({REAL(1.0), REAL(0.0), REAL(0.0)});
		}
	};

	/////////////////////////////////////////////////////////////////////////////////////
	///                     SCALAR FIELD ANALYZER TESTS                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ScalarFieldAnalyzer::GradientAt", "[field_analyzer][scalar]")
	{
			TEST_PRECISION_INFO();
		// f(x,y,z) = x² + y² + z²
		// ∇f = (2x, 2y, 2z)
		ParaboloidField field;
		ScalarFieldAnalyzer analyzer(field);

		Vec3Cart pos(REAL(1.0), REAL(2.0), REAL(3.0));
		Vec3Cart grad = analyzer.GradientAt(pos);

		REQUIRE_THAT(grad.X(), WithinRel(REAL(2.0), REAL(0.01)));
		REQUIRE_THAT(grad.Y(), WithinRel(REAL(4.0), REAL(0.01)));
		REQUIRE_THAT(grad.Z(), WithinRel(REAL(6.0), REAL(0.01)));
	}

	TEST_CASE("ScalarFieldAnalyzer::IsCriticalPoint", "[field_analyzer][scalar]")
	{
			TEST_PRECISION_INFO();
		ParaboloidField field;
		ScalarFieldAnalyzer analyzer(field);

		// Origin is a critical point (minimum)
		REQUIRE(analyzer.IsCriticalPoint(Vec3Cart(0, 0, 0), 1e-4));

		// (1, 1, 1) is NOT a critical point
		REQUIRE_FALSE(analyzer.IsCriticalPoint(Vec3Cart(1, 1, 1), 1e-4));
	}

	TEST_CASE("ScalarFieldAnalyzer::FindCriticalPoints", "[field_analyzer][scalar]")
	{
			TEST_PRECISION_INFO();
		ParaboloidField field;
		ScalarFieldAnalyzer analyzer(field);

		auto criticalPoints = analyzer.FindCriticalPoints(
			Vec3Cart(-1, -1, -1), Vec3Cart(1, 1, 1), 10, REAL(0.2));

		// Should find at least one critical point near origin
		REQUIRE(criticalPoints.size() >= 1);
		
		// Check that one of them is close to origin
		bool foundOrigin = false;
		for (const auto& pt : criticalPoints) {
			if (pt.NormL2() < REAL(0.3)) {
				foundOrigin = true;
				break;
			}
		}
		REQUIRE(foundOrigin);
	}

	TEST_CASE("ScalarFieldAnalyzer::LaplacianAt", "[field_analyzer][scalar]")
	{
			TEST_PRECISION_INFO();
		// f = x² + y² + z², Laplacian = 6
		NonHarmonicField field;
		ScalarFieldAnalyzer analyzer(field);

		Real laplacian = analyzer.LaplacianAt(Vec3Cart(1, 1, 1));
		REQUIRE_THAT(laplacian, WithinRel(REAL(6.0), REAL(0.05)));
	}

	TEST_CASE("ScalarFieldAnalyzer::IsHarmonic_true", "[field_analyzer][scalar][harmonic]")
	{
			TEST_PRECISION_INFO();
		// f = x² - y² is harmonic (Laplacian = 0)
		HarmonicField2D field;
		ScalarFieldAnalyzer analyzer(field);

		REQUIRE(analyzer.IsHarmonic(Vec3Cart(-1, -1, -1), Vec3Cart(1, 1, 1), 5, REAL(0.01)));
	}

	TEST_CASE("ScalarFieldAnalyzer::IsHarmonic_false", "[field_analyzer][scalar][harmonic]")
	{
			TEST_PRECISION_INFO();
		// f = x² + y² + z² is NOT harmonic (Laplacian = 6)
		NonHarmonicField field;
		ScalarFieldAnalyzer analyzer(field);

		REQUIRE_FALSE(analyzer.IsHarmonic(Vec3Cart(-1, -1, -1), Vec3Cart(1, 1, 1), 5, REAL(0.1)));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                     VECTOR FIELD ANALYZER TESTS                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorFieldAnalyzer::DivergenceAt", "[field_analyzer][vector]")
	{
			TEST_PRECISION_INFO();
		// F = (x, y, z), div = 3
		RadialOutwardField field;
		VectorFieldAnalyzer analyzer(field);

		Real div = analyzer.DivergenceAt(Vec3Cart(1, 2, 3));
		REQUIRE_THAT(div, WithinRel(REAL(3.0), REAL(0.01)));
	}

	TEST_CASE("VectorFieldAnalyzer::CurlAt", "[field_analyzer][vector]")
	{
			TEST_PRECISION_INFO();
		// F = (-y, x, 0), curl = (0, 0, 2)
		VortexField field;
		VectorFieldAnalyzer analyzer(field);

		Vec3Cart curl = analyzer.CurlAt(Vec3Cart(1, 2, 3));
		REQUIRE_THAT(curl.X(), WithinAbs(REAL(0.0), REAL(0.01)));
		REQUIRE_THAT(curl.Y(), WithinAbs(REAL(0.0), REAL(0.01)));
		REQUIRE_THAT(curl.Z(), WithinRel(REAL(2.0), REAL(0.01)));
	}

	TEST_CASE("VectorFieldAnalyzer::IsSolenoidal_true", "[field_analyzer][vector][solenoidal]")
	{
			TEST_PRECISION_INFO();
		// Rotation field is solenoidal
		RotationField field;
		VectorFieldAnalyzer analyzer(field);

		REQUIRE(analyzer.IsSolenoidal(Vec3Cart(-2, -2, -2), Vec3Cart(2, 2, 2), 5, REAL(0.01)));
	}

	TEST_CASE("VectorFieldAnalyzer::IsSolenoidal_false", "[field_analyzer][vector][solenoidal]")
	{
			TEST_PRECISION_INFO();
		// Radial field is NOT solenoidal (has sources)
		RadialOutwardField field;
		VectorFieldAnalyzer analyzer(field);

		REQUIRE_FALSE(analyzer.IsSolenoidal(Vec3Cart(-2, -2, -2), Vec3Cart(2, 2, 2), 5, REAL(0.01)));
	}

	TEST_CASE("VectorFieldAnalyzer::IsIrrotational_true", "[field_analyzer][vector][irrotational]")
	{
			TEST_PRECISION_INFO();
		// Gradient field is irrotational
		GradientField field;
		VectorFieldAnalyzer analyzer(field);

		REQUIRE(analyzer.IsIrrotational(Vec3Cart(-2, -2, -2), Vec3Cart(2, 2, 2), 5, REAL(0.01)));
	}

	TEST_CASE("VectorFieldAnalyzer::IsIrrotational_false", "[field_analyzer][vector][irrotational]")
	{
			TEST_PRECISION_INFO();
		// Vortex field is NOT irrotational
		VortexField field;
		VectorFieldAnalyzer analyzer(field);

		REQUIRE_FALSE(analyzer.IsIrrotational(Vec3Cart(-2, -2, -2), Vec3Cart(2, 2, 2), 5, REAL(0.01)));
	}

	TEST_CASE("VectorFieldAnalyzer::IsConservative", "[field_analyzer][vector][conservative]")
	{
			TEST_PRECISION_INFO();
		// Conservative field (gradient of potential)
		GradientField field;
		VectorFieldAnalyzer analyzer(field);

		REQUIRE(analyzer.IsConservative(Vec3Cart(-2, -2, -2), Vec3Cart(2, 2, 2), 5, REAL(0.01)));
	}

	TEST_CASE("VectorFieldAnalyzer::ConstantField_both_properties", "[field_analyzer][vector]")
	{
			TEST_PRECISION_INFO();
		// Constant field is BOTH solenoidal AND irrotational
		ConstantField field;
		VectorFieldAnalyzer analyzer(field);

		REQUIRE(analyzer.IsSolenoidal(Vec3Cart(-2, -2, -2), Vec3Cart(2, 2, 2), 5, REAL(0.01)));
		REQUIRE(analyzer.IsIrrotational(Vec3Cart(-2, -2, -2), Vec3Cart(2, 2, 2), 5, REAL(0.01)));
	}

	TEST_CASE("VectorFieldAnalyzer::MaxDivergence", "[field_analyzer][vector]")
	{
			TEST_PRECISION_INFO();
		RadialOutwardField field;
		VectorFieldAnalyzer analyzer(field);

		Real maxDiv = analyzer.MaxDivergence(Vec3Cart(-1, -1, -1), Vec3Cart(1, 1, 1), 5);
		// div(F) = 3 everywhere, so max should be ~3
		REQUIRE_THAT(maxDiv, WithinRel(REAL(3.0), REAL(0.05)));
	}

	TEST_CASE("VectorFieldAnalyzer::MaxCurlMagnitude", "[field_analyzer][vector]")
	{
			TEST_PRECISION_INFO();
		VortexField field;
		VectorFieldAnalyzer analyzer(field);

		Real maxCurl = analyzer.MaxCurlMagnitude(Vec3Cart(-1, -1, -1), Vec3Cart(1, 1, 1), 5);
		// curl = (0, 0, 2), magnitude = 2 everywhere
		REQUIRE_THAT(maxCurl, WithinRel(REAL(2.0), REAL(0.05)));
	}

	TEST_CASE("VectorFieldAnalyzer::CirculationAroundCircle_vortex", "[field_analyzer][vector][circulation]")
	{
			TEST_PRECISION_INFO();
		// Vortex field F = (-y, x, 0)
		// Circulation around circle of radius R in xy-plane = 2πR² (by Stokes: ∫curl·dA = 2 * πR²)
		VortexField field;
		VectorFieldAnalyzer analyzer(field);

		Real R = REAL(1.0);
		Real circulation = analyzer.CirculationAroundCircle(
			Vec3Cart(0, 0, 0), Vec3Cart(0, 0, 1), R, 200);

		Real expected = REAL(2.0) * Constants::PI * R * R;  // 2πR²
		REQUIRE_THAT(circulation, WithinRel(expected, REAL(0.02)));
	}

	TEST_CASE("VectorFieldAnalyzer::CirculationAroundCircle_conservative", "[field_analyzer][vector][circulation]")
	{
			TEST_PRECISION_INFO();
		// Conservative field should have zero circulation around any closed loop
		GradientField field;
		VectorFieldAnalyzer analyzer(field);

		Real circulation = analyzer.CirculationAroundCircle(
			Vec3Cart(0, 0, 0), Vec3Cart(0, 0, 1), REAL(1.0), 200);

		REQUIRE_THAT(circulation, WithinAbs(REAL(0.0), REAL(0.1)));
	}

} // namespace FieldAnalyzersTests
