///////////////////////////////////////////////////////////////////////////////////////////
// Moment of Inertia Tests - Rigid Body Dynamics
///////////////////////////////////////////////////////////////////////////////////////////
//
// Tests for moment of inertia calculations:
// - Discrete mass moment of inertia tensor
// - Continuous mass (numerical integration)
// - Known analytical results for standard shapes
// - Parallel axis theorem
// - Principal axes
//
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#include "base/Tensor.h"
#endif

#include "RigidBody/MomentOfInertiaCalculator.h"

using namespace MML;
using namespace MPL;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;

namespace MPL::Tests::RigidBody::MomentOfInertiaTests
{
	//===================================================================================
	// Single point mass
	//===================================================================================

	TEST_CASE("MomentOfInertia::Single_point_mass_at_origin", "[rigid_body][inertia]")
	{
		TEST_PRECISION_INFO();

		// Point mass at origin has zero moment of inertia about all axes
		DiscreteMass m(Vector3Cartesian(0.0, 0.0, 0.0), 5.0);
		DiscreteMassesConfig config({m});
		DiscreteMassMomentOfInertiaTensorCalculator calc(config);

		auto I = calc.calculate();

		REQUIRE_THAT(I(0, 0), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(1, 1), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(2, 2), WithinAbs(0.0, 1e-12));
	}

	TEST_CASE("MomentOfInertia::Single_point_mass_on_x_axis", "[rigid_body][inertia]")
	{
		TEST_PRECISION_INFO();

		// Point mass m at (r, 0, 0)
		// I_xx = m(y² + z²) = 0
		// I_yy = m(x² + z²) = mr²
		// I_zz = m(x² + y²) = mr²
		Real m = 2.0;
		Real r = 3.0;

		DiscreteMass mass(Vector3Cartesian(r, 0.0, 0.0), m);
		DiscreteMassesConfig config({mass});
		DiscreteMassMomentOfInertiaTensorCalculator calc(config);

		auto I = calc.calculate();

		REQUIRE_THAT(I(0, 0), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(1, 1), WithinAbs(m * r * r, 1e-12));
		REQUIRE_THAT(I(2, 2), WithinAbs(m * r * r, 1e-12));

		// Off-diagonal should be zero
		REQUIRE_THAT(I(0, 1), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(0, 2), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(1, 2), WithinAbs(0.0, 1e-12));
	}

	//===================================================================================
	// Two point masses (dumbbell)
	//===================================================================================

	TEST_CASE("MomentOfInertia::Dumbbell_symmetric", "[rigid_body][inertia]")
	{
		TEST_PRECISION_INFO();

		// Two equal masses at (±L, 0, 0)
		Real m = 1.0;
		Real L = 2.0;

		DiscreteMass m1(Vector3Cartesian(L, 0.0, 0.0), m);
		DiscreteMass m2(Vector3Cartesian(-L, 0.0, 0.0), m);
		DiscreteMassesConfig config({m1, m2});
		DiscreteMassMomentOfInertiaTensorCalculator calc(config);

		auto I = calc.calculate();

		// I_xx = 0 (both masses on x-axis)
		// I_yy = I_zz = 2mL²
		REQUIRE_THAT(I(0, 0), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(1, 1), WithinAbs(2.0 * m * L * L, 1e-12));
		REQUIRE_THAT(I(2, 2), WithinAbs(2.0 * m * L * L, 1e-12));

		// All products of inertia should be zero (symmetric about all planes)
		REQUIRE_THAT(I(0, 1), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(0, 2), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(1, 2), WithinAbs(0.0, 1e-12));
	}

	//===================================================================================
	// Four point masses (square)
	//===================================================================================

	TEST_CASE("MomentOfInertia::Square_arrangement", "[rigid_body][inertia]")
	{
		TEST_PRECISION_INFO();

		// Four equal masses at corners of square in xy-plane
		// (±a, ±a, 0)
		Real m = 1.0;
		Real a = 1.0;

		std::vector<DiscreteMass> masses = {
			DiscreteMass(Vector3Cartesian( a,  a, 0.0), m),
			DiscreteMass(Vector3Cartesian(-a,  a, 0.0), m),
			DiscreteMass(Vector3Cartesian(-a, -a, 0.0), m),
			DiscreteMass(Vector3Cartesian( a, -a, 0.0), m),
		};
		DiscreteMassesConfig config(masses);
		DiscreteMassMomentOfInertiaTensorCalculator calc(config);

		auto I = calc.calculate();

		// I_xx = 4m * a² (each mass contributes y²)
		// I_yy = 4m * a² (each mass contributes x²)
		// I_zz = 4m * 2a² = 8ma² (each contributes x² + y²)
		REQUIRE_THAT(I(0, 0), WithinAbs(4.0 * m * a * a, 1e-12));
		REQUIRE_THAT(I(1, 1), WithinAbs(4.0 * m * a * a, 1e-12));
		REQUIRE_THAT(I(2, 2), WithinAbs(8.0 * m * a * a, 1e-12));

		// Products of inertia: each pair cancels out due to symmetry
		REQUIRE_THAT(I(0, 1), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(0, 2), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(1, 2), WithinAbs(0.0, 1e-12));
	}

	//===================================================================================
	// Asymmetric configuration (non-zero products of inertia)
	//===================================================================================

	TEST_CASE("MomentOfInertia::Asymmetric_nonzero_products", "[rigid_body][inertia]")
	{
		TEST_PRECISION_INFO();

		// Single mass at (1, 1, 0) - should have non-zero I_xy
		Real m = 2.0;
		Real x = 1.0, y = 1.0, z = 0.0;

		DiscreteMass mass(Vector3Cartesian(x, y, z), m);
		DiscreteMassesConfig config({mass});
		DiscreteMassMomentOfInertiaTensorCalculator calc(config);

		auto I = calc.calculate();

		// I_xx = m(y² + z²) = m
		// I_yy = m(x² + z²) = m
		// I_zz = m(x² + y²) = 2m
		REQUIRE_THAT(I(0, 0), WithinAbs(m * 1.0, 1e-12));
		REQUIRE_THAT(I(1, 1), WithinAbs(m * 1.0, 1e-12));
		REQUIRE_THAT(I(2, 2), WithinAbs(m * 2.0, 1e-12));

		// I_xy = -mxy = -m
		REQUIRE_THAT(I(0, 1), WithinAbs(-m * x * y, 1e-12));
		REQUIRE_THAT(I(1, 0), WithinAbs(-m * x * y, 1e-12));

		// I_xz = I_yz = 0 (z = 0)
		REQUIRE_THAT(I(0, 2), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(I(1, 2), WithinAbs(0.0, 1e-12));
	}

	//===================================================================================
	// Tensor symmetry
	//===================================================================================

	TEST_CASE("MomentOfInertia::Tensor_is_symmetric", "[rigid_body][inertia]")
	{
		TEST_PRECISION_INFO();

		// Random-ish configuration
		std::vector<DiscreteMass> masses = {
			DiscreteMass(Vector3Cartesian(1.0, 2.0, 3.0), 1.0),
			DiscreteMass(Vector3Cartesian(-1.0, 0.5, -2.0), 2.0),
			DiscreteMass(Vector3Cartesian(0.5, -1.0, 1.5), 1.5),
		};
		DiscreteMassesConfig config(masses);
		DiscreteMassMomentOfInertiaTensorCalculator calc(config);

		auto I = calc.calculate();

		// Check symmetry: I_ij = I_ji
		REQUIRE_THAT(I(0, 1), WithinAbs(I(1, 0), 1e-12));
		REQUIRE_THAT(I(0, 2), WithinAbs(I(2, 0), 1e-12));
		REQUIRE_THAT(I(1, 2), WithinAbs(I(2, 1), 1e-12));
	}

	//===================================================================================
	// Trace relationship: I_xx + I_yy + I_zz = 2 * Σ m_i * r_i²
	//===================================================================================

	TEST_CASE("MomentOfInertia::Trace_relationship", "[rigid_body][inertia]")
	{
		TEST_PRECISION_INFO();

		std::vector<DiscreteMass> masses = {
			DiscreteMass(Vector3Cartesian(1.0, 0.0, 0.0), 1.0),
			DiscreteMass(Vector3Cartesian(0.0, 2.0, 0.0), 2.0),
			DiscreteMass(Vector3Cartesian(0.0, 0.0, 3.0), 3.0),
		};
		DiscreteMassesConfig config(masses);
		DiscreteMassMomentOfInertiaTensorCalculator calc(config);

		auto I = calc.calculate();

		// Trace = I_xx + I_yy + I_zz
		Real trace = I(0, 0) + I(1, 1) + I(2, 2);

		// Should equal 2 * Σ m_i * |r_i|²
		Real sum_mr2 = 0.0;
		for (const auto& m : masses) {
			Real r = m._position.NormL2();
			sum_mr2 += m._mass * r * r;
		}

		REQUIRE_THAT(trace, WithinAbs(2.0 * sum_mr2, 1e-10));
	}

} // namespace MPL::Tests::RigidBody::MomentOfInertiaTests
