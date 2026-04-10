#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector/VectorN.h"
#include "base/Vector/VectorTypes.h"
#include "base/Function.h"
#include "core/FieldOperations.h"
#include "core/MetricTensor.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

/**
 * COMPREHENSIVE TEST SUITE: Field Operations
 * 
 * Tests for differential operators on scalar and vector fields:
 * - Gradient: scalar field → vector field (∇f)
 * - Laplacian: scalar field → scalar field (∇²f)
 * - Divergence: vector field → scalar field (∇·F)
 * - Curl: vector field → vector field (∇×F)
 * 
 * Coordinate systems tested:
 * - Cartesian (arbitrary dimension for gradient/laplacian/divergence, 3D for curl)
 * - Spherical (r, θ, φ)
 * - Cylindrical (r, φ, z)
 * 
 * TEST STRATEGY:
 * - Use simple functions with known analytical derivatives
 * - Verify against analytically computed results
 * - Test fundamental vector calculus identities:
 *   - ∇×(∇f) = 0 (curl of gradient is zero)
 *   - ∇·(∇×F) = 0 (divergence of curl is zero)
 *   - For radial fields, verify spherical/cylindrical results
 */

namespace MML::Tests::Core::FieldOperationsTests
{
	///////////////////////////////////////////////////////////////////////////
	// TEST SCALAR FUNCTIONS WITH KNOWN ANALYTICAL DERIVATIVES
	///////////////////////////////////////////////////////////////////////////
	
	// f(x,y,z) = x² + y² + z²
	// ∇f = (2x, 2y, 2z)
	// ∇²f = 6
	class QuadraticScalarField : public IScalarFunction<3>
	{
	public:
		Real operator()(const VectorN<Real, 3>& pos) const override
		{
			return pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
		}
	};
	
	// f(x,y,z) = xyz
	// ∇f = (yz, xz, xy)
	// ∇²f = 0
	class ProductScalarField : public IScalarFunction<3>
	{
	public:
		Real operator()(const VectorN<Real, 3>& pos) const override
		{
			return pos[0] * pos[1] * pos[2];
		}
	};
	
	// f(x,y,z) = sin(x)*cos(y)*exp(z)
	// ∇f = (cos(x)*cos(y)*exp(z), -sin(x)*sin(y)*exp(z), sin(x)*cos(y)*exp(z))
	// ∇²f = -sin(x)*cos(y)*exp(z) - sin(x)*cos(y)*exp(z) + sin(x)*cos(y)*exp(z) = -sin(x)*cos(y)*exp(z)
	class TrigExpScalarField : public IScalarFunction<3>
	{
	public:
		Real operator()(const VectorN<Real, 3>& pos) const override
		{
			return sin(pos[0]) * cos(pos[1]) * exp(pos[2]);
		}
	};
	
	// 2D scalar field: f(x,y) = x² - y²
	// ∇f = (2x, -2y)
	// ∇²f = 2 - 2 = 0 (harmonic function!)
	class Harmonic2DField : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& pos) const override
		{
			return pos[0]*pos[0] - pos[1]*pos[1];
		}
	};
	
	///////////////////////////////////////////////////////////////////////////
	// TEST VECTOR FUNCTIONS WITH KNOWN ANALYTICAL PROPERTIES
	///////////////////////////////////////////////////////////////////////////
	
	// F(x,y,z) = (x, y, z) - radial outward field
	// ∇·F = 3
	// ∇×F = (0, 0, 0)
	class RadialVectorField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return pos;
		}
	};
	
	// F(x,y,z) = (y, -x, 0) - rotation about z-axis
	// ∇·F = 0 (incompressible)
	// ∇×F = (0, 0, -2)
	class RotationalVectorField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>{pos[1], -pos[0], REAL(0.0)};
		}
	};
	
	// F(x,y,z) = (yz, xz, xy) - this is ∇(xyz), so curl should be 0
	// ∇·F = z + z + y = 2z + y  Wait, let me recalculate:
	// ∂(yz)/∂x = 0, ∂(xz)/∂y = 0, ∂(xy)/∂z = 0
	// ∇·F = 0
	// ∇×F = (x-x, y-y, z-z) = (0, 0, 0)
	class GradientOfProductField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>{pos[1]*pos[2], pos[0]*pos[2], pos[0]*pos[1]};
		}
	};
	
	// F(x,y,z) = (x², y², z²)
	// ∇·F = 2x + 2y + 2z
	// ∇×F = (0, 0, 0) (each component depends only on its coordinate)
	class SquaredComponentsField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>{pos[0]*pos[0], pos[1]*pos[1], pos[2]*pos[2]};
		}
	};

	/*********************************************************************/
	/*****                 GRADIENT TESTS - CARTESIAN                *****/
	/*********************************************************************/

	TEST_CASE("Gradient_Cartesian_Quadratic_3D", "[field_operations][gradient][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		QuadraticScalarField f;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		auto grad = ScalarFieldOperations::GradientCart<3>(f, pos);
		
		// ∇f = (2x, 2y, 2z) = (2, 4, 6)
		REQUIRE_THAT(grad[0], WithinRel(REAL(2.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[1], WithinRel(REAL(4.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[2], WithinRel(REAL(6.0), TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("Gradient_Cartesian_Product_3D", "[field_operations][gradient][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		ProductScalarField f;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto grad = ScalarFieldOperations::GradientCart<3>(f, pos);
		
		// ∇f = (yz, xz, xy) = (12, 8, 6)
		REQUIRE_THAT(grad[0], WithinRel(REAL(12.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[1], WithinRel(REAL(8.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[2], WithinRel(REAL(6.0), TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("Gradient_Cartesian_2D", "[field_operations][gradient][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		Harmonic2DField f;
		VectorN<Real, 2> pos{REAL(3.0), REAL(4.0)};
		
		auto grad = ScalarFieldOperations::GradientCart<2>(f, pos);
		
		// ∇f = (2x, -2y) = (6, -8)
		REQUIRE_THAT(grad[0], WithinRel(REAL(6.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[1], WithinRel(REAL(-8.0), TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("Gradient_Cartesian_WithDerOrder", "[field_operations][gradient][cartesian][accuracy]")
	{
		TEST_PRECISION_INFO();
		
		TrigExpScalarField f;
		VectorN<Real, 3> pos{REAL(0.5), REAL(0.3), REAL(0.2)};
		
		// Analytical gradient
		Real x = pos[0], y = pos[1], z = pos[2];
		Real expected_x = cos(x) * cos(y) * exp(z);
		Real expected_y = -sin(x) * sin(y) * exp(z);
		Real expected_z = sin(x) * cos(y) * exp(z);
		
		// Test different derivative orders - higher should be more accurate
		auto grad4 = ScalarFieldOperations::GradientCart<3>(f, pos, 4);
		auto grad8 = ScalarFieldOperations::GradientCart<3>(f, pos, 8);
		
		// Both should be reasonably accurate
		REQUIRE_THAT(grad4[0], WithinRel(expected_x, REAL(1e-6)));
		REQUIRE_THAT(grad8[0], WithinRel(expected_x, TOL(1e-9, 1e-4)));
	}

	/*********************************************************************/
	/*****                 LAPLACIAN TESTS - CARTESIAN               *****/
	/*********************************************************************/

	TEST_CASE("Laplacian_Cartesian_Quadratic_3D", "[field_operations][laplacian][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		QuadraticScalarField f;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		Real lapl = ScalarFieldOperations::LaplacianCart<3>(f, pos);
		
		// ∇²f = 2 + 2 + 2 = 6 (constant everywhere)
		REQUIRE_THAT(lapl, WithinRel(REAL(6.0), TOL(1e-6, 5e-4)));
	}
	
	TEST_CASE("Laplacian_Cartesian_Harmonic_2D", "[field_operations][laplacian][cartesian][harmonic]")
	{
		TEST_PRECISION_INFO();
		
		Harmonic2DField f;
		VectorN<Real, 2> pos{REAL(3.0), REAL(4.0)};
		
		Real lapl = ScalarFieldOperations::LaplacianCart<2>(f, pos);
		
		// ∇²f = 2 - 2 = 0 (harmonic function!)
		// Note: numerical differentiation has ~TOL(1e-8, 1e-4) accuracy
		REQUIRE_THAT(lapl, WithinAbs(REAL(0.0), TOL(1e-7, 5e-3)));
	}
	
	TEST_CASE("Laplacian_Cartesian_Product_3D", "[field_operations][laplacian][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		ProductScalarField f;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Real lapl = ScalarFieldOperations::LaplacianCart<3>(f, pos);
		
		// f = xyz → ∂²f/∂x² = 0, ∂²f/∂y² = 0, ∂²f/∂z² = 0
		// ∇²f = 0
		// Note: numerical differentiation has ~TOL(1e-8, 1e-4) accuracy
		REQUIRE_THAT(lapl, WithinAbs(REAL(0.0), TOL(1e-8, 5e-4)));
	}
	
	TEST_CASE("Laplacian_Cartesian_TrigExp_3D", "[field_operations][laplacian][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		TrigExpScalarField f;
		VectorN<Real, 3> pos{REAL(0.5), REAL(0.3), REAL(0.2)};
		
		Real lapl = ScalarFieldOperations::LaplacianCart<3>(f, pos);
		
		// f = sin(x)*cos(y)*exp(z)
		// ∂²f/∂x² = -sin(x)*cos(y)*exp(z)
		// ∂²f/∂y² = -sin(x)*cos(y)*exp(z)
		// ∂²f/∂z² = sin(x)*cos(y)*exp(z)
		// ∇²f = -sin(x)*cos(y)*exp(z) - sin(x)*cos(y)*exp(z) + sin(x)*cos(y)*exp(z)
		//     = -sin(x)*cos(y)*exp(z)
		Real expected = -sin(pos[0]) * cos(pos[1]) * exp(pos[2]);
		
		REQUIRE_THAT(lapl, WithinRel(expected, TOL(1e-6, 5e-4)));
	}

	/*********************************************************************/
	/*****                DIVERGENCE TESTS - CARTESIAN               *****/
	/*********************************************************************/

	TEST_CASE("Divergence_Cartesian_Radial", "[field_operations][divergence][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		RadialVectorField F;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		Real div = VectorFieldOperations::DivCart<3>(F, pos);
		
		// F = (x, y, z) → ∇·F = 1 + 1 + 1 = 3
		REQUIRE_THAT(div, WithinRel(REAL(3.0), TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("Divergence_Cartesian_Rotational_IsZero", "[field_operations][divergence][cartesian][incompressible]")
	{
		TEST_PRECISION_INFO();
		
		RotationalVectorField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(1.0)};
		
		Real div = VectorFieldOperations::DivCart<3>(F, pos);
		
		// F = (y, -x, 0) → ∇·F = 0 + 0 + 0 = 0 (incompressible)
		REQUIRE_THAT(div, WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("Divergence_Cartesian_GradientField", "[field_operations][divergence][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		GradientOfProductField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Real div = VectorFieldOperations::DivCart<3>(F, pos);
		
		// F = (yz, xz, xy) → ∇·F = 0 + 0 + 0 = 0
		REQUIRE_THAT(div, WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("Divergence_Cartesian_SquaredComponents", "[field_operations][divergence][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		SquaredComponentsField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Real div = VectorFieldOperations::DivCart<3>(F, pos);
		
		// F = (x², y², z²) → ∇·F = 2x + 2y + 2z = 4 + 6 + 8 = 18
		REQUIRE_THAT(div, WithinRel(REAL(18.0), REAL(1e-6)));
	}

	/*********************************************************************/
	/*****                   CURL TESTS - CARTESIAN                  *****/
	/*********************************************************************/

	TEST_CASE("Curl_Cartesian_Rotational", "[field_operations][curl][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		RotationalVectorField F;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		auto curl = VectorFieldOperations::CurlCart(F, pos);
		
		// F = (y, -x, 0)
		// ∇×F = (∂0/∂y - ∂(-x)/∂z, ∂y/∂z - ∂0/∂x, ∂(-x)/∂x - ∂y/∂y)
		//     = (0 - 0, 0 - 0, -1 - 1) = (0, 0, -2)
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(curl[2], WithinRel(REAL(-2.0), TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("Curl_Cartesian_Irrotational_Radial", "[field_operations][curl][cartesian][irrotational]")
	{
		TEST_PRECISION_INFO();
		
		RadialVectorField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto curl = VectorFieldOperations::CurlCart(F, pos);
		
		// F = (x, y, z) → ∇×F = (0, 0, 0) (conservative field)
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("Curl_Cartesian_GradientOfScalar_IsZero", "[field_operations][curl][cartesian][identity]")
	{
		TEST_PRECISION_INFO();
		
		// F = ∇(xyz) = (yz, xz, xy)
		GradientOfProductField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto curl = VectorFieldOperations::CurlCart(F, pos);
		
		// ∇×(∇f) = 0 (fundamental vector identity)
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("Curl_Cartesian_SquaredComponents_IsZero", "[field_operations][curl][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		SquaredComponentsField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto curl = VectorFieldOperations::CurlCart(F, pos);
		
		// F = (x², y², z²) - each component depends only on its own coordinate
		// ∇×F = (0, 0, 0)
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****              VECTOR CALCULUS IDENTITIES                   *****/
	/*********************************************************************/

	TEST_CASE("Identity_CurlOfGradient_IsZero", "[field_operations][identity][fundamental]")
	{
		TEST_PRECISION_INFO();
		
		// For any scalar field f, ∇×(∇f) = 0
		TrigExpScalarField f;
		VectorN<Real, 3> pos{REAL(0.5), REAL(0.7), REAL(0.3)};
		
		// Compute gradient
		auto grad = ScalarFieldOperations::GradientCart<3>(f, pos);
		
		// Create a vector field from the gradient (wrapping in lambda)
		class GradientWrapper : public IVectorFunction<3>
		{
			const IScalarFunction<3>& _f;
		public:
			GradientWrapper(const IScalarFunction<3>& f) : _f(f) {}
			VectorN<Real, 3> operator()(const VectorN<Real, 3>& p) const override
			{
				return ScalarFieldOperations::GradientCart<3>(_f, p);
			}
		};
		
		GradientWrapper gradField(f);
		auto curl_of_grad = VectorFieldOperations::CurlCart(gradField, pos);
		
		// ∇×(∇f) should be zero
		Real magnitude = sqrt(curl_of_grad[0]*curl_of_grad[0] + 
		                      curl_of_grad[1]*curl_of_grad[1] + 
		                      curl_of_grad[2]*curl_of_grad[2]);
		REQUIRE_THAT(magnitude, WithinAbs(REAL(0.0), REAL(1e-6)));
	}
	
	TEST_CASE("Identity_DivergenceOfCurl_IsZero", "[field_operations][identity][fundamental]")
	{
		TEST_PRECISION_INFO();
		
		// For any vector field F, ∇·(∇×F) = 0
		RotationalVectorField F;
		VectorN<Real, 3> pos{REAL(1.5), REAL(2.5), REAL(3.5)};
		
		// Create a vector field from curl of F
		class CurlWrapper : public IVectorFunction<3>
		{
			const IVectorFunction<3>& _F;
		public:
			CurlWrapper(const IVectorFunction<3>& F) : _F(F) {}
			VectorN<Real, 3> operator()(const VectorN<Real, 3>& p) const override
			{
				return VectorFieldOperations::CurlCart(_F, p);
			}
		};
		
		CurlWrapper curlField(F);
		Real div_of_curl = VectorFieldOperations::DivCart<3>(curlField, pos);
		
		// ∇·(∇×F) should be zero
		REQUIRE_THAT(div_of_curl, WithinAbs(REAL(0.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****              SPHERICAL COORDINATE TESTS                   *****/
	/*********************************************************************/

	// f(r,θ,φ) = r² (depends only on r)
	// In spherical: ∇f = (2r, 0, 0)
	// ∇²f = (1/r²)∂/∂r(r²·2r) = (1/r²)∂/∂r(2r³) = (1/r²)·6r² = 6
	class SphericalRadialField : public IScalarFunction<3>
	{
	public:
		Real operator()(const VectorN<Real, 3>& pos) const override
		{
			return pos[0] * pos[0];  // r²
		}
	};

	TEST_CASE("Gradient_Spherical_RadialField", "[field_operations][gradient][spherical]")
	{
		TEST_PRECISION_INFO();
		
		SphericalRadialField f;
		Vec3Sph pos{REAL(2.0), REAL(0.5), REAL(1.0)};  // (r, θ, φ)
		
		auto grad = ScalarFieldOperations::GradientSpher(f, pos);
		
		// f = r² → ∇f = (∂f/∂r, (1/r)∂f/∂θ, (1/r·sinθ)∂f/∂φ) = (2r, 0, 0)
		REQUIRE_THAT(grad[0], WithinRel(REAL(4.0), REAL(1e-6)));  // 2r = 4
		REQUIRE_THAT(grad[1], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[2], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("Laplacian_Spherical_RadialField", "[field_operations][laplacian][spherical]")
	{
		TEST_PRECISION_INFO();
		
		SphericalRadialField f;
		Vec3Sph pos{REAL(2.0), REAL(0.8), REAL(1.2)};  // (r, θ, φ) - avoid poles
		
		Real lapl = ScalarFieldOperations::LaplacianSpher(f, pos);
		
		// f = r² → ∇²f = 6 (same as Cartesian x²+y²+z² since r² = x²+y²+z²)
		REQUIRE_THAT(lapl, WithinRel(REAL(6.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****             CYLINDRICAL COORDINATE TESTS                  *****/
	/*********************************************************************/

	// f(r,φ,z) = r² (depends only on cylindrical r)
	// ∇f = (2r, 0, 0)
	// ∇²f = (1/r)∂/∂r(r·2r) = (1/r)∂/∂r(2r²) = (1/r)·4r = 4
	class CylindricalRadialField : public IScalarFunction<3>
	{
	public:
		Real operator()(const VectorN<Real, 3>& pos) const override
		{
			return pos[0] * pos[0];  // r²
		}
	};

	TEST_CASE("Gradient_Cylindrical_RadialField", "[field_operations][gradient][cylindrical]")
	{
		TEST_PRECISION_INFO();
		
		CylindricalRadialField f;
		Vec3Cyl pos{REAL(3.0), REAL(1.0), REAL(2.0)};  // (r, φ, z)
		
		auto grad = ScalarFieldOperations::GradientCyl(f, pos);
		
		// f = r² → ∇f = (2r, 0, 0)
		REQUIRE_THAT(grad[0], WithinRel(REAL(6.0), REAL(1e-6)));  // 2r = 6
		REQUIRE_THAT(grad[1], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[2], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("Laplacian_Cylindrical_RadialField", "[field_operations][laplacian][cylindrical]")
	{
		TEST_PRECISION_INFO();
		
		CylindricalRadialField f;
		Vec3Cyl pos{REAL(3.0), REAL(1.5), REAL(2.0)};  // (r, φ, z)
		
		Real lapl = ScalarFieldOperations::LaplacianCyl(f, pos);
		
		// f = r² → ∇²f = (1/r)(r·2) + 2 = 2 + 2 = 4
		REQUIRE_THAT(lapl, WithinRel(REAL(4.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****             SPHERICAL DIVERGENCE TESTS                    *****/
	/*********************************************************************/

	// Radial field F = (r, 0, 0) in spherical coordinates
	// ∇·F = (1/r²)∂(r²·r)/∂r = (1/r²)∂(r³)/∂r = (1/r²)·3r² = 3
	class SphericalRadialVectorField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			// F = (Fᵣ, Fθ, Fφ) = (r, 0, 0)
			return VectorN<Real, 3>{pos[0], REAL(0.0), REAL(0.0)};
		}
	};

	TEST_CASE("Divergence_Spherical_RadialField", "[field_operations][divergence][spherical]")
	{
		TEST_PRECISION_INFO();
		
		SphericalRadialVectorField F;
		Vec3Sph pos{REAL(2.0), REAL(0.8), REAL(1.2)};  // (r, θ, φ) - avoid poles
		
		Real div = VectorFieldOperations::DivSpher(F, pos);
		
		// F = (r, 0, 0) → ∇·F = 3 (constant)
		REQUIRE_THAT(div, WithinRel(REAL(3.0), REAL(1e-5)));
	}

	// Purely θ-dependent field: F = (0, sinθ, 0)
	// ∇·F = (1/r·sinθ)∂(sinθ·sinθ)/∂θ = (1/r·sinθ)∂(sin²θ)/∂θ
	//     = (1/r·sinθ)·2sinθ·cosθ = (2cosθ)/r
	class SphericalThetaVectorField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			// F = (Fᵣ, Fθ, Fφ) = (0, sinθ, 0)
			return VectorN<Real, 3>{REAL(0.0), REAL(sin(pos[1])), REAL(0.0)};
		}
	};

	TEST_CASE("Divergence_Spherical_ThetaField", "[field_operations][divergence][spherical]")
	{
		TEST_PRECISION_INFO();
		
		SphericalThetaVectorField F;
		Real r = REAL(2.0);
		Real theta = REAL(0.8);  // θ in radians
		Vec3Sph pos{r, theta, REAL(1.2)};
		
		Real div = VectorFieldOperations::DivSpher(F, pos);
		
		// F = (0, sinθ, 0) → ∇·F = 2cosθ/r
		Real expected = REAL(2.0) * cos(theta) / r;
		REQUIRE_THAT(div, WithinRel(expected, REAL(1e-5)));
	}

	/*********************************************************************/
	/*****            CYLINDRICAL DIVERGENCE TESTS                   *****/
	/*********************************************************************/

	// Radial field F = (r, 0, 0) in cylindrical coordinates
	// ∇·F = (1/r)∂(r·r)/∂r + 0 + 0 = (1/r)·2r = 2
	class CylindricalRadialVectorField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			// F = (Fᵣ, Fφ, Fz) = (r, 0, 0)
			return VectorN<Real, 3>{pos[0], REAL(0.0), REAL(0.0)};
		}
	};

	TEST_CASE("Divergence_Cylindrical_RadialField", "[field_operations][divergence][cylindrical]")
	{
		TEST_PRECISION_INFO();
		
		CylindricalRadialVectorField F;
		Vec3Cyl pos{REAL(3.0), REAL(1.5), REAL(2.0)};  // (r, φ, z)
		
		Real div = VectorFieldOperations::DivCyl(F, pos);
		
		// F = (r, 0, 0) → ∇·F = 2 (constant)
		REQUIRE_THAT(div, WithinRel(REAL(2.0), REAL(1e-5)));
	}

	// Combined field: F = (r, 0, z)
	// ∇·F = (1/r)∂(r²)/∂r + 0 + ∂z/∂z = 2 + 1 = 3
	class CylindricalCombinedVectorField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			// F = (Fᵣ, Fφ, Fz) = (r, 0, z)
			return VectorN<Real, 3>{pos[0], REAL(0.0), pos[2]};
		}
	};

	TEST_CASE("Divergence_Cylindrical_CombinedField", "[field_operations][divergence][cylindrical]")
	{
		TEST_PRECISION_INFO();
		
		CylindricalCombinedVectorField F;
		Vec3Cyl pos{REAL(2.5), REAL(0.7), REAL(4.0)};  // (r, φ, z)
		
		Real div = VectorFieldOperations::DivCyl(F, pos);
		
		// F = (r, 0, z) → ∇·F = 2 + 1 = 3
		REQUIRE_THAT(div, WithinRel(REAL(3.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****               SPHERICAL CURL TESTS                        *****/
	/*********************************************************************/

	// Radial field F = (r², 0, 0) in spherical coordinates
	// ∇×F = 0 (radial fields are irrotational)
	class SphericalIrrotationalField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			Real r = pos[0];
			return VectorN<Real, 3>{r*r, REAL(0.0), REAL(0.0)};
		}
	};

	TEST_CASE("Curl_Spherical_Irrotational", "[field_operations][curl][spherical][irrotational]")
	{
		TEST_PRECISION_INFO();
		
		SphericalIrrotationalField F;
		Vec3Sph pos{REAL(2.0), REAL(0.8), REAL(1.2)};  // (r, θ, φ) - avoid poles
		
		auto curl = VectorFieldOperations::CurlSpher(F, pos);
		
		// Radial field → ∇×F = 0
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), REAL(1e-6)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), REAL(1e-6)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), REAL(1e-6)));
	}

	// Azimuthal field F = (0, 0, r·sinθ) in spherical coordinates
	// This is like a "φ-rotation" field
	class SphericalAzimuthalField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			Real r = pos[0];
			Real theta = pos[1];
			// F = (0, 0, r·sinθ)
			return VectorN<Real, 3>{REAL(0.0), REAL(0.0), REAL(r * sin(theta))};
		}
	};

	TEST_CASE("Curl_Spherical_AzimuthalField", "[field_operations][curl][spherical]")
	{
		TEST_PRECISION_INFO();
		
		SphericalAzimuthalField F;
		Real r = REAL(2.0);
		Real theta = REAL(0.8);
		Vec3Sph pos{r, theta, REAL(1.2)};
		
		auto curl = VectorFieldOperations::CurlSpher(F, pos);
		
		// For F = (0, 0, r·sinθ):
		// (∇×F)ᵣ = (1/r·sinθ)[∂(sinθ·r·sinθ)/∂θ - 0] = (1/r·sinθ)[∂(r·sin²θ)/∂θ]
		//        = (1/r·sinθ)·2r·sinθ·cosθ = 2cosθ
		// (∇×F)θ = (1/r)[0 - ∂(r·r·sinθ)/∂r] = (1/r)[-2r·sinθ] = -2sinθ
		// (∇×F)φ = (1/r)[0 - 0] = 0
		Real expected_r = REAL(2.0) * cos(theta);
		Real expected_theta = REAL(-2.0) * sin(theta);
		
		REQUIRE_THAT(curl[0], WithinRel(expected_r, REAL(1e-4)));
		REQUIRE_THAT(curl[1], WithinRel(expected_theta, REAL(1e-4)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), REAL(1e-6)));
	}

	/*********************************************************************/
	/*****               CYLINDRICAL CURL TESTS                      *****/
	/*********************************************************************/

	// Radial field F = (r, 0, 0) in cylindrical coordinates
	// ∇×F = 0 (radial fields are irrotational)
	class CylindricalIrrotationalField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>{pos[0], REAL(0.0), REAL(0.0)};
		}
	};

	TEST_CASE("Curl_Cylindrical_Irrotational", "[field_operations][curl][cylindrical][irrotational]")
	{
		TEST_PRECISION_INFO();
		
		CylindricalIrrotationalField F;
		Vec3Cyl pos{REAL(3.0), REAL(1.5), REAL(2.0)};  // (r, φ, z)
		
		auto curl = VectorFieldOperations::CurlCyl(F, pos);
		
		// Radial field → ∇×F = 0
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
	}

	// Azimuthal field F = (0, r, 0) in cylindrical coordinates
	// This represents solid body rotation about the z-axis
	// (∇×F)z = (1/r)[Fφ + r·∂Fφ/∂r - ∂Fᵣ/∂φ] = (1/r)[r + r·1 - 0] = (1/r)·2r = 2
	class CylindricalRotationalField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			// F = (Fᵣ, Fφ, Fz) = (0, r, 0) - solid body rotation
			return VectorN<Real, 3>{REAL(0.0), pos[0], REAL(0.0)};
		}
	};

	TEST_CASE("Curl_Cylindrical_Rotational", "[field_operations][curl][cylindrical]")
	{
		TEST_PRECISION_INFO();
		
		CylindricalRotationalField F;
		Vec3Cyl pos{REAL(2.0), REAL(1.0), REAL(3.0)};  // (r, φ, z)
		
		auto curl = VectorFieldOperations::CurlCyl(F, pos);
		
		// F = (0, r, 0) → (∇×F)z = 2 (uniform rotation)
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(curl[2], WithinRel(REAL(2.0), REAL(1e-5)));
	}

	// Field with z-dependence: F = (0, z, 0) in cylindrical coordinates
	// (∇×F)ᵣ = (1/r)∂Fz/∂φ - ∂Fφ/∂z = 0 - 1 = -1
	// (∇×F)φ = ∂Fᵣ/∂z - ∂Fz/∂r = 0 - 0 = 0
	// (∇×F)z = (1/r)[Fφ + r·∂Fφ/∂r - ∂Fᵣ/∂φ] = (1/r)[z + 0 - 0] = z/r
	class CylindricalZDependentField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			// F = (Fᵣ, Fφ, Fz) = (0, z, 0)
			return VectorN<Real, 3>{REAL(0.0), pos[2], REAL(0.0)};
		}
	};

	TEST_CASE("Curl_Cylindrical_ZDependentField", "[field_operations][curl][cylindrical]")
	{
		TEST_PRECISION_INFO();
		
		CylindricalZDependentField F;
		Real r = REAL(3.0);
		Real z = REAL(6.0);
		Vec3Cyl pos{r, REAL(1.0), z};  // (r, φ, z)
		
		auto curl = VectorFieldOperations::CurlCyl(F, pos);
		
		// F = (0, z, 0) → ∇×F = (-1, 0, z/r)
		REQUIRE_THAT(curl[0], WithinRel(REAL(-1.0), REAL(1e-6)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(curl[2], WithinRel(z / r, REAL(1e-5)));
	}

	/*********************************************************************/
	/*****        GENERAL COORDINATE GRADIENT/DIVERGENCE TESTS       *****/
	/*********************************************************************/

	TEST_CASE("Gradient_General_Cartesian_Matches_CartesianGradient", "[field_operations][gradient][general][metric]")
	{
		TEST_PRECISION_INFO();
		
		QuadraticScalarField f;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		// Using Cartesian metric tensor (identity matrix)
		MetricTensorCartesian3D cartMetric;
		
		auto grad_general = ScalarFieldOperations::Gradient<3>(f, pos, cartMetric);
		auto grad_cart = ScalarFieldOperations::GradientCart<3>(f, pos);
		
		// Both should give same result for Cartesian coordinates
		REQUIRE_THAT(grad_general[0], WithinRel(grad_cart[0], TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad_general[1], WithinRel(grad_cart[1], TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad_general[2], WithinRel(grad_cart[2], TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Divergence_General_Cartesian_Matches_CartesianDivergence", "[field_operations][divergence][general][metric]")
	{
		TEST_PRECISION_INFO();
		
		RadialVectorField F;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		// Using Cartesian metric tensor (identity matrix)
		MetricTensorCartesian3D cartMetric;
		
		Real div_general = VectorFieldOperations::Divergence<3>(F, pos, cartMetric);
		Real div_cart = VectorFieldOperations::DivCart<3>(F, pos);
		
		// Both should give same result for Cartesian coordinates
		REQUIRE_THAT(div_general, WithinRel(div_cart, REAL(1e-6)));
	}

	/*********************************************************************/
	/*****            EDGE CASES AND ERROR HANDLING                  *****/
	/*********************************************************************/

	TEST_CASE("Gradient_Cartesian_AtOrigin", "[field_operations][gradient][cartesian][edge]")
	{
		TEST_PRECISION_INFO();
		
		QuadraticScalarField f;
		VectorN<Real, 3> origin{REAL(0.0), REAL(0.0), REAL(0.0)};
		
		auto grad = ScalarFieldOperations::GradientCart<3>(f, origin);
		
		// ∇f = (2x, 2y, 2z) at origin = (0, 0, 0)
		REQUIRE_THAT(grad[0], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(grad[1], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(grad[2], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Laplacian_Cartesian_AtOrigin", "[field_operations][laplacian][cartesian][edge]")
	{
		TEST_PRECISION_INFO();
		
		QuadraticScalarField f;
		VectorN<Real, 3> origin{REAL(0.0), REAL(0.0), REAL(0.0)};
		
		Real lapl = ScalarFieldOperations::LaplacianCart<3>(f, origin);
		
		// ∇²f = 6 everywhere (constant)
		REQUIRE_THAT(lapl, WithinRel(REAL(6.0), REAL(1e-6)));
	}

	TEST_CASE("Divergence_Cartesian_AtOrigin", "[field_operations][divergence][cartesian][edge]")
	{
		TEST_PRECISION_INFO();
		
		RadialVectorField F;
		VectorN<Real, 3> origin{REAL(0.0), REAL(0.0), REAL(0.0)};
		
		Real div = VectorFieldOperations::DivCart<3>(F, origin);
		
		// ∇·F = 3 everywhere (constant)
		REQUIRE_THAT(div, WithinRel(REAL(3.0), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Curl_Cartesian_AtOrigin", "[field_operations][curl][cartesian][edge]")
	{
		TEST_PRECISION_INFO();
		
		RotationalVectorField F;
		VectorN<Real, 3> origin{REAL(0.0), REAL(0.0), REAL(0.0)};
		
		auto curl = VectorFieldOperations::CurlCart(F, origin);
		
		// F = (y, -x, 0) → ∇×F = (0, 0, -2) everywhere
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(curl[2], WithinRel(REAL(-2.0), TOL(1e-8, 1e-4)));
	}

	/*********************************************************************/
	/*****         HIGHER DIMENSIONAL TESTS (N=4, N=5)               *****/
	/*********************************************************************/

	// 4D scalar field: f(x,y,z,w) = x² + y² + z² + w²
	// ∇f = (2x, 2y, 2z, 2w)
	// ∇²f = 8
	class Quadratic4DField : public IScalarFunction<4>
	{
	public:
		Real operator()(const VectorN<Real, 4>& pos) const override
		{
			return pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2] + pos[3]*pos[3];
		}
	};

	TEST_CASE("Gradient_Cartesian_4D", "[field_operations][gradient][cartesian][4d]")
	{
		TEST_PRECISION_INFO();
		
		Quadratic4DField f;
		VectorN<Real, 4> pos{REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto grad = ScalarFieldOperations::GradientCart<4>(f, pos);
		
		// ∇f = (2x, 2y, 2z, 2w) = (2, 4, 6, 8)
		REQUIRE_THAT(grad[0], WithinRel(REAL(2.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[1], WithinRel(REAL(4.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[2], WithinRel(REAL(6.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(grad[3], WithinRel(REAL(8.0), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Laplacian_Cartesian_4D", "[field_operations][laplacian][cartesian][4d]")
	{
		TEST_PRECISION_INFO();
		
		Quadratic4DField f;
		VectorN<Real, 4> pos{REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Real lapl = ScalarFieldOperations::LaplacianCart<4>(f, pos);
		
		// ∇²f = 2 + 2 + 2 + 2 = 8
		REQUIRE_THAT(lapl, WithinRel(REAL(8.0), TOL(1e-6, 5e-4)));
	}

	// 4D radial vector field: F(x,y,z,w) = (x, y, z, w)
	// ∇·F = 1 + 1 + 1 + 1 = 4
	class Radial4DVectorField : public IVectorFunction<4>
	{
	public:
		VectorN<Real, 4> operator()(const VectorN<Real, 4>& pos) const override
		{
			return pos;
		}
	};

	TEST_CASE("Divergence_Cartesian_4D", "[field_operations][divergence][cartesian][4d]")
	{
		TEST_PRECISION_INFO();
		
		Radial4DVectorField F;
		VectorN<Real, 4> pos{REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Real div = VectorFieldOperations::DivCart<4>(F, pos);
		
		// ∇·F = 4
		REQUIRE_THAT(div, WithinRel(REAL(4.0), TOL(1e-8, 1e-4)));
	}

	/*********************************************************************/
	/*****                    PLACEHOLDER TEST                       *****/
	/*********************************************************************/

	TEST_CASE("Test_Field_Operations_Info", "[field_operations][info]") 
	{
		TEST_PRECISION_INFO();
		INFO("Field Operations test suite complete");
		SUCCEED();
	}

	/*********************************************************************/
	/*****          DETAILED API TESTS - SCALAR FIELD OPERATIONS     *****/
	/*********************************************************************/

	TEST_CASE("GradientCartDetailed_MatchesRawMethod", "[field_operations][gradient][cartesian][detailed]")
	{
		TEST_PRECISION_INFO();

		QuadraticScalarField f;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};

		auto raw = ScalarFieldOperations::GradientCart<3>(f, pos);
		auto result = ScalarFieldOperations::GradientCartDetailed<3>(f, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.algorithm_name == "GradientCart");
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE_THAT(result.value[0], WithinRel(raw[0], TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.value[1], WithinRel(raw[1], TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.value[2], WithinRel(raw[2], TOL(1e-12, 1e-5)));
	}

	TEST_CASE("GradientCartDetailed_WithErrorEstimate", "[field_operations][gradient][cartesian][detailed]")
	{
		TEST_PRECISION_INFO();

		QuadraticScalarField f;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};

		FieldOperationConfig config;
		config.estimate_error = true;

		auto result = ScalarFieldOperations::GradientCartDetailed<3>(f, pos, config);

		REQUIRE(result.IsSuccess());
		// ∇f = (2x, 2y, 2z) = (2, 4, 6)
		REQUIRE_THAT(result.value[0], WithinRel(REAL(2.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.value[1], WithinRel(REAL(4.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.value[2], WithinRel(REAL(6.0), TOL(1e-8, 1e-4)));
		// Error should be small and non-negative
		REQUIRE(result.error[0] >= 0.0);
		REQUIRE(result.error[1] >= 0.0);
		REQUIRE(result.error[2] >= 0.0);
		REQUIRE(result.error[0] < TOL(1e-5, 1e-3));
	}

	TEST_CASE("GradientSpherDetailed_MatchesRawMethod", "[field_operations][gradient][spherical][detailed]")
	{
		TEST_PRECISION_INFO();

		SphericalRadialField f;
		Vec3Sph pos{REAL(2.0), REAL(0.8), REAL(1.2)};

		auto raw = ScalarFieldOperations::GradientSpher(f, pos);
		auto result = ScalarFieldOperations::GradientSpherDetailed(f, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "GradientSpher");
		REQUIRE_THAT(result.value[0], WithinRel(raw[0], TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[1], WithinAbs(raw[1], TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.value[2], WithinAbs(raw[2], TOL(1e-8, 1e-4)));
	}

	TEST_CASE("GradientCylDetailed_MatchesRawMethod", "[field_operations][gradient][cylindrical][detailed]")
	{
		TEST_PRECISION_INFO();

		CylindricalRadialField f;
		Vec3Cyl pos{REAL(3.0), REAL(1.0), REAL(2.0)};

		auto raw = ScalarFieldOperations::GradientCyl(f, pos);
		auto result = ScalarFieldOperations::GradientCylDetailed(f, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "GradientCyl");
		REQUIRE_THAT(result.value[0], WithinRel(raw[0], TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[1], WithinAbs(raw[1], TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.value[2], WithinAbs(raw[2], TOL(1e-8, 1e-4)));
	}

	TEST_CASE("LaplacianCartDetailed_MatchesRawMethod", "[field_operations][laplacian][cartesian][detailed]")
	{
		TEST_PRECISION_INFO();

		QuadraticScalarField f;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};

		auto raw = ScalarFieldOperations::LaplacianCart<3>(f, pos);
		auto result = ScalarFieldOperations::LaplacianCartDetailed<3>(f, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "LaplacianCart");
		REQUIRE_THAT(result.value, WithinRel(raw, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("LaplacianCartDetailed_WithErrorEstimate", "[field_operations][laplacian][cartesian][detailed]")
	{
		TEST_PRECISION_INFO();

		QuadraticScalarField f;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};

		FieldOperationConfig config;
		config.estimate_error = true;

		auto result = ScalarFieldOperations::LaplacianCartDetailed<3>(f, pos, config);

		REQUIRE(result.IsSuccess());
		// ∇²f = 6
		REQUIRE_THAT(result.value, WithinRel(REAL(6.0), TOL(1e-6, 5e-4)));
		// Error should be small (scalar for Laplacian)
		REQUIRE(result.error >= 0.0);
		REQUIRE(result.error < TOL(1e-3, 5e-3));
	}

	TEST_CASE("LaplacianSpherDetailed_MatchesRawMethod", "[field_operations][laplacian][spherical][detailed]")
	{
		TEST_PRECISION_INFO();

		SphericalRadialField f;
		Vec3Sph pos{REAL(2.0), REAL(0.8), REAL(1.2)};

		auto raw = ScalarFieldOperations::LaplacianSpher(f, pos);
		auto result = ScalarFieldOperations::LaplacianSpherDetailed(f, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "LaplacianSpher");
		REQUIRE_THAT(result.value, WithinRel(raw, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("LaplacianCylDetailed_MatchesRawMethod", "[field_operations][laplacian][cylindrical][detailed]")
	{
		TEST_PRECISION_INFO();

		CylindricalRadialField f;
		Vec3Cyl pos{REAL(3.0), REAL(1.5), REAL(2.0)};

		auto raw = ScalarFieldOperations::LaplacianCyl(f, pos);
		auto result = ScalarFieldOperations::LaplacianCylDetailed(f, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "LaplacianCyl");
		REQUIRE_THAT(result.value, WithinRel(raw, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****          DETAILED API TESTS - VECTOR FIELD OPERATIONS     *****/
	/*********************************************************************/

	TEST_CASE("DivCartDetailed_MatchesRawMethod", "[field_operations][divergence][cartesian][detailed]")
	{
		TEST_PRECISION_INFO();

		RadialVectorField F;
		VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};

		auto raw = VectorFieldOperations::DivCart<3>(F, pos);
		auto result = VectorFieldOperations::DivCartDetailed<3>(F, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "DivCart");
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE_THAT(result.value, WithinRel(raw, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("DivCartDetailed_WithErrorEstimate", "[field_operations][divergence][cartesian][detailed]")
	{
		TEST_PRECISION_INFO();

		SquaredComponentsField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};

		FieldOperationConfig config;
		config.estimate_error = true;

		auto result = VectorFieldOperations::DivCartDetailed<3>(F, pos, config);

		REQUIRE(result.IsSuccess());
		// F = (x², y², z²) → ∇·F = 2x + 2y + 2z = 18
		REQUIRE_THAT(result.value, WithinRel(REAL(18.0), REAL(1e-6)));
		REQUIRE(result.error >= 0.0);
		REQUIRE(result.error < REAL(1e-3));
	}

	TEST_CASE("DivSpherDetailed_MatchesRawMethod", "[field_operations][divergence][spherical][detailed]")
	{
		TEST_PRECISION_INFO();

		SphericalRadialVectorField F;
		Vec3Sph pos{REAL(2.0), REAL(0.8), REAL(1.2)};

		auto raw = VectorFieldOperations::DivSpher(F, pos);
		auto result = VectorFieldOperations::DivSpherDetailed(F, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "DivSpher");
		REQUIRE_THAT(result.value, WithinRel(raw, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("DivCylDetailed_MatchesRawMethod", "[field_operations][divergence][cylindrical][detailed]")
	{
		TEST_PRECISION_INFO();

		CylindricalRadialVectorField F;
		Vec3Cyl pos{REAL(3.0), REAL(1.5), REAL(2.0)};

		auto raw = VectorFieldOperations::DivCyl(F, pos);
		auto result = VectorFieldOperations::DivCylDetailed(F, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "DivCyl");
		REQUIRE_THAT(result.value, WithinRel(raw, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("CurlCartDetailed_MatchesRawMethod", "[field_operations][curl][cartesian][detailed]")
	{
		TEST_PRECISION_INFO();

		RotationalVectorField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(1.0)};

		auto raw = VectorFieldOperations::CurlCart(F, pos);
		auto result = VectorFieldOperations::CurlCartDetailed(F, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "CurlCart");
		REQUIRE_THAT(result.value[0], WithinAbs(raw[0], TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[1], WithinAbs(raw[1], TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[2], WithinRel(raw[2], TOL(1e-10, 1e-5)));
	}

	TEST_CASE("CurlCartDetailed_WithErrorEstimate", "[field_operations][curl][cartesian][detailed]")
	{
		TEST_PRECISION_INFO();

		RotationalVectorField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(1.0)};

		FieldOperationConfig config;
		config.estimate_error = true;

		auto result = VectorFieldOperations::CurlCartDetailed(F, pos, config);

		REQUIRE(result.IsSuccess());
		// F = (y, -x, 0) → ∇×F = (0, 0, -2)
		REQUIRE_THAT(result.value[0], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[1], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[2], WithinRel(REAL(-2.0), TOL(1e-8, 1e-4)));
		// Per-component error should be small
		REQUIRE(result.error[0] >= 0.0);
		REQUIRE(result.error[1] >= 0.0);
		REQUIRE(result.error[2] >= 0.0);
	}

	TEST_CASE("CurlSpherDetailed_MatchesRawMethod", "[field_operations][curl][spherical][detailed]")
	{
		TEST_PRECISION_INFO();

		SphericalAzimuthalField F;
		Vec3Sph pos{REAL(2.0), REAL(0.8), REAL(1.2)};

		auto raw = VectorFieldOperations::CurlSpher(F, pos);
		auto result = VectorFieldOperations::CurlSpherDetailed(F, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "CurlSpher");
		REQUIRE_THAT(result.value[0], WithinRel(raw[0], TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[1], WithinRel(raw[1], TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[2], WithinAbs(raw[2], TOL(1e-8, 1e-4)));
	}

	TEST_CASE("CurlCylDetailed_MatchesRawMethod", "[field_operations][curl][cylindrical][detailed]")
	{
		TEST_PRECISION_INFO();

		CylindricalRotationalField F;
		Vec3Cyl pos{REAL(2.0), REAL(1.0), REAL(3.0)};

		auto raw = VectorFieldOperations::CurlCyl(F, pos);
		auto result = VectorFieldOperations::CurlCylDetailed(F, pos);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "CurlCyl");
		REQUIRE_THAT(result.value[0], WithinAbs(raw[0], TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[1], WithinAbs(raw[1], TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value[2], WithinRel(raw[2], TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Detailed_ExceptionPolicy_NoThrow", "[field_operations][detailed][exception_policy]")
	{
		TEST_PRECISION_INFO();

		// A singular point for spherical divergence (r=0)
		SphericalRadialVectorField F;
		Vec3Sph singular_pos{REAL(0.0), REAL(0.8), REAL(1.2)};

		FieldOperationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

		auto result = VectorFieldOperations::DivSpherDetailed(F, singular_pos, config,
		                                                       Singularity::DEFAULT_POLICY);

		// Should not throw, but status may indicate failure
		// (the Singularity policy inside DivSpher may still throw,
		//  which ExecuteFieldDetailed catches and converts to status)
		// We just verify we don't crash
		REQUIRE((result.IsSuccess() || !result.error_message.empty()));
	}
}
