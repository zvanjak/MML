#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/VectorTypes.h"
#include "base/Function.h"
#include "core/FieldOperations.h"
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
		REQUIRE_THAT(grad[0], WithinRel(REAL(2.0), REAL(1e-8)));
		REQUIRE_THAT(grad[1], WithinRel(REAL(4.0), REAL(1e-8)));
		REQUIRE_THAT(grad[2], WithinRel(REAL(6.0), REAL(1e-8)));
	}
	
	TEST_CASE("Gradient_Cartesian_Product_3D", "[field_operations][gradient][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		ProductScalarField f;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto grad = ScalarFieldOperations::GradientCart<3>(f, pos);
		
		// ∇f = (yz, xz, xy) = (12, 8, 6)
		REQUIRE_THAT(grad[0], WithinRel(REAL(12.0), REAL(1e-8)));
		REQUIRE_THAT(grad[1], WithinRel(REAL(8.0), REAL(1e-8)));
		REQUIRE_THAT(grad[2], WithinRel(REAL(6.0), REAL(1e-8)));
	}
	
	TEST_CASE("Gradient_Cartesian_2D", "[field_operations][gradient][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		Harmonic2DField f;
		VectorN<Real, 2> pos{REAL(3.0), REAL(4.0)};
		
		auto grad = ScalarFieldOperations::GradientCart<2>(f, pos);
		
		// ∇f = (2x, -2y) = (6, -8)
		REQUIRE_THAT(grad[0], WithinRel(REAL(6.0), REAL(1e-8)));
		REQUIRE_THAT(grad[1], WithinRel(REAL(-8.0), REAL(1e-8)));
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
		REQUIRE_THAT(grad8[0], WithinRel(expected_x, REAL(1e-9)));
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
		REQUIRE_THAT(lapl, WithinRel(REAL(6.0), REAL(1e-6)));
	}
	
	TEST_CASE("Laplacian_Cartesian_Harmonic_2D", "[field_operations][laplacian][cartesian][harmonic]")
	{
		TEST_PRECISION_INFO();
		
		Harmonic2DField f;
		VectorN<Real, 2> pos{REAL(3.0), REAL(4.0)};
		
		Real lapl = ScalarFieldOperations::LaplacianCart<2>(f, pos);
		
		// ∇²f = 2 - 2 = 0 (harmonic function!)
		// Note: numerical differentiation has ~1e-8 accuracy
		REQUIRE_THAT(lapl, WithinAbs(REAL(0.0), REAL(1e-7)));
	}
	
	TEST_CASE("Laplacian_Cartesian_Product_3D", "[field_operations][laplacian][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		ProductScalarField f;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Real lapl = ScalarFieldOperations::LaplacianCart<3>(f, pos);
		
		// f = xyz → ∂²f/∂x² = 0, ∂²f/∂y² = 0, ∂²f/∂z² = 0
		// ∇²f = 0
		// Note: numerical differentiation has ~1e-8 accuracy
		REQUIRE_THAT(lapl, WithinAbs(REAL(0.0), REAL(1e-8)));
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
		
		REQUIRE_THAT(lapl, WithinRel(expected, REAL(1e-6)));
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
		REQUIRE_THAT(div, WithinRel(REAL(3.0), REAL(1e-8)));
	}
	
	TEST_CASE("Divergence_Cartesian_Rotational_IsZero", "[field_operations][divergence][cartesian][incompressible]")
	{
		TEST_PRECISION_INFO();
		
		RotationalVectorField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(1.0)};
		
		Real div = VectorFieldOperations::DivCart<3>(F, pos);
		
		// F = (y, -x, 0) → ∇·F = 0 + 0 + 0 = 0 (incompressible)
		REQUIRE_THAT(div, WithinAbs(REAL(0.0), REAL(1e-10)));
	}
	
	TEST_CASE("Divergence_Cartesian_GradientField", "[field_operations][divergence][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		GradientOfProductField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Real div = VectorFieldOperations::DivCart<3>(F, pos);
		
		// F = (yz, xz, xy) → ∇·F = 0 + 0 + 0 = 0
		REQUIRE_THAT(div, WithinAbs(REAL(0.0), REAL(1e-10)));
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
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(curl[2], WithinRel(REAL(-2.0), REAL(1e-8)));
	}
	
	TEST_CASE("Curl_Cartesian_Irrotational_Radial", "[field_operations][curl][cartesian][irrotational]")
	{
		TEST_PRECISION_INFO();
		
		RadialVectorField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto curl = VectorFieldOperations::CurlCart(F, pos);
		
		// F = (x, y, z) → ∇×F = (0, 0, 0) (conservative field)
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), REAL(1e-10)));
	}
	
	TEST_CASE("Curl_Cartesian_GradientOfScalar_IsZero", "[field_operations][curl][cartesian][identity]")
	{
		TEST_PRECISION_INFO();
		
		// F = ∇(xyz) = (yz, xz, xy)
		GradientOfProductField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto curl = VectorFieldOperations::CurlCart(F, pos);
		
		// ∇×(∇f) = 0 (fundamental vector identity)
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), REAL(1e-8)));
	}
	
	TEST_CASE("Curl_Cartesian_SquaredComponents_IsZero", "[field_operations][curl][cartesian]")
	{
		TEST_PRECISION_INFO();
		
		SquaredComponentsField F;
		VectorN<Real, 3> pos{REAL(2.0), REAL(3.0), REAL(4.0)};
		
		auto curl = VectorFieldOperations::CurlCart(F, pos);
		
		// F = (x², y², z²) - each component depends only on its own coordinate
		// ∇×F = (0, 0, 0)
		REQUIRE_THAT(curl[0], WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(curl[1], WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(curl[2], WithinAbs(REAL(0.0), REAL(1e-10)));
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
		REQUIRE_THAT(grad[1], WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(grad[2], WithinAbs(REAL(0.0), REAL(1e-8)));
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
		REQUIRE_THAT(grad[1], WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(grad[2], WithinAbs(REAL(0.0), REAL(1e-8)));
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
	/*****                    PLACEHOLDER TEST                       *****/
	/*********************************************************************/

	TEST_CASE("Test_Field_Operations_Info", "[field_operations][info]") 
	{
		TEST_PRECISION_INFO();
		INFO("Field Operations test suite complete");
		SUCCEED();
	}
}
