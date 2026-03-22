#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Tensor.h"

#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/MetricTensor.h"

#include "mml/base/Geometry/Geometry3D.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Core::MetricTensorTests
{
	/********************************************************************************************************************/
	/********                           CARTESIAN METRIC TESTS                                                   ********/
	/********************************************************************************************************************/
	
	TEST_CASE("MetricTensorCartesian3D - Identity metric", "[MetricTensor][Cartesian]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		
		Vector3Cartesian pos(REAL(1.0), REAL(2.0), REAL(3.0));
		
		SECTION("Covariant metric is identity")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					Real expected = (i == j) ? REAL(1.0) : REAL(0.0);
					REQUIRE_THAT(g(i, j), WithinAbs(expected, REAL(1e-12)));
				}
			}
		}
		
		SECTION("Contravariant metric is identity")
		{
			auto g_inv = metric.GetContravariantMetric(pos);
			
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					Real expected = (i == j) ? REAL(1.0) : REAL(0.0);
					REQUIRE_THAT(g_inv(i, j), WithinAbs(expected, REAL(1e-10)));
				}
			}
		}
		
		SECTION("All Christoffel symbols vanish")
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						Real gamma = metric.GetChristoffelSymbolSecondKind(i, j, k, pos);
						INFO("Gamma^" << i << "_" << j << k);
						REQUIRE_THAT(gamma, WithinAbs(REAL(0.0), REAL(1e-8)));
					}
				}
			}
		}
	}
	
	/********************************************************************************************************************/
	/********                           SPHERICAL METRIC TESTS                                                   ********/
	/********************************************************************************************************************/
	
	TEST_CASE("MetricTensorSpherical - Metric components", "[MetricTensor][Spherical]")
	{
		TEST_PRECISION_INFO();
		MetricTensorSpherical metric;
		
		// Position: r=2, theta=pi/4, phi=pi/3
		Real r = REAL(2.0);
		Real theta = Constants::PI / REAL(4.0);
		Real phi = Constants::PI / REAL(3.0);
		VectorN<Real, 3> pos({r, theta, phi});
		
		SECTION("Diagonal covariant components")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			// g_rr = 1
			REQUIRE_THAT(g(0, 0), WithinAbs(REAL(1.0), REAL(1e-12)));
			
			// g_theta,theta = r^2 = 4
			REQUIRE_THAT(g(1, 1), WithinAbs(r * r, REAL(1e-12)));
			
			// g_phi,phi = r^2 * sin^2(theta)
			Real expected_gphi = r * r * std::sin(theta) * std::sin(theta);
			REQUIRE_THAT(g(2, 2), WithinAbs(expected_gphi, REAL(1e-12)));
		}
		
		SECTION("Off-diagonal components are zero")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			REQUIRE_THAT(g(0, 1), WithinAbs(REAL(0.0), REAL(1e-12)));
			REQUIRE_THAT(g(0, 2), WithinAbs(REAL(0.0), REAL(1e-12)));
			REQUIRE_THAT(g(1, 0), WithinAbs(REAL(0.0), REAL(1e-12)));
			REQUIRE_THAT(g(1, 2), WithinAbs(REAL(0.0), REAL(1e-12)));
			REQUIRE_THAT(g(2, 0), WithinAbs(REAL(0.0), REAL(1e-12)));
			REQUIRE_THAT(g(2, 1), WithinAbs(REAL(0.0), REAL(1e-12)));
		}
		
		SECTION("Contravariant metric is inverse of covariant")
		{
			auto g = metric.GetCovariantMetric(pos);
			auto g_inv = metric.GetContravariantMetric(pos);
			
			// g^rr = 1
			REQUIRE_THAT(g_inv(0, 0), WithinAbs(REAL(1.0), REAL(1e-10)));
			
			// g^theta,theta = 1/r^2
			REQUIRE_THAT(g_inv(1, 1), WithinAbs(REAL(1.0) / (r * r), REAL(1e-10)));
			
			// g^phi,phi = 1/(r^2 * sin^2(theta))
			Real sin_theta = std::sin(theta);
			REQUIRE_THAT(g_inv(2, 2), WithinAbs(REAL(1.0) / (r * r * sin_theta * sin_theta), REAL(1e-10)));
		}
	}
	
	TEST_CASE("MetricTensorSpherical - Christoffel symbols", "[MetricTensor][Spherical][Christoffel]")
	{
		TEST_PRECISION_INFO();
		MetricTensorSpherical metric;
		
		Real r = REAL(2.0);
		Real theta = Constants::PI / REAL(4.0);
		Real phi = Constants::PI / REAL(3.0);
		VectorN<Real, 3> pos({r, theta, phi});
		
		// Analytical Christoffel symbols for spherical coordinates:
		// Non-zero symbols:
		// Γ^r_θθ = -r
		// Γ^r_φφ = -r sin²θ
		// Γ^θ_rθ = Γ^θ_θr = 1/r
		// Γ^θ_φφ = -sinθ cosθ
		// Γ^φ_rφ = Γ^φ_φr = 1/r
		// Γ^φ_θφ = Γ^φ_φθ = cotθ
		
		SECTION("Γ^r_θθ = -r")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(0, 1, 1, pos);
			REQUIRE_THAT(gamma, WithinAbs(-r, REAL(1e-6)));
		}
		
		SECTION("Γ^r_φφ = -r sin²θ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(0, 2, 2, pos);
			Real expected = -r * std::sin(theta) * std::sin(theta);
			REQUIRE_THAT(gamma, WithinAbs(expected, REAL(1e-6)));
		}
		
		SECTION("Γ^θ_rθ = 1/r")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(1, 0, 1, pos);
			REQUIRE_THAT(gamma, WithinAbs(REAL(1.0) / r, REAL(1e-6)));
		}
		
		SECTION("Γ^θ_φφ = -sinθ cosθ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(1, 2, 2, pos);
			Real expected = -std::sin(theta) * std::cos(theta);
			REQUIRE_THAT(gamma, WithinAbs(expected, REAL(1e-6)));
		}
		
		SECTION("Γ^φ_rφ = 1/r")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(2, 0, 2, pos);
			REQUIRE_THAT(gamma, WithinAbs(REAL(1.0) / r, REAL(1e-6)));
		}
		
		SECTION("Γ^φ_θφ = cot θ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(2, 1, 2, pos);
			Real cot_theta = std::cos(theta) / std::sin(theta);
			REQUIRE_THAT(gamma, WithinAbs(cot_theta, REAL(1e-6)));
		}
	}
	
	/********************************************************************************************************************/
	/********                           CYLINDRICAL METRIC TESTS                                                 ********/
	/********************************************************************************************************************/
	
	TEST_CASE("MetricTensorCylindrical - Metric components", "[MetricTensor][Cylindrical]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCylindrical metric;
		
		// Position: rho=3, phi=pi/6, z=2
		Real rho = REAL(3.0);
		Real phi = Constants::PI / REAL(6.0);
		Real z = REAL(2.0);
		VectorN<Real, 3> pos({rho, phi, z});
		
		SECTION("Diagonal covariant components")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			// g_ρρ = 1
			REQUIRE_THAT(g(0, 0), WithinAbs(REAL(1.0), REAL(1e-12)));
			
			// g_φφ = ρ^2 = 9
			REQUIRE_THAT(g(1, 1), WithinAbs(rho * rho, REAL(1e-12)));
			
			// g_zz = 1
			REQUIRE_THAT(g(2, 2), WithinAbs(REAL(1.0), REAL(1e-12)));
		}
		
		SECTION("Off-diagonal components are zero")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			REQUIRE_THAT(g(0, 1), WithinAbs(REAL(0.0), REAL(1e-12)));
			REQUIRE_THAT(g(0, 2), WithinAbs(REAL(0.0), REAL(1e-12)));
			REQUIRE_THAT(g(1, 2), WithinAbs(REAL(0.0), REAL(1e-12)));
		}
	}
	
	TEST_CASE("MetricTensorCylindrical - Christoffel symbols", "[MetricTensor][Cylindrical][Christoffel]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCylindrical metric;
		
		Real rho = REAL(3.0);
		Real phi = Constants::PI / REAL(6.0);
		Real z = REAL(2.0);
		VectorN<Real, 3> pos({rho, phi, z});
		
		// Analytical Christoffel symbols for cylindrical coordinates:
		// Non-zero symbols:
		// Γ^ρ_φφ = -ρ
		// Γ^φ_ρφ = Γ^φ_φρ = 1/ρ
		
		SECTION("Γ^ρ_φφ = -ρ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(0, 1, 1, pos);
			REQUIRE_THAT(gamma, WithinAbs(-rho, REAL(1e-6)));
		}
		
		SECTION("Γ^φ_ρφ = 1/ρ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(1, 0, 1, pos);
			REQUIRE_THAT(gamma, WithinAbs(REAL(1.0) / rho, REAL(1e-6)));
		}
		
		SECTION("z-related Christoffel symbols vanish")
		{
			// All Christoffel symbols involving z should be zero
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					if (i == 2 || j == 2)  // involving z
					{
						Real gamma = metric.GetChristoffelSymbolSecondKind(i, j, 2, pos);
						INFO("Gamma^" << i << "_" << j << "2");
						REQUIRE_THAT(gamma, WithinAbs(REAL(0.0), REAL(1e-8)));
						
						gamma = metric.GetChristoffelSymbolSecondKind(2, i, j, pos);
						INFO("Gamma^2_" << i << j);
						REQUIRE_THAT(gamma, WithinAbs(REAL(0.0), REAL(1e-8)));
					}
				}
			}
		}
	}
	
	/********************************************************************************************************************/
	/********                           MINKOWSKI METRIC TESTS                                                   ********/
	/********************************************************************************************************************/
	
	TEST_CASE("MetricTensorMinkowski - Signature (-,+,+,+)", "[MetricTensor][Minkowski]")
	{
		TEST_PRECISION_INFO();
		MetricTensorMinkowski metric;
		
		VectorN<Real, 4> pos({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		
		SECTION("Diagonal components")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			// eta_00 = -1 (time component)
			REQUIRE_THAT(g(0, 0), WithinAbs(-REAL(1.0), REAL(1e-12)));
			
			// eta_11 = eta_22 = eta_33 = +1 (space components)
			REQUIRE_THAT(g(1, 1), WithinAbs(REAL(1.0), REAL(1e-12)));
			REQUIRE_THAT(g(2, 2), WithinAbs(REAL(1.0), REAL(1e-12)));
			REQUIRE_THAT(g(3, 3), WithinAbs(REAL(1.0), REAL(1e-12)));
		}
		
		SECTION("Off-diagonal components are zero")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (i != j)
					{
						INFO("g(" << i << "," << j << ")");
						REQUIRE_THAT(g(i, j), WithinAbs(REAL(0.0), REAL(1e-12)));
					}
				}
			}
		}
		
		SECTION("All Christoffel symbols vanish (flat spacetime)")
		{
			// Only test a representative sample (testing all 64 would be slow)
			for (int i = 0; i < 4; i++)
			{
				Real gamma = metric.GetChristoffelSymbolSecondKind(i, i, i, pos);
				INFO("Gamma^" << i << "_" << i << i);
				REQUIRE_THAT(gamma, WithinAbs(REAL(0.0), REAL(1e-8)));
			}
		}
	}
	
	/********************************************************************************************************************/
	/********                        METRIC FROM COORD TRANSFORMATION TESTS                                      ********/
	/********************************************************************************************************************/
	
	TEST_CASE("MetricTensorFromCoordTransf - Spherical from transformation", "[MetricTensor][FromTransf]")
	{
		TEST_PRECISION_INFO();
		
		CoordTransfSphericalToCartesian coordTransf;
		MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> metricFromTransf(coordTransf);
		MetricTensorSpherical metricDirect;
		
		Real r = REAL(2.0);
		Real theta = Constants::PI / REAL(4.0);
		Real phi = Constants::PI / REAL(3.0);
		VectorN<Real, 3> pos({r, theta, phi});
		
		SECTION("Metric from transformation matches direct definition")
		{
			auto g_transf = metricFromTransf.GetCovariantMetric(pos);
			auto g_direct = metricDirect.GetCovariantMetric(pos);
			
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					INFO("g(" << i << "," << j << ")");
					REQUIRE_THAT(g_transf(i, j), WithinAbs(g_direct(i, j), REAL(1e-8)));
				}
			}
		}
	}
	
	TEST_CASE("MetricTensorFromCoordTransf - Cylindrical from transformation", "[MetricTensor][FromTransf]")
	{
		TEST_PRECISION_INFO();
		
		CoordTransfCylindricalToCartesian coordTransf;
		MetricTensorFromCoordTransf<Vector3Cylindrical, Vector3Cartesian, 3> metricFromTransf(coordTransf);
		MetricTensorCylindrical metricDirect;
		
		Real rho = REAL(3.0);
		Real phi = Constants::PI / REAL(6.0);
		Real z = REAL(2.0);
		VectorN<Real, 3> pos({rho, phi, z});
		
		SECTION("Metric from transformation matches direct definition")
		{
			auto g_transf = metricFromTransf.GetCovariantMetric(pos);
			auto g_direct = metricDirect.GetCovariantMetric(pos);
			
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					INFO("g(" << i << "," << j << ")");
					REQUIRE_THAT(g_transf(i, j), WithinAbs(g_direct(i, j), REAL(1e-8)));
				}
			}
		}
	}
	
	/********************************************************************************************************************/
	/********                           COVARIANT DERIVATIVE TESTS                                               ********/
	/********************************************************************************************************************/
	
	TEST_CASE("Covariant derivative - Cartesian (reduces to partial)", "[MetricTensor][CovariantDerivative]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		
		// Vector field: V = (x^2, y^2, z^2)
		class QuadraticField : public IVectorFunction<3>
		{
		public:
			VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
			{
				return VectorN<Real, 3>({pos[0] * pos[0], pos[1] * pos[1], pos[2] * pos[2]});
			}
		};
		
		QuadraticField field;
		VectorN<Real, 3> pos({REAL(2.0), REAL(3.0), REAL(4.0)});
		
		SECTION("Covariant derivative equals partial derivative in flat space")
		{
			// ∇_j V^i = ∂_j V^i (since Christoffel symbols vanish)
			// ∂V^0/∂x^0 = 2x = 4
			// ∂V^1/∂x^1 = 2y = 6
			// ∂V^2/∂x^2 = 2z = 8
			
			auto nabla_0 = metric.CovariantDerivativeContravar(field, 0, pos);
			REQUIRE_THAT(nabla_0[0], WithinAbs(REAL(2.0) * pos[0], REAL(1e-6)));
			REQUIRE_THAT(nabla_0[1], WithinAbs(REAL(0.0), REAL(1e-6)));
			REQUIRE_THAT(nabla_0[2], WithinAbs(REAL(0.0), REAL(1e-6)));
			
			auto nabla_1 = metric.CovariantDerivativeContravar(field, 1, pos);
			REQUIRE_THAT(nabla_1[0], WithinAbs(REAL(0.0), REAL(1e-6)));
			REQUIRE_THAT(nabla_1[1], WithinAbs(REAL(2.0) * pos[1], REAL(1e-6)));
			REQUIRE_THAT(nabla_1[2], WithinAbs(REAL(0.0), REAL(1e-6)));
		}
	}
	
	/********************************************************************************************************************/
	/********                           CHRISTOFFEL SYMBOL FIRST KIND TESTS                                      ********/
	/********************************************************************************************************************/
	
	TEST_CASE("Christoffel first kind - Relation to second kind", "[MetricTensor][Christoffel]")
	{
		TEST_PRECISION_INFO();
		MetricTensorSpherical metric;
		
		Real r = REAL(2.0);
		Real theta = Constants::PI / REAL(4.0);
		Real phi = Constants::PI / REAL(3.0);
		VectorN<Real, 3> pos({r, theta, phi});
		
		SECTION("Γ_{ijk} = g_{km} Γ^m_{ij}")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			// Test for a few index combinations
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					for (int k = 0; k < 2; k++)
					{
						Real gamma_first = metric.GetChristoffelSymbolFirstKind(i, j, k, pos);
						
						// Compute from second kind
						Real gamma_computed = 0.0;
						for (int m = 0; m < 3; m++)
						{
							gamma_computed += g(m, k) * metric.GetChristoffelSymbolSecondKind(m, i, j, pos);
						}
						
						INFO("Gamma_" << i << j << k);
						REQUIRE_THAT(gamma_first, WithinAbs(gamma_computed, REAL(1e-6)));
					}
				}
			}
		}
	}
	
	/********************************************************************************************************************/
	/********                           ORIGINAL TEST (UPDATED)                                                  ********/
	/********************************************************************************************************************/
	
	TEST_CASE("Test_Metric_Tensors - Basic functionality", "[MetricTensor][basic]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metricCart;
		MetricTensorSpherical metricSpher;
		MetricTensorCylindrical metricCyl;

		CoordTransfSphericalToCartesian coordTransfSpherToCart;

		MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> metricSpherFromCart(coordTransfSpherToCart);
		MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> metricSpherFromCart2(CoordTransfSpherToCart);

		Vector3Cartesian pos(REAL(1.0), REAL(2.0), -REAL(1.0));
		Vector3Spherical posSpher = CoordTransfSpherToCart.transf(pos);
		Vector3Cylindrical posCyl = CoordTransfCylToCart.transf(pos);

		auto cart_metric = metricCart(pos);
		auto spher_metric = metricSpher(posSpher);
		auto cyl_metric = metricCyl(posCyl);
		
		// Verify Cartesian metric is identity
		REQUIRE_THAT(cart_metric(0, 0), WithinAbs(REAL(1.0), REAL(1e-12)));
		REQUIRE_THAT(cart_metric(1, 1), WithinAbs(REAL(1.0), REAL(1e-12)));
		REQUIRE_THAT(cart_metric(2, 2), WithinAbs(REAL(1.0), REAL(1e-12)));
	}

	TEST_CASE("MetricTensorSphericalContravar - Singularity at r=0 throws", "[MetricTensor][Spherical][Singularity]")
	{
		MetricTensorSphericalContravar metric;
		VectorN<Real, 3> at_origin{REAL(0.0), REAL(1.0), REAL(0.5)};

		// g^00 = 1 should be fine
		REQUIRE_THAT(metric.Component(0, 0, at_origin), WithinAbs(REAL(1.0), REAL(1e-12)));

		// g^11 = 1/r^2 should throw at r=0
		REQUIRE_THROWS(metric.Component(1, 1, at_origin));

		// g^22 = 1/(r^2*sin^2(theta)) should throw at r=0
		REQUIRE_THROWS(metric.Component(2, 2, at_origin));
	}

	TEST_CASE("MetricTensorSphericalContravar - Singularity at pole throws", "[MetricTensor][Spherical][Singularity]")
	{
		MetricTensorSphericalContravar metric;
		VectorN<Real, 3> at_pole{REAL(2.0), REAL(0.0), REAL(0.5)};  // theta=0 (north pole)

		// g^11 = 1/r^2 should be fine at r=2
		REQUIRE_THAT(metric.Component(1, 1, at_pole), WithinAbs(REAL(0.25), REAL(1e-12)));

		// g^22 = 1/(r^2*sin^2(0)) should throw at theta=0
		REQUIRE_THROWS(metric.Component(2, 2, at_pole));
	}
} // namespace MML::Tests::Core::MetricTensorTests
