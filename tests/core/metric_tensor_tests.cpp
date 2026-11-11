#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Function.h>
#include <mml/base/Tensor.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Vector/VectorTypes4D.h>

#include <mml/core/CoordTransf/CoordTransfBase.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/algorithms/Geodesic.h>
#include <mml/core/MetricTensor.h>

#include <mml/base/Geometry/Geometry3D.h>
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Core::MetricTensorTests
{
	class TrackingMetricTransform : public CoordTransfCylindricalToCartesian
	{
	public:
		mutable int jacobianCalls = 0;

		MatrixNM<Real, 3, 3> jacobian(const VectorN<Real, 3>&) const override
		{
			++jacobianCalls;
			MatrixNM<Real, 3, 3> jac;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					jac(i, j) = i == j ? static_cast<MML::Real>(2.0) : static_cast<MML::Real>(0.0);
			return jac;
		}
	};

	class NumericalJacobianTransform : public CoordTransf<Vector3Cartesian, Vector3Cartesian, 3>
	{
		static Real scaleX(const VectorN<Real, 3>& point) { return static_cast<MML::Real>(2.0) * point[0]; }
		static Real scaleY(const VectorN<Real, 3>& point) { return static_cast<MML::Real>(3.0) * point[1]; }
		static Real scaleZ(const VectorN<Real, 3>& point) { return static_cast<MML::Real>(4.0) * point[2]; }

		inline static ScalarFunction<3> _functions[3] = {
			ScalarFunction<3>{scaleX}, ScalarFunction<3>{scaleY}, ScalarFunction<3>{scaleZ}
		};

	public:
		Vector3Cartesian transf(const Vector3Cartesian& point) const override
		{
			return Vector3Cartesian{ scaleX(point), scaleY(point), scaleZ(point) };
		}

		const IScalarFunction<3>& coordTransfFunc(int i) const override { return _functions[i]; }
	};

	class UnitSphereMetric : public MetricTensorField<2>
	{
	public:
		UnitSphereMetric() : MetricTensorField<2>(0, 2) { }

		Real component(int i, int j, const VectorN<Real, 2>& pos) const override
		{
			if (i == 0 && j == 0)
				return static_cast<MML::Real>(1.0);
			if (i == 1 && j == 1)
				return std::sin(pos[0]) * std::sin(pos[0]);
			return static_cast<MML::Real>(0.0);
		}
	};

	class SchwarzschildMetric : public LorentzianMetric<4>
	{
	protected:
		Real _schwarzschildRadius;

	public:
		SchwarzschildMetric(Real schwarzschildRadius) : LorentzianMetric<4>(0, 2), _schwarzschildRadius(schwarzschildRadius) { }

		Real component(int i, int j, const VectorN<Real, 4>& pos) const override
		{
			Real r = pos[1];
			Real theta = pos[2];
			Real f = static_cast<MML::Real>(1.0) - _schwarzschildRadius / r;

			if (i == 0 && j == 0)
				return -f;
			if (i == 1 && j == 1)
				return static_cast<MML::Real>(1.0) / f;
			if (i == 2 && j == 2)
				return r * r;
			if (i == 3 && j == 3)
				return r * r * std::sin(theta) * std::sin(theta);
			return static_cast<MML::Real>(0.0);
		}
	};

	/********************************************************************************************************************/
	/********                           CARTESIAN METRIC TESTS                                                   ********/
	/********************************************************************************************************************/
	
	TEST_CASE("MetricTensorCartesian3D - Identity metric", "[MetricTensor][Cartesian]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		
		Vector3Cartesian pos(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		
		SECTION("Covariant metric is identity")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					Real expected = (i == j) ? static_cast<MML::Real>(1.0) : static_cast<MML::Real>(0.0);
					REQUIRE_THAT(g(i, j), WithinAbs(expected, TOL(1e-12, 1e-5)));
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
					Real expected = (i == j) ? static_cast<MML::Real>(1.0) : static_cast<MML::Real>(0.0);
					REQUIRE_THAT(g_inv(i, j), WithinAbs(expected, TOL(1e-10, 1e-5)));
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
						REQUIRE_THAT(gamma, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
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
		Real r = static_cast<MML::Real>(2.0);
		Real theta = Constants::PI / static_cast<MML::Real>(4.0);
		Real phi = Constants::PI / static_cast<MML::Real>(3.0);
		VectorN<Real, 3> pos({r, theta, phi});
		
		SECTION("Diagonal covariant components")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			// g_rr = 1
			REQUIRE_THAT(g(0, 0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
			
			// g_theta,theta = r^2 = 4
			REQUIRE_THAT(g(1, 1), WithinAbs(r * r, TOL(1e-12, 1e-5)));
			
			// g_phi,phi = r^2 * sin^2(theta)
			Real expected_gphi = r * r * std::sin(theta) * std::sin(theta);
			REQUIRE_THAT(g(2, 2), WithinAbs(expected_gphi, TOL(1e-12, 1e-5)));
		}
		
		SECTION("Off-diagonal components are zero")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			REQUIRE_THAT(g(0, 1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(0, 2), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(1, 0), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(1, 2), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(2, 0), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(2, 1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
		}
		
		SECTION("Contravariant metric is inverse of covariant")
		{
			auto g = metric.GetCovariantMetric(pos);
			auto g_inv = metric.GetContravariantMetric(pos);
			
			// g^rr = 1
			REQUIRE_THAT(g_inv(0, 0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
			
			// g^theta,theta = 1/r^2
			REQUIRE_THAT(g_inv(1, 1), WithinAbs(static_cast<MML::Real>(1.0) / (r * r), TOL(1e-10, 1e-5)));
			
			// g^phi,phi = 1/(r^2 * sin^2(theta))
			Real sin_theta = std::sin(theta);
			REQUIRE_THAT(g_inv(2, 2), WithinAbs(static_cast<MML::Real>(1.0) / (r * r * sin_theta * sin_theta), TOL(1e-10, 1e-5)));
		}
	}
	
	TEST_CASE("MetricTensorSpherical - Christoffel symbols", "[MetricTensor][Spherical][Christoffel]")
	{
		TEST_PRECISION_INFO();
		MetricTensorSpherical metric;
		
		Real r = static_cast<MML::Real>(2.0);
		Real theta = Constants::PI / static_cast<MML::Real>(4.0);
		Real phi = Constants::PI / static_cast<MML::Real>(3.0);
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
			REQUIRE_THAT(gamma, WithinAbs(-r, TOL(1e-6, 1e-3)));
		}
		
		SECTION("Γ^r_φφ = -r sin²θ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(0, 2, 2, pos);
			Real expected = -r * std::sin(theta) * std::sin(theta);
			REQUIRE_THAT(gamma, WithinAbs(expected, TOL(1e-6, 1e-3)));
		}
		
		SECTION("Γ^θ_rθ = 1/r")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(1, 0, 1, pos);
			REQUIRE_THAT(gamma, WithinAbs(static_cast<MML::Real>(1.0) / r, TOL(1e-6, 1e-3)));
		}
		
		SECTION("Γ^θ_φφ = -sinθ cosθ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(1, 2, 2, pos);
			Real expected = -std::sin(theta) * std::cos(theta);
			REQUIRE_THAT(gamma, WithinAbs(expected, TOL(1e-6, 1e-3)));
		}
		
		SECTION("Γ^φ_rφ = 1/r")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(2, 0, 2, pos);
			REQUIRE_THAT(gamma, WithinAbs(static_cast<MML::Real>(1.0) / r, TOL(1e-6, 1e-3)));
		}
		
		SECTION("Γ^φ_θφ = cot θ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(2, 1, 2, pos);
			Real cot_theta = std::cos(theta) / std::sin(theta);
			REQUIRE_THAT(gamma, WithinAbs(cot_theta, TOL(1e-6, 1e-3)));
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
		Real rho = static_cast<MML::Real>(3.0);
		Real phi = Constants::PI / static_cast<MML::Real>(6.0);
		Real z = static_cast<MML::Real>(2.0);
		VectorN<Real, 3> pos({rho, phi, z});
		
		SECTION("Diagonal covariant components")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			// g_ρρ = 1
			REQUIRE_THAT(g(0, 0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
			
			// g_φφ = ρ^2 = 9
			REQUIRE_THAT(g(1, 1), WithinAbs(rho * rho, TOL(1e-12, 1e-5)));
			
			// g_zz = 1
			REQUIRE_THAT(g(2, 2), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
		}
		
		SECTION("Off-diagonal components are zero")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			REQUIRE_THAT(g(0, 1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(0, 2), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(1, 2), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
		}
	}
	
	TEST_CASE("MetricTensorCylindrical - Christoffel symbols", "[MetricTensor][Cylindrical][Christoffel]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCylindrical metric;
		
		Real rho = static_cast<MML::Real>(3.0);
		Real phi = Constants::PI / static_cast<MML::Real>(6.0);
		Real z = static_cast<MML::Real>(2.0);
		VectorN<Real, 3> pos({rho, phi, z});
		
		// Analytical Christoffel symbols for cylindrical coordinates:
		// Non-zero symbols:
		// Γ^ρ_φφ = -ρ
		// Γ^φ_ρφ = Γ^φ_φρ = 1/ρ
		
		SECTION("Γ^ρ_φφ = -ρ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(0, 1, 1, pos);
			REQUIRE_THAT(gamma, WithinAbs(-rho, TOL(1e-6, 1e-3)));
		}
		
		SECTION("Γ^φ_ρφ = 1/ρ")
		{
			Real gamma = metric.GetChristoffelSymbolSecondKind(1, 0, 1, pos);
			REQUIRE_THAT(gamma, WithinAbs(static_cast<MML::Real>(1.0) / rho, TOL(1e-6, 1e-3)));
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
						REQUIRE_THAT(gamma, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
						
						gamma = metric.GetChristoffelSymbolSecondKind(2, i, j, pos);
						INFO("Gamma^2_" << i << j);
						REQUIRE_THAT(gamma, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
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
		
		VectorN<Real, 4> pos({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)});
		
		SECTION("Diagonal components")
		{
			auto g = metric.GetCovariantMetric(pos);
			
			// eta_00 = -1 (time component)
			REQUIRE_THAT(g(0, 0), WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
			
			// eta_11 = eta_22 = eta_33 = +1 (space components)
			REQUIRE_THAT(g(1, 1), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(2, 2), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(g(3, 3), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
		}

		SECTION("Metric signature type is Lorentzian")
		{
			REQUIRE(metric.IsLorentzian());
			REQUIRE_FALSE(metric.IsRiemannian());
			REQUIRE(MetricTensorCartesian3D().IsRiemannian());
			REQUIRE_FALSE(MetricTensorCartesian3D().IsLorentzian());
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
						REQUIRE_THAT(g(i, j), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
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
				REQUIRE_THAT(gamma, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
			}
		}

		SECTION("Metric contraction agrees with Vector4Minkowski scalar product")
		{
			Vector4Minkowski a{static_cast<MML::Real>(5.0), static_cast<MML::Real>(2.0), -static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0)};
			Vector4Minkowski b{static_cast<MML::Real>(4.0), -static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0)};
			auto eta = metric.GetCovariantMetric(pos);

			Real metric_product = static_cast<MML::Real>(0.0);
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					metric_product += eta(i, j) * a[i] * b[j];

			REQUIRE_THAT(ScalarProduct(a, b), WithinAbs(metric_product, TOL(1e-12, 1e-5)));
			REQUIRE_THAT(ScalarProduct(a, a), WithinAbs(-static_cast<MML::Real>(11.0), TOL(1e-12, 1e-5)));
		}

		SECTION("Eta contraction eta_munu eta^munu equals spacetime dimension")
		{
			auto eta = metric.GetCovariantMetric(pos);
			auto eta_inv = metric.GetContravariantMetric(pos);

			Real contraction = static_cast<MML::Real>(0.0);
			for (int mu = 0; mu < 4; mu++)
				for (int nu = 0; nu < 4; nu++)
					contraction += eta(mu, nu) * eta_inv(mu, nu);

			REQUIRE_THAT(contraction, WithinAbs(static_cast<MML::Real>(4.0), TOL(1e-12, 1e-5)));
		}

		SECTION("Classifies timelike spacelike and null intervals")
		{
			VectorN<Real, 4> timelike({static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)});
			VectorN<Real, 4> spacelike({static_cast<MML::Real>(3.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)});
			VectorN<Real, 4> null({static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(0.0)});

			REQUIRE_THAT(metric.IntervalSquared(timelike, pos), WithinAbs(-static_cast<MML::Real>(16.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(metric.IntervalSquared(spacelike, pos), WithinAbs(static_cast<MML::Real>(16.0), TOL(1e-12, 1e-5)));
			REQUIRE_THAT(metric.IntervalSquared(null, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));

			REQUIRE(metric.IsTimelike(timelike, pos));
			REQUIRE(metric.IsSpacelike(spacelike, pos));
			REQUIRE(metric.IsNull(null, pos));
			REQUIRE(metric.ClassifyInterval(timelike, pos) == LorentzianMetric<4>::IntervalType::Timelike);
			REQUIRE(metric.ClassifyInterval(spacelike, pos) == LorentzianMetric<4>::IntervalType::Spacelike);
			REQUIRE(metric.ClassifyInterval(null, pos) == LorentzianMetric<4>::IntervalType::Null);
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
		
		Real r = static_cast<MML::Real>(2.0);
		Real theta = Constants::PI / static_cast<MML::Real>(4.0);
		Real phi = Constants::PI / static_cast<MML::Real>(3.0);
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
					REQUIRE_THAT(g_transf(i, j), WithinAbs(g_direct(i, j), TOL(1e-8, 1e-4)));
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
		
		Real rho = static_cast<MML::Real>(3.0);
		Real phi = Constants::PI / static_cast<MML::Real>(6.0);
		Real z = static_cast<MML::Real>(2.0);
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
					REQUIRE_THAT(g_transf(i, j), WithinAbs(g_direct(i, j), TOL(1e-8, 1e-4)));
				}
			}
		}
	}

	TEST_CASE("MetricTensorFromCoordTransf - uses one virtual Jacobian evaluation", "[MetricTensor][FromTransf][Jacobian]")
	{
		TEST_PRECISION_INFO();
		TrackingMetricTransform transform;
		MetricTensorFromCoordTransf<Vector3Cylindrical, Vector3Cartesian, 3> metric(transform);
		const Vector3Cylindrical pos{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) };

		const auto covariantMetric = metric.GetCovariantMetric(pos);
		REQUIRE(transform.jacobianCalls == 1);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				REQUIRE(covariantMetric(i, j) == (i == j ? static_cast<MML::Real>(4.0) : static_cast<MML::Real>(0.0)));

		transform.jacobianCalls = 0;
		ScalarFunction<3> scalarField([](const VectorN<Real, 3>& point) {
			return point[0] * point[0] + point[1] * point[1] + point[2] * point[2];
		});
		const auto gradient = ScalarFieldOperations::Gradient<3>(scalarField, pos, metric);

		REQUIRE(transform.jacobianCalls == 1);
		REQUIRE(gradient.IsEqualTo(VectorN<Real, 3>{ static_cast<MML::Real>(0.5), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.5) }, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("MetricTensorFromCoordTransf - preserves numerical Jacobian fallback", "[MetricTensor][FromTransf][Jacobian]")
	{
		TEST_PRECISION_INFO();
		NumericalJacobianTransform transform;
		MetricTensorFromCoordTransf<Vector3Cartesian, Vector3Cartesian, 3> metric(transform);
		const auto covariantMetric = metric.GetCovariantMetric(Vector3Cartesian{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) });

		const Real expectedDiagonal[3] = { static_cast<MML::Real>(4.0), static_cast<MML::Real>(9.0), static_cast<MML::Real>(16.0) };
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				REQUIRE_THAT(covariantMetric(i, j), WithinAbs(i == j ? expectedDiagonal[i] : static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
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
		VectorN<Real, 3> pos({static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)});
		
		SECTION("Covariant derivative equals partial derivative in flat space")
		{
			// ∇_j V^i = ∂_j V^i (since Christoffel symbols vanish)
			// ∂V^0/∂x^0 = 2x = 4
			// ∂V^1/∂x^1 = 2y = 6
			// ∂V^2/∂x^2 = 2z = 8
			
			auto nabla_0 = metric.CovariantDerivativeContravar(field, 0, pos);
			REQUIRE_THAT(nabla_0[0], WithinAbs(static_cast<MML::Real>(2.0) * pos[0], TOL(1e-6, 1e-4)));
			REQUIRE_THAT(nabla_0[1], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-4)));
			REQUIRE_THAT(nabla_0[2], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-4)));
			
			auto nabla_1 = metric.CovariantDerivativeContravar(field, 1, pos);
			REQUIRE_THAT(nabla_1[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-4)));
			REQUIRE_THAT(nabla_1[1], WithinAbs(static_cast<MML::Real>(2.0) * pos[1], TOL(1e-6, 1e-4)));
			REQUIRE_THAT(nabla_1[2], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-4)));
		}
	}
	
	/********************************************************************************************************************/
	/********                           CHRISTOFFEL SYMBOL FIRST KIND TESTS                                      ********/
	/********************************************************************************************************************/
	
	TEST_CASE("Christoffel first kind - Relation to second kind", "[MetricTensor][Christoffel]")
	{
		TEST_PRECISION_INFO();
		MetricTensorSpherical metric;
		
		Real r = static_cast<MML::Real>(2.0);
		Real theta = Constants::PI / static_cast<MML::Real>(4.0);
		Real phi = Constants::PI / static_cast<MML::Real>(3.0);
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
						REQUIRE_THAT(gamma_first, WithinAbs(gamma_computed, static_cast<MML::Real>(1e-6)));
					}
				}
			}
		}
	}

	/********************************************************************************************************************/
	/********                           RIEMANN CURVATURE TESTS                                                  ********/
	/********************************************************************************************************************/

	TEST_CASE("Riemann curvature - Flat coordinate systems vanish", "[MetricTensor][Riemann]")
	{
		TEST_PRECISION_INFO();
		VectorN<Real, 3> cartPos({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		VectorN<Real, 3> sphericalPos({static_cast<MML::Real>(2.0), Constants::PI / static_cast<MML::Real>(4.0), Constants::PI / static_cast<MML::Real>(3.0)});
		VectorN<Real, 3> cylindricalPos({static_cast<MML::Real>(3.0), Constants::PI / static_cast<MML::Real>(6.0), static_cast<MML::Real>(2.0)});

		MetricTensorCartesian3D cartesian;
		MetricTensorSpherical spherical;
		MetricTensorCylindrical cylindrical;

		for (int rho = 0; rho < 3; rho++)
			for (int sigma = 0; sigma < 3; sigma++)
				for (int mu = 0; mu < 3; mu++)
					for (int nu = 0; nu < 3; nu++)
					{
						INFO("R^" << rho << "_" << sigma << mu << nu);
						REQUIRE_THAT(cartesian.GetRiemannCurvatureTensor(rho, sigma, mu, nu, cartPos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-7, 1e-3)));
						REQUIRE_THAT(spherical.GetRiemannCurvatureTensor(rho, sigma, mu, nu, sphericalPos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-2)));
						REQUIRE_THAT(cylindrical.GetRiemannCurvatureTensor(rho, sigma, mu, nu, cylindricalPos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-2)));
					}
	}

	TEST_CASE("Riemann curvature - Unit sphere has constant positive curvature", "[MetricTensor][Riemann]")
	{
		TEST_PRECISION_INFO();
		UnitSphereMetric metric;
		VectorN<Real, 2> pos({Constants::PI / static_cast<MML::Real>(4.0), Constants::PI / static_cast<MML::Real>(3.0)});

		Real sinTheta = std::sin(pos[0]);
		Real expected = sinTheta * sinTheta;
		Real rThetaPhiThetaPhi = metric.GetRiemannCurvatureTensor(0, 1, 0, 1, pos);
		Real rThetaPhiPhiTheta = metric.GetRiemannCurvatureTensor(0, 1, 1, 0, pos);

		REQUIRE_THAT(rThetaPhiThetaPhi, WithinAbs(expected, TOL(1e-5, 1e-2)));
		REQUIRE_THAT(rThetaPhiPhiTheta, WithinAbs(-expected, TOL(1e-5, 1e-2)));

		auto riemann = metric.GetRiemannCurvatureTensor(pos);
		REQUIRE(riemann.numContravar() == 1);
		REQUIRE(riemann.numCovar() == 3);
		REQUIRE(riemann.isContravar(3));
		REQUIRE(riemann.isCovar(0));
		REQUIRE_THAT(riemann(0, 1, 0, 1), WithinAbs(expected, TOL(1e-5, 1e-2)));
	}

	TEST_CASE("Ricci tensor and scalar - Flat metrics vanish", "[MetricTensor][Ricci]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		VectorN<Real, 3> pos({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});

		for (int mu = 0; mu < 3; mu++)
			for (int nu = 0; nu < 3; nu++)
			{
				INFO("R_" << mu << nu);
				REQUIRE_THAT(metric.GetRicciTensor(mu, nu, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-7, 1e-3)));
			}

		auto ricci = metric.GetRicciTensor(pos);
		REQUIRE(ricci.numContravar() == 0);
		REQUIRE(ricci.numCovar() == 2);
		REQUIRE_THAT(metric.GetRicciScalar(pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-7, 1e-3)));
	}

	TEST_CASE("Ricci tensor and scalar - Unit sphere", "[MetricTensor][Ricci]")
	{
		TEST_PRECISION_INFO();
		UnitSphereMetric metric;
		VectorN<Real, 2> pos({Constants::PI / static_cast<MML::Real>(4.0), Constants::PI / static_cast<MML::Real>(3.0)});

		Real sinTheta = std::sin(pos[0]);
		Real expectedPhiPhi = sinTheta * sinTheta;

		REQUIRE_THAT(metric.GetRicciTensor(0, 0, pos), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-2)));
		REQUIRE_THAT(metric.GetRicciTensor(1, 1, pos), WithinAbs(expectedPhiPhi, TOL(1e-5, 1e-2)));
		REQUIRE_THAT(metric.GetRicciTensor(0, 1, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-2)));
		REQUIRE_THAT(metric.GetRicciTensor(1, 0, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-2)));

		auto ricci = metric.GetRicciTensor(pos);
		REQUIRE_THAT(ricci(0, 0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-2)));
		REQUIRE_THAT(ricci(1, 1), WithinAbs(expectedPhiPhi, TOL(1e-5, 1e-2)));
		REQUIRE_THAT(metric.GetRicciScalar(pos), WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-5, 1e-2)));
	}

	TEST_CASE("Einstein tensor - Flat metrics vanish", "[MetricTensor][Einstein]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		VectorN<Real, 3> pos({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});

		for (int mu = 0; mu < 3; mu++)
			for (int nu = 0; nu < 3; nu++)
			{
				INFO("G_" << mu << nu);
				REQUIRE_THAT(metric.GetEinsteinTensor(mu, nu, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-7, 1e-3)));
			}

		auto einstein = metric.GetEinsteinTensor(pos);
		REQUIRE(einstein.numContravar() == 0);
		REQUIRE(einstein.numCovar() == 2);
		REQUIRE_THAT(einstein(0, 0), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-7, 1e-3)));
	}

	TEST_CASE("Einstein tensor - Unit sphere is zero in two dimensions", "[MetricTensor][Einstein]")
	{
		TEST_PRECISION_INFO();
		UnitSphereMetric metric;
		VectorN<Real, 2> pos({Constants::PI / static_cast<MML::Real>(4.0), Constants::PI / static_cast<MML::Real>(3.0)});

		for (int mu = 0; mu < 2; mu++)
			for (int nu = 0; nu < 2; nu++)
			{
				INFO("G_" << mu << nu);
				REQUIRE_THAT(metric.GetEinsteinTensor(mu, nu, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-2)));
			}

		auto einstein = metric.GetEinsteinTensor(pos);
		REQUIRE_THAT(einstein(0, 0), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-2)));
		REQUIRE_THAT(einstein(1, 1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-2)));
	}

	/********************************************************************************************************************/
	/********                           GEODESIC EQUATION TESTS                                                  ********/
	/********************************************************************************************************************/

	TEST_CASE("GeodesicEquationSystem - Cartesian derivatives are straight-line motion", "[MetricTensor][Geodesic]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		GeodesicEquationSystem<3> system(metric);

		Vector3Cartesian position(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Vector3Cartesian velocity(static_cast<MML::Real>(0.5), -static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
		Vector<Real> state = MakeGeodesicState<3>(position, velocity);
		Vector<Real> dstate(6);

		system.derivs(static_cast<MML::Real>(0.0), state, dstate);

		REQUIRE(system.getDim() == 6);
		REQUIRE(system.getVarName(0) == "q0");
		REQUIRE(system.getVarName(3) == "v0");
		for (int i = 0; i < 3; i++)
		{
			REQUIRE_THAT(dstate[i], WithinAbs(velocity[i], TOL(1e-12, 1e-5)));
			REQUIRE_THAT(dstate[3 + i], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
		}
	}

	TEST_CASE("IntegrateGeodesicFixedStep - Cartesian geodesic is a straight line", "[MetricTensor][Geodesic]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		Vector3Cartesian position(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Vector3Cartesian velocity(static_cast<MML::Real>(0.5), -static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
		Real lambdaEnd = static_cast<MML::Real>(2.0);

		auto solution = IntegrateGeodesicFixedStep<3>(metric, position, velocity, static_cast<MML::Real>(0.0), lambdaEnd, 40);
		Vector<Real> finalState = solution.getXValuesAtEnd();
		VectorN<Real, 3> finalPosition = GeodesicPositionFromState<3>(finalState);
		VectorN<Real, 3> finalVelocity = GeodesicVelocityFromState<3>(finalState);

		for (int i = 0; i < 3; i++)
		{
			REQUIRE_THAT(finalPosition[i], WithinAbs(position[i] + velocity[i] * lambdaEnd, TOL(1e-8, 1e-4)));
			REQUIRE_THAT(finalVelocity[i], WithinAbs(velocity[i], TOL(1e-8, 1e-4)));
		}
	}

	TEST_CASE("IntegrateGeodesicFixedStep - Unit sphere equator is a geodesic", "[MetricTensor][Geodesic]")
	{
		TEST_PRECISION_INFO();
		UnitSphereMetric metric;
		VectorN<Real, 2> position({Constants::PI / static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)});
		VectorN<Real, 2> velocity({static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});
		Real lambdaEnd = Constants::PI / static_cast<MML::Real>(2.0);

		auto solution = IntegrateGeodesicFixedStep<2>(metric, position, velocity, static_cast<MML::Real>(0.0), lambdaEnd, 120);
		Vector<Real> finalState = solution.getXValuesAtEnd();
		VectorN<Real, 2> finalPosition = GeodesicPositionFromState<2>(finalState);
		VectorN<Real, 2> finalVelocity = GeodesicVelocityFromState<2>(finalState);

		REQUIRE_THAT(finalPosition[0], WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), TOL(1e-6, 1e-3)));
		REQUIRE_THAT(finalPosition[1], WithinAbs(lambdaEnd, TOL(1e-6, 1e-3)));
		REQUIRE_THAT(finalVelocity[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-3)));
		REQUIRE_THAT(finalVelocity[1], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-6, 1e-3)));
	}

	/********************************************************************************************************************/
	/********                           CURVATURE VALIDATION SUITE                                              ********/
	/********************************************************************************************************************/

	TEST_CASE("Curvature validation - Schwarzschild vacuum is Ricci flat with nonzero Riemann", "[MetricTensor][CurvatureValidation][Schwarzschild]")
	{
		TEST_PRECISION_INFO();
		SchwarzschildMetric metric(static_cast<MML::Real>(2.0));
		VectorN<Real, 4> pos({static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), Constants::PI / static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.25)});

		REQUIRE(metric.IsLorentzian());
		REQUIRE_THAT(metric.component(0, 0, pos), WithinAbs(-static_cast<MML::Real>(0.8), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(metric.component(1, 1, pos), WithinAbs(static_cast<MML::Real>(1.25), TOL(1e-12, 1e-5)));

		for (int mu = 0; mu < 4; mu++)
			for (int nu = 0; nu < 4; nu++)
			{
				INFO("Ricci/Einstein component " << mu << "," << nu);
				REQUIRE_THAT(metric.GetRicciTensor(mu, nu, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-4, 1e-2)));
				REQUIRE_THAT(metric.GetEinsteinTensor(mu, nu, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-4, 1e-2)));
			}

		REQUIRE_THAT(metric.GetRicciScalar(pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-4, 1e-2)));
		REQUIRE(std::abs(metric.GetRiemannCurvatureTensor(1, 2, 1, 2, pos)) > TOL(1e-4, 1e-6));
	}

	/********************************************************************************************************************/
	/********                           PARALLEL TRANSPORT TESTS                                                ********/
	/********************************************************************************************************************/

	TEST_CASE("ParallelTransportEquationSystem - Cartesian transport keeps vector constant", "[MetricTensor][ParallelTransport]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		ParametricCurveFromStdFunc<3> curve([](Real lambda) {
			return Vector3Cartesian(lambda, static_cast<MML::Real>(2.0) * lambda, static_cast<MML::Real>(0.0));
		});
		ParallelTransportEquationSystem<3> system(metric, curve);

		Vector3Cartesian initialVector(static_cast<MML::Real>(1.0), -static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.5));
		Vector<Real> state = MakeParallelTransportState<3>(initialVector);
		Vector<Real> dstate(3);

		system.derivs(static_cast<MML::Real>(0.25), state, dstate);

		REQUIRE(system.getDim() == 3);
		REQUIRE(system.getVarName(0) == "V0");
		for (int i = 0; i < 3; i++)
			REQUIRE_THAT(dstate[i], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("IntegrateParallelTransportFixedStep - Cartesian transport keeps vector constant", "[MetricTensor][ParallelTransport]")
	{
		TEST_PRECISION_INFO();
		MetricTensorCartesian3D metric;
		ParametricCurveFromStdFunc<3> curve([](Real lambda) {
			return Vector3Cartesian(lambda, static_cast<MML::Real>(2.0) * lambda, static_cast<MML::Real>(0.0));
		});
		Vector3Cartesian initialVector(static_cast<MML::Real>(1.0), -static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.5));

		auto solution = IntegrateParallelTransportFixedStep<3>(metric, curve, initialVector, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 40);
		VectorN<Real, 3> finalVector = ParallelTransportVectorFromState<3>(solution.getXValuesAtEnd());

		for (int i = 0; i < 3; i++)
			REQUIRE_THAT(finalVector[i], WithinAbs(initialVector[i], TOL(1e-8, 1e-4)));
	}

	TEST_CASE("IntegrateParallelTransportFixedStep - Unit sphere equator keeps theta vector parallel", "[MetricTensor][ParallelTransport]")
	{
		TEST_PRECISION_INFO();
		UnitSphereMetric metric;
		ParametricCurveFromStdFunc<2> equator([](Real lambda) {
			return VectorN<Real, 2>({Constants::PI / static_cast<MML::Real>(2.0), lambda});
		});
		VectorN<Real, 2> initialVector({static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)});

		auto solution = IntegrateParallelTransportFixedStep<2>(metric, equator, initialVector, static_cast<MML::Real>(0.0), Constants::PI / static_cast<MML::Real>(2.0), 120);
		VectorN<Real, 2> finalVector = ParallelTransportVectorFromState<2>(solution.getXValuesAtEnd());

		REQUIRE_THAT(finalVector[0], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-6, 1e-3)));
		REQUIRE_THAT(finalVector[1], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-3)));
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

		Vector3Cartesian pos(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), -static_cast<MML::Real>(1.0));
		Vector3Spherical posSpher = CoordTransfSpherToCart.transf(pos);
		Vector3Cylindrical posCyl = CoordTransfCylToCart.transf(pos);

		auto cart_metric = metricCart(pos);
		auto spher_metric = metricSpher(posSpher);
		auto cyl_metric = metricCyl(posCyl);
		
		// Verify Cartesian metric is identity
		REQUIRE_THAT(cart_metric(0, 0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(cart_metric(1, 1), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(cart_metric(2, 2), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("MetricTensorSphericalContravar - Singularity at r=0 throws", "[MetricTensor][Spherical][Singularity]")
	{
		MetricTensorSphericalContravar metric;
		VectorN<Real, 3> at_origin{static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5)};

		// g^00 = 1 should be fine
		REQUIRE_THAT(metric.component(0, 0, at_origin), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));

		// g^11 = 1/r^2 should throw at r=0
		REQUIRE_THROWS(metric.component(1, 1, at_origin));

		// g^22 = 1/(r^2*sin^2(theta)) should throw at r=0
		REQUIRE_THROWS(metric.component(2, 2, at_origin));
	}

	TEST_CASE("MetricTensorSphericalContravar - Singularity at pole throws", "[MetricTensor][Spherical][Singularity]")
	{
		MetricTensorSphericalContravar metric;
		VectorN<Real, 3> at_pole{static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.5)};  // theta=0 (north pole)

		// g^11 = 1/r^2 should be fine at r=2
		REQUIRE_THAT(metric.component(1, 1, at_pole), WithinAbs(static_cast<MML::Real>(0.25), TOL(1e-12, 1e-5)));

		// g^22 = 1/(r^2*sin^2(0)) should throw at theta=0
		REQUIRE_THROWS(metric.component(2, 2, at_pole));
	}
	class AnalyticalSchwarzschildMetric : public SchwarzschildMetric {
	public:
		using SchwarzschildMetric::SchwarzschildMetric;

		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, 4>& pos) const override {
			ValidateChristoffelArguments(i, j, k, pos);
			const Real r = pos[1], theta = pos[2], rs = _schwarzschildRadius;
			const Real f = 1 - rs / r;
			const Real a = rs / (2 * r * (r - rs));
			if (i == 0 && ((j == 0 && k == 1) || (j == 1 && k == 0)))
				return a;
			if (i == 1 && j == 0 && k == 0)
				return f * rs / (2 * r * r);
			if (i == 1 && j == 1 && k == 1)
				return -a;
			if (i == 1 && j == 2 && k == 2)
				return -r * f;
			if (i == 1 && j == 3 && k == 3)
				return -r * f * std::sin(theta) * std::sin(theta);
			if (i == 2 && ((j == 1 && k == 2) || (j == 2 && k == 1)))
				return 1 / r;
			if (i == 2 && j == 3 && k == 3)
				return -std::sin(theta) * std::cos(theta);
			if (i == 3 && ((j == 1 && k == 3) || (j == 3 && k == 1)))
				return 1 / r;
			if (i == 3 && ((j == 2 && k == 3) || (j == 3 && k == 2)))
				return std::cos(theta) / std::sin(theta);
			return 0;
		}
	};

	template<class Builtin, int N>
	class WarpedBuiltinMetric : public Builtin {
	public:
		mutable int componentCalls = 0;

		Real component(int i, int j, const VectorN<Real, N>& pos) const override {
			++componentCalls;
			return Builtin::component(i, j, pos) * (i == 0 && j == 0 ? 1 + pos[0] * pos[0] : 1);
		}
	};

	template<int N>
	void CheckConnection(const MetricTensorField<N>& metric, const VectorN<Real, N>& pos, const Real (&expected)[N][N][N]) {
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				for (int k = 0; k < N; ++k) {
					INFO("indices " << i << ", " << j << ", " << k << "; position " << pos);
					const Real actual = metric.GetChristoffelSymbolSecondKind(i, j, k, pos);
					REQUIRE_THAT(actual, WithinAbs(expected[i][j][k], TOL(1e-12, 1e-5)));
					if (expected[i][j][k] == 0)
						REQUIRE(actual == 0);
					REQUIRE(actual == metric.GetChristoffelSymbolSecondKind(i, k, j, pos));
					REQUIRE_THAT(metric.GetChristoffelSymbolSecondKindNumerical(i, j, k, pos), WithinAbs(expected[i][j][k], TOL(1e-8, 1e-3)));
				}
	}

	TEST_CASE("Built-in analytical Christoffels cover every coefficient through base references", "[MetricTensor][Christoffel][wewi.8.2]") {
		MetricTensorCartesian3D cartesian;
		MetricTensorMinkowski minkowski;
		MetricTensorSpherical spherical;
		MetricTensorCylindrical cylindrical;
		const Real zero3[3][3][3]{};
		const Real zero4[4][4][4]{};
		for (const auto& pos :
				 {VectorN<Real, 3>{2, static_cast<MML::Real>(0.7), static_cast<MML::Real>(0.3)}, VectorN<Real, 3>{5, static_cast<MML::Real>(1.4), -1}, VectorN<Real, 3>{static_cast<MML::Real>(1.5), static_cast<MML::Real>(2.3), 2}}) {
			CheckConnection<3>(cartesian, pos, zero3);
			CheckConnection<4>(minkowski, VectorN<Real, 4>{pos[0], pos[1], pos[2], -3}, zero4);
			const Real r = pos[0], s = std::sin(pos[1]), c = std::cos(pos[1]);
			Real sphere[3][3][3]{};
			sphere[0][1][1] = -r;
			sphere[0][2][2] = -r * s * s;
			sphere[1][0][1] = sphere[1][1][0] = 1 / r;
			sphere[1][2][2] = -s * c;
			sphere[2][0][2] = sphere[2][2][0] = 1 / r;
			sphere[2][1][2] = sphere[2][2][1] = c / s;
			CheckConnection<3>(spherical, pos, sphere);
			Real cylinder[3][3][3]{};
			cylinder[0][1][1] = -r;
			cylinder[1][0][1] = cylinder[1][1][0] = 1 / r;
			CheckConnection<3>(cylindrical, pos, cylinder);
		}
	}

	template<class Builtin, int N>
	void CheckSubclassFallback(const VectorN<Real, N>& pos) {
		WarpedBuiltinMetric<Builtin, N> warped;
		const MetricTensorField<N>& base = warped;
		const Real expected = pos[0] / (1 + pos[0] * pos[0]);
		REQUIRE_THAT(base.GetChristoffelSymbolSecondKind(0, 0, 0, pos), WithinAbs(expected, TOL(1e-8, 1e-3)));
		REQUIRE(warped.componentCalls > 0);
		REQUIRE(base.GetChristoffelSymbolSecondKind(0, 0, 0, pos) == base.GetChristoffelSymbolSecondKindNumerical(0, 0, 0, pos));
	}

	TEST_CASE("Modified built-in subclasses retain numerical connection fallback", "[MetricTensor][Christoffel][wewi.8.2]") {
		const VectorN<Real, 3> pos{2, static_cast<MML::Real>(0.9), static_cast<MML::Real>(0.4)};
		CheckSubclassFallback<MetricTensorCartesian3D>(pos);
		CheckSubclassFallback<MetricTensorSpherical>(pos);
		CheckSubclassFallback<MetricTensorCylindrical>(pos);
		CheckSubclassFallback<MetricTensorMinkowski>(VectorN<Real, 4>{2, 3, 4, 5});

		MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> transformed(CoordTransfSpherToCart);
		const MetricTensorField<3>& base = transformed;
		MetricTensorSpherical spherical;
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k) {
					const Real actual = base.GetChristoffelSymbolSecondKind(i, j, k, pos);
					REQUIRE(actual == base.GetChristoffelSymbolSecondKindNumerical(i, j, k, pos));
					REQUIRE_THAT(actual, WithinAbs(spherical.GetChristoffelSymbolSecondKind(i, j, k, pos), TOL(1e-8, 1e-3)));
				}
	}

	TEST_CASE("User analytical Schwarzschild connection dispatches into curvature and geodesics", "[MetricTensor][Christoffel][wewi.8.2]") {
		const Real rs = 2;
		AnalyticalSchwarzschildMetric analytical(rs);
		SchwarzschildMetric numerical(rs);
		const MetricTensorField<4>& base = analytical;
		for (const auto& pos :
				 {VectorN<Real, 4>{0, 6, static_cast<MML::Real>(0.8), static_cast<MML::Real>(0.3)}, VectorN<Real, 4>{2, 10, static_cast<MML::Real>(1.4), -1}, VectorN<Real, 4>{-1, 20, static_cast<MML::Real>(2.1), 2}}) {
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					for (int k = 0; k < 4; ++k) {
						const Real gamma = base.GetChristoffelSymbolSecondKind(i, j, k, pos);
						REQUIRE(gamma == base.GetChristoffelSymbolSecondKind(i, k, j, pos));
						REQUIRE_THAT(gamma, WithinAbs(numerical.GetChristoffelSymbolSecondKind(i, j, k, pos), TOL(1e-8, 1e-3)));
						REQUIRE_THAT(gamma, WithinAbs(base.GetChristoffelSymbolSecondKindNumerical(i, j, k, pos), TOL(1e-8, 1e-3)));
					}
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					REQUIRE_THAT(base.GetRicciTensor(i, j, pos), WithinAbs(0, TOL(1e-8, 1e-2)));
			REQUIRE_THAT(base.GetRiemannCurvatureTensor(1, 2, 1, 2, pos), WithinAbs(-rs / (2 * pos[1]), TOL(1e-8, 1e-2)));

			GeodesicEquationSystem<4> system(base);
			Vector<Real> rhs(8);
			system.derivs(0, MakeGeodesicState<4>(pos, VectorN<Real, 4>{1, 0, 0, 0}), rhs);
			REQUIRE_THAT(rhs[5], WithinAbs(-(1 - rs / pos[1]) * rs / (2 * pos[1] * pos[1]), TOL(1e-12, 1e-5)));
		}
	}

	TEST_CASE("Spherical analytical connection integrates a Cartesian straight line", "[MetricTensor][Christoffel][Geodesic][wewi.8.2]") {
		MetricTensorSpherical metric;
		const VectorN<Real, 3> pos{2, static_cast<MML::Real>(1.1), static_cast<MML::Real>(0.3)}, velocity{static_cast<MML::Real>(0.2), static_cast<MML::Real>(-0.1), static_cast<MML::Real>(0.15)};
		const Real end = 0.4;
		const auto solution = IntegrateGeodesicFixedStep<3>(metric, pos, velocity, 0, end, 80);
		const auto state = solution.getXValuesAtEnd();
		const auto actual = CoordTransfSpherToCart.transf(GeodesicPositionFromState<3>(state));
		const auto expected = CoordTransfSpherToCart.transf(pos) + end * (CoordTransfSpherToCart.jacobian(pos) * velocity);
		for (int i = 0; i < 3; ++i)
			REQUIRE_THAT(actual[i], WithinAbs(expected[i], TOL(1e-8, 1e-3)));
	}

	TEST_CASE("Analytical connections reject singular positions even for structural zeros",
						"[MetricTensor][Christoffel][Singularity][wewi.8.2]") {
		MetricTensorSpherical spherical;
		MetricTensorCylindrical cylindrical;
		for (const auto& pos : {VectorN<Real, 3>{0, 1, 0}, VectorN<Real, 3>{2, 0, 0}, VectorN<Real, 3>{2, Constants::PI, 0}})
			REQUIRE_THROWS_AS(spherical.GetChristoffelSymbolSecondKind(0, 0, 0, pos), DomainError);
		REQUIRE_THROWS_AS(cylindrical.GetChristoffelSymbolSecondKind(2, 2, 2, VectorN<Real, 3>{0, 1, 2}), DomainError);
		REQUIRE_THROWS(spherical.GetChristoffelSymbolSecondKindNumerical(0, 0, 0, VectorN<Real, 3>{0, 1, 0}));
	}

	TEST_CASE("Spherical analytical connection parallel transports a constant Cartesian vector",
						"[MetricTensor][Christoffel][ParallelTransport][wewi.8.2]") {
		MetricTensorSpherical metric;
		ParametricCurveFromStdFunc<3> equator([](Real t) { return VectorN<Real, 3>{2, Constants::PI / 2, t}; });
		const Real end = 0.6;
		const auto solution = IntegrateParallelTransportFixedStep<3>(metric, equator, VectorN<Real, 3>{1, 0, 0}, 0, end, 80);
		const auto actual = ParallelTransportVectorFromState<3>(solution.getXValuesAtEnd());
		const VectorN<Real, 3> expected{std::cos(end), 0, -std::sin(end) / 2};
		for (int i = 0; i < 3; ++i)
			REQUIRE_THAT(actual[i], WithinAbs(expected[i], TOL(1e-8, 1e-3)));
	}

	TEST_CASE("Christoffel argument validation and inverse-storage rejection are explicit", "[MetricTensor][Christoffel][wewi.8.2]") {
		MetricTensorCartesian3D cartesian;
		MetricTensorSpherical spherical;
		MetricTensorCylindrical cylindrical;
		const VectorN<Real, 3> pos{2, 1, 3};
		const MetricTensorField<3>* metrics[]{&cartesian, &spherical, &cylindrical};
		for (const auto* metric : metrics) {
			REQUIRE_THROWS_AS(metric->GetChristoffelSymbolSecondKind(-1, 0, 0, pos), IndexError);
			REQUIRE_THROWS_AS(metric->GetChristoffelSymbolSecondKind(0, 3, 0, pos), IndexError);
			REQUIRE_THROWS_AS(metric->GetChristoffelSymbolSecondKindNumerical(0, 0, 3, pos), IndexError);
			REQUIRE_THROWS_AS(metric->GetChristoffelSymbolFirstKind(0, 0, 3, pos), IndexError);
			for (Real bad : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
				for (int axis = 0; axis < 3; ++axis) {
					auto invalid = pos;
					invalid[axis] = bad;
					REQUIRE_THROWS_AS(metric->GetChristoffelSymbolSecondKind(0, 0, 0, invalid), DomainError);
					REQUIRE_THROWS_AS(metric->GetChristoffelSymbolSecondKindNumerical(0, 0, 0, invalid), DomainError);
				}
		}
		MetricTensorSphericalContravar inverse;
		REQUIRE_THROWS_AS(inverse.GetChristoffelSymbolSecondKind(0, 1, 1, pos), DomainError);
		REQUIRE_THROWS_AS(inverse.GetChristoffelSymbolSecondKindNumerical(0, 1, 1, pos), DomainError);
		REQUIRE_THROWS_AS(inverse.GetChristoffelSymbolFirstKind(0, 1, 1, pos), DomainError);
		REQUIRE_THAT(inverse.component(1, 1, pos), WithinAbs(0.25, TOL(1e-12, 1e-5)));

		MetricTensorMinkowski minkowski;
		const VectorN<Real, 4> pos4{0, 1, 2, 3};
		REQUIRE_THROWS_AS(minkowski.GetChristoffelSymbolSecondKind(0, 0, 4, pos4), IndexError);
		auto invalid = pos4;
		invalid[3] = std::numeric_limits<Real>::quiet_NaN();
		REQUIRE_THROWS_AS(minkowski.GetChristoffelSymbolSecondKind(0, 0, 0, invalid), DomainError);
		REQUIRE_THROWS_AS(minkowski.GetChristoffelSymbolSecondKindNumerical(0, 0, 0, invalid), DomainError);
	}

	TEST_CASE("Numerical connection rejects nonfinite metric data instead of returning NaN", "[MetricTensor][Christoffel][wewi.8.2]") {
		class NonfiniteMetric : public MetricTensorCartesian3D {
		public:
			Real component(int, int, const VectorN<Real, 3>&) const override { return std::numeric_limits<Real>::quiet_NaN(); }
		};
		NonfiniteMetric metric;
		REQUIRE_THROWS_AS(metric.GetChristoffelSymbolSecondKind(0, 0, 0, VectorN<Real, 3>{1, 2, 3}), DomainError);
		REQUIRE_THROWS_AS(metric.GetChristoffelSymbolSecondKindNumerical(0, 0, 0, VectorN<Real, 3>{1, 2, 3}), DomainError);

		class NonfiniteStencilMetric : public NonfiniteMetric {
		public:
			Real component(int i, int j, const VectorN<Real, 3>& pos) const override {
				return pos[0] == 0 ? MetricTensorCartesian3D::component(i, j, pos) : NonfiniteMetric::component(i, j, pos);
			}
		};
		NonfiniteStencilMetric stencil;
		REQUIRE_THROWS_AS(stencil.GetChristoffelSymbolSecondKindNumerical(0, 0, 0, VectorN<Real, 3>{0, 1, 2}), DomainError);
	}
} // namespace MML::Tests::Core::MetricTensorTests
