#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#include <mml/core/Fields/FieldOperations.h>
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Core::SingularityPolicyTests {
	template<class Matrix>
	void CheckSentinel(const Matrix& matrix, SingularityPolicy policy) {
		for (int i = 0; i < matrix.rows(); ++i)
			for (int j = 0; j < matrix.cols(); ++j) {
				if (policy == SingularityPolicy::ReturnNaN)
					REQUIRE(std::isnan(matrix(i, j)));
				else if (policy == SingularityPolicy::ReturnInf)
					REQUIRE(matrix(i, j) == std::numeric_limits<Real>::infinity());
				else
					REQUIRE(matrix(i, j) == 0);
			}
	}

	template<class Transform>
	void CheckSingularTransform(const Transform& transform, const VectorN<Real, 3>& pos) {
		REQUIRE_THROWS_AS(transform.inverseJacobian(pos), SingularMatrixError);
		REQUIRE_THROWS_AS(transform.inverseJacobian(pos, SingularityPolicy::Throw), SingularMatrixError);
		REQUIRE_THROWS_AS(transform.inverseJacobian(pos, SingularityPolicy::Clamp), DomainError);
		for (auto policy : {SingularityPolicy::ReturnNaN, SingularityPolicy::ReturnInf, SingularityPolicy::ReturnZero})
			CheckSentinel(transform.inverseJacobian(pos, policy), policy);
	}

	TEST_CASE("Coordinate inverse Jacobians have explicit policies at chart singularities", "[coordinates][Singularity][wewi.8.3]") {
		CoordTransfSphericalToCartesian spherical;
		CoordTransfCartesianToSpherical inverseSpherical;
		CoordTransfCylindricalToCartesian cylindrical;
		CoordTransfCartesianToCylindrical inverseCylindrical;
		for (const auto& pos : {VectorN<Real, 3>{0, 0, 0}, VectorN<Real, 3>{0, 0, 2}, VectorN<Real, 3>{0, 0, -2}})
			CheckSingularTransform(spherical, pos);
		for (const auto& pos : {VectorN<Real, 3>{0, 1, 0}, VectorN<Real, 3>{2, 0, 0}, VectorN<Real, 3>{2, Constants::PI, 0}})
			CheckSingularTransform(inverseSpherical, pos);
		CheckSingularTransform(cylindrical, VectorN<Real, 3>{0, 0, 2});
		CheckSingularTransform(inverseCylindrical, VectorN<Real, 3>{0, 1, 2});
		CheckSingularTransform(spherical, VectorN<Real, 3>{Singularity::DEFAULT_SINGULARITY_TOL / 2, 0, 0});
		CheckSingularTransform(cylindrical, VectorN<Real, 3>{Singularity::DEFAULT_SINGULARITY_TOL / 2, 0, 1});
		REQUIRE_THROWS_AS(inverseSpherical.jacobian(VectorN<Real, 3>{0, 0, 2}), SingularMatrixError);
		REQUIRE_THROWS_AS(inverseCylindrical.jacobian(VectorN<Real, 3>{0, 0, 2}), SingularMatrixError);
	}

	TEST_CASE("Policy overload preserves virtual inverse Jacobian overrides and numerical fallback", "[coordinates][Singularity][wewi.8.3]") {
		class LegacyOverride : public CoordTransfCartesianToSpherical {
		public:
			mutable int calls = 0;
			MatrixNM<Real, 3, 3> inverseJacobian(const VectorN<Real, 3>&) const override {
				++calls;
				MatrixNM<Real, 3, 3> result;
				for (int i = 0; i < 3; ++i)
					result(i, i) = 2;
				return result;
			}
		};
		LegacyOverride legacy;
		const CoordTransfWithInverse<Vector3Cartesian, Vector3Spherical, 3>& base = legacy;
		REQUIRE(base.inverseJacobian(VectorN<Real, 3>{2, 1, 0}, SingularityPolicy::ReturnNaN)(0, 0) == 2);
		REQUIRE(legacy.calls == 1);

		class NumericalInverse : public CoordTransfCartesianToSpherical {
		public:
			using CoordTransfCartesianToSpherical::inverseJacobian;
			MatrixNM<Real, 3, 3> inverseJacobian(const VectorN<Real, 3>& pos) const override {
				return CoordTransfWithInverse<Vector3Cartesian, Vector3Spherical, 3>::inverseJacobian(pos);
			}
		};
		NumericalInverse numerical;
		CheckSingularTransform(numerical, VectorN<Real, 3>{0, 1, 0});
	}

	TEST_CASE("Regular inverse Jacobians remain correct for every policy and reject nonfinite input",
						"[coordinates][Singularity][wewi.8.3]") {
		const VectorN<Real, 3> pos{2, 1, static_cast<MML::Real>(0.4)};
		const auto cart = CoordTransfSpherToCart.transf(pos);
		const auto jac = CoordTransfSpherToCart.jacobian(pos);
		for (auto policy : {SingularityPolicy::Throw, SingularityPolicy::ReturnNaN, SingularityPolicy::ReturnInf, SingularityPolicy::ReturnZero,
												SingularityPolicy::Clamp}) {
			const auto product = CoordTransfSpherToCart.inverseJacobian(cart, policy) * jac;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					REQUIRE_THAT(product(i, j), WithinAbs(i == j ? 1 : 0, TOL(1e-12, 1e-5)));
			auto bad = cart;
			bad[0] = std::numeric_limits<Real>::quiet_NaN();
			REQUIRE_THROWS_AS(CoordTransfSpherToCart.inverseJacobian(bad, policy), DomainError);
		}
	}

	class ScaledMetric : public MetricTensorField<3> {
		Real _scale;
		bool _vary;

	public:
		ScaledMetric(Real scale, bool vary = false)
				: _scale(scale)
				, _vary(vary) {}
		Real component(int i, int j, const VectorN<Real, 3>& pos) const override {
			return i == j ? _scale * (_vary ? std::exp(2 * pos[0]) : 1) : 0;
		}
	};

	TEST_CASE("Metric inversion applies the shared matrix policy without regularization", "[MetricTensor][Singularity][wewi.8.3]") {
		MetricTensorSpherical metric;
		const VectorN<Real, 3> origin{0, 1, 0};
		REQUIRE_THROWS_AS(metric.GetContravariantMetric(origin), SingularMatrixError);
		REQUIRE_THROWS_AS(metric.GetContravariantMetric(origin, SingularityPolicy::Clamp), DomainError);
		for (auto policy : {SingularityPolicy::ReturnNaN, SingularityPolicy::ReturnInf, SingularityPolicy::ReturnZero})
			CheckSentinel(metric.GetContravariantMetric(origin, policy), policy);
		MetricTensorMinkowski minkowski;
		REQUIRE(minkowski.GetContravariantMetric(VectorN<Real, 4>{0, 0, 0, 0}, SingularityPolicy::ReturnNaN)(0, 0) == -1);
		for (Real scale : {static_cast<MML::Real>(1e-20), static_cast<MML::Real>(1e20)}) {
			ScaledMetric regular(scale);
			const auto product = regular.GetCovariantMetric(origin) * regular.GetContravariantMetric(origin);
			for (int i = 0; i < 3; ++i)
				REQUIRE_THAT(product(i, i), WithinAbs(1, TOL(1e-12, 1e-5)));
		}
		REQUIRE_THROWS_AS(metric.GetContravariantMetric(VectorN<Real, 3>{2, 1, 0}, static_cast<SingularityPolicy>(-1)), ArgumentError);
	}

	TEST_CASE("Degenerate determinant divergence throws or reports an unusable result", "[field_operations][Singularity][wewi.8.3]") {
		MetricTensorSpherical metric;
		VectorFunction<3> field([](const VectorN<Real, 3>& p) { return p; });
		FieldOperationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		for (const auto& pos : {VectorN<Real, 3>{0, 1, 0}, VectorN<Real, 3>{2, 0, 0}, VectorN<Real, 3>{2, Constants::PI, 0}}) {
			REQUIRE_THROWS_AS(VectorFieldOperations::DivergenceViaDet(field, pos, metric), DomainError);
			const auto result = VectorFieldOperations::DivergenceViaDetDetailed(field, pos, metric, config);
			REQUIRE_FALSE(result.IsSuccess());
			REQUIRE(result.status == AlgorithmStatus::InvalidInput);
			REQUIRE_FALSE(result.error_message.empty());
			REQUIRE(result.algorithm_name == "DivergenceViaDet");
		}
		config.exception_policy = EvaluationExceptionPolicy::Propagate;
		REQUIRE_THROWS_AS(VectorFieldOperations::DivergenceViaDetDetailed(field, VectorN<Real, 3>{0, 1, 0}, metric, config), DomainError);
	}

	TEST_CASE("Determinant divergence is invariant under tiny and huge constant metric scaling",
						"[field_operations][Singularity][wewi.8.3]") {
		VectorFunction<3> radial([](const VectorN<Real, 3>& p) { return p; });
		VectorFunction<3> constant([](const VectorN<Real, 3>&) { return VectorN<Real, 3>{1, 0, 0}; });
		const VectorN<Real, 3> pos{static_cast<MML::Real>(0.3), 1, 2};
		for (Real scale : {static_cast<MML::Real>(1e-20), static_cast<MML::Real>(1), static_cast<MML::Real>(1e20)}) {
			ScaledMetric flat(scale);
			ScaledMetric conformal(scale, true);
			REQUIRE_THAT(VectorFieldOperations::DivergenceViaDet(radial, pos, flat), WithinAbs(3, TOL(1e-8, 1e-3)));
			REQUIRE_THAT(VectorFieldOperations::DivergenceViaDet(constant, pos, conformal), WithinAbs(3, TOL(1e-7, 1e-2)));
		}
		for (Real scale : {std::sqrt(std::numeric_limits<Real>::min()), std::sqrt(std::numeric_limits<Real>::max())}) {
			ScaledMetric conformal(scale, true);
			REQUIRE_THAT(VectorFieldOperations::DivergenceViaDet(constant, pos, conformal), WithinAbs(3, TOL(1e-7, 1e-2)));
		}
		MetricTensorSpherical spherical;
		VectorFunction<3> radialSpherical([](const VectorN<Real, 3>& p) { return VectorN<Real, 3>{p[0], 0, 0}; });
		const Real smallRadius = 100 * std::sqrt(Precision::Eps);
		REQUIRE_THAT(VectorFieldOperations::DivergenceViaDet(radialSpherical, VectorN<Real, 3>{smallRadius, 1, 0}, spherical),
								 WithinAbs(3, TOL(1e-6, 1e-2)));
	}

	TEST_CASE("Detailed determinant divergence respects derivative configuration and counts actual evaluations",
						"[field_operations][Singularity][wewi.8.3]") {
		MetricTensorCartesian3D metric;
		int calls = 0;
		VectorFunctionFromStdFunc<3> field([&](const VectorN<Real, 3>& p) {
			++calls;
			return p;
		});
		for (int order : {0, 1, 2, 4, 6, 8}) {
			FieldOperationConfig config;
			config.derivative_order = order;
			config.estimate_error = false;
			calls = 0;
			const auto result = VectorFieldOperations::DivergenceViaDetDetailed(field, VectorN<Real, 3>{1, 2, 3}, metric, config);
			REQUIRE(result.IsSuccess());
			REQUIRE_THAT(result.value, WithinAbs(3, TOL(1e-6, 1e-2)));
			REQUIRE(result.function_evaluations == calls);
			REQUIRE(calls > 0);
		}
		FieldOperationConfig config;
		calls = 0;
		const auto estimated = VectorFieldOperations::DivergenceViaDetDetailed(field, VectorN<Real, 3>{1, 2, 3}, metric, config);
		REQUIRE(estimated.IsSuccess());
		REQUIRE(estimated.function_evaluations == calls);
		REQUIRE(std::isfinite(estimated.error));
		REQUIRE(estimated.error >= 0);
	}

	TEST_CASE("Signed and nondiagonal metrics use absolute determinant and rank checks",
						"[MetricTensor][field_operations][Singularity][wewi.8.3]") {
		class SwapMetric : public LorentzianMetric<2> {
		public:
			Real component(int i, int j, const VectorN<Real, 2>&) const override { return i == j ? Real(0) : Real(1); }
		};
		SwapMetric signedMetric;
		VectorFunction<2> radial([](const VectorN<Real, 2>& p) { return p; });
		const VectorN<Real, 2> pos{1, 2};
		REQUIRE_THAT(VectorFieldOperations::DivergenceViaDet(radial, pos, signedMetric), WithinAbs(2, TOL(1e-8, 1e-3)));
		class RankOneMetric : public MetricTensorField<2> {
		public:
			Real component(int i, int j, const VectorN<Real, 2>&) const override { return Real((i + 1) * (j + 1)); }
		};
		RankOneMetric rankOne;
		REQUIRE_THROWS_AS(rankOne.GetContravariantMetric(pos), SingularMatrixError);
		CheckSentinel(rankOne.GetContravariantMetric(pos, SingularityPolicy::ReturnNaN), SingularityPolicy::ReturnNaN);
		REQUIRE_THROWS_AS(VectorFieldOperations::DivergenceViaDet(radial, pos, rankOne), DomainError);
	}

	TEST_CASE("Invalid metric and stencil data never become successful determinant divergence",
						"[MetricTensor][field_operations][Singularity][wewi.8.3]") {
		class NonfiniteMetric : public MetricTensorField<3> {
		public:
			Real component(int, int, const VectorN<Real, 3>&) const override { return std::numeric_limits<Real>::quiet_NaN(); }
		};
		NonfiniteMetric bad;
		const VectorN<Real, 3> pos{0, 1, 2};
		VectorFunction<3> field([](const VectorN<Real, 3>& p) { return p; });
		FieldOperationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		REQUIRE_THROWS_AS(bad.GetContravariantMetric(pos, SingularityPolicy::ReturnZero), DomainError);
		const auto invalid = VectorFieldOperations::DivergenceViaDetDetailed(field, pos, bad, config);
		REQUIRE(invalid.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(invalid.IsSuccess());
		class DegenerateStencilMetric : public MetricTensorCartesian3D {
		public:
			Real component(int i, int j, const VectorN<Real, 3>& p) const override {
				return p[0] == 0 ? MetricTensorCartesian3D::component(i, j, p) : 0;
			}
		};
		DegenerateStencilMetric stencil;
		const auto failedStencil = VectorFieldOperations::DivergenceViaDetDetailed(field, pos, stencil, config);
		REQUIRE(failedStencil.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(failedStencil.IsSuccess());
		MetricTensorCartesian3D flat;
		VectorFunction<3> nanField([](const VectorN<Real, 3>&) { return VectorN<Real, 3>{std::numeric_limits<Real>::quiet_NaN(), 0, 0}; });
		for (bool check : {true, false}) {
			config.check_finite = check;
			const auto nonfinite = VectorFieldOperations::DivergenceViaDetDetailed(nanField, pos, flat, config);
			REQUIRE(nonfinite.status == AlgorithmStatus::NumericalInstability);
			REQUIRE_FALSE(nonfinite.IsSuccess());
		}
		class ExtremeStencilMetric : public ScaledMetric {
		public:
			using ScaledMetric::ScaledMetric;
			Real component(int i, int j, const VectorN<Real, 3>& p) const override {
				return p[0] == 0 ? (i == j ? 1 : 0) : ScaledMetric::component(i, j, p);
			}
		};
		for (Real scale : {std::numeric_limits<Real>::min(), std::numeric_limits<Real>::max()}) {
			ExtremeStencilMetric extreme(scale);
			const auto failure = VectorFieldOperations::DivergenceViaDetDetailed(field, pos, extreme, config);
			REQUIRE(failure.status == AlgorithmStatus::NumericalInstability);
			REQUIRE_FALSE(failure.IsSuccess());
		}
		config.derivative_order = 5;
		REQUIRE(VectorFieldOperations::DivergenceViaDetDetailed(field, pos, flat, config).status == AlgorithmStatus::InvalidInput);
	}

	TEST_CASE("Inverse Jacobians avoid squared-radius overflow and do not mask invalid coordinates", "[coordinates][Singularity][wewi.8.3]") {
		const Real magnitude = std::numeric_limits<Real>::max() / 10;
		const VectorN<Real, 3> large{magnitude, 2 * magnitude, 3 * magnitude};
		const auto spherical = CoordTransfSpherToCart.inverseJacobian(large, SingularityPolicy::Throw);
		const auto cylindrical = CoordTransfCylToCart.inverseJacobian(large, SingularityPolicy::Throw);
		REQUIRE_THAT(spherical(0, 0), WithinAbs(1 / std::sqrt(Real(14)), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(cylindrical(0, 0), WithinAbs(1 / std::sqrt(Real(5)), TOL(1e-12, 1e-5)));
		REQUIRE(spherical(2, 0) < 0);
		REQUIRE(cylindrical(1, 0) < 0);
		const Real max = std::numeric_limits<Real>::max();
		REQUIRE_THROWS_AS(CoordTransfSpherToCart.inverseJacobian(VectorN<Real, 3>{max, max, max}, SingularityPolicy::ReturnZero), DomainError);
		REQUIRE_THROWS_AS(CoordTransfCylToCart.inverseJacobian(VectorN<Real, 3>{max, max, 0}, SingularityPolicy::ReturnZero), DomainError);
		MetricTensorCartesian3D metric;
		REQUIRE_THROWS_AS(
			metric.GetContravariantMetric(VectorN<Real, 3>{max, std::numeric_limits<Real>::infinity(), 0}, SingularityPolicy::ReturnZero),
			DomainError);
	}
} // namespace MML::Tests::Core::SingularityPolicyTests
