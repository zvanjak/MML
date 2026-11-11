#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/ComplexFunction.h>
#include <mml/base/ChebyshevApproximation.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionBarycentric.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionPolynomial.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionRational.h>
#include <mml/core/Derivation.h>
#include <mml/core/Fields/ScalarFieldOperations.h>
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/algorithms/RootFinding.h>
#include <mml/algorithms/Optimization/Constraints/ProjectedGradient.h>
#endif

#include <complex>
#include <optional>
#include <string>
#include <type_traits>

using namespace MML;

namespace {
	template<class Result>
	void CheckUncomputed() {
		const Result result;
		CHECK_FALSE(result.IsSuccess());
		CHECK_FALSE(static_cast<bool>(result));
		CHECK_FALSE(IsUsable(result.status));
		CHECK(to_string(result.status) == "NotComputed");
		CHECK(result.algorithm_name.empty());
		CHECK(result.error_message.empty());
		CHECK(result.function_evaluations == 0);
	}

	template<class Scalar>
	void CheckScalarEstimate() {
		EvaluationResult<Scalar> result;
		CHECK((std::is_same_v<decltype(result.error), Scalar>));
		CHECK(std::isnan(result.error));
		CHECK_FALSE(HasEstimate(result.error));
		result.error = Scalar(0);
		CHECK(HasEstimate(result.error));
		result.error = std::numeric_limits<Scalar>::infinity();
		CHECK(HasEstimate(result.error));
		auto success = MakeEvaluationSuccessResult<EvaluationResult<Scalar>>("Computed");
		CHECK(success.IsSuccess());
		CHECK_FALSE(HasEstimate(success.error));
		auto failure = MakeEvaluationFailureResult<EvaluationResult<Scalar>>(
			AlgorithmStatus::NumericalInstability, "No evaluation", "Computed");
		CHECK_FALSE(failure.IsSuccess());
		CHECK_FALSE(HasEstimate(failure.error));
	}

	struct CustomError {
		int bound = 0;
		bool operator==(const CustomError& other) const { return bound == other.bound; }
	};

	template<class Error>
	void CheckOptionalEstimate() {
		EvaluationResult<Real, Error> result;
		CHECK((std::is_same_v<decltype(result.error), std::optional<Error>>));
		if constexpr (std::is_same_v<decltype(result.error), std::optional<Error>>) {
			CHECK_FALSE(result.error.has_value());
			CHECK_FALSE(HasEstimate(result.error));
			result.error = Error{};
			CHECK(result.error.has_value());
			CHECK(HasEstimate(result.error));
			CHECK(*result.error == Error{});
			result.error.reset();
			CHECK_FALSE(result.error.has_value());
		}
	}
}

TEST_CASE("Default shared result bases and inherited results cannot claim computed success",
	"[result-policy][uncomputed-result-contract]") {
	CheckUncomputed<IterativeResultBase>();
	CheckUncomputed<EvaluationResultBase>();
	CheckUncomputed<EvaluationResult<Real>>();
	CheckUncomputed<DerivativeResult<Real>>();
	CheckUncomputed<RootFinding::RootFindingResult>();
	CheckUncomputed<RootFinding::PolynomialRootSetResult>();
	CheckUncomputed<InterpolationResult>();
	CheckUncomputed<ChebyshevApproximationResult>();
	CheckUncomputed<LinearSolverResult<Real>>();
	CheckUncomputed<Optimization::ProjectedGradientResult<2>>();
	RootFinding::RootFindingResult root;
	CHECK_FALSE(HasEstimate(root.x_error));
	InterpolationResult interpolation;
	CHECK_FALSE(HasEstimate(interpolation.error_estimate));
	ChebyshevApproximationResult chebyshev;
	CHECK_FALSE(HasEstimate(chebyshev.coefficient_tail_error));
	CHECK_FALSE(HasEstimate(chebyshev.max_error_estimate));
	CHECK_FALSE(HasEstimate(chebyshev.rms_error_estimate));
	RootFinding::PolynomialRootSetResult polynomial;
	CHECK_FALSE(HasEstimate(polynomial.max_residual));
	CHECK_FALSE(HasEstimate(polynomial.max_relative_residual));
	Optimization::ProjectedGradientResult<2> projected;
	CHECK_FALSE(HasEstimate(projected.projected_gradient_norm));
	CHECK_FALSE(HasEstimate(projected.step_norm));
}

TEST_CASE("Generic scalar estimates retain their native type and distinguish unknown from computed zero",
	"[result-policy][uncomputed-result-contract]") {
	CheckScalarEstimate<float>();
	CheckScalarEstimate<double>();
	CheckScalarEstimate<long double>();
}

TEST_CASE("Generic fixed-shape and complex estimates do not invent zero components",
	"[result-policy][uncomputed-result-contract]") {
	EvaluationResult<Complex> complex;
	CHECK((std::is_same_v<decltype(complex.error), Complex>));
	CHECK(std::isnan(complex.error.real()));
	CHECK(std::isnan(complex.error.imag()));
	EvaluationResult<VectorN<float, 3>> vector;
	CHECK((std::is_same_v<decltype(vector.error), VectorN<float, 3>>));
	for (int i = 0; i < 3; ++i) CHECK(std::isnan(vector.error[i]));
	EvaluationResult<MatrixNM<double, 2, 3>> matrix;
	CHECK((std::is_same_v<decltype(matrix.error), MatrixNM<double, 2, 3>>));
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 3; ++j) CHECK(std::isnan(matrix.error(i, j)));
}

TEST_CASE("Non-NaN generic error types preserve absence separately from a computed zero or empty bound",
	"[result-policy][uncomputed-result-contract]") {
	CheckOptionalEstimate<int>();
	CheckOptionalEstimate<std::string>();
	CheckOptionalEstimate<CustomError>();
	CheckOptionalEstimate<VectorN<int, 2>>();
}

TEST_CASE("Explicit completed and checked-success factories remain usable without inventing estimates",
	"[result-policy][uncomputed-result-contract]") {
	auto checked = MakeSuccessResult<IterativeResultBase>(2, Real(0), "Checked");
	CHECK(checked.IsSuccess());
	CHECK(checked.converged);
	CHECK(HasEstimate(checked.achieved_tolerance));
	auto done = MakeCompletedResult<IterativeResultBase>("Direct");
	CHECK(done.IsSuccess());
	CHECK_FALSE(done.converged);
	CHECK_FALSE(HasEstimate(done.achieved_tolerance));
	auto evaluation = MakeEvaluationSuccessResult<EvaluationResult<Real>>("Evaluation", 3);
	CHECK(evaluation.IsSuccess());
	CHECK(evaluation.function_evaluations == 3);
	CHECK_FALSE(HasEstimate(evaluation.error));
}

TEST_CASE("Derivative disabled and failed estimates are unknown while computed estimates remain available",
	"[result-policy][uncomputed-result-contract]") {
	RealFunction function([](Real x) { return x * x; });
	DerivativeConfig config;
	config.estimate_error = false;
	auto disabled = Derivation::NDer4Detailed(function, 2, config);
	CHECK(disabled.IsSuccess());
	const Real roundingTolerance = std::max(Precision::DefaultToleranceStrict,
		Real(8) * std::numeric_limits<Real>::epsilon() * Real(4));
	CHECK_THAT(disabled.value, Catch::Matchers::WithinAbs(Real(4), roundingTolerance));
	CHECK_FALSE(HasEstimate(disabled.error));
	config.estimate_error = true;
	auto estimated = Derivation::NDer4Detailed(function, 2, config);
	CHECK(estimated.IsSuccess());
	CHECK(HasEstimate(estimated.error));
	CHECK(estimated.error >= 0);
	RealFunction zero([](Real) { return Real(0); });
	auto exact = Derivation::NDer4Detailed(zero, 2, config);
	CHECK(exact.error == 0);
	CHECK(HasEstimate(exact.error));
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	auto invalid = Derivation::NDer4Detailed(function, 2, Real(0), config);
	CHECK(invalid.status == AlgorithmStatus::InvalidInput);
	CHECK_FALSE(HasEstimate(invalid.error));
	CHECK_FALSE(invalid.error_message.empty());
	ComplexFunction complex([](Complex z) { return z * z; });
	config.estimate_error = false;
	auto complexDisabled = Derivation::NDer2ComplexDetailed(complex, Complex(2, 1), config);
	CHECK(complexDisabled.IsSuccess());
	CHECK_FALSE(HasEstimate(complexDisabled.error));
}

TEST_CASE("Chebyshev construction earns Completed or Success only after actually building its approximation",
	"[result-policy][uncomputed-result-contract]") {
	auto function = [](Real x) { return x * x; };
	auto fixed = ApproximateChebyshevDetailed(function, -1, 1, 8, 17);
	CHECK(fixed.IsSuccess());
	CHECK(fixed.status == AlgorithmStatus::Completed);
	CHECK_FALSE(fixed.converged);
	CHECK_THAT(fixed.approximation(Real(0.5)), Catch::Matchers::WithinAbs(Real(0.25), Precision::DefaultToleranceStrict));
	CHECK(HasEstimate(fixed.max_error_estimate));
	CHECK(fixed.max_error_estimate >= 0);
	ChebyshevApproximationConfig config;
	config.tolerance = Precision::DefaultToleranceStrict;
	auto adaptive = ApproximateChebyshevAdaptive(function, -1, 1, config);
	CHECK(adaptive.IsSuccess());
	CHECK(adaptive.converged);
	CHECK(adaptive.status == AlgorithmStatus::Success);
	CHECK(adaptive.max_error_estimate <= config.tolerance);
	CHECK(adaptive.coefficient_tail_error <= config.tolerance);
}

TEST_CASE("Interpolation completed evaluations stay usable without an invented default error estimate",
	"[result-policy][uncomputed-result-contract]") {
	Vector<Real> x{0, 1, 2}, y{0, 1, 4};
	BarycentricPolynomialInterp polynomial(x, y);
	auto result = polynomial.EvaluateDetailed(Real(0.5));
	CHECK(result.IsSuccess());
	CHECK(result.status == AlgorithmStatus::Completed);
	CHECK_FALSE(HasEstimate(result.error_estimate));
	CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(0.25), Precision::DefaultToleranceStrict));
	InterpolationConfig config;
	config.estimate_error = false;
	auto disabled = polynomial.EvaluateDetailed(Real(0.5), config);
	CHECK(disabled.IsSuccess());
	CHECK_FALSE(HasEstimate(disabled.error_estimate));
}

TEST_CASE("Composite estimate presence is type-correct and independent of result status",
	"[result-policy][uncomputed-result-contract]") {
	EvaluationResult<Vec3Cart> cartesian;
	EvaluationResult<Vec3Sph> spherical;
	EvaluationResult<Vec3Cyl> cylindrical;
	STATIC_REQUIRE((std::is_same_v<decltype(cartesian.error), Vec3Cart>));
	STATIC_REQUIRE((std::is_same_v<decltype(spherical.error), Vec3Sph>));
	STATIC_REQUIRE((std::is_same_v<decltype(cylindrical.error), Vec3Cyl>));
	CHECK_FALSE(HasEstimate(cartesian.error));
	CHECK_FALSE(HasEstimate(spherical.error));
	CHECK_FALSE(HasEstimate(cylindrical.error));
	for (int i = 0; i < 3; ++i) {
		CHECK(std::isnan(cartesian.error[i]));
		CHECK(std::isnan(spherical.error[i]));
		CHECK(std::isnan(cylindrical.error[i]));
	}
	cartesian.error = Vec3Cart{};
	CHECK(HasEstimate(cartesian.error));
	CHECK_FALSE(cartesian.IsSuccess());
	cartesian.error[1] = NoEstimate;
	CHECK_FALSE(HasEstimate(cartesian.error));
	EvaluationResult<Complex> complexResult;
	CHECK_FALSE(HasEstimate(complexResult.error));
	complexResult.error = Complex(0, 0);
	CHECK(HasEstimate(complexResult.error));
	complexResult.error = Complex(0, NoEstimate);
	CHECK_FALSE(HasEstimate(complexResult.error));
	EvaluationResult<MatrixNM<Real, 2, 3>> matrix;
	CHECK_FALSE(HasEstimate(matrix.error));
	matrix.error = MatrixNM<Real, 2, 3>{};
	CHECK(HasEstimate(matrix.error));
	matrix.error(1, 2) = NoEstimate;
	CHECK_FALSE(HasEstimate(matrix.error));
	const auto evaluated = MakeEvaluationSuccessResult<EvaluationResult<Vec3Cart>>("ExplicitEvaluation", 1);
	CHECK(evaluated.IsSuccess());
	CHECK_FALSE(HasEstimate(evaluated.error));
}

TEST_CASE("Optional estimates do not require default-constructible error payloads",
	"[result-policy][uncomputed-result-contract]") {
	struct Bound {
		int magnitude;
		explicit Bound(int value) : magnitude(value) {}
	};
	EvaluationResult<Real, Bound> result;
	CHECK_FALSE(HasEstimate(result.error));
	result.error.emplace(0);
	CHECK(HasEstimate(result.error));
	CHECK(result.error->magnitude == 0);
	CHECK_FALSE(result.IsSuccess());
	result.error.reset();
	CHECK_FALSE(HasEstimate(result.error));
}

TEST_CASE("Every interpolation populator reports completion and preserves available versus disabled errors",
	"[result-policy][uncomputed-result-contract]") {
	Vector<Real> x{0, 1, 2}, y{3, 5, 7};
	LinearInterpRealFunc linear(x, y);
	PolynomInterpRealFunc polynomial(x, y, 3);
	RationalInterpRealFunc rational(x, y, 3);
	BarycentricPolynomialInterp barycentric(x, y);
	BarycentricRationalInterp floaterHormann(x, y, 1);
	CHECK_FALSE(HasEstimate(polynomial.getLastErrorEst()));
	CHECK_FALSE(HasEstimate(rational.getLastErrorEst()));
	auto check = [&](const auto& function, bool hasEstimate) {
		InterpolationConfig config;
		const auto result = function.EvaluateDetailed(Real(0.5), config);
		CHECK(result.status == AlgorithmStatus::Completed);
		CHECK(result.IsSuccess());
		CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(4), Precision::DefaultToleranceStrict));
		CHECK(HasEstimate(result.error_estimate) == hasEstimate);
		CHECK(result.function_evaluations == 1);
		config.estimate_error = false;
		const auto disabled = function.EvaluateDetailed(Real(0.5), config);
		CHECK(disabled.IsSuccess());
		CHECK_FALSE(HasEstimate(disabled.error_estimate));
		CHECK(disabled.value == result.value);
		config.extrapolation_policy = ExtrapolationPolicy::Throw;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto rejected = function.EvaluateDetailed(3, config);
		CHECK(rejected.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(HasEstimate(rejected.error_estimate));
		CHECK_FALSE(rejected.error_message.empty());
	};
	check(linear, false);
	check(polynomial, true);
	check(rational, true);
	check(barycentric, false);
	check(floaterHormann, false);
	const auto polynomialNode = polynomial.EvaluateDetailed(x[0]);
	const auto rationalNode = rational.EvaluateDetailed(x[0]);
	CHECK(polynomialNode.error_estimate == 0);
	CHECK(rationalNode.error_estimate == 0);
	CHECK(HasEstimate(polynomialNode.error_estimate));
	CHECK(HasEstimate(rationalNode.error_estimate));
}

TEST_CASE("Polynomial root residual aggregation starts only after residuals have actually been evaluated",
	"[result-policy][uncomputed-result-contract]") {
	using namespace RootFinding;
	PolynomReal polynomial{-1, 0, 1};
	Vector<Complex> roots{Complex(-1, 0), Complex(1, 0)};
	const auto result = PolishPolynomialRoots(polynomial, roots);
	CHECK(result.IsSuccess());
	CHECK(HasEstimate(result.max_residual));
	CHECK(HasEstimate(result.max_relative_residual));
	Real expectedResidual = 0, expectedRelative = 0;
	for (const auto& root : result.roots) {
		const Real residual = std::abs(root.root * root.root - Real(1));
		CHECK_THAT(root.residual, Catch::Matchers::WithinAbs(residual, Precision::EqualityTolerance));
		expectedResidual = std::max(expectedResidual, root.residual);
		expectedRelative = std::max(expectedRelative, root.relative_residual);
	}
	CHECK(result.max_residual == expectedResidual);
	CHECK(result.max_relative_residual == expectedRelative);
	const auto empty = PolishPolynomialRoots(PolynomReal{Real(1)}, Vector<Complex>{});
	CHECK(empty.status == AlgorithmStatus::Completed);
	CHECK_FALSE(HasEstimate(empty.max_residual));
	CHECK_FALSE(HasEstimate(empty.max_relative_residual));
}

TEST_CASE("Projected gradient distinguishes unused steps from computed norms and iteration exhaustion",
	"[result-policy][uncomputed-result-contract][Optimization]") {
	using namespace Optimization;
	class Quadratic final : public IDifferentiableScalarFunction<1> {
	public:
		Real operator()(const VectorN<Real, 1>& point) const override {
			return point[0] * point[0];
		}
		void Gradient(const VectorN<Real, 1>& point, VectorN<Real, 1>& gradient) const override {
			gradient[0] = Real(2) * point[0];
		}
	} objective;
	BoundConstraints bounds(Vector<Real>{-10}, Vector<Real>{10});
	ProjectedGradientConfig config;
	config.max_iterations = 1;
	config.initial_step = Real(1) / 8;
	const auto initialOptimum = ProjectedGradientMinimize(objective, VectorN<Real, 1>{0}, bounds, config);
	CHECK(initialOptimum.IsSuccess());
	CHECK(initialOptimum.iterations_used == 0);
	CHECK(initialOptimum.projected_gradient_norm == 0);
	CHECK(HasEstimate(initialOptimum.projected_gradient_norm));
	CHECK_FALSE(HasEstimate(initialOptimum.step_norm));
	CHECK(initialOptimum.function_evaluations == 1);
	const auto exhausted = ProjectedGradientMinimize(objective, VectorN<Real, 1>{1}, bounds, config);
	CHECK(exhausted.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK_FALSE(exhausted.IsSuccess());
	CHECK_FALSE(exhausted.error_message.empty());
	CHECK(exhausted.iterations_used == 1);
	CHECK(exhausted.xmin[0] == Real(3) / 4);
	CHECK(exhausted.fmin == Real(9) / 16);
	CHECK(exhausted.step_norm == Real(1) / 4);
	CHECK(exhausted.projected_gradient_norm == Real(3) / 2);
	CHECK(HasEstimate(exhausted.step_norm));
	CHECK(HasEstimate(exhausted.projected_gradient_norm));
	CHECK(exhausted.function_evaluations == 2);
}
