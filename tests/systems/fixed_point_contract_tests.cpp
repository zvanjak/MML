#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/DynamicalSystem/FixedPointAnalysis.h>
#include <mml/systems/DynamicalSystem/DynamicalSystemBase.h>
#include <mml/systems/DynamicalSystemAnalyzer.h>
#endif
#include <catch2/catch_all.hpp>
#include <cmath>
#include <limits>
#include <functional>
#include <stdexcept>

using namespace MML;
using namespace MML::Systems;

namespace {
	class ScalarSystem : public DynamicalSystemBase<1, 0> {
	public:
		enum class Mode { Affine, NoRoot, Nonfinite };
		Mode mode;
		bool autonomous = true;
		int dimension = 1;
		mutable int derivativeCalls = 0;
		mutable int jacobianCalls = 0;
		std::function<void(const Vector<Real>&, Vector<Real>&)> derivativeOverride;
		std::function<void(const Vector<Real>&, Matrix<Real>&)> jacobianOverride;
		explicit ScalarSystem(Mode value = Mode::Affine) : mode(value) {}
		int getDim() const override { return dimension; }
		bool isAutonomous() const override { return autonomous; }
		void derivs(Real, const Vector<Real>& x, Vector<Real>& f) const override {
			++derivativeCalls;
			if (derivativeOverride) { derivativeOverride(x, f); return; }
			f[0] = mode == Mode::Affine ? x[0] - 1 :
				mode == Mode::NoRoot ? x[0] * x[0] + 1 : std::numeric_limits<Real>::quiet_NaN();
		}
		void jacobian(Real, const Vector<Real>& x, Matrix<Real>& matrix) const override {
			++jacobianCalls;
			if (jacobianOverride) { jacobianOverride(x, matrix); return; }
			matrix.resize(1, 1);
			matrix(0, 0) = mode == Mode::NoRoot ? 2 * x[0] : 1;
		}
	};
	template<int N>
	class LinearFlowFixture : public DynamicalSystemBase<N, 0> {
	public:
		Matrix<Real> matrix;
		explicit LinearFlowFixture(const Matrix<Real>& value) : matrix(value) {}
		void derivs(Real, const Vector<Real>& state, Vector<Real>& function) const override { function = matrix * state; }
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& value) const override { value = matrix; }
	};
	void CheckFailure(const FixedPoint<Real>& result, AlgorithmStatus status) {
		CHECK(result.status == status);
		CHECK_FALSE(result.error_message.empty());
		CHECK(result.type == FixedPointType::Unknown);
		CHECK_FALSE(result.isStable);
		CHECK(result.eigenvalues.empty());
	}
}

TEST_CASE("Fixed-point breakdown returns a coherent current location and specific status", "[FixedPoint][FixedPointContracts]") {
	ScalarSystem system(ScalarSystem::Mode::NoRoot);
	const auto result = FindFixedPointDetailed(system, Vector<Real>{1}, Precision::FixedPointTolerance, 5);
	CHECK(result.status == AlgorithmStatus::SingularMatrix);
	CHECK_FALSE(result.error_message.empty());
	CHECK(result.fixed_point.location[0] == 0);
	CHECK(result.fixed_point.convergenceResidual == result.fixed_point.location[0] * result.fixed_point.location[0] + 1);
	CHECK(result.fixed_point.iterations == 1);
	CHECK(result.fixed_point.type == FixedPointType::Unknown);
	CHECK_FALSE(result.fixed_point.isStable);
}

TEST_CASE("Fixed points reached on the final allowed update are fully classified", "[FixedPoint][FixedPointContracts]") {
	ScalarSystem system;
	const auto result = FindFixedPointDetailed(system, Vector<Real>{0}, Precision::FixedPointTolerance, 1);
	REQUIRE(result.IsSuccess());
	CHECK(result.fixed_point.location[0] == 1);
	CHECK(result.fixed_point.convergenceResidual == 0);
	CHECK(result.fixed_point.jacobian(0, 0) == 1);
	REQUIRE(result.fixed_point.eigenvalues.size() == 1);
	CHECK(result.fixed_point.eigenvalues[0] == std::complex<Real>(1, 0));
	CHECK(result.fixed_point.type == FixedPointType::UnstableNode);
	CHECK(result.fixed_point.iterations == 1);
}

TEST_CASE("Nonfinite residuals cannot satisfy zero-budget or multiple-start convergence", "[FixedPoint][FixedPointContracts]") {
	ScalarSystem system(ScalarSystem::Mode::Nonfinite);
	const auto result = FindFixedPointDetailed(system, Vector<Real>{0}, Precision::FixedPointTolerance, 0);
	CHECK_FALSE(result.IsSuccess());
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(result.error_message.empty());
	CHECK(FixedPointFinder::FindMultiple(system, {Vector<Real>{0}, Vector<Real>{1}}).empty());
}

TEST_CASE("Neutral rotation with a decaying axis is not an unstable focus", "[FixedPoint][FixedPointContracts]") {
	FixedPoint<Real> point;
	point.jacobian = Matrix<Real>(3, 3, {0, -1, 0, 1, 0, 0, 0, 0, -1});
	FixedPointFinder::ClassifyFixedPoint(point);
	CHECK(point.type == FixedPointType::Unknown);
	CHECK_FALSE(point.isStable);
	REQUIRE(point.eigenvalues.size() == 3);
	for (const auto& eigenvalue : point.eigenvalues)
		CHECK(eigenvalue.real() <= Precision::DefaultToleranceStrict);
}

TEST_CASE("Unusable finite eigensolutions do not certify fixed-point stability", "[FixedPoint][FixedPointContracts]") {
	FixedPoint<Real> point;
	const Real largest = std::numeric_limits<Real>::max();
	point.jacobian = Matrix<Real>(3, 3, {largest, largest / 2, 0, largest / 2, largest, 0, 0, 0, largest / 4});
	const auto eigen = EigenSolver::Solve(point.jacobian);
	REQUIRE_FALSE(eigen.converged);
	FixedPointFinder::ClassifyFixedPoint(point);
	CHECK(point.type == FixedPointType::Unknown);
	CHECK_FALSE(point.isStable);
	CHECK(point.eigenvalues.empty());
}

TEST_CASE("Large representable eigensolutions classify fixed points in original units", "[FixedPoint][FixedPointContracts]") {
	const Real scale = std::sqrt(std::numeric_limits<Real>::max()) * 1024;
	const Matrix<Real> base(3, 3, {2, 1, Real(0.5), 1, 3, 1, Real(0.5), 1, 4});
	for (Real sign : {Real{-1}, Real{1}}) {
		LinearFlowFixture<3> system(base * (sign * scale));
		const auto result = FindFixedPointDetailed(system, Vector<Real>{0, 0, 0});
		REQUIRE(result.IsSuccess());
		CHECK(result.fixed_point.isStable == (sign < 0));
		CHECK(result.fixed_point.type == (sign < 0 ? FixedPointType::StableNode : FixedPointType::UnstableNode));
		REQUIRE(result.fixed_point.eigenvalues.size() == 3);
		for (const auto& value : result.fixed_point.eigenvalues) {
			CHECK(std::isfinite(value.real()));
			CHECK(value.real() * sign > 0);
			CHECK(value.imag() == Real{0});
		}
		CHECK(FixedPointFinder::FindMultiple(system, {Vector<Real>{0, 0, 0}}).size() == 1);
	}
}

TEST_CASE("Fixed-point requests validate dimensions finite inputs tolerances budgets and autonomy before callbacks", "[FixedPointContracts]") {
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const Real inf = std::numeric_limits<Real>::infinity();
	for (int scenario = 0; scenario < 11; ++scenario) {
		ScalarSystem system;
		Vector<Real> guess{0};
		Real tolerance = Precision::FixedPointTolerance;
		int budget = 1;
		switch (scenario) {
		case 0: system.autonomous = false; break;
		case 1: system.dimension = 0; break;
		case 2: system.dimension = -1; break;
		case 3: guess = Vector<Real>{0, 1}; break;
		case 4: guess[0] = nan; break;
		case 5: guess[0] = inf; break;
		case 6: tolerance = 0; break;
		case 7: tolerance = -1; break;
		case 8: tolerance = nan; break;
		case 9: tolerance = inf; break;
		case 10: budget = -1; break;
		}
		const auto direct = FixedPointFinder::Find(system, guess, tolerance, budget);
		CheckFailure(direct, AlgorithmStatus::InvalidInput);
		CHECK_FALSE(HasEstimate(direct.convergenceResidual));
		CHECK(direct.function_evaluations == 0);
		const auto detailed = FindFixedPointDetailed(system, guess, tolerance, budget);
		CHECK(detailed.status == direct.status);
		CHECK(detailed.error_message == direct.error_message);
		CHECK(detailed.algorithm_name == "FixedPointFinder");
		CHECK(detailed.function_evaluations == 0);
		CHECK(system.derivativeCalls == 0);
		CHECK(system.jacobianCalls == 0);
	}
}

TEST_CASE("Zero-update budgets still evaluate initial convergence and real callback counts", "[FixedPointContracts]") {
	ScalarSystem system;
	auto root = FindFixedPointDetailed(system, Vector<Real>{1}, Precision::FixedPointTolerance, 0);
	REQUIRE(root.IsSuccess());
	CHECK(root.fixed_point.type == FixedPointType::UnstableNode);
	CHECK(root.fixed_point.iterations == 0);
	CHECK(root.function_evaluations == 1);
	CHECK(system.derivativeCalls == 1);
	CHECK(system.jacobianCalls == 1);
	auto stopped = FindFixedPointDetailed(system, Vector<Real>{0}, Precision::FixedPointTolerance, 0);
	CheckFailure(stopped.fixed_point, AlgorithmStatus::MaxIterationsExceeded);
	CHECK(stopped.fixed_point.location[0] == 0);
	CHECK(stopped.fixed_point.convergenceResidual == 1);
	CHECK(stopped.fixed_point.iterations == 0);
	CHECK(stopped.function_evaluations == 1);
	CHECK_FALSE(HasEstimate(stopped.fixed_point.jacobian(0, 0)));
	CHECK(system.jacobianCalls == 1);
	for (int budget : {1, 2, std::numeric_limits<int>::max()}) {
		const auto found = FindFixedPointDetailed(system, Vector<Real>{0}, Precision::FixedPointTolerance, budget);
		REQUIRE(found.IsSuccess());
		CHECK(found.function_evaluations == 2);
		CHECK(found.fixed_point.function_evaluations == 2);
		CHECK(found.fixed_point.iterations == 1);
	}
}

TEST_CASE("Final Newton classification evaluates the Jacobian at the returned root", "[FixedPointContracts]") {
	ScalarSystem system;
	system.derivativeOverride = [](const auto& x, auto& f) { f[0] = (x[0] - 1) * (1 + 6 * x[0] * x[0]); };
	system.jacobianOverride = [](const auto& x, auto& j) { j(0, 0) = 1 + 6 * x[0] * x[0] + (x[0] - 1) * 12 * x[0]; };
	const auto result = FindFixedPointDetailed(system, Vector<Real>{0}, Precision::FixedPointTolerance, 1);
	REQUIRE(result.IsSuccess());
	CHECK(result.fixed_point.location[0] == 1);
	CHECK(result.fixed_point.jacobian(0, 0) == 7);
	CHECK(result.fixed_point.eigenvalues[0] == std::complex<Real>(7, 0));
	CHECK(result.function_evaluations == 2);
}

TEST_CASE("Invalid derivative and Jacobian outputs cannot produce successful fixed points", "[FixedPointContracts]") {
	for (bool convergedInitially : {false, true}) {
		for (int scenario = 0; scenario < 7; ++scenario) {
			ScalarSystem system;
			switch (scenario) {
			case 0: system.derivativeOverride = [](const auto&, auto& f) { f.resize(0); }; break;
			case 1: system.derivativeOverride = [](const auto&, auto&) {}; break;
			case 2: system.derivativeOverride = [](const auto&, auto& f) { f[0] = std::numeric_limits<Real>::infinity(); }; break;
			case 3: system.jacobianOverride = [](const auto&, auto& j) { j.resize(1, 2); }; break;
			case 4: system.jacobianOverride = [](const auto&, auto& j) { j(0, 0) = std::numeric_limits<Real>::quiet_NaN(); }; break;
			case 5: system.jacobianOverride = [](const auto&, auto& j) { j(0, 0) = std::numeric_limits<Real>::infinity(); }; break;
			case 6: system.jacobianOverride = [](const auto&, auto&) {}; break;
			}
			const auto result = FindFixedPointDetailed(system, Vector<Real>{convergedInitially ? Real(1) : Real(0)});
			const auto status = scenario == 0 || scenario == 3 ? AlgorithmStatus::InvalidInput : AlgorithmStatus::NumericalInstability;
			CheckFailure(result.fixed_point, status);
			CHECK(result.status == status);
			CHECK(result.function_evaluations == 1);
			CHECK(result.fixed_point.iterations == 0);
			CHECK(result.fixed_point.location[0] == (convergedInitially ? Real(1) : Real(0)));
		}
	}
}

TEST_CASE("Nonfinite callback output after an update cannot reuse an older residual or classification", "[FixedPointContracts]") {
	ScalarSystem system;
	system.derivativeOverride = [](const auto& x, auto& f) {
		f[0] = x[0] == 0 ? Real(-1) : std::numeric_limits<Real>::quiet_NaN();
	};
	const auto result = FindFixedPointDetailed(system, Vector<Real>{0}, Precision::FixedPointTolerance, 1);
	CheckFailure(result.fixed_point, AlgorithmStatus::NumericalInstability);
	CHECK(result.fixed_point.location[0] == 1);
	CHECK(result.fixed_point.iterations == 1);
	CHECK(result.function_evaluations == 2);
	CHECK_FALSE(HasEstimate(result.fixed_point.convergenceResidual));
	CHECK_FALSE(HasEstimate(result.fixed_point.jacobian(0, 0)));
}

TEST_CASE("Newton step or state overflow preserves the last evaluated finite candidate", "[FixedPointContracts]") {
	for (bool stateOverflow : {false, true}) {
		ScalarSystem system;
		const Real maximum = std::numeric_limits<Real>::max();
		system.derivativeOverride = [=](const auto&, auto& f) { f[0] = stateOverflow ? -maximum / Real(1.5) : Real(1); };
		system.jacobianOverride = [=](const auto&, auto& j) { j(0, 0) = stateOverflow ? Real(1) : std::numeric_limits<Real>::denorm_min(); };
		const Real initial = stateOverflow ? maximum / 2 : Real(0);
		const auto result = FindFixedPointDetailed(system, Vector<Real>{initial});
		CheckFailure(result.fixed_point, AlgorithmStatus::NumericalInstability);
		CHECK(result.fixed_point.location[0] == initial);
		CHECK(std::isfinite(result.fixed_point.convergenceResidual));
		CHECK(result.fixed_point.iterations == 0);
		CHECK(result.function_evaluations == 1);
	}
	ScalarSystem system;
	system.dimension = 2;
	system.derivativeOverride = [](const auto&, auto& f) { f[0] = f[1] = std::numeric_limits<Real>::max(); };
	const auto normOverflow = FixedPointFinder::Find(system, Vector<Real>{0, 0});
	CheckFailure(normOverflow, AlgorithmStatus::NumericalInstability);
	CHECK(std::isinf(normOverflow.convergenceResidual));
}

TEST_CASE("Unexpected user callback exceptions preserve the existing Detailed exception policy", "[FixedPointContracts]") {
	for (bool fromJacobian : {false, true}) {
		ScalarSystem system;
		if (fromJacobian) system.jacobianOverride = [](const auto&, auto&) { throw std::runtime_error("Jacobian callback failure"); };
		else system.derivativeOverride = [](const auto&, auto&) { throw std::runtime_error("derivative callback failure"); };
		CHECK_THROWS_AS(FixedPointFinder::Find(system, Vector<Real>{0}), std::runtime_error);
		CHECK_THROWS_AS(FindFixedPointDetailed(system, Vector<Real>{0}), std::runtime_error);
		DynSysConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto result = FindFixedPointDetailed(system, Vector<Real>{0}, Precision::FixedPointTolerance, 1, config);
		CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
		CHECK(result.error_message.find("callback failure") != std::string::npos);
		CHECK_FALSE(HasEstimate(result.fixed_point.convergenceResidual));
	}
}

TEST_CASE("Repeated classification clears stale stability and rejects malformed matrices", "[FixedPointContracts]") {
	for (int scenario = 0; scenario < 4; ++scenario) {
		FixedPoint<Real> point;
		point.jacobian = Matrix<Real>(1, 1, Real(-1));
		FixedPointFinder::ClassifyFixedPoint(point);
		REQUIRE(point.isStable);
		if (scenario == 0) point.jacobian.resize(0, 0);
		if (scenario == 1) point.jacobian.resize(1, 2);
		if (scenario == 2) point.jacobian(0, 0) = std::numeric_limits<Real>::quiet_NaN();
		if (scenario == 3) point.jacobian(0, 0) = std::numeric_limits<Real>::infinity();
		FixedPointFinder::ClassifyFixedPoint(point);
		CheckFailure(point, scenario < 2 ? AlgorithmStatus::InvalidInput : AlgorithmStatus::NumericalInstability);
	}
}

TEST_CASE("Hyperbolic and neutral fixed-point types follow the complete finite spectrum", "[FixedPointContracts]") {
	struct Case { Matrix<Real> jacobian; FixedPointType type; bool stable; };
	const Real threshold = Precision::DefaultToleranceStrict;
	const std::vector<Case> cases{
		{Matrix<Real>(2, 2, {-1, 0, 0, -2}), FixedPointType::StableNode, true},
		{Matrix<Real>(2, 2, {1, 0, 0, 2}), FixedPointType::UnstableNode, false},
		{Matrix<Real>(2, 2, {-1, 0, 0, 2}), FixedPointType::Saddle, false},
		{Matrix<Real>(2, 2, {-1, -1, 1, -1}), FixedPointType::StableFocus, true},
		{Matrix<Real>(2, 2, {1, -1, 1, 1}), FixedPointType::UnstableFocus, false},
		{Matrix<Real>(3, 3, {-1, -1, 0, 1, -1, 0, 0, 0, 1}), FixedPointType::SaddleFocus, false},
		{Matrix<Real>(2, 2, {0, -1, 1, 0}), FixedPointType::Center, false},
		{Matrix<Real>(3, 3, {0, -1, 0, 1, 0, 0, 0, 0, -1}), FixedPointType::Unknown, false},
		{Matrix<Real>(3, 3, {0, -1, 0, 1, 0, 0, 0, 0, 1}), FixedPointType::Unknown, false},
		{Matrix<Real>(3, 3, {-1, 0, 0, 0, 0, 0, 0, 0, 1}), FixedPointType::Unknown, false},
		{Matrix<Real>(3, 3, {0, -1, 0, 1, 0, 0, 0, 0, 0}), FixedPointType::Unknown, false},
		{Matrix<Real>(1, 1, Real(0)), FixedPointType::Unknown, false},
		{Matrix<Real>(1, 1, threshold), FixedPointType::Unknown, false},
		{Matrix<Real>(1, 1, -threshold), FixedPointType::Unknown, false},
		{Matrix<Real>(1, 1, 2 * threshold), FixedPointType::UnstableNode, false},
		{Matrix<Real>(1, 1, -2 * threshold), FixedPointType::StableNode, true}
	};
	for (const auto& item : cases) {
		FixedPoint<Real> point;
		point.jacobian = item.jacobian;
		FixedPointFinder::ClassifyFixedPoint(point);
		REQUIRE(point.status == AlgorithmStatus::Success);
		CHECK(point.error_message.empty());
		CHECK(point.type == item.type);
		CHECK(point.isStable == item.stable);
		CHECK(point.eigenvalues.size() == static_cast<std::size_t>(item.jacobian.rows()));
	}
}

TEST_CASE("Classification failure retains the converged candidate but cannot enter multi-root results", "[FixedPointContracts]") {
	const Real largest = std::numeric_limits<Real>::max();
	LinearFlowFixture<3> system(Matrix<Real>(3, 3, {largest, largest / 2, 0, largest / 2, largest, 0, 0, 0, largest / 4}));
	const auto result = FindFixedPointDetailed(system, Vector<Real>{0, 0, 0});
	CHECK_FALSE(result.IsSuccess());
	CHECK(result.status == result.fixed_point.status);
	CHECK(result.error_message.find("eigensolve") != std::string::npos);
	CHECK(result.fixed_point.convergenceResidual == 0);
	CHECK(result.fixed_point.location.NormL2() == 0);
	CHECK(result.fixed_point.iterations == 0);
	CHECK(result.function_evaluations == 1);
	CHECK_FALSE(result.fixed_point.isStable);
	CHECK(result.fixed_point.eigenvalues.empty());
	CHECK(FixedPointFinder::FindMultiple(system, {Vector<Real>{0, 0, 0}}).empty());
}

TEST_CASE("Exhausted Newton budgets retain the final evaluated candidate without stale Jacobians", "[FixedPointContracts]") {
	ScalarSystem system;
	system.derivativeOverride = [](const auto& x, auto& f) { f[0] = x[0] * x[0] - 2; };
	system.jacobianOverride = [](const auto& x, auto& j) { j(0, 0) = 2 * x[0]; };
	const auto result = FindFixedPointDetailed(system, Vector<Real>{2}, Precision::FixedPointTolerance, 1);
	CheckFailure(result.fixed_point, AlgorithmStatus::MaxIterationsExceeded);
	CHECK(result.fixed_point.location[0] == Real(1.5));
	CHECK(result.fixed_point.convergenceResidual == Real(0.25));
	CHECK(result.fixed_point.iterations == 1);
	CHECK(result.function_evaluations == 2);
	CHECK_FALSE(HasEstimate(result.fixed_point.jacobian(0, 0)));
}

TEST_CASE("Multiple fixed points validate common parameters and exclude failed or duplicate candidates", "[FixedPointContracts]") {
	ScalarSystem system;
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	for (Real invalid : {Real(0), Real(-1), nan, std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(FixedPointFinder::FindMultiple(system, {}, invalid), ArgumentError);
		CHECK_THROWS_AS(FixedPointFinder::FindMultiple(system, {}, Precision::FixedPointTolerance, invalid), ArgumentError);
	}
	const auto points = FixedPointFinder::FindMultiple(system, {Vector<Real>{0}, Vector<Real>{nan}, Vector<Real>{3}, Vector<Real>{0, 1}});
	REQUIRE(points.size() == 1);
	CHECK(points[0].status == AlgorithmStatus::Success);
	CHECK(points[0].location[0] == 1);
	CHECK(FixedPointFinder::FindMultiple(system, {}).empty());
	system.autonomous = false;
	CHECK_THROWS_AS(FixedPointFinder::FindMultiple(system, {}), ArgumentError);
}

TEST_CASE("Analyzer reports distinguish neutral unknown unstable and saddle-focus equilibria", "[FixedPointContracts]") {
	LinearFlowFixture<2> center(Matrix<Real>(2, 2, {0, -1, 1, 0}));
	DynamicalSystemAnalyzer<> centerAnalyzer(center);
	const auto centerReport = centerAnalyzer.Analyze(Vector<Real>{0, 0}, {}, 1);
	CHECK(centerReport.numNeutralFixedPoints == 1);
	CHECK(centerReport.numUnstableFixedPoints == 0);
	CHECK(centerReport.numUnclassifiedFixedPoints == 0);
	CHECK(centerReport.summary.find("Neutral (linear): 1") != std::string::npos);
	LinearFlowFixture<3> neutral(Matrix<Real>(3, 3, {0, -1, 0, 1, 0, 0, 0, 0, -1}));
	DynamicalSystemAnalyzer<> neutralAnalyzer(neutral);
	const auto neutralReport = neutralAnalyzer.Analyze(Vector<Real>{0, 0, 0}, {}, 1);
	CHECK(neutralReport.numUnclassifiedFixedPoints == 1);
	CHECK(neutralReport.numUnstableFixedPoints == 0);
	CHECK(neutralReport.summary.find("Unclassified: 1") != std::string::npos);
	LinearFlowFixture<3> saddle(Matrix<Real>(3, 3, {-1, -1, 0, 1, -1, 0, 0, 0, 1}));
	DynamicalSystemAnalyzer<> saddleAnalyzer(saddle);
	const auto saddleReport = saddleAnalyzer.Analyze(Vector<Real>{0, 0, 0}, {}, 1);
	CHECK(saddleReport.numSaddlePoints == 1);
	CHECK(saddleReport.numUnstableFixedPoints == 0);
}

TEST_CASE("Adaptive Detailed template defaults preserve all configured and default call forms", "[ODETemplateContracts]") {
	ScalarSystem system;
	system.derivativeOverride = [](const auto&, auto& f) { f[0] = 1; };
	const Vector<Real> initial{0};
	const Real end = Real(0.25);
	const ODEAdaptiveConfig config;
	const auto defaultCall = ODEAdaptiveIntegrateDetailed(system, initial, 0, end);
	const auto explicitStepper = ODEAdaptiveIntegrateDetailed<DormandPrince5_Stepper>(system, initial, 0, end);
	const auto configured = ODEAdaptiveIntegrateDetailed(system, initial, 0, end, config);
	const auto explicitConfigured = ODEAdaptiveIntegrateDetailed<DormandPrince5_Stepper>(system, initial, 0, end, config);
	for (const auto* result : {&defaultCall, &explicitStepper, &configured, &explicitConfigured}) {
		REQUIRE(result->IsSuccess());
		REQUIRE(result->solution.has_value());
		CHECK(result->algorithm_name == "ODEAdaptiveIntegrator");
		CHECK(IsWithinAbs(result->solution->getXValuesAtEnd()[0], end, Precision::EqualityTolerance));
		CHECK(result->function_evaluations > 0);
		CHECK(result->accepted_steps == defaultCall.accepted_steps);
		CHECK(result->function_evaluations == defaultCall.function_evaluations);
	}
}
