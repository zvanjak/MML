#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/DynamicalSystem.h>
#endif
#include <catch2/catch_all.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

using namespace MML;
using namespace MML::Systems;

namespace {
	class EvaluationFlowFixture : public DynamicalSystemBase<2, 1> {
	public:
		bool rotation = false;
		int dimension = 2;
		mutable int calls = 0;
		int parameterWrites = 0;
		std::function<void(Real, const Vector<Real>&, Vector<Real>&)> callback;
		EvaluationFlowFixture() { _params[0] = 1; }
		int getDim() const override { return dimension; }
		void setParam(int index, Real value) override {
			++parameterWrites;
			DynamicalSystemBase<2, 1>::setParam(index, value);
		}
		void derivs(Real time, const Vector<Real>& state, Vector<Real>& derivative) const override {
			++calls;
			if (callback) { callback(time, state, derivative); return; }
			derivative[0] = rotation ? -state[1] : -_params[0] * state[0];
			derivative[1] = rotation ? state[0] : -2 * _params[0] * state[1];
		}
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
			matrix = rotation ? Matrix<Real>(2, 2, {0, -1, 1, 0}) :
				Matrix<Real>(2, 2, {-_params[0], 0, 0, -2 * _params[0]});
		}
	};
	DynSysConfig ConvertFailures() {
		DynSysConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		return config;
	}
}

TEST_CASE("Zero-time trajectories validate initial state before returning success", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto config = ConvertFailures();
	for (const auto& state : {Vector<Real>{}, Vector<Real>{0}, Vector<Real>{std::numeric_limits<Real>::quiet_NaN(), 0}}) {
		const auto result = IntegrateTrajectoryDetailed(flow, state, 0, Real(0.1), Real(0.01), config);
		CHECK(result.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(result.error_message.empty());
		CHECK(result.points.empty());
		CHECK(result.function_evaluations == 0);
	}
	CHECK(flow.calls == 0);
}

TEST_CASE("Valid zero-time trajectories return the initial state with zero RHS work", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto result = IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 2}, 0, Real(0.1));
	REQUIRE(result.IsSuccess());
	REQUIRE(result.points.size() == 1);
	CHECK(result.points[0] == Vector<Real>{1, 2});
	CHECK(result.function_evaluations == 0);
	CHECK(flow.calls == 0);
}

TEST_CASE("Zero-crossing requests validate state count and section without callbacks", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto config = ConvertFailures();
	for (const auto& result : {
		ComputePoincareSectionDetailed(flow, Vector<Real>{0}, PoincareSection<Real>{}, 0, Real(0.01), config),
		ComputePoincareSectionDetailed(flow, Vector<Real>{0, 0}, PoincareSection<Real>{}, -1, Real(0.01), config),
		ComputePoincareSectionDetailed(flow, Vector<Real>{0, 0}, PoincareSection<Real>{0, std::numeric_limits<Real>::quiet_NaN()}, 0, Real(0.01), config)}) {
		CHECK(result.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(result.error_message.empty());
		CHECK(result.function_evaluations == 0);
	}
	CHECK(flow.calls == 0);
}

TEST_CASE("Trajectory Detailed counts direct RHS callbacks rather than output points", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto result = IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 1, Real(0.25));
	REQUIRE(result.IsSuccess());
	CHECK(result.function_evaluations == flow.calls);
	CHECK(flow.calls > 0);
}

TEST_CASE("Poincare Detailed counts transient and event-search RHS callbacks", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	flow.rotation = true;
	const auto result = ComputePoincareSectionDetailed(flow, Vector<Real>{1, 0}, PoincareSection<Real>{0, 0, -1},
		1, Real(0.01), {}, 2, Real(0.1));
	REQUIRE(result.IsSuccess());
	REQUIRE(result.intersections.size() == 1);
	CHECK(result.function_evaluations == flow.calls);
	CHECK(flow.calls > 0);
}

TEST_CASE("Sweep Detailed counts RHS callbacks while preserving legacy nominal metadata", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto result = SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{1, 1}, 0, Real(0.1), Real(0.2), Real(0.01));
	REQUIRE(result.IsSuccess());
	CHECK(result.function_evaluations == flow.calls);
	CHECK(result.diagram.numTransientSteps == 10);
	CHECK(result.diagram.numRecordedPoints == 20);
	CHECK(flow.getParam(0) == 1);
}

TEST_CASE("Systems Detailed classifies invalid controls specifically", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto result = IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, -1, Real(0.1), Real(0.01), ConvertFailures());
	CHECK(result.status == AlgorithmStatus::InvalidInput);
	CHECK_FALSE(result.error_message.empty());
	CHECK(result.function_evaluations == 0);
	CHECK(flow.calls == 0);
}

TEST_CASE("Converted callback failures retain direct RHS work in every affected Detailed API", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	flow.callback = [](Real, const Vector<Real>&, Vector<Real>&) { throw std::runtime_error("RHS fixture failure"); };
	const auto trajectory = IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 1, Real(0.1), Real(0.01), ConvertFailures());
	CHECK(trajectory.status == AlgorithmStatus::AlgorithmSpecificFailure);
	CHECK(trajectory.function_evaluations == 1);
	const auto fixed = FindFixedPointDetailed(flow, Vector<Real>{1, 1}, Precision::FixedPointTolerance,
		Precision::FixedPointMaxIterations, ConvertFailures());
	CHECK(fixed.status == AlgorithmStatus::AlgorithmSpecificFailure);
	CHECK(fixed.function_evaluations == 1);
	CHECK(flow.calls == 2);
}

TEST_CASE("Trajectory validates all finite controls and dimensions before even zero work", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto config = ConvertFailures();
	for (Real invalid : {Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK(IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, invalid, 1, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	}
	for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK(IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 0, invalid, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
		CHECK(IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 0, 1, invalid, config).status == AlgorithmStatus::InvalidInput);
	}
	for (int invalid : {0, -1, std::numeric_limits<int>::max()}) {
		flow.dimension = invalid;
		CHECK(IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 0, 1, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	}
	CHECK(flow.calls == 0);
}

TEST_CASE("Zero Poincare crossings are a validated no-work success", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto result = ComputePoincareSectionDetailed(flow, Vector<Real>{1, 1}, {}, 0);
	CHECK(result.IsSuccess());
	CHECK(result.intersections.empty());
	CHECK(result.searched_time == 0);
	CHECK(result.function_evaluations == 0);
	CHECK(flow.calls == 0);
}

TEST_CASE("Poincare validates sections controls and representable search clocks before work", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto config = ConvertFailures();
	for (const auto& section : {PoincareSection<Real>{-1}, PoincareSection<Real>{2},
		PoincareSection<Real>{0, 0, 2}, PoincareSection<Real>{0, std::numeric_limits<Real>::infinity()}}) {
		CHECK(ComputePoincareSectionDetailed(flow, Vector<Real>{1, 1}, section, 0, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	}
	for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK(ComputePoincareSectionDetailed(flow, Vector<Real>{1, 1}, {}, 0, invalid, config).status == AlgorithmStatus::InvalidInput);
		CHECK(ComputePoincareSectionDetailed(flow, Vector<Real>{1, 1}, {}, 0, Real(0.01), config, invalid, 0).status == AlgorithmStatus::InvalidInput);
	}
	for (Real invalid : {Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK(ComputePoincareSectionDetailed(flow, Vector<Real>{1, 1}, {}, 0, Real(0.01), config, 1, invalid).status == AlgorithmStatus::InvalidInput);
	}
	const Real huge = std::numeric_limits<Real>::max();
	CHECK(ComputePoincareSectionDetailed(flow, Vector<Real>{1, 1}, {}, 0, Real(0.01), config, huge, huge).status == AlgorithmStatus::InvalidInput);
	CHECK(ComputePoincareSectionDetailed(flow, Vector<Real>{1, 1}, {}, 0, Real(0.01), config, 1, huge).status == AlgorithmStatus::InvalidInput);
	CHECK(flow.calls == 0);
}

TEST_CASE("Sweep rejects malformed requests and unsafe nominal estimates before mutation", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto config = ConvertFailures();
	for (Real invalid : {Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK(SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{1, 1}, 0, invalid, 0, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
		CHECK(SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{1, 1}, 0, 0, invalid, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	}
	for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK(SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{1, 1}, 0, 0, 0, invalid, config).status == AlgorithmStatus::InvalidInput);
	}
	for (const auto& state : {Vector<Real>{}, Vector<Real>{std::numeric_limits<Real>::quiet_NaN(), 0}}) {
		CHECK(SweepBifurcationDetailed(flow, 0, 1, 2, 2, state, 0, 0, 0, 1, config).status == AlgorithmStatus::InvalidInput);
	}
	const auto excessive = SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{1, 1}, 0,
		0, Real(std::numeric_limits<int>::max()) + 1, 1, config);
	CHECK(excessive.status == AlgorithmStatus::ResourceLimitExceeded);
	CHECK(excessive.function_evaluations == 0);
	const Real huge = std::numeric_limits<Real>::max();
	const auto extreme = SweepBifurcationDetailed(flow, 0, -huge, huge, 2, Vector<Real>{1, 1}, 0, 0, 0, 1, config);
	REQUIRE(extreme.IsSuccess());
	CHECK(extreme.diagram.parameterValues == std::vector<Real>{-huge, huge});
	CHECK(extreme.function_evaluations == 0);
	CHECK(flow.getParam(0) == 1);
	CHECK(flow.parameterWrites == 3);
	CHECK(flow.calls == 0);
}

TEST_CASE("Legacy sweep estimates round and range-check without narrowing overflow", "[SystemsEvaluationContracts]") {
	CHECK(DynSysDetail::NominalIntervalCount(Real(0.5), 1) == 1);
	CHECK(DynSysDetail::NominalIntervalCount(0, 1) == 0);
	const Real roundedLimit = Real(std::numeric_limits<int>::max());
	if (static_cast<long double>(roundedLimit) > std::numeric_limits<int>::max()) {
		const Real belowLimit = std::nextafter(roundedLimit, Real(0));
		REQUIRE(static_cast<long double>(belowLimit) <= std::numeric_limits<int>::max());
		CHECK(DynSysDetail::NominalIntervalCount(belowLimit, 1) == static_cast<int>(belowLimit));
		CHECK_THROWS_AS(DynSysDetail::NominalIntervalCount(roundedLimit, 1), ResourceLimitError);
	} else {
		CHECK(DynSysDetail::NominalIntervalCount(roundedLimit, 1) == std::numeric_limits<int>::max());
	}
	CHECK_THROWS_AS(DynSysDetail::NominalIntervalCount(Real(std::numeric_limits<int>::max()) + 1, 1), ResourceLimitError);
	CHECK_THROWS_AS(DynSysDetail::NominalIntervalCount(std::numeric_limits<Real>::max(), std::numeric_limits<Real>::denorm_min()), ResourceLimitError);
}

TEST_CASE("Zero-duration sweeps retain legacy metadata and report zero RHS work", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto result = SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{1, 1}, 0, 0, 0);
	REQUIRE(result.IsSuccess());
	CHECK(result.diagram.parameterValues.size() == 2);
	CHECK(result.diagram.attractorValues.size() == 2);
	CHECK(result.diagram.numTransientSteps == 0);
	CHECK(result.diagram.numRecordedPoints == 0);
	CHECK(result.function_evaluations == 0);
	CHECK(flow.calls == 0);
	CHECK(flow.getParam(0) == 1);
}

TEST_CASE("Poincare budget exhaustion retains partial intersections and actual work", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	flow.rotation = true;
	const auto result = ComputePoincareSectionDetailed(flow, Vector<Real>{1, 0}, {0, 0, -1}, 2,
		Real(0.01), {}, 2, 0);
	CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK_FALSE(result.error_message.empty());
	REQUIRE(result.intersections.size() == 1);
	CHECK(IsNearlyZero(result.intersections[0][0], Precision::CoarseTolerance));
	CHECK(result.searched_time == 2);
	CHECK(result.function_evaluations == flow.calls);
	CHECK(flow.calls > 0);
	CHECK_THROWS_AS(PhaseSpaceAnalyzer::ComputePoincareSection(flow, Vector<Real>{1, 0}, {0, 0, -1},
		2, Real(0.01), 2, 0), ConvergenceError);
}

TEST_CASE("Every integration API rejects malformed and unwritten RHS outputs safely", "[SystemsEvaluationContracts]") {
	for (bool finiteChecks : {true, false}) {
		for (bool wrongShape : {true, false}) {
			EvaluationFlowFixture flow;
			auto config = ConvertFailures();
			config.check_finite = finiteChecks;
			flow.callback = [=](Real, const Vector<Real>&, Vector<Real>& derivative) {
				if (wrongShape) derivative.resize(1);
			};
			const AlgorithmStatus expected = wrongShape ? AlgorithmStatus::InvalidInput : AlgorithmStatus::NumericalInstability;
			const auto trajectory = IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 1, Real(0.1), Real(0.01), config);
			const auto section = ComputePoincareSectionDetailed(flow, Vector<Real>{1, 1}, {}, 1, Real(0.01), config, 1, 0);
			const auto sweep = SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{1, 1}, 0, 0, 1, Real(0.01), config);
			for (const EvaluationResultBase* result : {static_cast<const EvaluationResultBase*>(&trajectory),
				static_cast<const EvaluationResultBase*>(&section), static_cast<const EvaluationResultBase*>(&sweep)}) {
				CHECK(result->status == expected);
				CHECK_FALSE(result->error_message.empty());
				CHECK(result->function_evaluations == 1);
			}
			CHECK(flow.calls == 3);
			CHECK(trajectory.points.empty());
			CHECK(section.intersections.empty());
			CHECK(sweep.diagram.parameterValues.empty());
			CHECK(flow.getParam(0) == 1);
		}
	}
}

TEST_CASE("Optional finite-stage checks cannot relax mandatory input and callback safety", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	flow.callback = [](Real, const Vector<Real>&, Vector<Real>& derivative) { derivative = Vector<Real>{1, 2}; };
	const Vector<Real> nonfinite{std::numeric_limits<Real>::quiet_NaN(), 0};
	int evaluations = 0;
	Vector<Real> derivative(2);
	DynamicalAnalysisDetail::CheckedODESystem strict(flow, evaluations, true);
	CHECK_THROWS_AS(strict.derivs(0, nonfinite, derivative), NumericalMethodError);
	CHECK(evaluations == 0);
	DynamicalAnalysisDetail::CheckedODESystem optional(flow, evaluations, false);
	CHECK_NOTHROW(optional.derivs(0, nonfinite, derivative));
	CHECK(evaluations == 1);
	CHECK(flow.calls == 1);
	CHECK(derivative == Vector<Real>{1, 2});
	auto config = ConvertFailures();
	config.check_finite = false;
	CHECK(IntegrateTrajectoryDetailed(flow, nonfinite, 0, 1, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	CHECK(ComputePoincareSectionDetailed(flow, nonfinite, {}, 0, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	CHECK(flow.calls == 1);
}

TEST_CASE("Unrepresentable trajectories fail with actual work when optional checks are disabled", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	auto config = ConvertFailures();
	config.check_finite = false;
	const Real huge = std::numeric_limits<Real>::max();
	flow.callback = [=](Real, const Vector<Real>&, Vector<Real>& derivative) { derivative = Vector<Real>{huge, huge}; };
	const auto result = IntegrateTrajectoryDetailed(flow, Vector<Real>{huge, huge}, 1, 1, Real(0.1), config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(result.error_message.empty());
	CHECK(result.points.empty());
	CHECK(result.function_evaluations == flow.calls);
	CHECK(flow.calls > 0);
}

TEST_CASE("Detailed failure categories preserve callback work and propagation policy", "[SystemsEvaluationContracts]") {
	const auto config = ConvertFailures();
	for (AlgorithmStatus status : {AlgorithmStatus::InvalidInput, AlgorithmStatus::SingularMatrix,
		AlgorithmStatus::NumericalInstability, AlgorithmStatus::ResourceLimitExceeded, AlgorithmStatus::MaxIterationsExceeded}) {
		EvaluationFlowFixture flow;
		flow.callback = [=](Real, const Vector<Real>&, Vector<Real>&) {
			switch (status) {
				case AlgorithmStatus::InvalidInput: throw ArgumentError("fixture category");
				case AlgorithmStatus::SingularMatrix: throw SingularMatrixError("fixture category");
				case AlgorithmStatus::NumericalInstability: throw MatrixNumericalError("fixture category");
				case AlgorithmStatus::ResourceLimitExceeded: throw ResourceLimitError("fixture category");
				case AlgorithmStatus::MaxIterationsExceeded: throw ConvergenceError("fixture category");
				default: throw std::logic_error("invalid test category");
			}
		};
		const auto result = IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 1, 1, Real(0.01), config);
		CHECK(result.status == status);
		CHECK(result.error_message == "fixture category");
		CHECK(result.function_evaluations == 1);
		CHECK(flow.calls == 1);
		CHECK_THROWS_WITH(IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 1, 1), "fixture category");
		CHECK(flow.calls == 2);
	}
}

TEST_CASE("Trajectory output retains the initial state final endpoint and decay accuracy", "[SystemsEvaluationContracts]") {
	EvaluationFlowFixture flow;
	const auto result = IntegrateTrajectoryDetailed(flow, Vector<Real>{1, 1}, 1, Real(0.3));
	REQUIRE(result.IsSuccess());
	REQUIRE(result.points.size() == 5);
	CHECK(result.points.front() == Vector<Real>{1, 1});
	CHECK(IsWithinAbs(result.points.back()[0], std::exp(Real(-1)), Precision::CoarseTolerance));
	CHECK(IsWithinAbs(result.points.back()[1], std::exp(Real(-2)), Precision::CoarseTolerance));
	CHECK(result.function_evaluations == flow.calls);
}
