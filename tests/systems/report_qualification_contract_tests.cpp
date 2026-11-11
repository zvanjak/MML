#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/DynamicalSystem.h>
#endif
#include <catch2/catch_all.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

using namespace MML;
using namespace MML::Systems;

namespace {
	template<int N>
	class ReportFlow : public DynamicalSystemBase<N, 0> {
		Matrix<Real> _matrix;
		Vector<Real> _offset;
	public:
		bool failOnDerivative = false;
		explicit ReportFlow(Matrix<Real> matrix, Vector<Real> offset = Vector<Real>(N, Real(0)))
			: _matrix(std::move(matrix)), _offset(std::move(offset)) {}
		void derivs(Real, const Vector<Real>& state, Vector<Real>& derivative) const override {
			if (failOnDerivative) throw std::runtime_error("report callback sentinel");
			derivative = _matrix * state + _offset;
		}
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override { matrix = _matrix; }
	};

	class TwoEquilibria : public DynamicalSystemBase<1, 0> {
	public:
		void derivs(Real, const Vector<Real>& state, Vector<Real>& derivative) const override {
			derivative[0] = state[0] * state[0] - 1;
		}
		void jacobian(Real, const Vector<Real>& state, Matrix<Real>& matrix) const override {
			matrix = Matrix<Real>(1, 1, {2 * state[0]});
		}
	};
}

TEST_CASE("Zero-exponent drift reports do not claim bounded or convergent trajectories", "[ReportQualificationContracts][ReportBaseline]") {
	ReportFlow<1> drift(Matrix<Real>(1, 1, {0}), Vector<Real>{1});
	const auto report = DynamicalSystemAnalyzer<>(drift).Analyze(Vector<Real>{0}, {}, 1);
	CHECK(report.lyapunov.maxExponent == 0);
	CHECK_FALSE(report.isChaotic);
	CHECK(report.summary.find("Trajectories converge or remain bounded") == std::string::npos);
	CHECK(report.summary.find("No boundedness or convergence conclusion") != std::string::npos);
}

TEST_CASE("Positive-exponent reports present a growth heuristic rather than certified chaos", "[ReportQualificationContracts][ReportBaseline]") {
	ReportFlow<1> expanding(Matrix<Real>(1, 1, {1}));
	const auto report = DynamicalSystemAnalyzer<>(expanding).Analyze(Vector<Real>{1}, {}, 1);
	CHECK(report.isChaotic);
	CHECK(report.summary.find("Finite-time positive growth: detected") != std::string::npos);
	CHECK(report.summary.find("CHAOTIC DYNAMICS") == std::string::npos);
	CHECK(report.summary.find("not a proof of chaos") != std::string::npos);
}

TEST_CASE("Negative-sum reports do not infer a bounded attractor or measured fractal dimension", "[ReportQualificationContracts][ReportBaseline]") {
	ReportFlow<2> expanding(Matrix<Real>(2, 2, {1, 0, 0, -2}));
	const auto report = DynamicalSystemAnalyzer<>(expanding).Analyze(Vector<Real>{1, 0}, {}, 1);
	CHECK(report.isChaotic);
	CHECK(report.isDissipative);
	CHECK(IsWithinAbs(report.kaplanYorkeDimension, Real(1.5), Precision::CoarseTolerance));
	CHECK(report.summary.find("Volume contraction indicator: detected") != std::string::npos);
	CHECK(report.summary.find("Attractor has dimension") == std::string::npos);
	CHECK(report.summary.find("Kaplan-Yorke dimension estimate") != std::string::npos);
}

TEST_CASE("Near-zero exponent sums do not establish Hamiltonian structure or exact volume preservation", "[ReportQualificationContracts][ReportBaseline]") {
	ReportFlow<1> drift(Matrix<Real>(1, 1, {0}), Vector<Real>{1});
	const auto report = DynamicalSystemAnalyzer<>(drift).Analyze(Vector<Real>{0}, {}, 1);
	CHECK(report.isConservative);
	CHECK(report.summary.find("Near-zero-sum indicator: detected") != std::string::npos);
	CHECK(report.summary.find("not proof of Hamiltonian structure") != std::string::npos);
	CHECK(report.summary.find("CONSERVATIVE (") == std::string::npos);
}

TEST_CASE("Default reports deterministically represent uncomputed analysis", "[ReportQualificationContracts]") {
	DynamicalSystemReport<> defaultInitialized;
	DynamicalSystemReport<> valueInitialized{};
	for (const auto* report : {&defaultInitialized, &valueInitialized}) {
		CHECK(report->dimension == 0);
		CHECK(report->numParameters == 0);
		CHECK(report->numStableFixedPoints == 0);
		CHECK(report->numUnstableFixedPoints == 0);
		CHECK(report->numSaddlePoints == 0);
		CHECK(report->numNeutralFixedPoints == 0);
		CHECK(report->numUnclassifiedFixedPoints == 0);
		CHECK_FALSE(report->isChaotic);
		CHECK_FALSE(report->isDissipative);
		CHECK_FALSE(report->isConservative);
		CHECK(std::isnan(report->kaplanYorkeDimension));
		CHECK(std::isnan(report->lyapunovTime));
		CHECK(report->stateNames.empty());
		CHECK(report->paramNames.empty());
		CHECK(report->currentParams.size() == 0);
		CHECK(report->fixedPoints.empty());
		CHECK(report->lyapunov.exponents.size() == 0);
		CHECK(report->trajectory.empty());
		CHECK(report->poincareSection.empty());
		CHECK(report->summary.empty());
	}
}

TEST_CASE("Report metadata distinguishes computed estimates from uncomputed defaults", "[ReportQualificationContracts]") {
	ReportFlow<1> decay(Matrix<Real>(1, 1, {-1}));
	const auto report = DynamicalSystemAnalyzer<>(decay).Analyze(Vector<Real>{1}, {}, Real(0.25));
	CHECK(report.dimension == 1);
	CHECK(report.numParameters == 0);
	CHECK(report.stateNames.size() == 1);
	CHECK(report.currentParams.size() == 0);
	CHECK(report.lyapunov.totalTime == Real(0.25));
	CHECK(report.summary.find("Finite-time duration: 0.25") != std::string::npos);
	CHECK(report.kaplanYorkeDimension == 0);
	CHECK(std::isinf(report.lyapunovTime));
	CHECK(report.summary.find("No positive-growth time scale reported") != std::string::npos);
	CHECK(report.summary.find("not a certified predictability horizon") != std::string::npos);
	CHECK(report.trajectory.empty());
	CHECK(report.poincareSection.empty());
}

TEST_CASE("Report flags and quick-query numeric behavior retain their existing thresholds", "[ReportQualificationContracts]") {
	for (Real rate : {Real(-1), Real(0), Real(1), Precision::LyapunovChaosThreshold / 2}) {
		ReportFlow<1> flow(Matrix<Real>(1, 1, {rate}));
		DynamicalSystemAnalyzer<> analyzer(flow);
		const auto report = analyzer.Analyze(Vector<Real>{1}, {}, 1);
		CHECK(report.isChaotic == (report.lyapunov.maxExponent > Precision::LyapunovChaosThreshold));
		CHECK(report.isDissipative == (report.lyapunov.sum < -Precision::LooseTolerance));
		CHECK(report.isConservative == IsNearlyZero(report.lyapunov.sum, Precision::CoarseTolerance));
		CHECK(report.kaplanYorkeDimension == report.lyapunov.kaplanYorkeDimension);
		CHECK(analyzer.IsChaotic(Vector<Real>{1}, 1) == report.isChaotic);
		CHECK(analyzer.IsDissipative(Vector<Real>{1}, 1) == report.isDissipative);
		CHECK(analyzer.GetFractalDimension(Vector<Real>{1}, 1) == report.kaplanYorkeDimension);
		if (report.isChaotic) {
			CHECK(report.lyapunovTime == Real(1) / report.lyapunov.maxExponent);
			CHECK(analyzer.GetLyapunovTime(Vector<Real>{1}, 1) == report.lyapunovTime);
		} else {
			CHECK(std::isinf(report.lyapunovTime));
			CHECK(std::isinf(analyzer.GetLyapunovTime(Vector<Real>{1}, 1)));
		}
	}
}

TEST_CASE("Independent volume tolerances can overlap without hiding either report indicator", "[ReportQualificationContracts]") {
	const Real rate = -std::sqrt(Precision::LooseTolerance * Precision::CoarseTolerance);
	ReportFlow<1> flow(Matrix<Real>(1, 1, {rate}));
	const auto report = DynamicalSystemAnalyzer<>(flow).Analyze(Vector<Real>{1}, {}, 1);
	REQUIRE(report.isDissipative);
	REQUIRE(report.isConservative);
	CHECK(report.summary.find("Volume contraction indicator: detected") != std::string::npos);
	CHECK(report.summary.find("Near-zero-sum indicator: detected") != std::string::npos);
	CHECK(report.summary.find("independent tolerances and can overlap") != std::string::npos);
}

TEST_CASE("Subthreshold positive growth is not relabeled as a nonpositive exponent or regular dynamics", "[ReportQualificationContracts]") {
	ReportFlow<1> flow(Matrix<Real>(1, 1, {Precision::LyapunovChaosThreshold / 2}));
	const auto report = DynamicalSystemAnalyzer<>(flow).Analyze(Vector<Real>{1}, {}, 1);
	REQUIRE(report.lyapunov.maxExponent > 0);
	CHECK_FALSE(report.isChaotic);
	CHECK(report.summary.find("Finite-time positive growth: not detected") != std::string::npos);
	CHECK(report.summary.find("REGULAR DYNAMICS") == std::string::npos);
	CHECK(report.summary.find("Trajectories converge or remain bounded") == std::string::npos);
}

TEST_CASE("Report fixed-point searches are explicitly nonexhaustive and local", "[ReportQualificationContracts]") {
	TwoEquilibria flow;
	const auto report = DynamicalSystemAnalyzer<>(flow).Analyze(Vector<Real>{-1}, {Vector<Real>{1}}, 1);
	REQUIRE(report.fixedPoints.size() == 1);
	CHECK(report.fixedPoints.front().location[0] == 1);
	CHECK(report.numUnstableFixedPoints == 1);
	CHECK_FALSE(report.isChaotic);
	CHECK(report.summary.find("Search is non-exhaustive") != std::string::npos);
	CHECK(report.summary.find("not global stability") != std::string::npos);
}

TEST_CASE("Callback errors cannot become successful report classifications", "[ReportQualificationContracts]") {
	ReportFlow<1> flow(Matrix<Real>(1, 1, {-1}));
	flow.failOnDerivative = true;
	CHECK_THROWS_WITH(DynamicalSystemAnalyzer<>(flow).Analyze(Vector<Real>{1}, {}, 1), "report callback sentinel");
}

TEST_CASE("Report copies preserve unavailable scalar estimates without shared aliases", "[ReportQualificationContracts]") {
	DynamicalSystemReport<> original;
	auto copy = original;
	copy.dimension = 1;
	copy.isChaotic = true;
	copy.kaplanYorkeDimension = 1;
	CHECK(original.dimension == 0);
	CHECK_FALSE(original.isChaotic);
	CHECK(std::isnan(original.kaplanYorkeDimension));
	CHECK(std::isnan(original.lyapunovTime));
	CHECK(copy.kaplanYorkeDimension == 1);
}
