#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/RootFinding.h>
#include <mml/algorithms/Optimization/Constraints/AugmentedLagrangian.h>
#include <mml/algorithms/Optimization/Constraints/BoxConstrainedPowell.h>
#include <mml/algorithms/Optimization/Multidim/MultidimSolvers.h>
#endif

using namespace MML;

namespace {
	class ScriptedObjective final : public IScalarFunction<2> {
	public:
		explicit ScriptedObjective(std::function<Real()> callback) : _callback(std::move(callback)) {}
		Real operator()(const VectorN<Real, 2>&) const override { return _callback(); }
	private:
		std::function<Real()> _callback;
	};

	MinimumBracket AccountingMinimumBracket()
	{
		RealFunction function([](Real x) {
			const Real offset = x - Real(1) / 3;
			return 1 + offset * offset;
		});
		return Minimization::BracketMinimum(function, -1, 0);
	}

	MinimumBracket UnitCenteredMinimumBracket(bool narrow = false)
	{
		MinimumBracket bracket;
		bracket.ax = narrow ? 1 - Precision::Eps : 0;
		bracket.bx = 1;
		bracket.cx = narrow ? 1 + Precision::Eps : 2;
		bracket.fa = (bracket.ax - 1) * (bracket.ax - 1);
		bracket.fb = 0;
		bracket.fc = (bracket.cx - 1) * (bracket.cx - 1);
		bracket.valid = true;
		return bracket;
	}

	void CheckBrentOutcome(const Minimization1DResult& result, int calls, bool converged, Real value = 0)
	{
		CHECK(result.converged == converged);
		CHECK(result.status == (converged ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded));
		CHECK(result.IsSuccess() == converged);
		CHECK(static_cast<bool>(result) == converged);
		CHECK(result.error_message.empty() == converged);
		CHECK(result.function_evaluations == calls);
		CHECK(result.iterations_used == result.iterations);
		CHECK(result.iterations > 0);
		CHECK_FALSE(result.algorithm_name.empty());
		CHECK(result.elapsed_time_ms >= 0);
		CHECK(result.xmin == 1);
		CHECK(result.fmin == value);
	}

	template<int N>
	void CheckPowellShapeRejection(Optimization::Powell& optimizer, const IScalarFunction<N>& function,
		const VectorN<Real, N>& start, Matrix<Real>& directions)
	{
		const auto before = directions;
		bool rejected = false;
		try {
			const auto unexpected = optimizer.Minimize(function, start, directions);
			(void)unexpected;
		} catch (const Optimization::MultidimOptimizationInputError& error) {
			rejected = true;
			const std::string message(error.what());
			CHECK(message.find("Powell::Minimize") != std::string::npos);
			CHECK(message.find("direction matrix") != std::string::npos);
			CHECK(message.find("expected " + std::to_string(N) + " x " + std::to_string(N)) != std::string::npos);
			CHECK(message.find("got " + std::to_string(before.rows()) + " x " + std::to_string(before.cols())) != std::string::npos);
		}
		CHECK(rejected);
		CHECK(optimizer.getNumFuncEvals() == 0);
		REQUIRE(directions.rows() == before.rows());
		REQUIRE(directions.cols() == before.cols());
		for (int i = 0; i < before.rows(); ++i) {
			for (int j = 0; j < before.cols(); ++j) {
				if (std::isnan(before(i, j))) CHECK(std::isnan(directions(i, j)));
				else CHECK(directions(i, j) == before(i, j));
			}
		}
	}
}

TEST_CASE("Callback accounting: scalar roots count endpoints samples and final checks",
	"[callback-accounting][RootFinding]")
{
	using namespace RootFinding;
	using Solver = RootFindingResult (*)(const IRealFunction&, Real, Real, const RootFindingConfig&);
	const Solver solvers[] = {&FindRootBisectionDetailed, &FindRootFalsePositionDetailed,
		&FindRootSecantDetailed, &FindRootRiddersDetailed, &FindRootBrentDetailed,
		&FindRootNewtonDetailed};
	for (const auto solve : solvers) {
		for (const int limit : {1, 100}) {
			RootFindingConfig config;
			config.max_iterations = limit;
			int calls = 0;
			RealFunctionFromStdFunc function([&](Real x) { ++calls; return x * x - 2; });
			const auto result = solve(function, 1, 2, config);
			CAPTURE(result.algorithm_name, limit);
			CHECK(result.function_evaluations == calls);
			CHECK(result.iterations_used <= limit);

			calls = 0;
			const auto invalid = solve(function, 1, 1, config);
			CHECK(invalid.status == AlgorithmStatus::InvalidInput);
			CHECK(invalid.function_evaluations == 0);
			CHECK(calls == 0);

			calls = 0;
			RealFunctionFromStdFunc endpoint([&](Real x) { ++calls; return x - 1; });
			const auto completed = solve(endpoint, 1, 2, config);
			CHECK(completed.IsSuccess());
			CHECK(completed.iterations_used == 0);
			CHECK(completed.function_evaluations == calls);
			CHECK(calls == 2);

			calls = 0;
			RealFunctionFromStdFunc nonfinite([&](Real) {
				++calls;
				return std::numeric_limits<Real>::quiet_NaN();
			});
			const auto failed = solve(nonfinite, 1, 2, config);
			CHECK(failed.status == AlgorithmStatus::NumericalInstability);
			CHECK(failed.function_evaluations == calls);
			CHECK(calls == 2);
		}
	}
}

TEST_CASE("Callback accounting: Newton separates analytic derivatives from objective samples",
	"[callback-accounting][RootFinding]")
{
	RootFinding::RootFindingConfig config;
	config.max_iterations = 1;
	int calls = 0, derivativeCalls = 0;
	RealFunctionFromStdFunc function([&](Real x) { ++calls; return x * x - 2; });
	RealFunctionFromStdFunc derivative([&](Real x) { ++derivativeCalls; return 2 * x; });
	const auto analytic = RootFinding::FindRootNewtonDetailed(function, derivative, 1, 2, config);
	CHECK(analytic.function_evaluations == calls);
	CHECK(calls == 4);
	CHECK(derivativeCalls == 1);
	calls = 0;
	const auto numerical = RootFinding::FindRootNewtonDetailed(function, 1, 2, config);
	CHECK(numerical.function_evaluations == calls);
	CHECK(calls == 8);
	CHECK(numerical.iterations_used == analytic.iterations_used);
}

TEST_CASE("Callback accounting: rejected root refinements include throwing callbacks",
	"[callback-accounting][RootFinding]")
{
	using namespace RootFinding;
	for (const bool exact : {false, true}) {
		for (const int failAt : {1, 2, 3}) {
			if (exact && failAt != 1) continue;
			FindAllRealRootsConfig config;
			config.isolation.num_intervals = exact ? 2 : 1;
			const int isolationCalls = config.isolation.num_intervals + 1;
			int calls = 0;
			auto function = [&](Real x) {
				if (++calls == isolationCalls + failAt)
					throw std::runtime_error("controlled root callback");
				return exact ? x : x * x - 2;
			};
			const auto result = FindAllRealRootsInInterval(function, exact ? Real(-1) : Real(1), exact ? Real(1) : Real(2), config);
			REQUIRE(result.rejected.size() == 1);
			CHECK(result.rejected.front().refinement.function_evaluations == failAt);
			CHECK(result.rejected.front().refinement.error_message == "controlled root callback");
			CHECK(result.function_evaluations == calls);
			CHECK(calls == isolationCalls + failAt);
		}
	}
}

TEST_CASE("Callback accounting: tangent refinements count actual failed and final calls",
	"[callback-accounting][RootFinding]")
{
	using namespace RootFinding;
	FindAllRealRootsConfig config;
	config.isolation.num_intervals = 2;
	config.isolation.tangent_refinement_iterations = 8;
	config.isolation.tangent_tolerance = 1;
	config.refinement.f_tolerance = 1;
	int isolationCalls = 0;
	RealFunctionFromStdFunc oracle([&](Real x) {
		++isolationCalls;
		const Real offset = x - Real(1) / 3;
		return offset * offset;
	});
	const auto candidates = FindRootIntervals(oracle, 0, 1, config.isolation);
	REQUIRE(candidates.size() == 1);
	REQUIRE(candidates.front().type == RootCandidateType::TangentCandidate);
	for (const bool throwing : {false, true}) {
		for (int failAt = 0; failAt <= 13; ++failAt) {
			CAPTURE(throwing, failAt);
			int calls = 0;
			auto function = [&](Real x) {
				if (++calls == isolationCalls + failAt && failAt > 0) {
					if (throwing) throw std::runtime_error("controlled tangent callback");
					return std::numeric_limits<Real>::quiet_NaN();
				}
				const Real offset = x - Real(1) / 3;
				return offset * offset;
			};
			const auto result = FindAllRealRootsInInterval(function, 0, 1, config);
			CHECK(result.function_evaluations == calls);
			if (failAt == 0) {
				REQUIRE(result.accepted.size() == 1);
				CHECK(result.accepted.front().refinement.function_evaluations == calls - isolationCalls);
				CHECK(result.accepted.front().refinement.iterations_used == 8);
			} else {
				REQUIRE(result.rejected.size() == 1);
				const auto& refinement = result.rejected.front().refinement;
				CHECK(refinement.status == AlgorithmStatus::NumericalInstability);
				CHECK_FALSE(refinement.error_message.empty());
				if (throwing) CHECK(refinement.error_message == "controlled tangent callback");
				CHECK(refinement.function_evaluations == calls - isolationCalls);
			}
		}
	}
}

TEST_CASE("Callback accounting: Nelder Mead counts attempts at every callback boundary",
	"[callback-accounting][NelderMead]")
{
	using namespace Optimization;
	for (const bool expansion : {false, true}) {
		for (int failAt = 1; failAt <= (expansion ? 5 : 7); ++failAt) {
			CAPTURE(expansion, failAt);
			int calls = 0;
			NelderMead optimizer;
			ScriptedObjective function([&]() {
				++calls;
				CHECK(optimizer.getNumFuncEvals() == calls);
				if (calls == failAt) throw std::runtime_error("controlled simplex callback");
				if (calls <= 3) return Real(calls - 1);
				return expansion ? Real(-1) : Real(3);
			});
			CHECK_THROWS_AS(optimizer.Minimize(function, VectorN<Real, 2>{0, 0}), std::runtime_error);
			CHECK(calls == failAt);
			CHECK(optimizer.getNumFuncEvals() == calls);
		}
	}
}

TEST_CASE("Callback accounting: Nelder Mead resets counts on reuse and validation failures",
	"[callback-accounting][NelderMead]")
{
	using namespace Optimization;
	int calls = 0;
	NelderMead optimizer;
	ScriptedObjective constant([&]() { ++calls; return Real(1); });
	const auto first = optimizer.Minimize(constant, VectorN<Real, 2>{0, 0});
	CHECK(first.converged);
	CHECK(first.function_evaluations == calls);
	CHECK(first.iterations == 0);
	CHECK(first.iterations_used == 0);
	CHECK(calls == 3);

	calls = 0;
	ScriptedObjective failed([&]() {
		if (++calls == 2) return std::numeric_limits<Real>::quiet_NaN();
		return Real(1);
	});
	CHECK_THROWS_AS(optimizer.Minimize(failed, VectorN<Real, 2>{0, 0}), MultidimOptimizationInputError);
	CHECK(optimizer.getNumFuncEvals() == calls);
	CHECK(calls == 2);

	calls = 0;
	CHECK_THROWS_AS(optimizer.Minimize(constant, VectorN<Real, 2>{0, 0}, Real(0)),
		MultidimOptimizationInputError);
	CHECK(calls == 0);
	CHECK(optimizer.getNumFuncEvals() == 0);
	const auto reused = optimizer.Minimize(constant, VectorN<Real, 2>{0, 0});
	CHECK(reused.function_evaluations == calls);
	CHECK(calls == 3);
}

TEST_CASE("Callback accounting: Nelder Mead batch budget semantics remain unchanged",
	"[callback-accounting][NelderMead]")
{
	using namespace Optimization;
	for (const int budget : {1, 3, 4}) {
		int calls = 0;
		NelderMead optimizer(Precision::OptimizationGradientTolerance, budget);
		ScriptedObjective function([&]() {
			++calls;
			return calls <= 3 ? Real(calls - 1) : Real(3);
		});
		const auto result = optimizer.Minimize(function, VectorN<Real, 2>{0, 0});
		CHECK_FALSE(result.converged);
		CHECK(result.function_evaluations == calls);
		CHECK(optimizer.getNumFuncEvals() == calls);
		CHECK(result.iterations == (budget <= 3 ? 0 : 1));
		CHECK(result.iterations_used == result.iterations);
		CHECK(calls == (budget <= 3 ? 3 : 7));
	}
}

TEST_CASE("Callback accounting: ALM includes exhausted inner work and final objectives",
	"[callback-accounting][ALM]")
{
	using namespace Optimization;
	for (const bool bounded : {false, true}) {
		ALMConfig config;
		config.max_inner_iterations = 1;
		AugmentedLagrangian<2> optimizer(config);
		if (bounded) optimizer.SetBounds(-2, 2);
		int calls = 0;
		auto objective = [&](const VectorN<Real, 2>& x) {
			++calls;
			return 1 + x[0] * x[0] + x[1] * x[1];
		};
		const auto result = optimizer.Optimize(objective, VectorN<Real, 2>{1, 1});
		CHECK_FALSE(result.converged);
		CHECK(result.inner_iterations == 0);
		CHECK(result.outer_iterations == 1);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == (bounded ? 5 : 4));
	}
}

TEST_CASE("Callback accounting: scalar counter overflow is rejected before another invocation",
	"[callback-accounting][RootFinding]")
{
	int evaluations = std::numeric_limits<int>::max() - 1;
	int calls = 0;
	RealFunctionFromStdFunc function([&](Real) { ++calls; return Real(1); });
	CHECK(RootFinding::Detail::EvaluateFunction(function, 0, evaluations) == 1);
	CHECK(evaluations == std::numeric_limits<int>::max());
	CHECK(calls == 1);
	CHECK_THROWS_AS(RootFinding::Detail::EvaluateFunction(function, 0, evaluations), RootFindingError);
	CHECK(evaluations == std::numeric_limits<int>::max());
	CHECK(calls == 1);
}

TEST_CASE("Callback accounting: every Nelder Mead public input form resets failed work",
	"[callback-accounting][NelderMead]")
{
	using namespace Optimization;
	for (int form = 0; form < 5; ++form) {
		CAPTURE(form);
		int calls = 0;
		NelderMead optimizer;
		ScriptedObjective function([&]() { ++calls; return Real(1); });
		const auto completed = optimizer.Minimize(function, VectorN<Real, 2>{0, 0});
		CHECK(completed.function_evaluations == 3);
		calls = 0;
		if (form == 0) {
			CHECK_THROWS(optimizer.Minimize(function, VectorN<Real, 2>{0, 0}, Vector<Real>(1, 1)));
		} else if (form == 1) {
			CHECK_THROWS(optimizer.Minimize(function, VectorN<Real, 2>{0, 0}, Vector<Real>(2, Real(0))));
		} else if (form == 2) {
			CHECK_THROWS(optimizer.Minimize(function, Matrix<Real>(2, 2)));
		} else if (form == 3) {
			Matrix<Real> simplex(3, 2);
			simplex(0, 0) = std::numeric_limits<Real>::infinity();
			CHECK_THROWS(optimizer.Minimize(function, simplex));
		} else {
			optimizer.setMaxIter(0);
			CHECK_THROWS(optimizer.Minimize(function, VectorN<Real, 2>{0, 0}));
		}
		CHECK(calls == 0);
		CHECK(optimizer.getNumFuncEvals() == 0);
	}
}

TEST_CASE("Callback accounting: Nelder Mead convenience and maximizing APIs retain actual counts",
	"[callback-accounting][NelderMead]")
{
	using namespace Optimization;
	int calls = 0;
	ScriptedObjective function([&]() { ++calls; return Real(1); });
	MultidimOptimizationConfig config;
	config.max_iterations = 1;
	const VectorN<Real, 2> start{0, 0};
	for (int form = 0; form < 4; ++form) {
		calls = 0;
		const auto result = form == 0 ? NelderMeadMinimize(function, start)
			: form == 1 ? NelderMeadMaximize(function, start)
			: form == 2 ? NelderMeadMinimize(function, start, config)
			: NelderMeadMaximize(function, start, config);
		CHECK(result.converged);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == 3);
		CHECK(result.iterations == 0);
	}
}

TEST_CASE("Callback accounting: discarding rejected diagnostics preserves aggregate work",
	"[callback-accounting][RootFinding]")
{
	RootFinding::FindAllRealRootsConfig config;
	config.isolation.num_intervals = 1;
	config.keep_rejected_candidates = false;
	int calls = 0;
	auto function = [&](Real x) {
		if (++calls == 3) throw std::runtime_error("controlled refinement callback");
		return x * x - 2;
	};
	const auto result = RootFinding::FindAllRealRootsInInterval(function, 1, 2, config);
	CHECK(result.rejected.empty());
	CHECK(result.candidates_rejected == 1);
	CHECK_FALSE(result.complete);
	CHECK(result.function_evaluations == calls);
	CHECK(calls == 3);
}

TEST_CASE("Callback accounting: nonfinite Nelder Mead trial samples remain counted",
	"[callback-accounting][NelderMead]")
{
	using namespace Optimization;
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity()}) {
		int calls = 0;
		NelderMead optimizer(Precision::OptimizationGradientTolerance, 4);
		ScriptedObjective function([&]() {
			if (++calls > 3) return invalid;
			return Real(calls - 1);
		});
		const auto result = optimizer.Minimize(function, VectorN<Real, 2>{0, 0});
		CHECK_FALSE(result.converged);
		CHECK(result.function_evaluations == calls);
		CHECK(optimizer.getNumFuncEvals() == calls);
		CHECK(calls >= 4);
	}
}

TEST_CASE("Callback accounting: bounded Nelder Mead counts its final projected objective",
	"[callback-accounting][BoxNelderMead]")
{
	using namespace Optimization;
	for (const bool flat : {false, true}) {
		int calls = 0;
		ScalarFunctionFromStdFunc<2> objective([&](const VectorN<Real, 2>& x) {
			++calls;
			return flat ? Real(1) : 1 + x[0] * x[0] + x[1] * x[1];
		});
		BoundConstraints bounds(Vector<Real>{-2, -2}, Vector<Real>{2, 2});
		MultidimOptimizationConfig config;
		config.max_iterations = 1;
		const auto result = BoxConstrainedNelderMeadMinimize(objective, VectorN<Real, 2>{1, 1}, bounds, config);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == 4);
		CHECK(result.iterations == 0);
		CHECK(result.iterations_used == 0);
		CHECK(result.converged == flat);
		CHECK(result.fmin == (flat ? 1 : 3));
	}
}

TEST_CASE("Callback accounting: 1D bracket overloads count fresh objective samples",
	"[callback-accounting][optimization-accounting][Optimization]")
{
	const auto bracket = AccountingMinimumBracket();
	for (int family = 0; family < 3; ++family) {
		for (const int limit : {1, 100}) {
			CAPTURE(family, limit);
			int calls = 0, derivativeCalls = 0;
			RealFunctionFromStdFunc function([&](Real x) {
				++calls;
				const Real offset = x - Real(1) / 3;
				return 1 + offset * offset;
			});
			RealFunctionFromStdFunc derivative([&](Real x) {
				++derivativeCalls;
				return 2 * (x - Real(1) / 3);
			});
			const auto result = family == 0
				? Minimization::GoldenSectionSearch(function, bracket, Precision::OptimizationTolerance, limit)
				: family == 1
				? Minimization::BrentMinimize(function, bracket, Precision::OptimizationTolerance, limit)
				: Minimization::BrentMinimizeWithDeriv(function, derivative, bracket, Precision::OptimizationTolerance, limit);
			CHECK(result.function_evaluations == calls);
			CHECK(result.iterations_used == result.iterations);
			CHECK(result.iterations <= limit);
			CHECK(result.converged == (limit == 100));
			if (limit == 1) CHECK(calls == (family == 0 ? 3 : 2));
			if (family == 2) {
				CHECK(derivativeCalls > 0);
				CHECK(result.function_evaluations != calls + derivativeCalls);
			} else CHECK(derivativeCalls == 0);
		}
	}
}

TEST_CASE("Callback accounting: 1D interval overloads include bracketing without double counting",
	"[callback-accounting][optimization-accounting][Optimization]")
{
	for (int family = 0; family < 3; ++family) {
		for (const int limit : {1, 100}) {
			CAPTURE(family, limit);
			int calls = 0;
			RealFunctionFromStdFunc function([&](Real x) {
				++calls;
				const Real offset = x - Real(1) / 3;
				return 1 + offset * offset;
			});
			RealFunction derivative([](Real x) { return 2 * (x - Real(1) / 3); });
			const auto bracket = Minimization::BracketMinimum(function, -1, 0);
			REQUIRE(bracket.valid);
			const int bracketingCalls = calls;
			calls = 0;
			const auto refinement = family == 0
				? Minimization::GoldenSectionSearch(function, bracket, Precision::OptimizationTolerance, limit)
				: family == 1
				? Minimization::BrentMinimize(function, bracket, Precision::OptimizationTolerance, limit)
				: Minimization::BrentMinimizeWithDeriv(function, derivative, bracket, Precision::OptimizationTolerance, limit);
			const int refinementCalls = calls;
			calls = 0;
			const auto combined = family == 0
				? Minimization::GoldenSectionSearch(function, -1, 0, Precision::OptimizationTolerance, limit)
				: family == 1
				? Minimization::BrentMinimize(function, -1, 0, Precision::OptimizationTolerance, limit)
				: Minimization::BrentMinimizeWithDeriv(function, derivative, -1, 0, Precision::OptimizationTolerance, limit);
			CHECK(refinement.function_evaluations == refinementCalls);
			CHECK(combined.function_evaluations == calls);
			CHECK(calls == bracketingCalls + refinementCalls);
			CHECK(combined.iterations == refinement.iterations);
			CHECK(combined.xmin == refinement.xmin);
			CHECK(combined.fmin == refinement.fmin);
		}
	}
}

TEST_CASE("Callback accounting: 1D initial completion still counts its samples",
	"[callback-accounting][optimization-accounting][Optimization]")
{
	MinimumBracket bracket;
	bracket.ax = 1;
	bracket.bx = 1 + Precision::Eps;
	bracket.cx = 1 + 2 * Precision::Eps;
	bracket.valid = true;
	for (int family = 0; family < 3; ++family) {
		int calls = 0, derivativeCalls = 0;
		RealFunctionFromStdFunc function([&](Real) { ++calls; return Real(1); });
		RealFunctionFromStdFunc derivative([&](Real) { ++derivativeCalls; return Real(0); });
		const auto result = family == 0 ? Minimization::GoldenSectionSearch(function, bracket)
			: family == 1 ? Minimization::BrentMinimize(function, bracket)
			: Minimization::BrentMinimizeWithDeriv(function, derivative, bracket);
		CHECK(result.converged);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == (family == 0 ? 2 : 1));
		CHECK(result.iterations == (family == 0 ? 0 : 1));
		CHECK(derivativeCalls == (family == 2 ? 1 : 0));
	}
}

TEST_CASE("Callback accounting: failed 1D bracketing retains its work",
	"[callback-accounting][optimization-accounting][Optimization]")
{
	int calls = 0, derivativeCalls = 0;
	RealFunctionFromStdFunc decreasing([&](Real x) { ++calls; return -std::asinh(x); });
	RealFunctionFromStdFunc derivative([&](Real x) {
		++derivativeCalls;
		return -1 / std::hypot(Real(1), x);
	});
	const auto bracket = Minimization::BracketMinimum(decreasing, 0, 1);
	REQUIRE_FALSE(bracket.valid);
	const int expectedCalls = calls;
	RealFunctionFromStdFunc increasing([&](Real x) { ++calls; return std::asinh(x); });
	Minimization1DConfig config;
	for (int family = 0; family < 9; ++family) {
		calls = 0;
		const auto result = family == 0 ? Minimization::GoldenSectionSearch(decreasing, 0, 1)
			: family == 1 ? Minimization::BrentMinimize(decreasing, 0, 1)
			: family == 2 ? Minimization::BrentMinimizeWithDeriv(decreasing, derivative, 0, 1)
			: family == 3 ? Minimization::GoldenSectionSearch(decreasing, 0, 1, config)
			: family == 4 ? Minimization::BrentMinimize(decreasing, 0, 1, config)
			: family == 5 ? Minimization::BrentMinimizeWithDeriv(decreasing, derivative, 0, 1, config)
			: family == 6 ? Minimization::GoldenSectionMaximize(increasing, 0, 1)
			: family == 7 ? Minimization::BrentMaximize(increasing, 0, 1)
			: Minimization::BrentMaximize(increasing, 0, 1, config);
		CHECK_FALSE(result.converged);
		CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
		CHECK_FALSE(result.error_message.empty());
		CHECK(result.iterations == 0);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == expectedCalls);
		CHECK(derivativeCalls == 0);
	}
}

TEST_CASE("Scalar bracketing cannot certify overflowed extrapolation coordinates",
	"[callback-accounting][optimization-accounting][optimization-contracts]") {
	int calls = 0;
	RealFunctionFromStdFunc decreasing([&](Real x) {
		++calls;
		CHECK(std::isfinite(x));
		return -x;
	});
	const Real max = std::numeric_limits<Real>::max();
	const auto bracket = Minimization::BracketMinimum(decreasing, max / Real(4), max / Real(2));
	CHECK_FALSE(bracket.valid);
	CHECK(calls >= 2);
	CHECK(calls <= 4);
}

TEST_CASE("Callback accounting: 1D config callable and maximization APIs preserve objective totals",
	"[callback-accounting][optimization-accounting][Optimization]")
{
	for (int form = 0; form < 6; ++form) {
		for (const int limit : {1, 100}) {
			int calls = 0, derivativeCalls = 0;
			auto function = [&](Real x) {
				++calls;
				const Real offset = x - Real(1) / 3;
				const Real value = 1 + offset * offset;
				return form >= 3 ? -value : value;
			};
			auto derivative = [&](Real x) { ++derivativeCalls; return 2 * (x - Real(1) / 3); };
			Minimization1DConfig config;
			config.max_iterations = limit;
			const auto result = form == 0 ? Minimization::GoldenSectionSearch(function, -1, 0, config)
				: form == 1 ? Minimization::BrentMinimize(function, -1, 0, config)
				: form == 2 ? Minimization::BrentMinimizeWithDeriv(function, derivative, -1, 0, config)
				: form == 3 ? Minimization::GoldenSectionMaximize(function, -1, 0, config.tolerance, limit)
				: form == 4 ? Minimization::BrentMaximize(function, -1, 0, config.tolerance, limit)
				: Minimization::BrentMaximize(function, -1, 0, config);
			CHECK(result.function_evaluations == calls);
			CHECK(result.iterations <= limit);
			CHECK(result.converged == (limit == 100));
			CHECK((form >= 3 ? result.fmin < 0 : result.fmin > 0));
			if (form == 2) CHECK(derivativeCalls > 0);
			else CHECK(derivativeCalls == 0);
		}
	}
}

TEST_CASE("Callback accounting: Powell includes extrapolation and replacement line searches",
	"[callback-accounting][optimization-accounting][Powell]")
{
	using namespace Optimization;
	auto value = [](const VectorN<Real, 2>& x) {
		const Real sum = x[0] + x[1];
		return 1 + sum * sum + (x[0] * x[0] + x[1] * x[1]) / 100;
	};
	int referenceCalls = 0;
	ScalarFunctionFromStdFunc<2> reference([&](const VectorN<Real, 2>& x) { ++referenceCalls; return value(x); });
	const VectorN<Real, 2> start{-5, 5};
	VectorN<Real, 2> point = start;
	const Real initialValue = reference(point);
	CHECK(initialValue > 1);
	VectorN<Real, 2> firstDirection{1, 0}, secondDirection{0, 1};
	const Real firstMinimum = LineMinimizer::Minimize(reference, point, firstDirection);
	const Real secondMinimum = LineMinimizer::Minimize(reference, point, secondDirection);
	CHECK(secondMinimum < firstMinimum);
	const VectorN<Real, 2> extrapolated = Real(2) * point - start;
	const int sweepCalls = referenceCalls;

	std::vector<VectorN<Real, 2>> points;
	ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>& x) { points.push_back(x); return value(x); });
	Powell optimizer(Precision::OptimizationTolerance, 1);
	const auto result = optimizer.Minimize(function, start);
	REQUIRE(points.size() > static_cast<std::size_t>(sweepCalls + 1));
	CHECK(points[sweepCalls][0] == extrapolated[0]);
	CHECK(points[sweepCalls][1] == extrapolated[1]);
	CHECK(result.function_evaluations == static_cast<int>(points.size()));
	CHECK(result.iterations == 1);
	CHECK(result.iterations_used == 1);
	CHECK_FALSE(result.converged);
}

TEST_CASE("Callback accounting: Powell direct custom-direction and config APIs count callbacks",
	"[callback-accounting][optimization-accounting][Powell]")
{
	using namespace Optimization;
	for (const int limit : {1, 200}) {
		for (int form = 0; form < 3; ++form) {
			int calls = 0;
			ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>& x) {
				++calls;
				return 1 + x[0] * x[0] + 2 * x[1] * x[1];
			});
			MultidimOptimizationConfig config;
			config.max_iterations = limit;
			config.tolerance = Precision::OptimizationTolerance;
			Powell optimizer(config.tolerance, limit);
			Matrix<Real> directions(2, 2, Real(0));
			directions(0, 0) = directions(1, 1) = 1;
			const VectorN<Real, 2> start{-5, 5};
			const auto result = form == 0 ? optimizer.Minimize(function, start)
				: form == 1 ? optimizer.Minimize(function, start, directions)
				: PowellMinimize(function, start, config);
			CHECK(result.function_evaluations == calls);
			if (form != 2) CHECK(optimizer.getNumFuncEvals() == calls);
			CHECK(calls > result.iterations);
			CHECK(result.iterations <= limit);
			CHECK(result.iterations_used == result.iterations);
			CHECK(result.converged == (limit == 200));
		}
	}
}

TEST_CASE("Callback accounting: bounded Powell includes the final projected check",
	"[callback-accounting][optimization-accounting][BoxPowell]")
{
	using namespace Optimization;
	for (const int limit : {1, 200}) {
		int calls = 0;
		bool feasible = true;
		VectorN<Real, 2> last;
		ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>& x) {
			++calls;
			last = x;
			feasible = feasible && x[0] >= -2 && x[0] <= 2 && x[1] >= -2 && x[1] <= 2;
			const Real sum = x[0] + x[1];
			return 1 + sum * sum + (x[0] * x[0] + x[1] * x[1]) / 100;
		});
		BoundConstraints bounds(Vector<Real>{-2, -2}, Vector<Real>{2, 2});
		MultidimOptimizationConfig config;
		config.max_iterations = limit;
		config.tolerance = Precision::OptimizationTolerance;
		const auto result = BoxConstrainedPowellMinimize(function, VectorN<Real, 2>{-1, 1}, bounds, config);
		CHECK(result.function_evaluations == calls);
		CHECK(feasible);
		CHECK(result.xmin[0] == last[0]);
		CHECK(result.xmin[1] == last[1]);
		CHECK(result.iterations_used == result.iterations);
		CHECK(result.iterations <= limit);
	}
}

TEST_CASE("Callback accounting: 1D validation and callback exceptions preserve zero and partial work",
	"[callback-accounting][optimization-accounting][Optimization]")
{
	const auto bracket = AccountingMinimumBracket();
	for (int family = 0; family < 3; ++family) {
		int calls = 0;
		RealFunctionFromStdFunc function([&](Real x) { ++calls; return 1 + x * x; });
		RealFunction derivative([](Real x) { return 2 * x; });
		const Real invalidBound = std::numeric_limits<Real>::quiet_NaN();
		CHECK_THROWS_AS((family == 0 ? Minimization::GoldenSectionSearch(function, invalidBound, 1)
			: family == 1 ? Minimization::BrentMinimize(function, invalidBound, 1)
			: Minimization::BrentMinimizeWithDeriv(function, derivative, invalidBound, 1)), NumericInputError);
		CHECK(calls == 0);
		CHECK_THROWS_AS((family == 0 ? Minimization::GoldenSectionSearch(function, bracket, Real(0))
			: family == 1 ? Minimization::BrentMinimize(function, bracket, Real(0))
			: Minimization::BrentMinimizeWithDeriv(function, derivative, bracket, Real(0))), NumericInputError);
		CHECK(calls == 0);

		for (const bool throwing : {false, true}) {
			for (const int failAt : {1, 3}) {
				calls = 0;
				RealFunctionFromStdFunc failed([&](Real x) {
					if (++calls == failAt) {
						if (throwing) throw std::logic_error("controlled 1D callback");
						return std::numeric_limits<Real>::quiet_NaN();
					}
					const Real offset = x - Real(1) / 3;
					return 1 + offset * offset;
				});
				auto solve = [&]() {
					return family == 0 ? Minimization::GoldenSectionSearch(failed, bracket)
						: family == 1 ? Minimization::BrentMinimize(failed, bracket)
						: Minimization::BrentMinimizeWithDeriv(failed, derivative, bracket);
				};
				if (throwing) CHECK_THROWS_AS(solve(), std::logic_error);
				else CHECK_THROWS_AS(solve(), NumericInputError);
				CHECK(calls == failAt);
			}
		}
	}
}

TEST_CASE("Callback accounting: Powell exposes throwing counts in initial nested and extrapolated calls",
	"[callback-accounting][optimization-accounting][Powell]")
{
	using namespace Optimization;
	const auto value = [](const VectorN<Real, 2>& x) {
		const Real sum = x[0] + x[1];
		return 1 + sum * sum + (x[0] * x[0] + x[1] * x[1]) / 100;
	};
	const VectorN<Real, 2> start{-5, 5};
	int referenceCalls = 0;
	ScalarFunctionFromStdFunc<2> reference([&](const VectorN<Real, 2>& x) { ++referenceCalls; return value(x); });
	VectorN<Real, 2> point = start, first{1, 0}, second{0, 1};
	const Real initial = reference(point);
	const Real firstMinimum = LineMinimizer::Minimize(reference, point, first);
	const Real secondMinimum = LineMinimizer::Minimize(reference, point, second);
	CHECK(initial > firstMinimum);
	CHECK(firstMinimum > secondMinimum);
	const int extrapolationCall = referenceCalls + 1;
	for (const int failAt : {1, 2, 3, extrapolationCall, extrapolationCall + 1}) {
		int calls = 0;
		Powell optimizer(Precision::OptimizationTolerance, 1);
		ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>& x) {
			++calls;
			CHECK(optimizer.getNumFuncEvals() == calls);
			if (calls == failAt) throw std::logic_error("controlled Powell callback");
			return value(x);
		});
		CHECK_THROWS_AS(optimizer.Minimize(function, start), std::logic_error);
		CHECK(calls == failAt);
		CHECK(optimizer.getNumFuncEvals() == calls);
	}
}

TEST_CASE("Callback accounting: Powell resets counts for reuse invalid inputs and nonfinite samples",
	"[callback-accounting][optimization-accounting][Powell]")
{
	using namespace Optimization;
	Powell optimizer;
	int calls = 0;
	ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>& x) {
		++calls;
		return 1 + x[0] * x[0] + x[1] * x[1];
	});
	for (int form = 0; form < 4; ++form) {
		optimizer.setMaxIter(200);
		optimizer.setFtol(Precision::OptimizationTolerance);
		calls = 0;
		const auto completed = optimizer.Minimize(function, VectorN<Real, 2>{1, 1});
		CHECK(completed.converged);
		CHECK(completed.function_evaluations == calls);
		CHECK(optimizer.getNumFuncEvals() == calls);
		calls = 0;
		if (form == 0)
			CHECK_THROWS(optimizer.Minimize(function, VectorN<Real, 2>{std::numeric_limits<Real>::infinity(), 1}));
		else if (form == 1) {
			Matrix<Real> directions(2, 2);
			directions(0, 0) = std::numeric_limits<Real>::infinity();
			CHECK_THROWS(optimizer.Minimize(function, VectorN<Real, 2>{1, 1}, directions));
		} else if (form == 2) {
			optimizer.setFtol(0);
			CHECK_THROWS(optimizer.Minimize(function, VectorN<Real, 2>{1, 1}));
		} else {
			optimizer.setMaxIter(0);
			CHECK_THROWS(optimizer.Minimize(function, VectorN<Real, 2>{1, 1}));
		}
		CHECK(calls == 0);
		CHECK(optimizer.getNumFuncEvals() == 0);
	}
	optimizer.setMaxIter(200);
	for (const int failAt : {1, 2, 3}) {
		calls = 0;
		ScalarFunctionFromStdFunc<2> failed([&](const VectorN<Real, 2>& x) {
			if (++calls == failAt) return std::numeric_limits<Real>::quiet_NaN();
			return 1 + x[0] * x[0] + x[1] * x[1];
		});
		if (failAt == 1)
			CHECK_THROWS_AS(optimizer.Minimize(failed, VectorN<Real, 2>{1, 1}), MultidimOptimizationInputError);
		else
			CHECK_THROWS_AS(optimizer.Minimize(failed, VectorN<Real, 2>{1, 1}), NumericInputError);
		CHECK(calls == failAt);
		CHECK(optimizer.getNumFuncEvals() == calls);
	}
}

TEST_CASE("Callback accounting: bounded Powell final callback failures propagate unchanged",
	"[callback-accounting][optimization-accounting][BoxPowell]")
{
	using namespace Optimization;
	const BoundConstraints bounds(Vector<Real>{-2, -2}, Vector<Real>{2, 2});
	const VectorN<Real, 2> start{1, 1};
	int calls = 0;
	ScalarFunctionFromStdFunc<2> reference([&](const VectorN<Real, 2>& x) {
		++calls;
		return 1 + x[0] * x[0] + x[1] * x[1];
	});
	const auto completed = BoxConstrainedPowellMinimize(reference, start, bounds);
	REQUIRE(completed.converged);
	const int finalCall = calls;
	for (const bool throwing : {false, true}) {
		calls = 0;
		ScalarFunctionFromStdFunc<2> failed([&](const VectorN<Real, 2>& x) {
			if (++calls == finalCall) {
				if (throwing) throw std::logic_error("controlled final Powell callback");
				return std::numeric_limits<Real>::quiet_NaN();
			}
			return 1 + x[0] * x[0] + x[1] * x[1];
		});
		if (throwing) CHECK_THROWS_AS(BoxConstrainedPowellMinimize(failed, start, bounds), std::logic_error);
		else CHECK_THROWS_AS(BoxConstrainedPowellMinimize(failed, start, bounds), MultidimOptimizationInputError);
		CHECK(calls == finalCall);
	}
}

TEST_CASE("Callback accounting: 1D counter overflow rejects the next callback before invocation",
	"[callback-accounting][optimization-accounting][Optimization]")
{
	int calls = 0, evaluations = std::numeric_limits<int>::max() - 1;
	RealFunctionFromStdFunc function([&](Real) {
		++calls;
		CHECK(evaluations == std::numeric_limits<int>::max());
		return Real(1);
	});
	Minimization::Detail::CountingRealFunction counted(function, evaluations);
	CHECK(counted(0) == 1);
	CHECK(calls == 1);
	CHECK_THROWS_AS(counted(0), OptimizationError);
	CHECK(calls == 1);
	CHECK(evaluations == std::numeric_limits<int>::max());
}

TEST_CASE("Callback accounting: Powell counts both failed bracketing attempts in line-search retry",
	"[callback-accounting][optimization-accounting][Powell]")
{
	using namespace Optimization;
	int referenceCalls = 0;
	RealFunctionFromStdFunc reference([&](Real x) { ++referenceCalls; return -std::asinh(x); });
	const auto first = Minimization::BracketMinimum(reference, 0, 1);
	REQUIRE_FALSE(first.valid);
	const auto second = Minimization::BracketMinimum(reference, 0, Real(1) / 10);
	REQUIRE_FALSE(second.valid);
	int calls = 0;
	ScalarFunctionFromStdFunc<1> function([&](const VectorN<Real, 1>& x) {
		++calls;
		return -std::asinh(x[0]);
	});
	Powell optimizer;
	CHECK_THROWS_AS(optimizer.Minimize(function, VectorN<Real, 1>{0}), MultidimOptimizationError);
	CHECK(optimizer.getNumFuncEvals() == calls);
	CHECK(calls == 1 + referenceCalls);
}

TEST_CASE("Optimization contracts: scalar Brent success agrees with status and usability",
	"[optimization-contracts][optimization-accounting][Optimization]")
{
	for (int family = 0; family < 2; ++family) {
		for (bool narrow : {false, true}) {
			for (bool reversed : {false, true}) {
				auto bracket = UnitCenteredMinimumBracket(narrow);
				if (reversed) {
					std::swap(bracket.ax, bracket.cx);
					std::swap(bracket.fa, bracket.fc);
				}
				int calls = 0, derivativeCalls = 0;
				RealFunctionFromStdFunc function([&](Real x) { ++calls; return (x - 1) * (x - 1); });
				RealFunctionFromStdFunc derivative([&](Real x) { ++derivativeCalls; return 2 * (x - 1); });
				const auto result = family == 0 ? Minimization::BrentMinimize(function, bracket)
					: Minimization::BrentMinimizeWithDeriv(function, derivative, bracket);
				CheckBrentOutcome(result, calls, true);
				if (narrow) {
					CHECK(calls == 1);
					CHECK(result.iterations == 1);
					CHECK(derivativeCalls == family);
				}
			}
		}
	}
}

TEST_CASE("Optimization contracts: derivative Brent short-step success retains the accepted minimum",
	"[optimization-contracts][optimization-accounting][Optimization]")
{
	std::vector<Real> points;
	int derivativeCalls = 0;
	RealFunctionFromStdFunc function([&](Real x) { points.push_back(x); return (x - 1) * (x - 1); });
	RealFunctionFromStdFunc derivative([&](Real x) { ++derivativeCalls; return 2 * (x - 1); });
	const auto result = Minimization::BrentMinimizeWithDeriv(function, derivative,
		UnitCenteredMinimumBracket(), Precision::OptimizationTolerance, 2);
	CheckBrentOutcome(result, static_cast<int>(points.size()), true);
	REQUIRE(points.size() == 3);
	CHECK(points[0] == 1);
	CHECK(points[1] == Real(1) / 2);
	CHECK(points[2] != 1);
	CHECK((points[2] - 1) * (points[2] - 1) > result.fmin);
	CHECK(derivativeCalls == 2);
	CHECK(result.iterations == 2);
}

TEST_CASE("Optimization contracts: Brent interval callable and config wrappers preserve outcomes",
	"[optimization-contracts][optimization-accounting][Optimization]")
{
	for (int family = 0; family < 2; ++family) {
		for (int limit : {1, 100}) {
			for (int form = 0; form < 6; ++form) {
				int calls = 0;
				auto value = [&](Real x) { ++calls; return (x - 1) * (x - 1); };
				auto slope = [](Real x) { return 2 * (x - 1); };
				RealFunctionFromStdFunc function(value), derivative(slope);
				Minimization1DConfig config;
				config.max_iterations = limit;
				Minimization1DResult result;
				if (form == 0)
					result = family == 0 ? Minimization::BrentMinimize(function, UnitCenteredMinimumBracket(), config.tolerance, limit)
						: Minimization::BrentMinimizeWithDeriv(function, derivative, UnitCenteredMinimumBracket(), config.tolerance, limit);
				else if (form == 1)
					result = family == 0 ? Minimization::BrentMinimize(value, UnitCenteredMinimumBracket(), config.tolerance, limit)
						: Minimization::BrentMinimizeWithDeriv(value, slope, UnitCenteredMinimumBracket(), config.tolerance, limit);
				else if (form == 2)
					result = family == 0 ? Minimization::BrentMinimize(function, 0, 1, config.tolerance, limit)
						: Minimization::BrentMinimizeWithDeriv(function, derivative, 0, 1, config.tolerance, limit);
				else if (form == 3)
					result = family == 0 ? Minimization::BrentMinimize(value, 0, 1, config.tolerance, limit)
						: Minimization::BrentMinimizeWithDeriv(value, slope, 0, 1, config.tolerance, limit);
				else if (form == 4)
					result = family == 0 ? Minimization::BrentMinimize(function, 0, 1, config)
						: Minimization::BrentMinimizeWithDeriv(function, derivative, 0, 1, config);
				else
					result = family == 0 ? Minimization::BrentMinimize(value, 0, 1, config)
						: Minimization::BrentMinimizeWithDeriv(value, slope, 0, 1, config);
				CheckBrentOutcome(result, calls, limit != 1);
				CHECK(result.iterations <= limit);
				if (limit == 1) CHECK(result.iterations == 1);
			}
		}
	}
}

TEST_CASE("Optimization contracts: Brent maximization preserves successful and exhausted statuses",
	"[optimization-contracts][optimization-accounting][Optimization]")
{
	for (int limit : {1, 100}) {
		for (int form = 0; form < 4; ++form) {
			int calls = 0;
			auto value = [&](Real x) { ++calls; return 2 - (x - 1) * (x - 1); };
			RealFunctionFromStdFunc function(value);
			Minimization1DConfig config;
			config.max_iterations = limit;
			const auto result = form == 0 ? Minimization::BrentMaximize(function, 0, 1, config.tolerance, limit)
				: form == 1 ? Minimization::BrentMaximize(value, 0, 1, config.tolerance, limit)
				: form == 2 ? Minimization::BrentMaximize(function, 0, 1, config)
				: Minimization::BrentMaximize(value, 0, 1, config);
			CheckBrentOutcome(result, calls, limit != 1, 2);
			CHECK(result.iterations <= limit);
		}
	}
}

TEST_CASE("Optimization contracts: scalar defaults and failed brackets remain unusable",
	"[optimization-contracts][optimization-accounting][Optimization]")
{
	const Minimization1DResult notComputed;
	CHECK(notComputed.status == AlgorithmStatus::NotComputed);
	CHECK_FALSE(notComputed.converged);
	CHECK_FALSE(notComputed.IsSuccess());
	CHECK_FALSE(static_cast<bool>(notComputed));
	CHECK(notComputed.function_evaluations == 0);
	for (int family = 0; family < 2; ++family) {
		for (int form = 0; form < 4; ++form) {
			int calls = 0, derivativeCalls = 0;
			auto value = [&](Real x) { ++calls; return -std::asinh(x); };
			auto slope = [&](Real) { ++derivativeCalls; return Real(0); };
			RealFunctionFromStdFunc function(value), derivative(slope);
			Minimization1DConfig config;
			Minimization1DResult result;
			if (form == 0)
				result = family == 0 ? Minimization::BrentMinimize(function, 0, 1)
					: Minimization::BrentMinimizeWithDeriv(function, derivative, 0, 1);
			else if (form == 1)
				result = family == 0 ? Minimization::BrentMinimize(value, 0, 1)
					: Minimization::BrentMinimizeWithDeriv(value, slope, 0, 1);
			else if (form == 2)
				result = family == 0 ? Minimization::BrentMinimize(function, 0, 1, config)
					: Minimization::BrentMinimizeWithDeriv(function, derivative, 0, 1, config);
			else
				result = family == 0 ? Minimization::BrentMinimize(value, 0, 1, config)
					: Minimization::BrentMinimizeWithDeriv(value, slope, 0, 1, config);
			CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
			CHECK_FALSE(result.converged);
			CHECK_FALSE(result.IsSuccess());
			CHECK_FALSE(static_cast<bool>(result));
			CHECK_FALSE(result.error_message.empty());
			CHECK(result.function_evaluations == calls);
			CHECK(calls > 0);
			CHECK(result.iterations == 0);
			CHECK(derivativeCalls == 0);
		}
	}
}

TEST_CASE("Optimization contracts: Powell rejects safely oversized and rectangular directions",
	"[optimization-contracts][optimization-accounting][Powell]")
{
	using namespace Optimization;
	int calls = 0;
	ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>& x) {
		++calls;
		return 1 + x[0] * x[0] + 2 * x[1] * x[1];
	});
	Powell optimizer(Precision::OptimizationTolerance, 1);
	const VectorN<Real, 2> start{-5, 5};
	const int shapes[][2] = {{3, 3}, {3, 2}, {2, 3}};
	for (const auto& shape : shapes) {
		calls = 0;
		const auto previous = optimizer.Minimize(function, start);
		CHECK(previous.function_evaluations == calls);
		CHECK(optimizer.getNumFuncEvals() > 0);
		calls = 0;
		Matrix<Real> directions(shape[0], shape[1], Real(0));
		directions(0, 0) = directions(1, 1) = 1;
		CheckPowellShapeRejection(optimizer, function, start, directions);
		CHECK(calls == 0);
	}
}

TEST_CASE("Optimization contracts: Powell checks shape before inspecting nonfinite entries",
	"[optimization-contracts][optimization-accounting][Powell]")
{
	using namespace Optimization;
	int calls = 0;
	ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>&) { ++calls; return Real(1); });
	Powell optimizer;
	Matrix<Real> directions(3, 3, Real(0));
	directions(0, 0) = std::numeric_limits<Real>::quiet_NaN();
	CheckPowellShapeRejection(optimizer, function, VectorN<Real, 2>{1, 1}, directions);
	CHECK(calls == 0);
}

TEST_CASE("Optimization contracts: Powell rejects undersized and empty directions without entry access",
	"[optimization-contracts][optimization-accounting][Powell]")
{
	using namespace Optimization;
	int calls = 0;
	ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>&) { ++calls; return Real(1); });
	Powell optimizer;
	const int shapes[][2] = {{1, 1}, {1, 2}, {2, 1}, {1, 3}, {3, 1}, {0, 0}, {0, 2}, {2, 0}, {0, 3}, {3, 0}};
	for (const auto& shape : shapes) {
		Matrix<Real> directions(shape[0], shape[1], Real(7));
		CheckPowellShapeRejection(optimizer, function, VectorN<Real, 2>{1, 1}, directions);
		CHECK(calls == 0);
	}
	Matrix<Real> empty;
	CheckPowellShapeRejection(optimizer, function, VectorN<Real, 2>{1, 1}, empty);
	CHECK(calls == 0);
}

TEST_CASE("Optimization contracts: Powell retains finite-entry validation and recovers with custom directions",
	"[optimization-contracts][optimization-accounting][Powell]")
{
	using namespace Optimization;
	int calls = 0;
	ScalarFunctionFromStdFunc<2> function([&](const VectorN<Real, 2>& x) {
		++calls;
		return 1 + x[0] * x[0] + 2 * x[1] * x[1];
	});
	Powell optimizer(Precision::OptimizationTolerance, 1);
	const VectorN<Real, 2> start{-5, 5};
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		calls = 0;
		const auto previous = optimizer.Minimize(function, start);
		CHECK(previous.function_evaluations == calls);
		CHECK(optimizer.getNumFuncEvals() > 0);
		calls = 0;
		Matrix<Real> directions(2, 2, Real(0));
		directions(0, 0) = bad;
		bool rejected = false;
		try {
			const auto unexpected = optimizer.Minimize(function, start, directions);
			(void)unexpected;
		} catch (const MultidimOptimizationInputError& error) {
			rejected = true;
			CHECK(std::string(error.what()).find("direction matrix contains NaN/Inf") != std::string::npos);
		}
		CHECK(rejected);
		CHECK(calls == 0);
		CHECK(optimizer.getNumFuncEvals() == 0);
		if (std::isnan(bad)) CHECK(std::isnan(directions(0, 0)));
		else CHECK(directions(0, 0) == bad);
	}
	Matrix<Real> directions(2, 2, Real(0));
	directions(0, 1) = directions(1, 0) = 1;
	optimizer.setMaxIter(200);
	calls = 0;
	const auto recovered = optimizer.Minimize(function, start, directions);
	CHECK(recovered.converged);
	CHECK(recovered.status == AlgorithmStatus::Success);
	CHECK(recovered.function_evaluations == calls);
	CHECK(optimizer.getNumFuncEvals() == calls);
	CHECK(IsNearlyEqual(recovered.xmin[0], Real(0), Precision::OptimizationTolerance));
	CHECK(IsNearlyEqual(recovered.xmin[1], Real(0), Precision::OptimizationTolerance));
	CHECK(recovered.fmin == 1);
}

TEST_CASE("Optimization contracts: Powell direction-shape checks follow the function dimension",
	"[optimization-contracts][optimization-accounting][Powell]")
{
	using namespace Optimization;
	SECTION("One dimension") {
		int calls = 0;
		ScalarFunctionFromStdFunc<1> function([&](const VectorN<Real, 1>&) { ++calls; return Real(1); });
		Powell optimizer;
		const int shapes[][2] = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 2}};
		for (const auto& shape : shapes) {
			Matrix<Real> directions(shape[0], shape[1], Real(1));
			CheckPowellShapeRejection(optimizer, function, VectorN<Real, 1>{1}, directions);
			CHECK(calls == 0);
		}
	}
	SECTION("Three dimensions") {
		int calls = 0;
		ScalarFunctionFromStdFunc<3> function([&](const VectorN<Real, 3>&) { ++calls; return Real(1); });
		Powell optimizer;
		const int shapes[][2] = {{2, 2}, {2, 3}, {3, 2}, {3, 4}, {4, 3}, {0, 0}};
		for (const auto& shape : shapes) {
			Matrix<Real> directions(shape[0], shape[1], Real(1));
			CheckPowellShapeRejection(optimizer, function, VectorN<Real, 3>{1, -2, 3}, directions);
			CHECK(calls == 0);
		}
	}
}
