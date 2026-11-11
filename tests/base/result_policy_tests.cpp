///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        result_policy_tests.cpp                                             ///
///  Description: Contract tests for docs/RESULT_POLICY.md vocabulary: AlgorithmStatus ///
///               Completed, NoEstimate / HasEstimate, IsUsable, result helpers, and  ///
///               the explicit-conversion rule.                                       ///
///                                                                                   ///
///  MinimalMathLibrary-wewi.3.1                                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml/base/AlgorithmTypes.h>
#include <mml/core/Integration.h>
#include <mml/core/Integration/MonteCarloIntegration.h>
#include <mml/core/ComplexAnalysis.h>
#include <mml/core/SparseSolvers/IterativeSolverBase.h>
#include <mml/base/Function.h>
#include <mml/algorithms/Optimization/OptimizationMultidim.h>
#include <mml/algorithms/Optimization/LinearProgramming.h>
#include <mml/core/SparseSolvers/IterativeSolvers.h>
#include <mml/tools/DataLoader.h>
#include <mml/tools/Persistence.h>
#include <filesystem>
#include <fstream>

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include <type_traits>

using namespace MML;

namespace MML::Tests::Base::ResultPolicyTests {

	TEST_CASE("IOResult is shared, explicit and preserves file and offset context", "[result-policy][io][wewi.3.5]") {
		STATIC_REQUIRE_FALSE(std::is_convertible_v<IOResult, bool>);
		STATIC_REQUIRE(std::is_constructible_v<bool, IOResult>);
		STATIC_REQUIRE(
			std::is_same_v<decltype(Persistence::SaveJson(std::declval<std::ostream&>(), std::declval<const Vector<Real>&>())), IOResult>);
		STATIC_REQUIRE(
			std::is_same_v<decltype(Data::LoadCSVSafe(std::declval<const std::string&>(), std::declval<Data::Dataset&>())), IOResult>);
		const IOResult unopened;
		REQUIRE_FALSE(unopened);
		REQUIRE(unopened.offset == -1);
		REQUIRE(unopened.offset_kind == IOOffsetKind::Unknown);
		const auto failed = SerializeFailure(IOError::READ_FAILED, "truncated").At(17).WithFile("sample.mmlb");
		REQUIRE_FALSE(failed);
		REQUIRE(failed.filename == "sample.mmlb");
		REQUIRE(failed.offset == 17);
		REQUIRE(failed.offset_kind == IOOffsetKind::Byte);

		const auto json = Persistence::ParseJson("{\"value\":}");
		REQUIRE_FALSE(json.result);
		REQUIRE(json.result.error == IOError::MALFORMED_INPUT);
		REQUIRE(json.result.offset == 9);
		REQUIRE(json.result.offset_kind == IOOffsetKind::Byte);

		Vector<Real> vector;
		const std::string missing = (std::filesystem::temp_directory_path() / "mml_missing_result_contract" / "vector.mmlb").string();
		const auto load = Persistence::LoadBinary(missing, vector);
		REQUIRE(load.error == IOError::FILE_NOT_OPENED);
		REQUIRE(load.filename == missing);
		REQUIRE(load.offset == -1);
	}

	TEST_CASE("Dataset and persistence failures carry actual file locations", "[result-policy][io][wewi.3.5]") {
		const auto path = std::filesystem::temp_directory_path() /
											("mml_io_contract_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".csv");
		const std::string filename = path.string();
		{
			std::ofstream file(path);
			file << "a,b\n\n1,2\n3\n";
			REQUIRE(file.good());
		}
		Data::Dataset dataset;
		dataset.name = "unchanged";
		Data::CSVLoadOptions options;
		options.strictSchema = true;
		const auto schema = Data::LoadCSV(filename, dataset, options);
		REQUIRE_FALSE(schema);
		REQUIRE(schema.error == IOError::SCHEMA_MISMATCH);
		REQUIRE(schema.filename == filename);
		REQUIRE(schema.offset == 4);
		REQUIRE(schema.offset_kind == IOOffsetKind::Line);
		REQUIRE(dataset.name == "unchanged");
		{
			std::ofstream file(path);
			file << "a,b\n1,2\n3,4\n";
		}
		const auto csv = Data::LoadData(filename, dataset);
		REQUIRE(csv);
		REQUIRE(csv.filename == filename);
		REQUIRE(dataset.NumRows() == 2);
		REQUIRE(dataset.NumColumns() == 2);
		{
			std::ofstream file(path);
			file << "[{\"a\":1},{\"a\":2}]";
		}
		const auto jsonSuccess = Data::LoadJSONSafe(filename, dataset);
		REQUIRE(jsonSuccess);
		REQUIRE(jsonSuccess.filename == filename);
		REQUIRE(dataset.NumRows() == 2);
		REQUIRE(dataset.NumColumns() == 1);
		dataset.name = "unchanged";
		{
			std::ofstream file(path);
			file << "[{\"a\":}]";
		}
		const auto json = Data::LoadJSONSafe(filename, dataset);
		REQUIRE_FALSE(json);
		REQUIRE(json.filename == filename);
		REQUIRE(json.offset == 6);
		REQUIRE(json.offset_kind == IOOffsetKind::Byte);
		REQUIRE(dataset.name == "unchanged");
		{
			std::ofstream file(path);
			file << "2 2\n1 2\n3";
		}
		Matrix<Real> matrix;
		const auto truncated = Persistence::LoadMatrixFromFile(filename, matrix);
		REQUIRE_FALSE(truncated);
		REQUIRE(truncated.filename == filename);
		REQUIRE(truncated.offset == 11);
		REQUIRE(truncated.offset_kind == IOOffsetKind::Byte);
		const auto saved = Persistence::SaveJson(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)}, filename);
		REQUIRE(saved);
		REQUIRE(saved.filename == filename);
		Vector<Real> loaded;
		const auto load = Persistence::LoadJson(filename, loaded);
		REQUIRE(load);
		REQUIRE(load.filename == filename);
		REQUIRE(loaded == Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)});
		const auto binarySaved = Persistence::SaveBinary(loaded, filename);
		REQUIRE(binarySaved);
		const auto truncatedSize = std::filesystem::file_size(path) - 2;
		std::filesystem::resize_file(path, truncatedSize);
		const auto binaryFailure = Persistence::LoadBinary(filename, loaded);
		REQUIRE_FALSE(binaryFailure);
		REQUIRE(binaryFailure.error == IOError::TRUNCATED_INPUT);
		REQUIRE(binaryFailure.filename == filename);
		REQUIRE(binaryFailure.offset == static_cast<std::int64_t>(truncatedSize));
		REQUIRE(binaryFailure.offset_kind == IOOffsetKind::Byte);
		{
			std::ofstream file(path);
			file.setstate(std::ios::badbit);
			const auto writeFailure = FinishOutputFile(file, filename, IOResult::Success());
			REQUIRE_FALSE(writeFailure);
			REQUIRE(writeFailure.error == IOError::WRITE_FAILED);
			REQUIRE(writeFailure.filename == filename);
		}
		REQUIRE(std::filesystem::remove(path));
	}

	TEST_CASE("Optimization and Krylov share the iterative result contract", "[result-policy][wewi.3.4]") {
		using namespace Optimization;
		STATIC_REQUIRE(std::is_base_of_v<IterativeResultBase, Minimization1DResult>);
		STATIC_REQUIRE(std::is_base_of_v<IterativeResultBase, MultidimOptimizationResult>);
		STATIC_REQUIRE(std::is_base_of_v<IterativeResultBase, LinearProgrammingResult>);
		STATIC_REQUIRE(std::is_base_of_v<IterativeResultBase, SparseSolvers::KrylovResult<Real>>);
		REQUIRE_FALSE(static_cast<bool>(SparseSolvers::KrylovResult<Real>{}));
		STATIC_REQUIRE_FALSE(std::is_convertible_v<LinearProgrammingResult, bool>);
		STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, OptimizationError>);
		STATIC_REQUIRE(ToAlgorithmStatus(LPStatus::MaxIterations) == AlgorithmStatus::MaxIterationsExceeded);
		STATIC_REQUIRE(ToAlgorithmStatus(LPStatus::Unbounded) == AlgorithmStatus::AlgorithmSpecificFailure);

		RealFunction f([](Real x) { return (x - static_cast<MML::Real>(2.0)) * (x - static_cast<MML::Real>(2.0)); });
		const auto minimum = Minimization::BrentMinimize(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(4.0));
		REQUIRE(minimum.converged);
		REQUIRE(minimum.iterations_used == minimum.iterations);
		REQUIRE(minimum.iterations_used > 0);
		REQUIRE_FALSE(minimum.algorithm_name.empty());

		struct CountedConstant : IScalarFunction<2> {
			mutable int calls = 0;
			Real operator()(const VectorN<Real, 2>&) const override {
				++calls;
				return static_cast<MML::Real>(3.0);
			}
		} constant;
		const auto simplex = NelderMead().Minimize(constant, VectorN<Real, 2>{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)});
		REQUIRE(simplex.converged);
		REQUIRE(simplex.iterations_used == 0);
		REQUIRE(simplex.iterations == simplex.iterations_used);
		REQUIRE(simplex.function_evaluations == constant.calls);
		REQUIRE(simplex.function_evaluations == 3);
		REQUIRE(simplex.achieved_tolerance == static_cast<MML::Real>(0.0));

		LinearProgram lp(2);
		lp.SetObjective({static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)}, LPObjective::Minimize);
		lp.AddConstraint({static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)}, LPConstraintType::GreaterEqual, static_cast<MML::Real>(1.0));
		LPConfig config;
		config.maxIterations = 0;
		SimplexSolver solver(config);
		for (const auto result : {solver.Solve(lp), solver.SolveWithSensitivity(lp)}) {
			REQUIRE(result.lp_status == LPStatus::MaxIterations);
			REQUIRE(result.status == AlgorithmStatus::MaxIterationsExceeded);
			REQUIRE_FALSE(result.converged);
			REQUIRE_FALSE(static_cast<bool>(result));
			REQUIRE(result.x.size() == 2);
			REQUIRE(result.achieved_tolerance == Catch::Approx(static_cast<MML::Real>(1.0)));
			REQUIRE_FALSE(result.error_message.empty());
			REQUIRE_FALSE(result.algorithm_name.empty());
			REQUIRE_FALSE(HasEstimate(result.basisConditionNumber));
		}
		const auto optimum = SimplexSolver().Solve(lp);
		REQUIRE(optimum.status == AlgorithmStatus::Success);
		REQUIRE(optimum.converged);
		REQUIRE(optimum.error_message.empty());
		REQUIRE(optimum.iterations_used == optimum.iterations);
		REQUIRE(optimum.achieved_tolerance <= LPConfig().tolerance);
	}

	TEST_CASE("Krylov diagnostics describe the returned iterate on every outcome", "[result-policy][SparseSolvers][wewi.3.4]") {
		auto A = SparseMatrix::Laplacian1D<Real>(9);
		Vector<Real> b(9, static_cast<MML::Real>(1.0)), x(9, static_cast<MML::Real>(0.0));
		SparseSolvers::SolverConfig<Real> config;
		config.maxIterations = 1;
		const auto limited = SparseSolvers::solveCG(A, b, x, config);
		REQUIRE(limited.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE_FALSE(limited.converged);
		REQUIRE(limited.iterations_used == 1);
		REQUIRE(limited.residualNorm == Catch::Approx((b - A * x).NormL2()));
		REQUIRE(limited.achieved_tolerance == limited.residualNorm);
		REQUIRE_FALSE(limited.algorithm_name.empty());
		REQUIRE_FALSE(limited.error_message.empty());

		const Vector<Real> wrongSize(8);
		const auto invalid = SparseSolvers::solveCG(A, wrongSize, x);
		REQUIRE(invalid.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(HasEstimate(invalid.achieved_tolerance));
		REQUIRE(std::isnan(invalid.residualNorm));

		SparseSolvers::GMRES<Real> gmres(0);
		REQUIRE(gmres.solve(A, b, x).status == AlgorithmStatus::InvalidInput);
		b = Vector<Real>(9, config.absTolerance / static_cast<MML::Real>(4.0));
		const auto tiny = SparseSolvers::solveCG(A, b, x);
		REQUIRE(tiny.converged);
		REQUIRE(tiny.residualNorm == Catch::Approx(b.NormL2()).margin(static_cast<MML::Real>(0.0)));
		REQUIRE(tiny.residualNorm > static_cast<MML::Real>(0.0));
		config.absTolerance = static_cast<MML::Real>(0.0);
		b = Vector<Real>(9, static_cast<MML::Real>(0.0));
		x = Vector<Real>(9, static_cast<MML::Real>(1.0));
		for (const auto zero : {SparseSolvers::solveCG(A, b, x, config),
				 SparseSolvers::solveBiCGSTAB(A, b, x, config),
				 SparseSolvers::solveGMRES(A, b, x, 3, config)}) {
			REQUIRE(zero.converged);
			REQUIRE(zero.residualNorm == static_cast<MML::Real>(0.0));
			REQUIRE(zero.relativeResidual == static_cast<MML::Real>(0.0));
		}
	}

TEST_CASE("AlgorithmStatus::Completed is a distinct, usable, non-converged status", "[result-policy][wewi.3.1]")
{
    REQUIRE(to_string(AlgorithmStatus::Completed) == "Completed");
    REQUIRE(static_cast<int>(AlgorithmStatus::Completed) != static_cast<int>(AlgorithmStatus::Success));

    STATIC_REQUIRE(IsUsable(AlgorithmStatus::Success));
    STATIC_REQUIRE(IsUsable(AlgorithmStatus::Completed));
    STATIC_REQUIRE_FALSE(IsUsable(AlgorithmStatus::MaxIterationsExceeded));
    STATIC_REQUIRE_FALSE(IsUsable(AlgorithmStatus::Stalled));
    STATIC_REQUIRE_FALSE(IsUsable(AlgorithmStatus::ToleranceUnachievable));
    STATIC_REQUIRE_FALSE(IsUsable(AlgorithmStatus::NumericalInstability));
    STATIC_REQUIRE_FALSE(IsUsable(AlgorithmStatus::SingularMatrix));
    STATIC_REQUIRE_FALSE(IsUsable(AlgorithmStatus::InvalidInput));
    STATIC_REQUIRE_FALSE(IsUsable(AlgorithmStatus::AlgorithmSpecificFailure));
}

TEST_CASE("NoEstimate is NaN and HasEstimate distinguishes it from any computed value", "[result-policy][wewi.3.1]")
{
    REQUIRE(std::isnan(NoEstimate));
    REQUIRE_FALSE(HasEstimate(NoEstimate));
    REQUIRE(HasEstimate(static_cast<MML::Real>(0.0)));            // a computed zero is a real estimate
    REQUIRE(HasEstimate(static_cast<MML::Real>(1e-300)));
    REQUIRE(HasEstimate(std::numeric_limits<Real>::infinity()));   // computed, if useless
    REQUIRE_FALSE(NoEstimate == NoEstimate);    // cannot be mistaken for a number by accident
}

TEST_CASE("Result helpers populate status, converged and estimates consistently", "[result-policy][wewi.3.1]")
{
    struct R : IterativeResultBase { Real payload = 0.0; };

    auto ok = MakeSuccessResult<R>(5, static_cast<MML::Real>(1e-12), "Algo");
    REQUIRE(ok.status == AlgorithmStatus::Success);
    REQUIRE(ok.converged);
    REQUIRE(ok.IsSuccess());
    REQUIRE(static_cast<bool>(ok));
    REQUIRE(HasEstimate(ok.achieved_tolerance));
    REQUIRE(ok.error_message.empty());
    REQUIRE(ok.iterations_used == 5);

    auto done = MakeCompletedResult<R>("Direct", 1);
    REQUIRE(done.status == AlgorithmStatus::Completed);
    REQUIRE_FALSE(done.converged);
    REQUIRE(done.IsSuccess());                 // usable
    REQUIRE(static_cast<bool>(done));
    REQUIRE_FALSE(HasEstimate(done.achieved_tolerance));
    REQUIRE(done.error_message.empty());
    REQUIRE(done.iterations_used == 1);

    auto bad = MakeFailureResult<R>(AlgorithmStatus::MaxIterationsExceeded, "Failed after 100 iterations", 100, "Algo");
    REQUIRE_FALSE(bad.converged);
    REQUIRE_FALSE(bad.IsSuccess());
    REQUIRE_FALSE(static_cast<bool>(bad));
    REQUIRE_FALSE(bad.error_message.empty());

    // operator bool is explicit: `bool b = result;` must not compile
    STATIC_REQUIRE_FALSE(std::is_convertible_v<R, bool>);
    STATIC_REQUIRE(std::is_constructible_v<bool, R>);

    struct E : EvaluationResultBase { Real value = 0.0; };
    E completed;
    completed.status = AlgorithmStatus::Completed;
    REQUIRE(completed.IsSuccess());
    E failed = MakeEvaluationFailureResult<E>(AlgorithmStatus::NumericalInstability, "NaN encountered", "Eval");
    REQUIRE_FALSE(failed.IsSuccess());
}

TEST_CASE("Integration results follow the policy: no implicit conversion, fixed rules are Completed", "[result-policy][wewi.3.1]")
{
    STATIC_REQUIRE_FALSE(std::is_convertible_v<IntegrationResult, Real>);
    STATIC_REQUIRE(std::is_constructible_v<Real, IntegrationResult>);

    auto fixed = IntegrationResult::FixedRule(static_cast<MML::Real>(2.5));
    REQUIRE(fixed.value == static_cast<MML::Real>(2.5));
    REQUIRE_FALSE(fixed.converged);
    REQUIRE_FALSE(fixed.error_estimated);
    REQUIRE_FALSE(HasEstimate(fixed.error_estimate));

    IntegrationDetailedResult d;
    d.value = static_cast<MML::Real>(1.0);
    d.MarkFixedRule(10);
    REQUIRE(d.status == AlgorithmStatus::Completed);
    REQUIRE(d.IsSuccess());
    REQUIRE_FALSE(d.converged);
    REQUIRE_FALSE(HasEstimate(d.error_estimate));
    REQUIRE(d.function_evaluations == 10);
}

TEST_CASE("Plain integration returns the value and throws when the result is not usable (wewi.3.3)", "[result-policy][wewi.3.3]")
{
    RealFunctionFromStdFunc f([](Real x) { return x * x; });

    // plain: a Real, equal to the Detailed value
    Real a = IntegrateTrap(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0));
    STATIC_REQUIRE(std::is_same_v<decltype(IntegrateTrap(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0))), Real>);
    auto d = IntegrateTrapDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0));
    REQUIRE(a == d.value);
    REQUIRE(d.converged);
    REQUIRE(d.IsUsable());

    // Detailed with the plain argument list (eps) and with the config-bearing list both exist
    auto d2 = IntegrateTrapDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), static_cast<MML::Real>(1e-6));
    auto d3 = IntegrateTrapDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), IntegrationConfig{}, static_cast<MML::Real>(1e-6));
    REQUIRE(d2.value == d3.value);

    // unreachable tolerance: Detailed says so, plain throws instead of returning a partial value
    auto bad = IntegrateTrapDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-8), static_cast<MML::Real>(1e-12));
    REQUIRE_FALSE(bad.converged);
    REQUIRE_FALSE(bad.IsUsable());
    REQUIRE(bad.status == AlgorithmStatus::MaxIterationsExceeded);
    REQUIRE_THROWS_AS(IntegrateTrap(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-8), static_cast<MML::Real>(1e-12)), ConvergenceError);
    REQUIRE_THROWS_AS(bad.ValueOrThrow("test"), ConvergenceError);

    // fixed rule: Completed is usable, so the plain call returns
    REQUIRE(IntegrateGauss10(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)) == Catch::Approx(1.0 / 3.0));
    REQUIRE(IntegrateGauss10Detailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)).IsUsable());

    // exactly-zero integrand and zero-width interval converge (criterion is <=)
    RealFunctionFromStdFunc zero([](Real) { return static_cast<MML::Real>(0.0); });
    REQUIRE(IntegrateTrap(zero, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)) == static_cast<MML::Real>(0.0));
    REQUIRE(IntegrateTrapDetailed(zero, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)).converged);
    REQUIRE(IntegrateTrap(f, static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)) == static_cast<MML::Real>(0.0));
}

TEST_CASE("All integration result types are one BasicIntegrationResult family (wewi.3.3)", "[result-policy][wewi.3.3]")
{
    STATIC_REQUIRE(std::is_same_v<IntegrationResult, BasicIntegrationResult<Real>>);
    STATIC_REQUIRE(std::is_same_v<IntegrationDetailedResult, IntegrationResult>);
    STATIC_REQUIRE(std::is_same_v<ComplexAnalysis::ContourIntegrationResult, BasicIntegrationResult<Complex>>);
    STATIC_REQUIRE(std::is_base_of_v<IntegrationResult, Integration::GKResult>);
    STATIC_REQUIRE(std::is_base_of_v<IntegrationResult, Integration::AdaptiveResult2D>);
    STATIC_REQUIRE(std::is_base_of_v<IntegrationResult, Integration::AdaptiveResult3D>);
    STATIC_REQUIRE(std::is_base_of_v<IntegrationResult, MonteCarloResult>);
    STATIC_REQUIRE(std::is_base_of_v<EvaluationResultBase, IntegrationResult>);

    // conversions are explicit everywhere
    STATIC_REQUIRE_FALSE(std::is_convertible_v<IntegrationResult, Real>);
    STATIC_REQUIRE_FALSE(std::is_convertible_v<ComplexAnalysis::ContourIntegrationResult, Complex>);
		STATIC_REQUIRE_FALSE(std::is_convertible_v<SparseSolvers::KrylovResult<Real>, bool>);
		STATIC_REQUIRE(std::is_constructible_v<bool, SparseSolvers::KrylovResult<Real>>);

		// the field-wise constructor derives the status
    REQUIRE(IntegrationResult(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-9), 5, true).status == AlgorithmStatus::Success);
    REQUIRE(IntegrationResult(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-3), 20, false).status == AlgorithmStatus::MaxIterationsExceeded);
    REQUIRE(IntegrationResult(static_cast<MML::Real>(1.0), NoEstimate, 1, false, false).status == AlgorithmStatus::Completed);
    REQUIRE_FALSE(HasEstimate(IntegrationResult(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), 1, false, false).error_estimate));   // no estimate -> NoEstimate regardless
}

} // namespace MML::Tests::Base::ResultPolicyTests
