// TOLERANCE_POLICY in mml/core (task MinimalMathLibrary-wewi.4.3): literals became named Precision::
// constants, and GMRES's lucky-breakdown test is relative to the Arnoldi start norm (it used to be an
// absolute 1000*eps, so any system scaled down far enough "broke down" on the first iteration).
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/SparseMatrix/SparseMatrix.h>
#include <mml/core/SparseSolvers/GMRES.h>
#include <mml/core/SparseSolvers/ConjugateGradient.h>
#include <mml/core/SparseSolvers/IterativeSolverBase.h>
#include <mml/core/Integration/GaussKronrod.h>
#include <mml/core/Integration/MonteCarloIntegration.h>
#include <mml/core/Integration/Integration1D.h>
#endif

using namespace MML;

namespace MML::Tests::Core::ToleranceCoreTests
{
	TEST_CASE("ToleranceCore_GMRES_breakdown_is_scale_invariant", "[tolerance][core][GMRES]")
	{
		const int n = 30;
		auto A = SparseMatrix::Laplacian1D<Real>(n);
		Vector<Real> b(n, 1.0), x0(n, 0.0);

		SparseSolvers::SolverConfig<Real> cfg;
		cfg.tolerance = TOL(1e-10, Precision::CoarseTolerance);
		cfg.absTolerance = 0;
		cfg.maxIterations = 200;

		SparseSolvers::GMRES<Real> solver(n);
		solver.setConfig(cfg);
		Vector<Real> xPlain = x0;
		auto plain = solver.solve(A, b, xPlain);
		INFO(plain.error_message);
		REQUIRE(plain.converged);
		REQUIRE(plain.iterations_used > 1);
		REQUIRE(plain.relativeResidual <= cfg.tolerance);

		// scale the whole system by 1e-20: identical solution, must take the same path (no fake breakdown)
		const Real s = MML::Real(1e-20);
		SparseMatrix::SparseMatrixCOO<Real> coo(n, n);
		for (int i = 0; i < n; ++i)
			for (int k = A.rowPointers()[i]; k < A.rowPointers()[i + 1]; ++k)
				coo.addEntry(i, A.colIndices()[k], s * A.values()[k]);
		SparseMatrix::SparseMatrixCSR<Real> As(coo);
		Vector<Real> bs(n, s);
		cfg.absTolerance = 0.0;                      // the absolute floor must not short-circuit the scaled system either
		SparseSolvers::GMRES<Real> solverS(n);
		solverS.setConfig(cfg);
		Vector<Real> xScaled = x0;
		auto scaled = solverS.solve(As, bs, xScaled);
		INFO(scaled.error_message);
		REQUIRE(scaled.converged);
		REQUIRE(scaled.relativeResidual <= cfg.tolerance);
		REQUIRE(scaled.iterations_used == plain.iterations_used);
		REQUIRE(SparseSolvers::norm2(A * xPlain - b) / SparseSolvers::norm2(b) <= cfg.tolerance);
		REQUIRE(SparseSolvers::norm2(As * xScaled - bs) / SparseSolvers::norm2(bs) <= cfg.tolerance);
		for (int i = 0; i < n; ++i) {
			const Real expected = Real((i + 1) * (n - i)) / 2;
			const Real budget = TOL(1e-8, Precision::RoundoffZeroThreshold) * std::abs(expected);
			REQUIRE(IsWithinAbs(xPlain[i], expected, budget));
			REQUIRE(IsWithinAbs(xScaled[i], expected, budget));
			REQUIRE(IsWithinAbs(xScaled[i], xPlain[i], 2 * budget));
		}
	}

	TEST_CASE("ToleranceCore_named_defaults", "[tolerance][core]")
	{
		SparseSolvers::SolverConfig<Real> cfg;
		REQUIRE(cfg.tolerance == Precision::IterativeSolverTolerance);
		REQUIRE(cfg.absTolerance == Precision::RoundoffZeroThreshold);

		MonteCarloConfig mc;
		REQUIRE(mc.target_error == Precision::MonteCarloTargetRelativeError);

		static_assert(Precision::IntegrationRoundoffFloorScale == 50 * Precision::Eps);
		static_assert(Precision::KrylovBreakdownThreshold == 1000 * Precision::Eps);
		static_assert(Precision::SVDCompletionTolerancePerDimension == 32 * Precision::Eps);

		// Gauss-Kronrod default tolerance is the integration family constant
		Real v = Integration::IntegrateGaussKronrod([](Real x) { return std::sin(x); }, 0.0, Constants::PI);
		REQUIRE(v == Catch::Approx(2.0).epsilon(10 * Precision::IntegrationTolerance));
	}
}
