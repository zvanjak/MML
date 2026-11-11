// TOLERANCE_POLICY in mml/core (task MinimalMathLibrary-wewi.4.3): literals became named Precision::
// constants, and GMRES's lucky-breakdown test is relative to the Arnoldi start norm (it used to be an
// absolute 1000*eps, so any system scaled down far enough "broke down" on the first iteration).
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
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
		cfg.tolerance = 1e-10;
		cfg.maxIterations = 200;

		SparseSolvers::GMRES<Real> solver(n);
		solver.setConfig(cfg);
		Vector<Real> xPlain = x0;
		auto plain = solver.solve(A, b, xPlain);
		REQUIRE(plain.converged);

		// scale the whole system by 1e-20: identical solution, must take the same path (no fake breakdown)
		const Real s = 1e-20;
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
		REQUIRE(scaled.converged);
		REQUIRE(scaled.iterations_used == plain.iterations_used);
		for (int i = 0; i < n; ++i)
			REQUIRE(xScaled[i] == Catch::Approx(xPlain[i]).epsilon(1e-8));
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
