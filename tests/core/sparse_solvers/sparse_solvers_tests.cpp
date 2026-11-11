#include <catch2/catch_all.hpp>

#include <mml/base/SparseMatrix/SparseMatrix.h>
#include <mml/core/SparseSolvers/IterativeSolvers.h>

#include <vector>

namespace MML::Tests::Core::SparseSolversTests
{
    TEST_CASE("Sparse matrix Laplacian solves with conjugate gradient", "[SparseMatrix][SparseSolvers]")
    {
        auto A = SparseMatrix::Laplacian1D<double>(3);
				Vector<double> b{1.0, 2.0, 1.0};
				Vector<double> x(3, 0.0);

				SparseSolvers::SolverConfig<double> config;
        config.setTolerance(1e-12).setMaxIterations(50).setValidateSPD(true);

        auto result = SparseSolvers::solveCG(A, b, x, config);

				REQUIRE(result.converged);
				REQUIRE(x[0] == Catch::Approx(2.0).epsilon(1e-12));
        REQUIRE(x[1] == Catch::Approx(3.0).epsilon(1e-12));
        REQUIRE(x[2] == Catch::Approx(2.0).epsilon(1e-12));
    }

    TEST_CASE("GMRES applies its preconditioner in the Arnoldi iteration", "[SparseMatrix][SparseSolvers][GMRES]")
    {
        auto A = SparseMatrix::DiagonalCSR(std::vector<double>{1.0, 10.0, 100.0});
				Vector<double> b{1.0, 10.0, 100.0};
				Vector<double> x(3, 0.0);

				auto preconditioner = std::make_shared<SparseSolvers::JacobiPreconditioner<double>>();
        preconditioner->setup(A);

        SparseSolvers::SolverConfig<double> config;
        config.setTolerance(1e-12).setMaxIterations(1);

        auto result = SparseSolvers::solveGMRES(A, b, x, preconditioner, 1, config);

				REQUIRE(result.converged);
				REQUIRE(result.iterations_used == 1);
				REQUIRE(x[0] == Catch::Approx(1.0).epsilon(1e-12));
        REQUIRE(x[1] == Catch::Approx(1.0).epsilon(1e-12));
        REQUIRE(x[2] == Catch::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Stateful preconditioners require setup before apply", "[SparseMatrix][SparseSolvers][Preconditioner]")
    {
			const Vector<double> residual{1.0, 2.0};
			Vector<double> result;

			SECTION("Jacobi") {
				SparseSolvers::JacobiPreconditioner<double> preconditioner;
				REQUIRE_THROWS_AS(preconditioner.apply(residual, result), std::logic_error);
        }
        SECTION("SSOR")
        {
            SparseSolvers::SSORPreconditioner<double> preconditioner;
            REQUIRE_THROWS_AS(preconditioner.apply(residual, result), std::logic_error);
        }
        SECTION("ILU0")
        {
            SparseSolvers::ILU0Preconditioner<double> preconditioner;
            REQUIRE_THROWS_AS(preconditioner.apply(residual, result), std::logic_error);
        }
        SECTION("Block Jacobi")
        {
            SparseSolvers::BlockJacobiPreconditioner<double> preconditioner;
            REQUIRE_THROWS_AS(preconditioner.apply(residual, result), std::logic_error);
        }
    }

    TEST_CASE("Preconditioners reject missing or zero pivots", "[SparseMatrix][SparseSolvers][Preconditioner]")
    {
        SparseMatrix::SparseMatrixCSR<double> missingDiagonal(2, 2,
            {1.0}, {0}, {0, 1, 1});
        SparseMatrix::SparseMatrixCSR<double> zeroDiagonal(2, 2,
            {1.0, 0.0}, {0, 1}, {0, 1, 2});

        SECTION("SSOR")
        {
            SparseSolvers::SSORPreconditioner<double> preconditioner;
            REQUIRE_THROWS_AS(preconditioner.setup(missingDiagonal), SingularMatrixError);
            REQUIRE_THROWS_AS(preconditioner.setup(zeroDiagonal), SingularMatrixError);
        }
        SECTION("ILU0")
        {
            SparseSolvers::ILU0Preconditioner<double> preconditioner;
            REQUIRE_THROWS_AS(preconditioner.setup(missingDiagonal), SingularMatrixError);
            REQUIRE_THROWS_AS(preconditioner.setup(zeroDiagonal), SingularMatrixError);
        }
    }

    TEST_CASE("Conjugate gradient optionally rejects non-SPD matrices", "[SparseMatrix][SparseSolvers][CG]")
    {
        SparseMatrix::SparseMatrixCOO<double> entries(2, 2);
        entries.addEntry(0, 0, 1.0);
        entries.addEntry(0, 1, 2.0);
        entries.addEntry(1, 0, 2.0);
        entries.addEntry(1, 1, 1.0);
        SparseMatrix::SparseMatrixCSR<double> A(entries);
				Vector<double> b{1.0, 1.0};
				Vector<double> x(2, 0.0);

				SparseSolvers::SolverConfig<double> config;
        config.setValidateSPD(true);

        const auto result = SparseSolvers::solveCG(A, b, x, config);

				REQUIRE(result.status == AlgorithmStatus::InvalidInput);
				REQUIRE(result.error_message == "CG requires a symmetric positive-definite matrix");
		}

    TEST_CASE("AutoSolver routes only symmetric diagonally-dominant matrices to CG", "[SparseMatrix][SparseSolvers][AutoSolver][wewi.1.12]")
    {
        // Old heuristic: diag > 0 on the first ten rows => CG. This symmetric INDEFINITE matrix
        // (eigenvalues 3 and -1) passed it and was handed to CG (analysis v4, CORE C4).
        SparseMatrix::SparseMatrixCOO<double> indefinite(2, 2);
        indefinite.addEntry(0, 0, 1.0);
        indefinite.addEntry(0, 1, 2.0);
        indefinite.addEntry(1, 0, 2.0);
        indefinite.addEntry(1, 1, 1.0);
        SparseMatrix::SparseMatrixCSR<double> A(indefinite);
        REQUIRE_FALSE(SparseSolvers::AutoSolver<double>::isDiagonallyDominantWithPositiveDiagonal(A));
        auto solverA = SparseSolvers::AutoSolver<double>::create(A);
        REQUIRE(dynamic_cast<SparseSolvers::BiCGSTAB<double>*>(solverA.get()) != nullptr);

        // Symmetric, positive diagonal, but the LAST row breaks dominance (beyond the old 10-row window)
        const int n = 12;
        SparseMatrix::SparseMatrixCOO<double> lateBreak(n, n);
        for (int i = 0; i < n; ++i) {
            lateBreak.addEntry(i, i, 4.0);
            if (i + 1 < n) { lateBreak.addEntry(i, i + 1, -1.0); lateBreak.addEntry(i + 1, i, -1.0); }
        }
        lateBreak.addEntry(0, n - 1, 10.0);
        lateBreak.addEntry(n - 1, 0, 10.0);
        SparseMatrix::SparseMatrixCSR<double> B(lateBreak);
        REQUIRE(B.IsSymmetric());
        REQUIRE_FALSE(SparseSolvers::AutoSolver<double>::isDiagonallyDominantWithPositiveDiagonal(B));
        REQUIRE(dynamic_cast<SparseSolvers::BiCGSTAB<double>*>(SparseSolvers::AutoSolver<double>::create(B).get()) != nullptr);

        // Genuinely SPD (1D Laplacian + shift): CG
        SparseMatrix::SparseMatrixCOO<double> spd(n, n);
        for (int i = 0; i < n; ++i) {
            spd.addEntry(i, i, 2.5);
            if (i + 1 < n) { spd.addEntry(i, i + 1, -1.0); spd.addEntry(i + 1, i, -1.0); }
        }
        SparseMatrix::SparseMatrixCSR<double> C(spd);
        REQUIRE(SparseSolvers::AutoSolver<double>::isDiagonallyDominantWithPositiveDiagonal(C));
        REQUIRE(dynamic_cast<SparseSolvers::ConjugateGradient<double>*>(SparseSolvers::AutoSolver<double>::create(C).get()) != nullptr);

        // Weak dominance everywhere with no strict row (singular Laplacian with Neumann ends) is NOT accepted
        SparseMatrix::SparseMatrixCOO<double> singular(3, 3);
        singular.addEntry(0, 0, 1.0);  singular.addEntry(0, 1, -1.0);
        singular.addEntry(1, 0, -1.0); singular.addEntry(1, 1, 2.0); singular.addEntry(1, 2, -1.0);
        singular.addEntry(2, 1, -1.0); singular.addEntry(2, 2, 1.0);
        SparseMatrix::SparseMatrixCSR<double> D(singular);
        REQUIRE_FALSE(SparseSolvers::AutoSolver<double>::isDiagonallyDominantWithPositiveDiagonal(D));

        // and autoSolve still solves the indefinite system (via BiCGSTAB)
				Vector<double> b{1.0, 1.0}, x(2, 0.0);
				const auto result = SparseSolvers::autoSolve(A, b, x);
				REQUIRE(result.converged);
				REQUIRE(x[0] == Catch::Approx(1.0 / 3.0).margin(1e-8));
        REQUIRE(x[1] == Catch::Approx(1.0 / 3.0).margin(1e-8));
    }
}