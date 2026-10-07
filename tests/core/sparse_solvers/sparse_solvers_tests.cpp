#include <catch2/catch_all.hpp>

#include <mml/base/SparseMatrix/SparseMatrix.h>
#include <mml/core/SparseSolvers/IterativeSolvers.h>

#include <vector>

namespace MML::Tests::Core::SparseSolversTests
{
    TEST_CASE("Sparse matrix Laplacian solves with conjugate gradient", "[SparseMatrix][SparseSolvers]")
    {
        auto A = SparseMatrix::laplacian1D<double>(3);
        std::vector<double> b{1.0, 2.0, 1.0};
        std::vector<double> x(3, 0.0);

        SparseSolvers::SolverConfig<double> config;
        config.setTolerance(1e-12).setMaxIterations(50).setValidateSPD(true);

        auto result = SparseSolvers::solveCG(A, b, x, config);

        REQUIRE(result.converged());
        REQUIRE(x[0] == Catch::Approx(2.0).epsilon(1e-12));
        REQUIRE(x[1] == Catch::Approx(3.0).epsilon(1e-12));
        REQUIRE(x[2] == Catch::Approx(2.0).epsilon(1e-12));
    }

    TEST_CASE("GMRES applies its preconditioner in the Arnoldi iteration", "[SparseMatrix][SparseSolvers][GMRES]")
    {
        auto A = SparseMatrix::diagonalCSR(std::vector<double>{1.0, 10.0, 100.0});
        std::vector<double> b{1.0, 10.0, 100.0};
        std::vector<double> x(3, 0.0);

        auto preconditioner = std::make_shared<SparseSolvers::JacobiPreconditioner<double>>();
        preconditioner->setup(A);

        SparseSolvers::SolverConfig<double> config;
        config.setTolerance(1e-12).setMaxIterations(1);

        auto result = SparseSolvers::solveGMRES(A, b, x, preconditioner, 1, config);

        REQUIRE(result.converged());
        REQUIRE(result.iterations == 1);
        REQUIRE(x[0] == Catch::Approx(1.0).epsilon(1e-12));
        REQUIRE(x[1] == Catch::Approx(1.0).epsilon(1e-12));
        REQUIRE(x[2] == Catch::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Stateful preconditioners require setup before apply", "[SparseMatrix][SparseSolvers][Preconditioner]")
    {
        const std::vector<double> residual{1.0, 2.0};
        std::vector<double> result;

        SECTION("Jacobi")
        {
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
        std::vector<double> b{1.0, 1.0};
        std::vector<double> x(2, 0.0);

        SparseSolvers::SolverConfig<double> config;
        config.setValidateSPD(true);

        const auto result = SparseSolvers::solveCG(A, b, x, config);

        REQUIRE(result.status == SparseSolvers::SolverStatus::InvalidInput);
        REQUIRE(result.message == "CG requires a symmetric positive-definite matrix");
    }
}
