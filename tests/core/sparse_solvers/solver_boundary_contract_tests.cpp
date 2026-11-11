#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/SparseSolvers/IterativeSolvers.h>
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/core/LinAlgEqSolvers/LinAlgSVD.h>
#endif

#include <limits>
#include <type_traits>

namespace MML::Tests::Core::SolverBoundaryContracts
{
	using namespace SparseSolvers;
	using SparseMatrix::SparseMatrixCSR;

	template<class T>
	concept HasKrylovSurface = requires {
		typename SolverConfig<T>;
		typename KrylovResult<T>;
		typename IterativeSolverBase<T>;
		typename ConjugateGradient<T>;
		typename BiCGSTAB<T>;
		typename GMRES<T>;
		typename AutoSolver<T>;
	};

	template<class Result>
	concept HasOwnedSolution = requires(Result result) { result.solution; };

	template<template<class> class Surface, class T>
	concept NamesScalarSurface = requires { typename Surface<T>; };

	template<template<class> class Surface>
	void CheckScalarConstraint()
	{
		STATIC_REQUIRE((NamesScalarSurface<Surface, float>));
		STATIC_REQUIRE((NamesScalarSurface<Surface, double>));
		STATIC_REQUIRE((NamesScalarSurface<Surface, long double>));
		STATIC_REQUIRE_FALSE((NamesScalarSurface<Surface, Complex>));
		STATIC_REQUIRE_FALSE((NamesScalarSurface<Surface, int>));
	}

	SparseMatrixCSR<Real> CoupledMatrix(Real scale = 1)
	{
		return {2, 2, {4 * scale, scale, scale, 3 * scale}, {0, 1, 0, 1}, {0, 2, 4}};
	}

	void CheckVector(const Vector<Real>& actual, const Vector<Real>& expected)
	{
		REQUIRE(actual.size() == expected.size());
		for (int i = 0; i < actual.size(); ++i)
			CHECK(std::abs(actual[i] - expected[i]) <= Precision::DefaultToleranceStrict * std::max(Real(1), std::abs(expected[i])));
	}

	template<class Solver>
	void CheckKrylovAlias(Solver& solver)
	{
		const auto matrix = CoupledMatrix();
		const Vector<Real> rhs{1, 2};
		Vector<Real> expected = rhs;
		const auto ordinary = solver.solve(matrix, rhs, expected);
		REQUIRE(ordinary.converged);
		Vector<Real> alias = rhs;
		const auto result = solver.solve(matrix, alias, alias);
		CHECK(result.converged);
		CheckVector(alias, expected);
		const auto residual = matrix * alias - rhs;
		CHECK(norm2(residual) < Precision::DefaultToleranceStrict);
		CHECK(std::abs(result.residualNorm - norm2(residual)) < Precision::DefaultToleranceStrict);
	}

	template<class Solver>
	void CheckOwnedSolve(Solver& solver)
	{
		const auto matrix = CoupledMatrix();
		const Vector<Real> rhs{1, 2};
		const auto result = solver.solve(matrix, rhs);
		REQUIRE(result.converged);
		CHECK(HasOwnedSolution<decltype(result)>);
		if constexpr (HasOwnedSolution<decltype(result)>) {
			REQUIRE(result.solution.size() == rhs.size());
			CHECK(norm2(matrix * result.solution - rhs) < Precision::DefaultToleranceStrict);
		}
		Vector<Real> external(rhs.size(), Real(0));
		const auto ordinary = solver.solve(matrix, rhs, external);
		REQUIRE(ordinary.converged);
		CHECK(ordinary.solution.empty());
		solver.config().maxIterations = 0;
		const auto partial = solver.solve(matrix, rhs);
		CHECK(partial.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(partial.solution.size() == rhs.size());
		CHECK(norm2(partial.solution) == 0);
		CHECK(std::abs(partial.residualNorm - norm2(rhs)) < Precision::DefaultToleranceStrict);
	}

	TEST_CASE("Preconditioner parameters reject invalid domains before setup", "[SparseSolvers][solver-boundary]")
	{
		for (int size : {0, -1, std::numeric_limits<int>::min()})
			CHECK_THROWS_AS(BlockJacobiPreconditioner<Real>(size), std::invalid_argument);
		BlockJacobiPreconditioner<Real> block;
		CHECK_THROWS_AS(block.setBlockSize(0), std::invalid_argument);
		for (Real omega : {Real(0), Real(-1), Real(2), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
			CHECK_THROWS_AS(SSORPreconditioner<Real>(omega), std::invalid_argument);
			SSORPreconditioner<Real> ssor;
			CHECK_THROWS_AS(ssor.setOmega(omega), std::invalid_argument);
		}
	}

	TEST_CASE("Block Jacobi supports meaningful in-place application", "[SparseSolvers][solver-boundary]")
	{
		BlockJacobiPreconditioner<Real> preconditioner(2);
		preconditioner.setup(CoupledMatrix());
		const Vector<Real> rhs{1, 3};
		Vector<Real> separate;
		preconditioner.apply(rhs, separate);
		CheckVector(separate, Vector<Real>{0, 1});
		Vector<Real> alias = rhs;
		preconditioner.apply(alias, alias);
		CheckVector(alias, separate);
	}

	TEST_CASE("SSOR setup installs an immutable snapshot", "[SparseSolvers][solver-boundary]")
	{
		auto matrix = CoupledMatrix();
		SSORPreconditioner<Real> preconditioner;
		preconditioner.setup(matrix);
		const Vector<Real> rhs{1, 2};
		Vector<Real> original, after;
		preconditioner.apply(rhs, original);
		matrix.values()[1] = 2;
		matrix.values()[2] = 2;
		preconditioner.apply(rhs, after);
		CheckVector(after, original);
		auto copy = preconditioner;
		copy.apply(rhs, after);
		CheckVector(after, original);
	}

	TEST_CASE("Krylov RHS aliases retain the original mathematical problem", "[SparseSolvers][solver-boundary]")
	{
		ConjugateGradient<Real> cg;
		BiCGSTAB<Real> bicg;
		GMRES<Real> gmres;
		CheckKrylovAlias(cg);
		CheckKrylovAlias(bicg);
		CheckKrylovAlias(gmres);
	}

	TEST_CASE("Zero-guess Krylov convenience retains its computed solution", "[SparseSolvers][solver-boundary]")
	{
		ConjugateGradient<Real> cg;
		BiCGSTAB<Real> bicg;
		GMRES<Real> gmres;
		CheckOwnedSolve(cg);
		CheckOwnedSolve(bicg);
		CheckOwnedSolve(gmres);
	}

	TEST_CASE("Owned Krylov failure retains the nonzero last iterate", "[SparseSolvers][solver-boundary]")
	{
		ConjugateGradient<Real> solver;
		solver.config().maxIterations = 1;
		const auto matrix = CoupledMatrix();
		const Vector<Real> rhs{1, 2};
		const auto result = solver.solve(matrix, rhs);
		REQUIRE(result.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(result.iterations_used == 1);
		CheckVector(result.solution, Vector<Real>{Real(1) / 4, Real(1) / 2});
		CHECK(std::abs(result.residualNorm - norm2(matrix * result.solution - rhs)) < Precision::DefaultToleranceStrict);
	}

	TEST_CASE("Gauss Jordan rejects incompatible simultaneous output aliases", "[solver-boundary][GaussJordan]")
	{
		Matrix<Real> matrix(2, 2, {4, 1, 1, 3});
		const Matrix<Real> original = matrix;
		CHECK_THROWS_AS(GaussJordanSolver<Real>::SolveInPlace(matrix, matrix), std::invalid_argument);
		CHECK(matrix == original);
	}

	TEST_CASE("AutoSolver does not certify a reducible singular PSD matrix", "[SparseSolvers][solver-boundary]")
	{
		const SparseMatrixCSR<Real> singular(3, 3, {1, -1, -1, 1, 1}, {0, 1, 0, 1, 2}, {0, 2, 4, 5});
		const Vector<Real> nullVector{1, 1, 0};
		CHECK(norm2(singular * nullVector) == 0);
		CHECK_FALSE(AutoSolver<Real>::isDiagonallyDominantWithPositiveDiagonal(singular));
		const auto solver = AutoSolver<Real>::create(singular);
		CHECK(dynamic_cast<BiCGSTAB<Real>*>(solver.get()) != nullptr);
	}

	TEST_CASE("Krylov supported scalar domain is visible without member instantiation", "[SparseSolvers][solver-boundary]")
	{
		CHECK(HasKrylovSurface<float>);
		CHECK(HasKrylovSurface<double>);
		CHECK(HasKrylovSurface<long double>);
		CHECK_FALSE(HasKrylovSurface<Complex>);
		CHECK_FALSE(HasKrylovSurface<int>);
		CheckScalarConstraint<SolverConfig>();
		CheckScalarConstraint<KrylovResult>();
		CheckScalarConstraint<IterativeSolverBase>();
		CheckScalarConstraint<ConjugateGradient>();
		CheckScalarConstraint<BiCGSTAB>();
		CheckScalarConstraint<GMRES>();
		CheckScalarConstraint<AutoSolver>();
		CheckScalarConstraint<SSORPreconditioner>();
		CheckScalarConstraint<BlockJacobiPreconditioner>();
	}

	TEST_CASE("Block partition arithmetic is bounded and retains remainder blocks", "[SparseSolvers][solver-boundary]")
	{
		for (int size : {1, 2, 3, std::numeric_limits<int>::max()}) {
			CAPTURE(size);
			BlockJacobiPreconditioner<Real> preconditioner(size);
			const auto matrix = SparseMatrix::DiagonalCSR(std::vector<Real>{2, 3, 4});
			REQUIRE_NOTHROW(preconditioner.setup(matrix));
			Vector<Real> alias{2, 6, 12};
			preconditioner.apply(alias, alias);
			CheckVector(alias, Vector<Real>{1, 2, 3});
		}
		BlockJacobiPreconditioner<Real> empty;
		CHECK_THROWS_AS(empty.setup(SparseMatrixCSR<Real>()), MatrixDimensionError);
		CHECK_THROWS_AS(empty.setup(SparseMatrixCSR<Real>(2, 3)), MatrixDimensionError);
	}

	TEST_CASE("Preconditioner snapshots accept temporaries and survive source destruction", "[SparseSolvers][solver-boundary]")
	{
		SSORPreconditioner<Real> temporary;
		temporary.setup(CoupledMatrix());
		SSORPreconditioner<Real> destroyed;
		{
			const auto local = CoupledMatrix();
			destroyed.setup(local);
		}
		auto copied = destroyed;
		auto moved = std::move(destroyed);
		const Vector<Real> rhs{1, 2};
		Vector<Real> expected, actual;
		temporary.apply(rhs, expected);
		copied.apply(rhs, actual);
		CheckVector(actual, expected);
		moved.apply(rhs, actual);
		CheckVector(actual, expected);
		Vector<Real> alias = rhs;
		moved.apply(alias, alias);
		CheckVector(alias, expected);
	}

	TEST_CASE("Existing preconditioner in-place actions remain consistent", "[SparseSolvers][solver-boundary]")
	{
		const auto matrix = CoupledMatrix();
		const Vector<Real> rhs{1, 2};
		auto check = [&](auto preconditioner) {
			preconditioner.setup(matrix);
			Vector<Real> expected, alias = rhs;
			preconditioner.apply(rhs, expected);
			preconditioner.apply(alias, alias);
			CheckVector(alias, expected);
		};
		check(IdentityPreconditioner<Real>());
		check(JacobiPreconditioner<Real>());
		check(ILU0Preconditioner<Real>());
		check(SSORPreconditioner<Real>());
	}

	TEST_CASE("Preconditioner setup and setters preserve installed state on rejection", "[SparseSolvers][solver-boundary]")
	{
		SSORPreconditioner<Real> ssor;
		BlockJacobiPreconditioner<Real> block(2);
		const auto matrix = CoupledMatrix();
		ssor.setup(matrix);
		block.setup(matrix);
		const Vector<Real> rhs{1, 2};
		Vector<Real> ssorBefore, blockBefore, actual;
		ssor.apply(rhs, ssorBefore);
		block.apply(rhs, blockBefore);
		CHECK_THROWS_AS(ssor.setOmega(std::numeric_limits<Real>::quiet_NaN()), ArgumentError);
		CHECK_THROWS_AS(block.setBlockSize(0), ArgumentError);
		const SparseMatrixCSR<Real> singular(2, 2, {1, 1, 1, 1}, {0, 1, 0, 1}, {0, 2, 4});
		CHECK_THROWS_AS(block.setup(singular), SingularMatrixError);
		const auto zeroDiagonal = SparseMatrix::DiagonalCSR(std::vector<Real>{1, 0});
		CHECK_THROWS_AS(ssor.setup(zeroDiagonal), SingularMatrixError);
		CHECK_THROWS_AS(ssor.setup(SparseMatrixCSR<Real>()), MatrixDimensionError);
		CHECK_THROWS_AS(ssor.setup(SparseMatrixCSR<Real>(2, 3)), MatrixDimensionError);
		auto nonfinite = CoupledMatrix();
		nonfinite.values()[1] = std::numeric_limits<Real>::infinity();
		CHECK_THROWS_AS(ssor.setup(nonfinite), ArgumentError);
		CHECK_THROWS_AS(block.setup(nonfinite), ArgumentError);
		ssor.apply(rhs, actual);
		CheckVector(actual, ssorBefore);
		block.apply(rhs, actual);
		CheckVector(actual, blockBefore);
		block.setBlockSize(1);
		block.apply(rhs, actual);
		CheckVector(actual, blockBefore);
		block.setup(matrix);
		block.apply(rhs, actual);
		CheckVector(actual, Vector<Real>{Real(1) / 4, Real(2) / 3});
		CHECK(ssor.omega() == 1);
		ssor.setOmega(Real(3) / 2);
		ssor.apply(rhs, actual);
		CHECK(norm2(actual - ssorBefore) > Precision::DefaultToleranceStrict);
	}

	TEST_CASE("Preconditioner inverse actions preserve scale and fail explicitly outside range", "[SparseSolvers][solver-boundary]")
	{
		for (Real scale : {Testing::Tol(Real(1e-100), Real(1e-15)), Real(1), Testing::Tol(Real(1e100), Real(1e15))}) {
			CAPTURE(scale);
			const auto matrix = CoupledMatrix(scale);
			const Vector<Real> exact{2, -1};
			const Vector<Real> rhs = matrix * exact;
			BlockJacobiPreconditioner<Real> block(2);
			block.setup(matrix);
			Vector<Real> result;
			block.apply(rhs, result);
			CheckVector(result, exact);
			CHECK(norm2(matrix * result - rhs) / norm2(rhs) < Precision::DefaultToleranceStrict);
			SSORPreconditioner<Real> ssor;
			ssor.setup(SparseMatrix::DiagonalCSR(std::vector<Real>{2 * scale, 4 * scale}));
			const Vector<Real> diagonalRhs{2 * scale, -4 * scale};
			ssor.apply(diagonalRhs, result);
			CheckVector(result, Vector<Real>{1, -1});
		}
		BlockJacobiPreconditioner<Real> block(2);
		block.setup(CoupledMatrix());
		Vector<Real> result{7, 8};
		const Vector<Real> wrongSize{1};
		CHECK_THROWS_AS(block.apply(wrongSize, result), std::invalid_argument);
		CheckVector(result, Vector<Real>{7, 8});
		const Vector<Real> invalid{std::numeric_limits<Real>::quiet_NaN(), 1};
		CHECK_THROWS_AS(block.apply(invalid, result), ArgumentError);
		SSORPreconditioner<Real> ssor;
		ssor.setup(CoupledMatrix());
		CHECK_THROWS_AS(ssor.apply(invalid, result), ArgumentError);
	}

	TEST_CASE("Tall SVD consumes the complete RHS before resizing an aliased output", "[solver-boundary][SVD]")
	{
		for (Real scale : {Testing::Tol(Real(1e-100), Real(1e-15)), Real(1), Testing::Tol(Real(1e100), Real(1e15))}) {
			CAPTURE(scale);
			const Matrix<Real> matrix(3, 2, {scale, 0, 0, scale, scale, scale});
			const Vector<Real> rhs{2 * scale, 3 * scale, 5 * scale};
			SVDecompositionSolver<Real> svd(matrix);
			Vector<Real> ordinary;
			svd.Solve(rhs, ordinary);
			Vector<Real> alias = rhs;
			svd.Solve(alias, alias);
			CheckVector(alias, Vector<Real>{2, 3});
			CheckVector(alias, ordinary);
			CHECK(norm2(matrix * alias - rhs) / norm2(rhs) < Precision::DefaultToleranceStrict);
			Vector<Real> wrongSize{1};
			CHECK_THROWS_AS(svd.Solve(wrongSize, wrongSize), VectorDimensionError);
			CheckVector(wrongSize, Vector<Real>{1});
		}
	}

	TEST_CASE("Gauss Jordan distinct outputs retain inverse and multi-RHS solutions", "[solver-boundary][GaussJordan]")
	{
		const Matrix<Real> original(2, 2, {4, 1, 1, 3});
		const Matrix<Real> rhs(2, 2, {1, 2, 3, 4});
		auto matrix = original, solved = rhs;
		GaussJordanSolver<Real>::SolveInPlace(matrix, solved);
		const auto reconstructed = original * solved;
		const auto identity = original * matrix;
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j) {
				CHECK(std::abs(reconstructed(i, j) - rhs(i, j)) < Precision::DefaultToleranceStrict);
				CHECK(std::abs(identity(i, j) - Real(i == j)) < Precision::DefaultToleranceStrict);
			}
	}

	TEST_CASE("Square SVD matrix RHS aliases preserve independent columns", "[solver-boundary][SVD]")
	{
		const Matrix<Real> matrix(2, 2, {4, 1, 1, 3});
		// Unit-scale RHS values isolate aliasing from platform-dependent decomposition roundoff.
		const Matrix<Real> rhs(2, 2, {Real(0.25), Real(0.5), Real(0.75), Real(1)});
		SVDecompositionSolver<Real> solver(matrix);
		auto ordinary = rhs;
		solver.Solve(rhs, ordinary);
		auto alias = rhs;
		solver.Solve(alias, alias);
		const auto reconstructed = matrix * alias;
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j) {
				CHECK(alias(i, j) == ordinary(i, j));
				CHECK(std::abs(reconstructed(i, j) - rhs(i, j)) < Precision::DefaultToleranceStrict);
			}
	}

	TEST_CASE("AutoSolver requires strict dominance per connected component at every scale", "[SparseSolvers][solver-boundary]")
	{
		for (Real scale : {Testing::Tol(Real(1e-100), Real(1e-15)), Real(1), Testing::Tol(Real(1e100), Real(1e15))}) {
			const SparseMatrixCSR<Real> anchored(4, 4, {scale, -scale, -scale, 2 * scale, scale, -scale, -scale, 2 * scale},
				{0, 1, 0, 1, 2, 3, 2, 3}, {0, 2, 4, 6, 8});
			REQUIRE(AutoSolver<Real>::isDiagonallyDominantWithPositiveDiagonal(anchored));
			CHECK(dynamic_cast<ConjugateGradient<Real>*>(AutoSolver<Real>::create(anchored).get()) != nullptr);
			const SparseMatrixCSR<Real> unanchored(4, 4, {scale, -scale, -scale, scale, scale, -scale, -scale, 2 * scale},
				{0, 1, 0, 1, 2, 3, 2, 3}, {0, 2, 4, 6, 8});
			CHECK_FALSE(AutoSolver<Real>::isDiagonallyDominantWithPositiveDiagonal(unanchored));
		}
		const auto laplacian = SparseMatrix::Laplacian1D<Real>(5);
		CHECK(AutoSolver<Real>::isDiagonallyDominantWithPositiveDiagonal(laplacian));
		CHECK_FALSE(AutoSolver<Real>::isDiagonallyDominantWithPositiveDiagonal(SparseMatrixCSR<Real>()));
		CHECK_FALSE(AutoSolver<Real>::isDiagonallyDominantWithPositiveDiagonal(SparseMatrixCSR<Real>(2, 3)));
		auto nonfinite = CoupledMatrix();
		nonfinite.values()[0] = std::numeric_limits<Real>::infinity();
		CHECK_FALSE(AutoSolver<Real>::isDiagonallyDominantWithPositiveDiagonal(nonfinite));
	}
}
