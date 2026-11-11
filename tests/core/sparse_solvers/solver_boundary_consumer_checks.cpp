#include <mml/core/SparseSolvers/IterativeSolvers.h>
#include <mml/core/LinAlgEqSolvers/LinAlgSVD.h>

#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>

namespace
{
	thread_local bool trackAllocations = false;
	thread_local std::size_t allocationCount = 0;

	struct AllocationScope {
		AllocationScope() { allocationCount = 0; trackAllocations = true; }
		~AllocationScope() { trackAllocations = false; }
	};

	template<class Operation>
	std::size_t Allocations(Operation&& operation)
	{
		AllocationScope scope;
		operation();
		return allocationCount;
	}

	void Require(bool condition, const char* message)
	{
		if (!condition) throw std::runtime_error(message);
	}
}

// Isolated executable: allocation instrumentation must not replace new/delete in the main test suite.
void* operator new(std::size_t size)
{
	for (;;) {
		if (void* storage = std::malloc(size == 0 ? 1 : size)) {
			if (trackAllocations) ++allocationCount;
			return storage;
		}
		const auto handler = std::get_new_handler();
		if (!handler) throw std::bad_alloc();
		handler();
	}
}

void operator delete(void* storage) noexcept { std::free(storage); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* storage) noexcept { ::operator delete(storage); }
void operator delete(void* storage, std::size_t) noexcept { ::operator delete(storage); }
void operator delete[](void* storage, std::size_t) noexcept { ::operator delete(storage); }

int main()
{
	try {
		using namespace MML;
		using namespace MML::SparseSolvers;
		const SparseMatrix::SparseMatrixCSR<Real> matrix(2, 2, {4, 1, 1, 3}, {0, 1, 0, 1}, {0, 2, 4});
		const Vector<Real> rhs{1, 3};

		BlockJacobiPreconditioner<Real> block(2);
		block.setup(matrix);
		Vector<Real> output(2, Real(0)), alias = rhs;
		const auto ordinaryBlock = Allocations([&] { block.apply(rhs, output); });
		const auto aliasedBlock = Allocations([&] { block.apply(alias, alias); });
		Require(ordinaryBlock == 0, "Nonaliased block application allocated");
		Require(aliasedBlock == 1, "Aliased block application must allocate exactly one RHS snapshot");
		Require(norm2(alias - output) < Precision::DefaultToleranceStrict, "Block allocation probe lost correctness");

		SSORPreconditioner<Real> ssor;
		ssor.setup(matrix);
		const auto ssorCalls = Allocations([&] { ssor.apply(rhs, output); });
		Require(ssorCalls == 1, "SSOR apply must allocate only its existing triangular workspace, not a matrix copy");

		const Matrix<Real> tall(3, 2, {1, 0, 0, 1, 1, 1});
		SVDecompositionSolver<Real> svd(tall);
		const Vector<Real> tallRhs{2, 3, 5};
		Vector<Real> tallAlias = tallRhs;
		const auto ordinarySvd = Allocations([&] { svd.Solve(tallRhs, output); });
		const auto aliasedSvd = Allocations([&] { svd.Solve(tallAlias, tallAlias); });
		Require(ordinarySvd == 1 && aliasedSvd == ordinarySvd, "SVD alias must use its existing projection workspace only");
		Require(norm2(tallAlias - Vector<Real>{2, 3}) < Precision::DefaultToleranceStrict, "SVD allocation probe lost correctness");

		ConjugateGradient<Real> solver;
		Vector<Real> external(2, Real(0));
		KrylovResult<Real> ordinaryResult, ownedResult;
		const auto ordinarySolve = Allocations([&] { ordinaryResult = solver.solve(matrix, rhs, external); });
		const auto ownedSolve = Allocations([&] { ownedResult = solver.solve(matrix, rhs); });
		Require(ordinaryResult.converged && ownedResult.converged, "Krylov allocation probe did not converge");
		Require(ownedSolve == ordinarySolve + 1, "Zero-guess result copied instead of moving its one owned initial vector");
		Require(ordinaryResult.solution.empty(), "Explicit-output solve unexpectedly owns a solution copy");
		Require(norm2(matrix * ownedResult.solution - rhs) < Precision::DefaultToleranceStrict, "Owned solution failed its true residual check");

		external = rhs;
		alias = rhs;
		const auto separateRhs = Allocations([&] { ordinaryResult = solver.solve(matrix, rhs, external); });
		const auto aliasedRhs = Allocations([&] { ownedResult = solver.solve(matrix, alias, alias); });
		Require(ordinaryResult.converged && ownedResult.converged, "Aliased Krylov allocation probe did not converge");
		Require(aliasedRhs == separateRhs + 1, "Krylov RHS alias must add exactly one immutable snapshot");
		Require(norm2(alias - external) < Precision::DefaultToleranceStrict, "Krylov alias changed the solved problem");

		std::printf("Allocation checks passed: block %zu/%zu; SSOR %zu; SVD %zu/%zu; Krylov explicit/owned %zu/%zu; RHS separate/alias %zu/%zu.\n",
			ordinaryBlock, aliasedBlock, ssorCalls, ordinarySvd, aliasedSvd, ordinarySolve, ownedSolve, separateRhs, aliasedRhs);
		return 0;
	}
	catch (const std::exception& error) {
		std::fprintf(stderr, "Solver boundary consumer check failed: %s\n", error.what());
		return 1;
	}
}
