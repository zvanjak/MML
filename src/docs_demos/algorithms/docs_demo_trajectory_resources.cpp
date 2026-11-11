#include <mml/algorithms/ODESolvers/ODESolverStiff.h>
#include <iostream>
#include <stdexcept>

namespace {
	class ConstantTrajectoryODE : public MML::IODESystemWithJacobian {
	public:
		int getDim() const override { return 1; }
		void derivs(MML::Real, const MML::Vector<MML::Real>&, MML::Vector<MML::Real>& f) const override { f[0] = 0; }
		void jacobian(MML::Real, const MML::Vector<MML::Real>&, MML::Matrix<MML::Real>& j) const override { j(0, 0) = 0; }
	};
}

void Docs_Demo_TrajectoryResources()
{
	using namespace MML;
	Matrix<Real> narrow(1, 100001);
	narrow(0, 100000) = 7;
	ConstantTrajectoryODE system;
	StiffSolverConfig config;
	config.step_size = 1;
	config.max_steps = 100001;
	auto result = SolveBackwardEuler(system, 0, Vector<Real>{3}, 100001, config);
	if (result.status != AlgorithmStatus::Success ||
		result.solution.size() != 100002 ||
		result.total_steps + 1 != result.solution.size() ||
		result.solution.getXValuesAtEnd()[0] != Real(3))
		throw std::runtime_error("Element-budget trajectory contract failed");
	std::cout << to_string(result.status) << ": long/narrow matrix has " << narrow.cols()
	          << " columns and last value " << narrow(0, 100000)
	          << "\nRetained points: " << result.solution.size()
	          << ", accepted intervals: " << result.total_steps
	          << ", last saved time: " << result.solution.getTValue(result.solution.size() - 1) << '\n';
	bool rejected = false;
	try {
		result.solution.reserveForPoint(result.solution.max_size());
	}
	catch (const ResourceLimitError& error) {
		rejected = true;
		std::cout << "Over-budget reservation rejected: " << error.what() << '\n';
	}
	if (!rejected || result.solution.size() != 100002 ||
		result.solution.getXValuesAtEnd()[0] != Real(3))
		throw std::runtime_error("Failed reservation did not preserve the trajectory");
	bool matrixRejected = false;
	try {
		Matrix<double> huge(1000000, 1000000);
	}
	catch (const MatrixAllocationError& error) {
		matrixRejected = true;
		std::cout << "Allocation rejected: " << error.what() << '\n';
	}
	if (!matrixRejected)
		throw std::runtime_error("Matrix element budget was not enforced");
}
