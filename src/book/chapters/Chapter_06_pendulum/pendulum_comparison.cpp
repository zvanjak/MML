#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Vector.h"

#include "mml/tools/ConsolePrinter.h"

#include "mml/algorithms/ODESystemSolver.h"
#include "mml/algorithms/ODESystemStepCalculators.h"
#include "mml/algorithms/ODESystemSteppers.h"

#include "mml/algorithms/FunctionsAnalyzer.h"

#include "mpl/Mechanics/Pendulums.h"
#endif

using namespace MML;
using namespace MPL;

void ComparePendulumPeriods(
	PendulumODE& sys,
	Vector<Real>& initCond,
	Real t1, Real t2, int expectNumSteps, Real minSaveInterval)
{
	ODESystemFixedStepSolver	fixedSolver(sys, StepCalculators::RK4_Basic);
	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);

	Vector<Real> initAnglesDeg{ 1, 2, 5, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 };
	Vector<Real> initAngles(initAnglesDeg.size());
	for (int i = 0; i < initAngles.size(); i++)
		initAngles[i] = Utils::DegToRad(initAnglesDeg[i]);

	TablePrinter<double, double> data("Angle", 5, 1,
		{ "Linear ", "Exact  ", "Diff %", "Fix.step.sim", "Diff %", "Adapt.step.sim", "Diff % " },
		{ {10,6,'F'}, {10,6,'F'}, {7,3,'F'}, {13,7,'F'}, {8,5,'F'}, {13,7,'F'}, {8,5,'F'} });

	for (int i = 0; i < initAngles.size(); i++)
	{
		initCond[0] = initAngles[i];

		// solving with fixed and adaptive solver
		ODESystemSolution solF = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);
		ODESystemSolution solA = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-06, 0.01);

		// calculating period from solutions
		PolynomInterpRealFunc solFInterp = solF.getSolAsPolyInterp(0, 3);
		PolynomInterpRealFunc solAInterp = solA.getSolAsPolyInterp(0, 3);
		RealFunctionAnalyzer fa3(solFInterp);
		RealFunctionAnalyzer fa4(solAInterp);
		Real periodSimulFixed = fa3.calcRootsPeriod(t1, t2, 1000);
		Real periodSimulAdapt = fa4.calcRootsPeriod(t1, t2, 1000);

		// analytical formulas
		Real periodLin = sys.calcPeriodLinearized();
		Real periodExact = sys.calcExactPeriod(initAngles[i]);
		Real diffPercent = Abs(periodLin - periodExact) / periodExact * 100;
		Real diffPercentFixed = Abs(2 * periodSimulFixed - periodExact) / periodExact * 100;
		Real diffPercentAdapt = Abs(2 * periodSimulAdapt - periodExact) / periodExact * 100;

		data.addRow(initAnglesDeg[i], { double(periodLin), double(periodExact), double(diffPercent),
				double(2 * periodSimulFixed), double(diffPercentFixed), double(2 * periodSimulAdapt), double(diffPercentAdapt) });
	}

	data.Print();
}

void CompareAdaptiveStepCounts(
	PendulumODE& sys,
	Vector<Real> initCond,
	Real t1, Real t2, Real minSaveInterval)
{
	Vector<Real> acc{ 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10 };
	Vector<Real> numSteps(acc.size());

	std::cout << "\nAccuracy   Num steps   OK steps   Bad steps" << std::endl;
	for (int i = 0; i < acc.size(); i++)
	{
		ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);
		ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, acc[i], 0.01);
		numSteps[i] = solAdapt.getTotalNumSteps();

		std::cout << std::setw(7) << acc[i] << "        " << std::setw(3) << numSteps[i] << "        "
			<< std::setw(3) << solAdapt.getNumStepsOK() << "        "
			<< std::setw(3) << solAdapt.getNumStepsBad() << std::endl;
	}
}
