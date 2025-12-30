#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "tools/Visualizer.h"
#include "tools/Serializer.h"

#include "interfaces/IODESystem.h"
#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"
#include "algorithms/ODEAdaptiveIntegrator.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                       IODESystem INTERFACE DEMOS                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Creating ODE systems using the IODESystem interface
void Docs_Demo_IODESystem_Interface()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: IODESystem Interface\n";
	std::cout << "==========================================================================\n";
	
	// The IODESystem interface defines:
	//   - getDim(): returns system dimension
	//   - derivs(t, x, dxdt): computes derivatives dx/dt = f(t, x)
	//   - getVarName(ind): optional, returns variable names
	
	std::cout << "\nIODESystem interface methods:\n";
	std::cout << "  int getDim() const                              - system dimension\n";
	std::cout << "  void derivs(t, x, dxdt) const                   - compute derivatives\n";
	std::cout << "  std::string getVarName(int ind) const           - variable names (optional)\n";
	
	// Example 1: Using a predefined ODE system
	std::cout << "\n--- Example 1: Predefined Simple Harmonic Oscillator ---\n";
	TestBeds::SimpleHarmonicOscillatorODE sho(2.0);  // omega = 2.0
	
	std::cout << "System dimension: " << sho.getDim() << std::endl;
	std::cout << "Variable names: ";
	for (int i = 0; i < sho.getDim(); i++)
		std::cout << sho.getVarName(i) << " ";
	std::cout << std::endl;
	
	// Evaluate derivatives at a point
	Vector<Real> x({1.0, 0.5});    // x=1, v=0.5
	Vector<Real> dxdt(2);
	sho.derivs(0.0, x, dxdt);
	std::cout << "At (x,v) = (1.0, 0.5), derivatives (dx/dt, dv/dt) = ";
	dxdt.Print(std::cout, 6, 3);
	std::cout << std::endl;
}

// Demo: Creating custom ODE systems with ODESystem class
void Docs_Demo_ODESystem_Creation()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Creating ODE Systems\n";
	std::cout << "==========================================================================\n";
	
	// Method 1: ODESystem with function pointer
	std::cout << "\n--- Method 1: ODESystem with function pointer ---\n";
	
	// Define a simple decay ODE: dy/dt = -k*y, k=1
	auto decayFunc = [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
		dydt[0] = -1.0 * y[0];
	};
	
	ODESystem decayODE(1, decayFunc);
	std::cout << "Created decay ODE: dy/dt = -y\n";
	std::cout << "Dimension: " << decayODE.getDim() << std::endl;
	
	Vector<Real> y({10.0}), dydt(1);
	decayODE.derivs(0.0, y, dydt);
	std::cout << "At y=10: dy/dt = " << dydt[0] << std::endl;
	
	// Method 2: Use predefined systems from TestBeds
	std::cout << "\n--- Method 2: Predefined systems from TestBeds ---\n";
	
	// Lorenz system (chaotic attractor)
	TestBeds::LorenzSystemODE lorenz(10.0, 28.0, 8.0/3.0);  // sigma, rho, beta
	std::cout << "Lorenz system (3D chaotic): dim = " << lorenz.getDim() << std::endl;
	
	// Van der Pol oscillator (nonlinear)
	TestBeds::VanDerPolODE vanderpol(0.3);  // mu parameter
	std::cout << "Van der Pol oscillator: dim = " << vanderpol.getDim() << std::endl;
	
	// Lotka-Volterra (predator-prey dynamics)
	TestBeds::LotkaVolterraODE lotkavolterra;
	std::cout << "Lotka-Volterra predator-prey: dim = " << lotkavolterra.getDim() << std::endl;
	
	// Simple pendulum
	TestBeds::SimplePendulumODE pendulum(9.81);
	std::cout << "Simple pendulum: dim = " << pendulum.getDim() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       STEP CALCULATORS (FIXED-STEP)                                ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Available step calculators for fixed-step integration
void Docs_Demo_Step_Calculators()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Fixed-Step Calculators\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nAvailable step calculators:\n";
	std::cout << "  EulerStep_Calculator           - 1st order (simple, low accuracy)\n";
	std::cout << "  EulerCromer_StepCalculator     - 1st order (better for oscillatory)\n";
	std::cout << "  VelocityVerlet_StepCalculator  - 2nd order (symplectic, Hamiltonian)\n";
	std::cout << "  Leapfrog_StepCalculator        - 2nd order (symplectic, energy-conserving)\n";
	std::cout << "  RungeKutta4_StepCalculator     - 4th order (general purpose)\n";
	std::cout << "  RK5_CashKarp_StepCalculator    - 5th order (with error estimate)\n";
	
	// Compare methods on simple harmonic oscillator
	std::cout << "\n--- Comparing methods on Simple Harmonic Oscillator ---\n";
	TestBeds::SimpleHarmonicOscillatorODE sho(1.0);
	Vector<Real> y0({1.0, 0.0});  // x=1, v=0
	Real t_end = 2.0 * Constants::PI;  // One complete period
	int steps = 100;
	
	// Euler method
	EulerStep_Calculator euler;
	ODESystemFixedStepSolver eulerSolver(sho, euler);
	auto solEuler = eulerSolver.integrate(y0, 0.0, t_end, steps);
	
	// RK4 method
	RungeKutta4_StepCalculator rk4;
	ODESystemFixedStepSolver rk4Solver(sho, rk4);
	auto solRK4 = rk4Solver.integrate(y0, 0.0, t_end, steps);
	
	// Leapfrog method (symplectic)
	Leapfrog_StepCalculator leapfrog;
	ODESystemFixedStepSolver leapfrogSolver(sho, leapfrog);
	auto solLeapfrog = leapfrogSolver.integrate(y0, 0.0, t_end, steps);
	
	// Exact solution at t = 2*pi: x = 1, v = 0 (returns to start)
	std::cout << "After one period (t=2π), expected: x=1.0, v=0.0\n";
	std::cout << "  Euler:    x=" << solEuler.getXValue(steps, 0) 
	          << ", v=" << solEuler.getXValue(steps, 1) << std::endl;
	std::cout << "  RK4:      x=" << solRK4.getXValue(steps, 0)
	          << ", v=" << solRK4.getXValue(steps, 1) << std::endl;
	std::cout << "  Leapfrog: x=" << solLeapfrog.getXValue(steps, 0)
	          << ", v=" << solLeapfrog.getXValue(steps, 1) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       ADAPTIVE INTEGRATORS                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Modern adaptive integrators
void Docs_Demo_Adaptive_Integrators()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Adaptive Step-Size Integrators\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nAvailable adaptive integrators:\n";
	std::cout << "  ODEAdaptiveIntegrator<CashKarp_Stepper>     - 5th order Cash-Karp\n";
	std::cout << "  ODEAdaptiveIntegrator<DormandPrince5_Stepper>  - 5th order DP5 (FSAL)\n";
	std::cout << "  ODEAdaptiveIntegrator<DormandPrince8_Stepper>  - 8th order DP8 (high accuracy)\n";
	std::cout << "\nConvenience aliases:\n";
	std::cout << "  CashKarpIntegrator, DormandPrince5Integrator, DormandPrince8Integrator\n";
	
	// Compare on Lorenz system
	std::cout << "\n--- Adaptive integration of Lorenz system ---\n";
	TestBeds::LorenzSystemODE lorenz(10.0, 28.0, 8.0/3.0);
	Vector<Real> y0({1.0, 1.0, 1.0});
	
	Real eps = 1e-8;  // Error tolerance
	Real saveInterval = 0.1;
	
	// Cash-Karp integrator
	ODEAdaptiveIntegrator<CashKarp_Stepper> ckIntegrator(lorenz);
	auto solCK = ckIntegrator.integrate(y0, 0.0, 10.0, saveInterval, eps);
	
	// Dormand-Prince 5th order (uses FSAL optimization)
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> dp5Integrator(lorenz);
	auto solDP5 = dp5Integrator.integrate(y0, 0.0, 10.0, saveInterval, eps);
	
	std::cout << "Cash-Karp:       " << solCK.getNumStepsOK() << " OK steps, "
	          << solCK.getNumStepsBad() << " rejected" << std::endl;
	std::cout << "Dormand-Prince5: " << solDP5.getNumStepsOK() << " OK steps, "
	          << solDP5.getNumStepsBad() << " rejected" << std::endl;
	
	std::cout << "\nFinal state (t=10):\n";
	std::cout << "  CK:  "; solCK.getXValuesAtEnd().Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  DP5: "; solDP5.getXValuesAtEnd().Print(std::cout, 8, 4); std::cout << std::endl;
}

// Demo: Legacy solver interface (backward compatibility)
void Docs_Demo_Legacy_Solver_Interface()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Legacy ODESystemSolver Interface\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nLegacy interface (still supported):\n";
	std::cout << "  ODESystemSolver<RK5_CashKarp_Stepper> solver(system);\n";
	std::cout << "  solver.integrate(y0, t1, t2, saveInterval, eps, h1, hmin);\n";
	std::cout << "\nPreferred modern interface:\n";
	std::cout << "  CashKarpIntegrator integrator(system);\n";
	std::cout << "  integrator.integrate(y0, t1, t2, saveInterval, eps);\n";
	
	// Example using legacy interface
	TestBeds::SimpleHarmonicOscillatorODE sho(1.0);
	Vector<Real> y0({1.0, 0.0});
	
	ODESystemSolver<RK5_CashKarp_Stepper> legacySolver(sho);
	auto sol = legacySolver.integrate(y0, 0.0, 10.0, 0.1, 1e-6, 0.01, 0.0);
	
	std::cout << "\nLegacy solver result:\n";
	std::cout << "  Saved points: " << sol.size() << std::endl;
	std::cout << "  Final state: "; sol.getXValuesAtEnd().Print(std::cout, 8, 4);
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       ODESystemSolution POST-PROCESSING                            ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Working with ODESystemSolution
void Docs_Demo_ODESystemSolution_Basics()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: ODESystemSolution Basics\n";
	std::cout << "==========================================================================\n";
	
	// Integrate Van der Pol oscillator
	TestBeds::VanDerPolODE vanderpol(1.0);
	Vector<Real> y0({2.0, 0.0});
	
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(vanderpol);
	ODESystemSolution sol = integrator.integrate(y0, 0.0, 20.0, 0.1, 1e-6);
	
	std::cout << "\nODESystemSolution properties:\n";
	std::cout << "  System dimension: " << sol.getSysDim() << std::endl;
	std::cout << "  Time range: [" << sol.getT1() << ", " << sol.getT2() << "]" << std::endl;
	std::cout << "  Saved points: " << sol.size() << std::endl;
	std::cout << "  Steps OK/rejected: " << sol.getNumStepsOK() << "/" << sol.getNumStepsBad() << std::endl;
	
	std::cout << "\nAccessing data:\n";
	std::cout << "  getTValue(i)         - time at point i\n";
	std::cout << "  getXValue(i, comp)   - component value at point i\n";
	std::cout << "  getXValuesAtEnd()    - final state vector\n";
	std::cout << "  getTValues()         - all time points as Vector\n";
	std::cout << "  getXValues()         - all values as Matrix (dim x points)\n";
	
	std::cout << "\nFirst 5 points:\n";
	for (int i = 0; i < 5; i++) {
		std::cout << "  t=" << sol.getTValue(i) << ": x=" << sol.getXValue(i, 0)
		          << ", v=" << sol.getXValue(i, 1) << std::endl;
	}
}

// Demo: Solution interpolation
void Docs_Demo_ODESystemSolution_Interpolation()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Solution Interpolation\n";
	std::cout << "==========================================================================\n";
	
	// Integrate with coarse output
	TestBeds::SimpleHarmonicOscillatorODE sho(1.0);
	Vector<Real> y0({1.0, 0.0});
	
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(sho);
	ODESystemSolution sol = integrator.integrate(y0, 0.0, 2*Constants::PI, 0.5, 1e-8);
	
	std::cout << "\nInterpolation methods:\n";
	std::cout << "  getSolAsLinInterp(comp)         - linear interpolation\n";
	std::cout << "  getSolAsPolyInterp(comp, order) - polynomial interpolation\n";
	std::cout << "  getSolAsSplineInterp(comp)      - cubic spline interpolation\n";
	
	// Create interpolators
	auto linInterp = sol.getSolAsLinInterp(0);       // Linear for x
	auto splineInterp = sol.getSolAsSplineInterp(0); // Spline for x
	
	std::cout << "\nCompare interpolation at t=π (exact x=cos(π)=-1):\n";
	Real t_test = Constants::PI;
	std::cout << "  Linear:  x(" << t_test << ") = " << linInterp(t_test) << std::endl;
	std::cout << "  Spline:  x(" << t_test << ") = " << splineInterp(t_test) << std::endl;
	std::cout << "  Exact:   x(" << t_test << ") = " << std::cos(t_test) << std::endl;
}

// Demo: Solution as parametric curve
void Docs_Demo_ODESystemSolution_ParametricCurve()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Solution as Parametric Curve\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nParametric curve methods:\n";
	std::cout << "  getSolAsParamCurve2D(ind1, ind2)       - 2D phase space curve\n";
	std::cout << "  getSolAsParamCurve3D(ind1, ind2, ind3) - 3D phase space curve\n";
	
	// 2D example: Van der Pol limit cycle
	std::cout << "\n--- Van der Pol Limit Cycle (2D phase space) ---\n";
	TestBeds::VanDerPolODE vanderpol(1.0);
	Vector<Real> y0_vdp({0.1, 0.0});  // Start near origin
	
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator2d(vanderpol);
	ODESystemSolution sol2d = integrator2d.integrate(y0_vdp, 0.0, 30.0, 0.1, 1e-6);
	
	auto curve2d = sol2d.getSolAsParamCurve2D(0, 1);  // x vs v
	std::cout << "Created 2D parametric curve (x, v) from t=0 to t=30\n";
	std::cout << "Curve at t=15: "; curve2d(15.0).Print(std::cout, 6, 3); std::cout << std::endl;
	
	// 3D example: Lorenz attractor
	std::cout << "\n--- Lorenz Attractor (3D phase space) ---\n";
	TestBeds::LorenzSystemODE lorenz(10.0, 28.0, 8.0/3.0);
	Vector<Real> y0_lorenz({1.0, 1.0, 1.0});
	
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator3d(lorenz);
	ODESystemSolution sol3d = integrator3d.integrate(y0_lorenz, 0.0, 50.0, 0.05, 1e-8);
	
	auto curve3d = sol3d.getSolAsParamCurve3D(0, 1, 2);  // x, y, z
	std::cout << "Created 3D parametric curve (x, y, z) from t=0 to t=50\n";
	std::cout << "Curve at t=25: "; curve3d(25.0).Print(std::cout, 6, 3); std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       PHYSICS EXAMPLES                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Simple harmonic oscillator with analytical verification
void Docs_Demo_Physics_SHO()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Simple Harmonic Oscillator\n";
	std::cout << "==========================================================================\n";
	
	// x'' + ω²x = 0, written as system: x' = v, v' = -ω²x
	Real omega = 2.0;
	TestBeds::SimpleHarmonicOscillatorODE sho(omega);
	
	// Initial conditions: x(0) = 1, v(0) = 0
	Vector<Real> y0({1.0, 0.0});
	Real T = 2.0 * Constants::PI / omega;  // Period
	
	std::cout << "ω = " << omega << ", Period T = " << T << std::endl;
	std::cout << "Initial: x(0) = 1, v(0) = 0\n";
	std::cout << "Analytical: x(t) = cos(ωt), v(t) = -ω·sin(ωt)\n";
	
	// Integrate for 3 periods
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(sho);
	auto sol = integrator.integrate(y0, 0.0, 3*T, 0.1, 1e-10);
	
	std::cout << "\nVerification at t = T, 2T, 3T (should return to x=1, v=0):\n";
	for (int n = 1; n <= 3; n++) {
		Real t = n * T;
		auto interp_x = sol.getSolAsSplineInterp(0);
		auto interp_v = sol.getSolAsSplineInterp(1);
		Real x_num = interp_x(t), v_num = interp_v(t);
		Real x_exact = std::cos(omega * t), v_exact = -omega * std::sin(omega * t);
		std::cout << "  t=" << n << "T: x=" << x_num << " (err=" << std::abs(x_num - x_exact) << ")"
		          << ", v=" << v_num << " (err=" << std::abs(v_num - v_exact) << ")" << std::endl;
	}
}

// Demo: Damped harmonic oscillator
void Docs_Demo_Physics_Damped_Oscillator()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Damped Harmonic Oscillator\n";
	std::cout << "==========================================================================\n";
	
	// x'' + 2ζω₀x' + ω₀²x = 0
	Real omega0 = 1.0, zeta = 0.1;  // Underdamped
	TestBeds::DampedHarmonicOscillatorODE damped(omega0, zeta);
	
	Vector<Real> y0({1.0, 0.0});
	
	std::cout << "ω₀ = " << omega0 << ", ζ = " << zeta << " (underdamped)\n";
	std::cout << "Initial: x(0) = 1, v(0) = 0\n";
	
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(damped);
	auto sol = integrator.integrate(y0, 0.0, 30.0, 0.2, 1e-8);
	
	std::cout << "\nAmplitude decay over time:\n";
	for (Real t = 0; t <= 30.0; t += 5.0) {
		auto interp = sol.getSolAsSplineInterp(0);
		std::cout << "  t=" << t << ": x=" << interp(t) << std::endl;
	}
}

// Demo: Lorenz chaotic system
void Docs_Demo_Physics_Lorenz()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Lorenz Chaotic System\n";
	std::cout << "==========================================================================\n";
	
	// Lorenz equations: dx/dt = σ(y-x), dy/dt = x(ρ-z)-y, dz/dt = xy-βz
	Real sigma = 10.0, rho = 28.0, beta = 8.0/3.0;
	TestBeds::LorenzSystemODE lorenz(sigma, rho, beta);
	
	std::cout << "Parameters: σ=" << sigma << ", ρ=" << rho << ", β=" << beta << std::endl;
	std::cout << "These values produce chaotic behavior.\n";
	
	Vector<Real> y0({1.0, 1.0, 1.0});
	
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(lorenz);
	auto sol = integrator.integrate(y0, 0.0, 50.0, 0.02, 1e-8);
	
	std::cout << "\nIntegration statistics:\n";
	std::cout << "  Total saved points: " << sol.size() << std::endl;
	std::cout << "  Steps OK/rejected: " << sol.getNumStepsOK() << "/" << sol.getNumStepsBad() << std::endl;
	
	// Demonstrate sensitive dependence on initial conditions
	std::cout << "\n--- Sensitive Dependence on Initial Conditions ---\n";
	Vector<Real> y0_perturbed({1.0 + 1e-10, 1.0, 1.0});  // Tiny perturbation
	auto sol_perturbed = integrator.integrate(y0_perturbed, 0.0, 50.0, 0.02, 1e-8);
	
	auto interp1 = sol.getSolAsSplineInterp(0);
	auto interp2 = sol_perturbed.getSolAsSplineInterp(0);
	
	std::cout << "Perturbation: Δx(0) = 1e-10\n";
	std::cout << "Divergence of x-component over time:\n";
	for (Real t : {10.0, 20.0, 30.0, 40.0}) {
		Real diff = std::abs(interp1(t) - interp2(t));
		std::cout << "  t=" << t << ": |Δx| = " << diff << std::endl;
	}
}

// Demo: Predator-prey dynamics (Lotka-Volterra)
void Docs_Demo_Physics_Predator_Prey()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Predator-Prey Dynamics (Lotka-Volterra)\n";
	std::cout << "==========================================================================\n";
	
	// dx/dt = αx - βxy (prey)
	// dy/dt = δxy - γy  (predator)
	TestBeds::LotkaVolterraODE lv;  // Default parameters
	
	std::cout << "Lotka-Volterra equations model predator-prey interaction.\n";
	std::cout << "Prey grows exponentially but is eaten by predators.\n";
	std::cout << "Predators grow when prey is abundant but die without food.\n";
	
	Vector<Real> y0({10.0, 5.0});  // 10 prey, 5 predators
	
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(lv);
	auto sol = integrator.integrate(y0, 0.0, 50.0, 0.1, 1e-6);
	
	std::cout << "\nPopulation dynamics:\n";
	std::cout << "  t=0:  prey=" << y0[0] << ", predators=" << y0[1] << std::endl;
	
	auto interp_prey = sol.getSolAsSplineInterp(0);
	auto interp_pred = sol.getSolAsSplineInterp(1);
	
	for (Real t : {10.0, 20.0, 30.0, 40.0, 50.0}) {
		std::cout << "  t=" << t << ": prey=" << interp_prey(t) 
		          << ", predators=" << interp_pred(t) << std::endl;
	}
	
	std::cout << "\nThe populations oscillate - classic predator-prey cycles.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       ORIGINAL DEMO (Runge-Kutta 4th order)                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_ODE_solvers_RungeKutta_4th_order()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Runge-Kutta 4th Order Methods\n";
	std::cout << "==========================================================================\n";
	
	// get it from predefined test-bed 
	TestBeds::LorenzSystemODE sys0(10.0, 28.0, 8.0 / 3.0);

	const double h1 = 0.01, hmin = 0.0;
	Vector<Real> ystart0({ 2.0, 1.0, 1.0 });

	std::cout << "\n--- Fixed-step RK4 ---\n";
	RungeKutta4_StepCalculator   rk4FixedStepCalc;
	ODESystemFixedStepSolver solver(sys0, rk4FixedStepCalc);
	ODESystemSolution sol = solver.integrate(ystart0, 0.0, 2.0, 200);

	std::cout << "Integrated from t=0 to t=2 with 200 fixed steps\n";
	std::cout << "Final state: "; sol.getXValuesAtEnd().Print(std::cout, 8, 4);
	std::cout << std::endl;

	std::cout << "\n--- Adaptive Cash-Karp (5th order) ---\n";

	ystart0[0] = 2.0;
	ystart0[1] = 1.0;
	ystart0[2] = 1.0;
	ODESystemSolver<RK5_CashKarp_Stepper> solver2(sys0);
	ODESystemSolution sol2 = solver2.integrate(ystart0, 0.0, 2.0, 0.1, 1e-06, h1, hmin);

	std::cout << "Integrated from t=0 to t=2 with adaptive steps\n";
	std::cout << "Steps OK/rejected: " << sol2.getNumStepsOK() << "/" << sol2.getNumStepsBad() << std::endl;
	std::cout << "Final state: "; sol2.getXValuesAtEnd().Print(std::cout, 8, 4);
	std::cout << std::endl;

	Visualizer::VisualizeODESysSolAsMultiFunc(sol2, "Lorenz system solution - RK 4th order", std::vector<std::string>{"x", "y", "z"}, "lorenz_RK_4th.txt");
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       MASTER DEMO FUNCTION                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_ODE_solvers()
{
	std::cout << "\n\n";
	std::cout << "##########################################################################\n";
	std::cout << "#                    ODE SYSTEM SOLVER DEMOS                            #\n";
	std::cout << "##########################################################################\n";
	
	// IODESystem Interface
	Docs_Demo_IODESystem_Interface();
	Docs_Demo_ODESystem_Creation();
	
	// Step Calculators and Solvers
	Docs_Demo_Step_Calculators();
	Docs_Demo_Adaptive_Integrators();
	Docs_Demo_Legacy_Solver_Interface();
	
	// ODESystemSolution post-processing
	Docs_Demo_ODESystemSolution_Basics();
	Docs_Demo_ODESystemSolution_Interpolation();
	Docs_Demo_ODESystemSolution_ParametricCurve();
	
	// Physics examples
	Docs_Demo_Physics_SHO();
	Docs_Demo_Physics_Damped_Oscillator();
	Docs_Demo_Physics_Lorenz();
	Docs_Demo_Physics_Predator_Prey();
	
	// Original demo
	Docs_Demo_ODE_solvers_RungeKutta_4th_order();
}