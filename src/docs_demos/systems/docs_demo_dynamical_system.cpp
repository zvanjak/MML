#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/systems/DynamicalSystem.h>

#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>
#include <mml/algorithms/ODESolvers/ODESolverFixedStep.h>
#include <mml/algorithms/ODESolvers/ODEStepCalculators.h>
#endif

#include <set>
#include <limits>
#include <stdexcept>

using namespace MML;
using namespace MML::Systems;

///////////////////////////////////////////////////////////////////////////////////////////
///                         BASIC CONTINUOUS SYSTEMS                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_BasicSystems()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Built-in Dynamical Systems\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nMML provides classic dynamical systems ready to use.\n";
	
	// Lorenz System
	std::cout << "\n--- Lorenz System ---\n";
	std::cout << "dx/dt = σ(y - x)\n";
	std::cout << "dy/dt = x(ρ - z) - y\n";
	std::cout << "dz/dt = xy - βz\n";
	
	LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);  // σ, ρ, β
	
	std::cout << "\nSystem properties:\n";
	std::cout << "  Dimension: " << lorenz.getDim() << "\n";
	std::cout << "  Parameters: " << lorenz.getNumParam() << "\n";
	std::cout << "  σ = " << lorenz.getParam(0) << " (\"" << lorenz.getParamName(0) << "\")\n";
	std::cout << "  ρ = " << lorenz.getParam(1) << " (\"" << lorenz.getParamName(1) << "\")\n";
	std::cout << "  β = " << lorenz.getParam(2) << " (\"" << lorenz.getParamName(2) << "\")\n";
	std::cout << "  Dissipative: " << (lorenz.isDissipative() ? "yes" : "no") << "\n";
	std::cout << "  Has analytical Jacobian: " << (lorenz.hasAnalyticalJacobian() ? "yes" : "no") << "\n";
	
	Vector<Real> x0 = lorenz.getDefaultInitialCondition();
	std::cout << "  Default IC: "; x0.Print(std::cout, 6, 2); std::cout << "\n";
	
	// Van der Pol
	std::cout << "\n--- Van der Pol Oscillator ---\n";
	std::cout << "dx/dt = y\n";
	std::cout << "dy/dt = μ(1 - x²)y - x\n";
	
	VanDerPolSystem vdp(1.0);
	std::cout << "  Dimension: " << vdp.getDim() << "\n";
	std::cout << "  μ = " << vdp.getParam(0) << "\n";
	std::cout << "  State vars: " << vdp.getStateName(0) << ", " << vdp.getStateName(1) << "\n";
	
	// Rössler
	std::cout << "\n--- Rössler System ---\n";
	std::cout << "dx/dt = -y - z\n";
	std::cout << "dy/dt = x + ay\n";
	std::cout << "dz/dt = b + z(x - c)\n";
	
	RosslerSystem rossler(0.2, 0.2, 5.7);  // a, b, c
	std::cout << "  Parameters: a=" << rossler.getParam(0) 
	          << ", b=" << rossler.getParam(1) 
	          << ", c=" << rossler.getParam(2) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         INTEGRATING TRAJECTORIES                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_Trajectories()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Integrating Trajectories\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nAll continuous systems work with MML's ODE solvers.\n";
	
	LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
	
	// Use fixed-step RK4 solver
	RungeKutta4_StepCalculator rk4;
	ODESystemFixedStepSolver solver(lorenz, rk4);
	
	Vector<Real> x0({1.0, 1.0, 1.0});
	
	std::cout << "\n--- Lorenz trajectory from (1,1,1) ---\n";
	auto sol = solver.integrate(x0, 0.0, 10.0, 1000);  // 1000 steps
	
	std::cout << "Integrated " << sol.getNumSteps() << " steps from t=0 to t=10\n";
	std::cout << "\nSample points:\n";
	std::cout << "  t=0.0:  (" << sol.getXValue(0, 0) << ", " << sol.getXValue(0, 1) << ", " << sol.getXValue(0, 2) << ")\n";
	int mid = sol.getNumSteps() / 2;
	std::cout << "  t=5.0:  (" << sol.getXValue(mid, 0) << ", " << sol.getXValue(mid, 1) << ", " << sol.getXValue(mid, 2) << ")\n";
	int last = sol.getNumSteps();
	std::cout << "  t=10.0: (" << sol.getXValue(last, 0) << ", " 
	          << sol.getXValue(last, 1) << ", " << sol.getXValue(last, 2) << ")\n";
	
	// Van der Pol limit cycle
	std::cout << "\n--- Van der Pol limit cycle ---\n";
	VanDerPolSystem vdp(1.0);
	ODESystemFixedStepSolver vdpSolver(vdp, rk4);
	
	Vector<Real> vdp_x0({2.0, 0.0});
	auto vdp_sol = vdpSolver.integrate(vdp_x0, 0.0, 20.0, 2000);
	
	std::cout << "Integrated Van der Pol for t=0 to t=20\n";
	int lastIdx = vdp_sol.getNumSteps();
	std::cout << "Final state: (" << vdp_sol.getXValue(lastIdx, 0) << ", " << vdp_sol.getXValue(lastIdx, 1) << ")\n";
	std::cout << "(Approaches limit cycle)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         FIXED POINT ANALYSIS                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_FixedPointContracts()
{
	class DemoSystem : public DynamicalSystemBase<1, 0> {
	public:
		bool noRoot = false;
		bool autonomous = true;
		bool isAutonomous() const override { return autonomous; }
		void derivs(Real, const Vector<Real>& x, Vector<Real>& f) const override {
			f[0] = noRoot ? x[0] * x[0] + 1 : x[0] - 1;
		}
		void jacobian(Real, const Vector<Real>& x, Matrix<Real>& matrix) const override {
			matrix(0, 0) = noRoot ? 2 * x[0] : 1;
		}
	};
	DemoSystem system;
	const auto finalUpdate = FindFixedPointDetailed(system, Vector<Real>{0}, Precision::FixedPointTolerance, 1);
	if (!finalUpdate.IsSuccess() || finalUpdate.fixed_point.location[0] != 1
		|| finalUpdate.fixed_point.eigenvalues.size() != 1 || finalUpdate.function_evaluations != 2)
		throw std::runtime_error("The final permitted Newton update was not finalized");
	system.noRoot = true;
	const auto failed = FindFixedPointDetailed(system, Vector<Real>{1}, Precision::FixedPointTolerance, 5);
	if (failed.status != AlgorithmStatus::SingularMatrix || failed.fixed_point.location[0] != 0
		|| failed.fixed_point.convergenceResidual != 1 || failed.error_message.empty())
		throw std::runtime_error("Singular Newton breakdown did not preserve a coherent candidate");
	system.autonomous = false;
	const auto invalid = FixedPointFinder::Find(system, Vector<Real>{1});
	if (invalid.status != AlgorithmStatus::InvalidInput || invalid.function_evaluations != 0)
		throw std::runtime_error("A time-dependent system was accepted as an autonomous equilibrium");

	FixedPoint<Real> linear;
	linear.jacobian = Matrix<Real>(2, 2, {0, -1, 1, 0});
	FixedPointFinder::ClassifyFixedPoint(linear);
	if (linear.status != AlgorithmStatus::Success || linear.type != FixedPointType::Center || linear.isStable)
		throw std::runtime_error("The linear center label was misclassified");
	linear.jacobian = Matrix<Real>(3, 3, {0, -1, 0, 1, 0, 0, 0, 0, -1});
	FixedPointFinder::ClassifyFixedPoint(linear);
	if (linear.status != AlgorithmStatus::Success || linear.type != FixedPointType::Unknown || linear.isStable)
		throw std::runtime_error("A neutral/decaying spectrum was incorrectly called unstable");
	VanDerPolSystem vdp(1);
	const auto roots = FixedPointFinder::FindMultiple(vdp,
		{Vector<Real>{0, 0}, Vector<Real>{1, 1}, Vector<Real>{-1, -1}, Vector<Real>{1, -1}});
	if (roots.size() != 1 || roots.front().status != AlgorithmStatus::Success)
		throw std::runtime_error("Multiple-start fixed-point validation or deduplication failed");
	std::cout << "Fixed-point contracts passed: final-update classification, coherent singular result, "
		"autonomous-only inputs, neutral linear labels and finite unique roots.\n";
}

void Docs_Demo_DynamicalSystem_FixedPoints()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Fixed Point Analysis\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFixed points are equilibria where dx/dt = 0.\n";
	std::cout << "FixedPointFinder uses Newton's method to locate them.\n";
	
	// Van der Pol has one fixed point at origin
	std::cout << "\n--- Van der Pol Fixed Point ---\n";
	VanDerPolSystem vdp(1.0);
	
	Vector<Real> guess({0.1, 0.1});
	auto fp = FixedPointFinder::Find(vdp, guess);
	if (fp.status != AlgorithmStatus::Success)
		throw std::runtime_error(fp.error_message);
	
	std::cout << "Initial guess: "; guess.Print(std::cout, 6, 2); std::cout << "\n";
	std::cout << "Fixed point found: "; fp.location.Print(std::cout, 10, 6); std::cout << "\n";
	std::cout << "Convergence residual: " << std::scientific << fp.convergenceResidual << "\n";
	std::cout << "Newton iterations: " << fp.iterations << "\n";
	std::cout << "Type: " << to_string(fp.type) << "\n";
	std::cout << "Asymptotic stability certified by linearization: " << (fp.isStable ? "yes" : "no") << "\n";
	std::cout << "Direct function evaluations: " << fp.function_evaluations << "\n";
	
	std::cout << "\nEigenvalues of Jacobian:\n";
	for (size_t i = 0; i < fp.eigenvalues.size(); i++) {
		std::cout << "  λ" << i << " = " << std::fixed << std::setprecision(6) 
		          << fp.eigenvalues[i].real();
		if (std::abs(fp.eigenvalues[i].imag()) > 1e-10)
			std::cout << " + " << fp.eigenvalues[i].imag() << "i";
		std::cout << "\n";
	}
	
	// Lorenz fixed points
	std::cout << "\n--- Lorenz Fixed Points ---\n";
	LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
	
	// Origin
	auto fp_origin = FixedPointFinder::Find(lorenz, Vector<Real>({0.0, 0.0, 0.0}));
	std::cout << "Origin: "; fp_origin.location.Print(std::cout, 8, 4); 
	std::cout << " - " << to_string(fp_origin.type) << "\n";
	
	// C+ fixed point (around sqrt(β(ρ-1)), sqrt(β(ρ-1)), ρ-1)
	Real val = std::sqrt(8.0/3.0 * (28.0 - 1.0));
	auto fp_cplus = FixedPointFinder::Find(lorenz, Vector<Real>({val, val, 27.0}));
	std::cout << "C+ point: "; fp_cplus.location.Print(std::cout, 8, 4);
	std::cout << " - " << to_string(fp_cplus.type) << "\n";
	
	// C- fixed point
	auto fp_cminus = FixedPointFinder::Find(lorenz, Vector<Real>({-val, -val, 27.0}));
	std::cout << "C- point: "; fp_cminus.location.Print(std::cout, 8, 4);
	std::cout << " - " << to_string(fp_cminus.type) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         FIXED POINT TYPES                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_FixedPointTypes()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Fixed Point Classification\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFixed points are classified by their Jacobian eigenvalues:\n";
	std::cout << "- StableNode: All eigenvalues have negative real parts\n";
	std::cout << "- UnstableNode: All eigenvalues have positive real parts\n";
	std::cout << "- Saddle: Mixed signs (some positive, some negative)\n";
	std::cout << "- StableFocus: Complex eigenvalues with negative real part\n";
	std::cout << "- UnstableFocus: Complex eigenvalues with positive real part\n";
	std::cout << "- Center: Purely imaginary eigenvalues\n";
	
	// Different Van der Pol μ values give different stability
	std::cout << "\n--- Van der Pol: Effect of μ on Stability ---\n";
	
	std::vector<Real> mus = {0.1, 0.5, 1.0, 2.0, 5.0};
	for (Real mu : mus) {
		VanDerPolSystem vdp(mu);
		auto fp = FixedPointFinder::Find(vdp, Vector<Real>({0.0, 0.0}));
		
		std::cout << "μ = " << std::fixed << std::setprecision(1) << mu 
		          << ": " << std::setw(16) << to_string(fp.type);
		
		// Show eigenvalues
		std::cout << "  eigenvalues: ";
		for (const auto& ev : fp.eigenvalues) {
			std::cout << std::setprecision(3) << ev.real();
			if (std::abs(ev.imag()) > 1e-6)
				std::cout << "±" << std::abs(ev.imag()) << "i";
			std::cout << " ";
		}
		std::cout << "\n";
	}
	
	std::cout << "\nNote: For μ > 0 the origin is unstable; the focus/node label depends on μ.\n";
	std::cout << "The system has a stable limit cycle, not a stable fixed point.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         LYAPUNOV EXPONENTS                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_LyapunovExponents()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Lyapunov Exponents\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nLyapunov estimates measure finite-time tangent growth along a trajectory.\n";
	std::cout << "- Positive growth alone does not certify chaos or an attractor.\n";
	std::cout << "- Near-zero growth alone does not classify the dynamics.\n";
	std::cout << "- Negative growth measures local contraction, not global stability.\n";
	std::cout << "- Kaplan-Yorke is a spectrum-based dimension estimate.\n";
	
	// Lorenz system - famous chaotic attractor
	std::cout << "\n--- Lorenz System (σ=10, ρ=28, β=8/3) ---\n";
	LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
	
	std::cout << "Computing Lyapunov exponents (this may take a moment)...\n";
	auto lyap_lorenz = LyapunovAnalyzer::Compute(lorenz,
		Vector<Real>({1.0, 1.0, 1.0}),  // Initial condition
		500.0,   // Total time
		1.0,     // Orthonormalization interval
		0.01     // Initial adaptive IVP step
	);
	
	std::cout << "\nLyapunov exponents: ";
	lyap_lorenz.exponents.Print(std::cout, 10, 4);
	std::cout << "\n";
	std::cout << "Max exponent (λ₁): " << std::fixed << std::setprecision(4) 
	          << lyap_lorenz.maxExponent << "\n";
	std::cout << "Sum of exponents: " << lyap_lorenz.sum << "\n";
	std::cout << "Kaplan-Yorke dimension: " << lyap_lorenz.kaplanYorkeDimension << "\n";
	std::cout << "Positive-growth heuristic: " << (lyap_lorenz.isChaotic ? "detected" : "not detected") << "\n";
	std::cout << "Number of QR orthonormalizations: " << lyap_lorenz.numOrthonormalizations << "\n";
	
	std::cout << "\nExpected for Lorenz: λ₁ ≈ 0.9, λ₂ ≈ 0, λ₃ ≈ -14.6\n";
	std::cout << "Sum ≈ -(σ + 1 + β) = -" << (10 + 1 + 8.0/3.0) << " (dissipative)\n";
	
	// Rössler - simpler chaotic system
	std::cout << "\n--- Rössler System (a=0.2, b=0.2, c=5.7) ---\n";
	RosslerSystem rossler(0.2, 0.2, 5.7);
	
	auto lyap_rossler = LyapunovAnalyzer::Compute(rossler,
		Vector<Real>({1.0, 1.0, 0.0}),
		500.0, 1.0, 0.01
	);
	
	std::cout << "Lyapunov exponents: ";
	lyap_rossler.exponents.Print(std::cout, 10, 4);
	std::cout << "\n";
	std::cout << "Max exponent: " << lyap_rossler.maxExponent << "\n";
	std::cout << "Positive-growth heuristic: " << (lyap_rossler.isChaotic ? "detected" : "not detected") << "\n";
	std::cout << "Kaplan-Yorke dimension: " << lyap_rossler.kaplanYorkeDimension << "\n";
}

void Docs_Demo_SystemsEvaluationContracts()
{
	class CountedDecay : public DynamicalSystemBase<1, 1> {
	public:
		mutable int calls = 0;
		CountedDecay() { setParam(0, 1); }
		void derivs(Real, const Vector<Real>& x, Vector<Real>& derivative) const override {
			++calls;
			derivative[0] = -_params[0] * x[0];
		}
	};
	CountedDecay decay;
	auto idle = IntegrateTrajectoryDetailed(decay, Vector<Real>{1}, 0, 0.25);
	auto trajectory = IntegrateTrajectoryDetailed(decay, Vector<Real>{1}, 1, 0.25);
	std::cout << "Zero-time points and RHS calls: "
		<< idle.points.size() << ", " << idle.function_evaluations << "\n";
	std::cout << "Trajectory RHS calls: " << trajectory.function_evaluations
		<< ", observed: " << decay.calls << "\n";
	if (!idle.IsSuccess() || idle.points.size() != 1 || idle.function_evaluations != 0 ||
		!trajectory.IsSuccess() || trajectory.function_evaluations != decay.calls ||
		!IsWithinAbs(trajectory.points.back()[0], std::exp(Real(-1)), Precision::CoarseTolerance))
		throw std::logic_error("Trajectory evaluation contract demo failed");

	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	auto invalid = IntegrateTrajectoryDetailed(decay, Vector<Real>{}, 0, 0.25, 0.01, config);
	std::cout << "Invalid input: " << invalid.error_message << "\n";
	if (invalid.status != AlgorithmStatus::InvalidInput || invalid.function_evaluations != 0)
		throw std::logic_error("Zero-work input validation demo failed");

	const int beforeSweep = decay.calls;
	auto sweep = SweepBifurcationDetailed(decay, 0, 1, 2, 2, Vector<Real>{1}, 0, 0, 0.1, 0.01);
	std::cout << "Sweep RHS calls: " << sweep.function_evaluations
		<< ", legacy record estimate: " << sweep.diagram.numRecordedPoints << "\n";
	if (!sweep.IsSuccess() || sweep.function_evaluations != decay.calls - beforeSweep ||
		sweep.diagram.numRecordedPoints != 10 || decay.getParam(0) != 1)
		throw std::logic_error("Sweep evaluation contract demo failed");

	const int beforeSection = decay.calls;
	auto section = ComputePoincareSectionDetailed(decay, Vector<Real>{1}, {}, 1, 0.01, config, 0.1, 0);
	std::cout << "Bounded noncrossing search RHS calls: " << section.function_evaluations << "\n";
	if (section.status != AlgorithmStatus::MaxIterationsExceeded ||
		section.function_evaluations != decay.calls - beforeSection || !section.intersections.empty())
		throw std::logic_error("Poincare evaluation contract demo failed");

	const Real bound = std::numeric_limits<Real>::max();
	auto wide = BifurcationAnalyzer::Sweep(
		decay, 0, -bound, bound, 3, Vector<Real>{0}, 0, 0, 0);
	DynamicalSystemAnalyzer<> analyzer(decay);
	auto scan = analyzer.ScanChaosVsParameter(0, 1, 2, 1, Vector<Real>{0}, 0.02);
	auto roots = analyzer.FindFixedPointsInBox(Vector<Real>{-2}, Vector<Real>{2}, 3, 3);
	std::cout << "Extreme-span midpoint: " << wide.parameterValues[1] << "\n";
	std::cout << "Single scan parameter: " << scan[0].first << "\n";
	std::cout << "Grid roots: " << roots.size() << "\n";
	if (wide.parameterValues != std::vector<Real>{-bound, 0, bound} || scan.size() != 1 ||
		scan[0].first != 1 || !IsWithinAbs(scan[0].second, Real(-1), Precision::CoarseTolerance) ||
		roots.size() != 1 || !IsNearlyZero(roots[0].location[0]) || decay.getParam(0) != 1)
		throw std::logic_error("Parameter scan and grid contract demo failed");

	bool rejected = false;
	try { analyzer.FindFixedPointsInBox(Vector<Real>{-2}, Vector<Real>{2}, 3, 2); }
	catch (const ResourceLimitError& error) {
		rejected = true;
		std::cout << "Grid rejected before work: " << error.what() << "\n";
	}
	if (!rejected) throw std::logic_error("Grid budget contract demo failed");

	class DrivenClock : public DynamicalSystemBase<1, 1> {
	public:
		DrivenClock() { setParam(0, 7); }
		bool isAutonomous() const override { return false; }
		void derivs(Real time, const Vector<Real>&, Vector<Real>& derivative) const override {
			derivative[0] = std::cos(time);
		}
	};
	DrivenClock driven;
	auto afterTransient = BifurcationAnalyzer::Sweep(
		driven, 0, 1, 2, 1, Vector<Real>{0}, 0, Constants::PI / 2, Constants::PI, 0.02);
	auto continued = BifurcationAnalyzer::Sweep(
		driven, 0, 1, 3, 3, Vector<Real>{0}, 0, 0, Constants::PI, 0.02);
	std::cout << "Maxima after transient: " << afterTransient.attractorValues[0].size() << "\n";
	std::cout << "Continued maxima: " << continued.attractorValues[0].size()
		<< ", " << continued.attractorValues[1].size()
		<< ", " << continued.attractorValues[2].size() << "\n";
	if (!afterTransient.attractorValues[0].empty() ||
		continued.attractorValues[0].size() != 1 || !continued.attractorValues[1].empty() ||
		continued.attractorValues[2].size() != 1 || driven.getParam(0) != 7 ||
		!IsWithinAbs(continued.attractorValues[0][0], Real(1), Precision::CoarseTolerance) ||
		!IsWithinAbs(continued.attractorValues[2][0], Real(1), Precision::CoarseTolerance))
		throw std::logic_error("Bifurcation physical-clock demo failed");
}

void Docs_Demo_LyapunovContracts()
{
	LorenzSystem lorenz(10.0, 28.0, 8.0 / 3.0);
	auto result = LyapunovAnalyzer::Compute(lorenz, Vector<Real>({1.0, 1.0, 1.0}), 20.0, 0.25, 0.01);
	std::cout << "Lyapunov exponents: " << result.exponents << "\n";
	std::cout << "Max exponent: " << result.maxExponent << "\n";
	std::cout << "Sum: " << result.sum << "\n";
	std::cout << "Kaplan-Yorke dimension: " << result.kaplanYorkeDimension << "\n";
	std::cout << "Positive-growth heuristic: " << result.isChaotic << "\n";

	class StrongDecay : public DynamicalSystemBase<1, 0> {
	public:
		void derivs(Real, const Vector<Real>& x, Vector<Real>& derivative) const override {
			derivative[0] = -100 * x[0];
		}
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
			matrix.resize(1, 1);
			matrix(0, 0) = -100;
		}
	};
	StrongDecay decay;
	auto contraction = LyapunovAnalyzer::Compute(decay, Vector<Real>{0}, 1.0, 1.0, 0.01);
	std::cout << "Expected rate -100, estimated: " << contraction.maxExponent << "\n";
	std::cout << "Accepted QR blocks: " << contraction.numOrthonormalizations << "\n";
	std::cout << "QR refinements: " << contraction.qr_refinements << "\n";
	if (!IsWithinAbs(contraction.maxExponent, Real(-100), Precision::CoarseTolerance))
		throw std::logic_error("Strong-contraction Lyapunov demo is outside tolerance");

	TentMap tent(2.0);
	Vector<Real> x({0.3});
	for (int iteration = 0; iteration < 100; ++iteration)
		x = tent.iterate(x);
	std::cout << "Analytical lambda: " << tent.analyticalLyapunov() << "\n";
	auto lyap = DiscreteMapLyapunov<1>::Compute(tent, x, 10000, 1000);
	std::cout << "Numerical lambda: " << lyap.maxExponent << "\n";
	if (!IsNearlyEqual(lyap.maxExponent, tent.analyticalLyapunov()))
		throw std::logic_error("Tent-map Lyapunov demo is outside tolerance");

	LogisticMap singular(2);
	bool rejected = false;
	try { DiscreteMapLyapunov<1>::Compute(singular, Vector<Real>{Real(0.5)}, 8, 0); }
	catch (const SingularMatrixError& error) {
		rejected = true;
		std::cout << "Expected singular-frame diagnostic: " << error.what() << "\n";
	}
	if (!rejected)
		throw std::logic_error("Zero-stretch Lyapunov demo did not report singularity");
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DISCRETE MAPS - LOGISTIC                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_LogisticMap()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Logistic Map\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nThe logistic map: x_{n+1} = r * x_n * (1 - x_n)\n";
	std::cout << "is the simplest system exhibiting chaos.\n";
	
	std::cout << "\n--- Behavior for Different r Values ---\n";
	
	std::vector<std::pair<Real, std::string>> r_values = {
		{2.5, "fixed point"},
		{3.2, "period-2 cycle"},
		{3.5, "period-4 cycle"},
		{3.9, "chaos"},
		{4.0, "full chaos"}
	};
	
	for (const auto& [r, desc] : r_values) {
		LogisticMap logistic(r);
		
		Vector<Real> x({0.5});
		
		// Iterate past transient
		for (int i = 0; i < 1000; i++)
			x = logistic.iterate(x);
		
		// Collect attractor points
		std::vector<Real> attractor;
		for (int i = 0; i < 20; i++) {
			x = logistic.iterate(x);
			// Add unique values (within tolerance)
			bool isNew = true;
			for (Real v : attractor) {
				if (std::abs(v - x[0]) < 1e-6) {
					isNew = false;
					break;
				}
			}
			if (isNew && attractor.size() < 10)
				attractor.push_back(x[0]);
		}
		
		std::cout << "r = " << std::fixed << std::setprecision(1) << r 
		          << " (" << desc << "): ";
		std::cout << std::setprecision(4);
		for (size_t i = 0; i < std::min(attractor.size(), size_t(5)); i++)
			std::cout << attractor[i] << " ";
		if (attractor.size() > 5) std::cout << "...";
		std::cout << "\n";
	}
	
	// Lyapunov exponent for chaotic case
	std::cout << "\n--- Lyapunov Exponent ---\n";
	LogisticMap logistic4(4.0);
	auto lyap = DiscreteMapLyapunov<1>::Compute(logistic4, Vector<Real>({0.3}), 10000, 1000);
	
	std::cout << "r = 4.0: λ = " << std::setprecision(6) << lyap.maxExponent << "\n";
	std::cout << "Analytical: λ = log(2) = " << std::log(2.0) << "\n";
	std::cout << "Positive-growth heuristic: " << (lyap.isChaotic ? "detected" : "not detected") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DISCRETE MAPS - HÉNON                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_HenonMap()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Hénon Map\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nThe Hénon map is a classic 2D chaotic map:\n";
	std::cout << "x_{n+1} = 1 - a*x_n² + y_n\n";
	std::cout << "y_{n+1} = b*x_n\n";
	
	HenonMap henon(1.4, 0.3);  // Classic parameters
	
	std::cout << "\nClassic parameters: a = " << henon.getA() << ", b = " << henon.getB() << "\n";

	// IDiscreteMap is an IFunction with the IParametrized mixin (docs: Discrete Maps intro)
	std::cout << "getNumParam() = " << henon.getNumParam() << "\n";
	henon.setParams(Vector<Real>({1.2, 0.3}));          // VectorDimensionError if the size is not 2
	Vector<Real> next = henon(Vector<Real>({0.1, 0.1})); // == henon.iterate(...)
	std::cout << "henon(0.1, 0.1) with a = 1.2: (" << next[0] << ", " << next[1] << ")\n";
	henon.setParams(Vector<Real>({1.4, 0.3}));
	
	// Generate trajectory
	std::cout << "\n--- Trajectory on Strange Attractor ---\n";
	Vector<Real> x({0.0, 0.0});
	
	// Skip transient
	for (int i = 0; i < 1000; i++)
		x = henon.iterate(x);
	
	// Sample attractor points
	std::cout << "Points on attractor (after transient):\n";
	for (int i = 0; i < 5; i++) {
		x = henon.iterate(x);
		std::cout << "  (" << std::fixed << std::setprecision(6) 
		          << x[0] << ", " << x[1] << ")\n";
	}
	
	// Compute Lyapunov exponents
	std::cout << "\n--- Lyapunov Spectrum ---\n";
	auto lyap = DiscreteMapLyapunov<2>::Compute(henon, Vector<Real>({0.0, 0.0}), 50000, 1000);
	
	std::cout << "Lyapunov exponents: λ₁ = " << std::setprecision(4) << lyap.exponents[0]
	          << ", λ₂ = " << lyap.exponents[1] << "\n";
	std::cout << "Sum: " << (lyap.exponents[0] + lyap.exponents[1]) << "\n";
	std::cout << "Positive-growth heuristic: " << (lyap.isChaotic ? "detected" : "not detected") << "\n";
	
	// Jacobian determinant
	std::cout << "\nJacobian determinant: " << henon.jacobianDeterminant() << " (= -b, constant)\n";
	std::cout << "This means the map contracts area by factor |b| = 0.3\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DISCRETE MAPS - STANDARD MAP                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_StandardMap()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Standard Map (Chirikov-Taylor)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nThe standard map models a kicked rotor:\n";
	std::cout << "p_{n+1} = p_n + K*sin(θ_n)\n";
	std::cout << "θ_{n+1} = θ_n + p_{n+1}   (mod 2π)\n";
	
	std::cout << "\n--- Different K Values ---\n";
	std::cout << "K = 0: Integrable (all circles)\n";
	std::cout << "Near K = 0.97: familiar torus-breakup regime, not a global chaos certificate\n";
	std::cout << "K > 1: Widespread chaos with islands\n";
	
	std::vector<Real> K_values = {0.5, 0.9716, 2.0};
	
	for (Real K : K_values) {
		StandardMap stdmap(K);
		
		std::cout << "\nK = " << std::fixed << std::setprecision(4) << K << ":\n";
		std::cout << "  Area preserving: " << (stdmap.isAreaPreserving() ? "yes" : "no") << "\n";
		
		// Compute Lyapunov exponent
		auto lyap = DiscreteMapLyapunov<2>::Compute(stdmap, 
			Vector<Real>({0.5, 1.0}), 10000, 500);
		
		std::cout << "  Max Lyapunov exponent: " << std::setprecision(4) << lyap.maxExponent << "\n";
		std::cout << "  Positive-growth heuristic: " << (lyap.isChaotic ? "detected" : "not detected") << "\n";
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         TENT MAP                                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_TentMap()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Tent Map\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nThe tent map is piecewise-linear:\n";
	std::cout << "x_{n+1} = μ*x_n        if x_n < 0.5\n";
	std::cout << "x_{n+1} = μ*(1 - x_n)  if x_n ≥ 0.5\n";
	
	std::cout << "\nExact Lyapunov exponent: λ = log(μ) for μ ∈ (1, 2]\n";
	
	std::cout << "\n--- Comparing Numerical vs Analytical ---\n";
	
	std::vector<Real> mus = {1.5, 1.8, 2.0};
	
	for (Real mu : mus) {
		TentMap tent(mu);
		
		// Numerical Lyapunov
		auto lyap = DiscreteMapLyapunov<1>::Compute(tent, Vector<Real>({0.3}), 10000, 500);
		
		std::cout << "μ = " << std::fixed << std::setprecision(1) << mu << ": ";
		std::cout << "numerical λ = " << std::setprecision(6) << lyap.maxExponent;
		std::cout << ", analytical λ = " << tent.analyticalLyapunov();
		std::cout << ", error = " << std::abs(lyap.maxExponent - tent.analyticalLyapunov()) << "\n";
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         BIFURCATION ANALYSIS                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_Bifurcation()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Bifurcation Analysis\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFinite post-transient samples illustrate parameter dependence.\n";
	std::cout << "Rounded distinct-value counts do not certify periodicity or chaos.\n";
	
	// Simple demonstration with Logistic map
	std::cout << "\n--- Logistic Map Period-Doubling Route to Chaos ---\n";
	std::cout << "\nSampling attractor at different r values:\n";
	
	std::cout << std::setw(6) << "r" << std::setw(15) << "Rounded sample" << "  Values\n";
	std::cout << std::string(50, '-') << "\n";
	
	std::vector<Real> r_samples = {2.8, 3.0, 3.2, 3.45, 3.55, 3.7, 3.9};
	
	for (Real r : r_samples) {
		LogisticMap logistic(r);
		
		Vector<Real> x({0.5});
		
		// Transient
		for (int i = 0; i < 2000; i++)
			x = logistic.iterate(x);
		
		// Collect unique attractor points
		std::set<Real> unique_vals;
		for (int i = 0; i < 100; i++) {
			x = logistic.iterate(x);
			// Round to detect periodicity
			Real rounded = std::round(x[0] * 10000) / 10000.0;
			unique_vals.insert(rounded);
		}
		
		std::string type;
		if (unique_vals.size() <= 1) type = "<=1 value";
		else if (unique_vals.size() <= 2) type = "<=2 values";
		else if (unique_vals.size() <= 4) type = "<=4 values";
		else if (unique_vals.size() <= 8) type = "<=8 values";
		else type = "many values";
		
		std::cout << std::fixed << std::setprecision(2) << std::setw(6) << r 
		          << std::setw(15) << type << "  ";
		
		int count = 0;
		for (Real v : unique_vals) {
			if (count++ >= 4) { std::cout << "..."; break; }
			std::cout << std::setprecision(4) << v << " ";
		}
		std::cout << "\n";
	}
	
	std::cout << "\nFeigenbaum constant δ ≈ 4.669 controls period-doubling spacing.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         PARAMETER MANIPULATION                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_ModelDomains()
{
	LorenzSystem lorenz;
	lorenz.setParam(1, 150); // Suggested rho range is [0, 100], not an admissibility limit.
	std::cout << "Accepted rho outside suggested range: " << lorenz.getParam(1) << "\n";

	DoublePendulumSystem pendulum(1, 1, 1, 1, -9.80665);
	pendulum.setParam(4, 0); // Zero or reversed gravity is supported.
	const auto original = pendulum.getParams();
	bool rejected = false;
	try {
		pendulum.setParams(Vector<Real>{2, 1, 0, 1, 0});
	}
	catch (const ArgumentError& error) {
		rejected = true;
		std::cout << "Rejected zero length: " << error.what() << "\n";
	}
	if (!rejected || pendulum.getParams() != original)
		throw std::runtime_error("Pendulum domain rejection must preserve all parameters");
	std::cout << "Parameters unchanged: " << (pendulum.getParams() == original) << "\n";

	HenonMap map;
	const auto mapOriginal = map.getParams();
	rejected = false;
	try {
		map.setParams(Vector<Real>{2, std::numeric_limits<Real>::quiet_NaN()});
	}
	catch (const ArgumentError& error) {
		rejected = true;
		std::cout << "Rejected nonfinite map parameter: " << error.what() << "\n";
	}
	if (!rejected || map.getParams() != mapOriginal)
		throw std::runtime_error("Map domain rejection must preserve all parameters");

	rejected = false;
	try {
		lorenz.getStateName(-1);
	}
	catch (const IndexError& error) {
		rejected = true;
		std::cout << "Rejected invalid metadata index: " << error.what() << "\n";
	}
	if (!rejected)
		throw std::runtime_error("Invalid metadata indices must be rejected");
}

void Docs_Demo_DynamicalSystem_Parameters()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Parameter Manipulation\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nAll systems inherit parameter management from DynamicalSystemBase.\n";
	
	LorenzSystem lorenz;
	
	std::cout << "\n--- Lorenz System Parameters ---\n";
	
	// Get parameter info
	for (int i = 0; i < lorenz.getNumParam(); i++) {
		auto [minVal, maxVal] = lorenz.getParamRange(i);
		std::cout << "  " << lorenz.getParamName(i) << " = " << lorenz.getParam(i)
		          << " (suggested range: [" << minVal << ", " << maxVal << "])\n";
	}
	
	// Modify parameters
	std::cout << "\nModifying parameters:\n";
	lorenz.setParam(1, 24.0);
	std::cout << "  Set rho = " << lorenz.getParam(1) << " (parameter change, not a dynamics classification)\n";
	
	// Get/set all parameters at once
	Vector<Real> params = lorenz.getParams();
	std::cout << "\nAll parameters: "; params.Print(std::cout, 8, 2); std::cout << "\n";
	
	// Restore the default parameters.
	lorenz.setParams(Vector<Real>({10.0, 28.0, 8.0/3.0}));
	std::cout << "Restored default parameters: sigma=" << lorenz.getParam(0)
	          << ", ρ=" << lorenz.getParam(1) << ", β=" << lorenz.getParam(2) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         JACOBIAN COMPUTATION                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_Jacobian()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Jacobian Computation\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nThe Jacobian matrix J = ∂f/∂x determines local stability.\n";
	std::cout << "Systems can provide analytical Jacobians for efficiency.\n";
	
	LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
	
	std::cout << "\n--- Lorenz Jacobian at Different Points ---\n";
	std::cout << "Has analytical Jacobian: " << (lorenz.hasAnalyticalJacobian() ? "yes" : "no") << "\n";
	
	Matrix<Real> J(3, 3);
	
	// At origin
	Vector<Real> x_origin({0.0, 0.0, 0.0});
	lorenz.jacobian(0.0, x_origin, J);
	std::cout << "\nAt origin (0,0,0):\n";
	J.Print(std::cout, 10, 4);
	
	// At a point on the attractor
	Vector<Real> x_attractor({-8.5, -8.5, 27.0});  // Near C- fixed point
	lorenz.jacobian(0.0, x_attractor, J);
	std::cout << "\nAt (-8.5, -8.5, 27):\n";
	J.Print(std::cout, 10, 4);
	
	// Van der Pol Jacobian
	std::cout << "\n--- Van der Pol Jacobian at Origin ---\n";
	VanDerPolSystem vdp(1.0);
	Matrix<Real> J_vdp(2, 2);
	vdp.jacobian(0.0, Vector<Real>({0.0, 0.0}), J_vdp);
	J_vdp.Print(std::cout, 10, 4);
	
	std::cout << "\nNote: J = [[0, 1], [-1, μ]] at origin\n";
	std::cout << "Eigenvalues have positive real part for μ > 0 → unstable.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         FINITE-TIME INDICATOR WORKFLOW                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_ChaosDetection()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Qualified Finite-Time Indicators\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFinite-time indicators are not a chaos certificate.\n";
	std::cout << "Check duration, transients and numerical accuracy separately.\n";
	std::cout << "Boundedness, invariant sets and aperiodicity require independent evidence.\n";
	
	LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
	
	std::cout << "\n--- Testing Lorenz System ---\n";
	
	// 1. Check fixed points
	std::cout << "\n1. Fixed Point Analysis:\n";
	auto fp = FixedPointFinder::Find(lorenz, Vector<Real>({0.0, 0.0, 0.0}));
	std::cout << "   Origin: " << to_string(fp.type) << " (unstable)\n";
	
	Real val = std::sqrt(8.0/3.0 * 27.0);
	auto fp2 = FixedPointFinder::Find(lorenz, Vector<Real>({val, val, 27.0}));
	std::cout << "   C+ point: " << to_string(fp2.type) << "\n";
	
	// 2. Compute Lyapunov exponents
	std::cout << "\n2. Lyapunov Exponents:\n";
	auto lyap = LyapunovAnalyzer::Compute(lorenz,
		Vector<Real>({1.0, 1.0, 1.0}), 20.0, 0.25, 0.01);
	if (lyap.isChaotic && lyap.sum < -Precision::LooseTolerance) {
		std::cout << "Positive tangent growth with a negative finite-time sum\n";
		std::cout << "Kaplan-Yorke estimate: " << lyap.kaplanYorkeDimension << "\n";
		std::cout << "Boundedness and chaotic-attractor assumptions require independent evidence\n";
	}
	
	std::cout << "   λ₁ = " << std::fixed << std::setprecision(4) << lyap.exponents[0] 
	          << (lyap.exponents[0] > 0 ? " > 0 ✓" : " ≤ 0") << "\n";
	std::cout << "   λ₂ = " << lyap.exponents[1] << "\n";
	std::cout << "   λ₃ = " << lyap.exponents[2] << "\n";
	
	std::cout << "\n3. Finite-Time Volume Indicator:\n";
	std::cout << "   Sum of exponents: " << lyap.sum << "\n";
	std::cout << "   Negative-sum heuristic: " << (lyap.sum < -Precision::LooseTolerance ? "detected" : "not detected") << "\n";
	std::cout << "   A negative sum is not a boundedness proof.\n";
	
	std::cout << "\n4. Spectrum-Based Dimension Estimate:\n";
	std::cout << "   Kaplan-Yorke estimate: " << lyap.kaplanYorkeDimension << "\n";
	std::cout << "   Not a measured fractal dimension or evidence of a strange attractor.\n";
	
	std::cout << "\n--- Indicator ---\n";
	std::cout << "Positive-growth heuristic: " << (lyap.isChaotic ? "detected" : "not detected") << "\n";
}

void Docs_Demo_ReportQualifications()
{
	DynamicalSystemReport<> empty;
	if (empty.dimension != 0 || empty.isChaotic || !std::isnan(empty.kaplanYorkeDimension) ||
		!std::isnan(empty.lyapunovTime) || !empty.summary.empty())
		throw std::logic_error("Uncomputed report contract demo failed");

	class ConstantDrift : public DynamicalSystemBase<1, 0> {
	public:
		void derivs(Real, const Vector<Real>&, Vector<Real>& derivative) const override {
			derivative[0] = 1;
		}
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
			matrix = Matrix<Real>(1, 1, {0});
		}
	};
	ConstantDrift drift;
	auto driftReport = DynamicalSystemAnalyzer<>(drift).Analyze(Vector<Real>{0}, {}, 1);
	std::cout << driftReport.summary;
	if (driftReport.isChaotic || !driftReport.isConservative ||
		driftReport.summary.find("No boundedness or convergence conclusion") == std::string::npos ||
		driftReport.summary.find("Trajectories converge or remain bounded") != std::string::npos)
		throw std::logic_error("Drift report qualification demo failed");

	class ExpandingContraction : public DynamicalSystemBase<2, 0> {
	public:
		void derivs(Real, const Vector<Real>& state, Vector<Real>& derivative) const override {
			derivative[0] = state[0];
			derivative[1] = -2 * state[1];
		}
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
			matrix = Matrix<Real>(2, 2, {1, 0, 0, -2});
		}
	};
	ExpandingContraction expanding;
	auto expandingReport = DynamicalSystemAnalyzer<>(expanding).Analyze(Vector<Real>{1, 0}, {}, 1);
	std::cout << expandingReport.summary;
	if (!expandingReport.isChaotic || !expandingReport.isDissipative ||
		!IsWithinAbs(expandingReport.kaplanYorkeDimension, Real(1.5), Precision::CoarseTolerance) ||
		expandingReport.summary.find("Attractor has dimension") != std::string::npos)
		throw std::logic_error("Expanding negative-sum report demo failed");

	LorenzSystem lorenz;
	DynamicalSystemAnalyzer<Real> analyzer(lorenz);
	DynamicalSystemReport<Real> report = analyzer.Analyze(Vector<Real>({1, 1, 1}), {}, 20);
	std::cout << report.summary; // local fixed points and qualified finite-time indicators
	if (report.lyapunov.totalTime != 20 || report.summary.empty())
		throw std::logic_error("Bounded Lorenz report demo failed");
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN ENTRY POINT                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_SystemsApiCustomModels()
{
	class DampedOscillator : public DynamicalSystemBase<2, 2> {
	public:
		DampedOscillator(Real frequency = 1, Real damping = 0.2) {
			setParams(Vector<Real>{frequency, damping});
			_stateNames = {"position", "velocity"};
			_paramNames = {"frequency", "damping"};
		}
		void derivs(Real, const Vector<Real>& x, Vector<Real>& derivative) const override {
			derivative[0] = x[1];
			derivative[1] = -_params[0] * _params[0] * x[0] - _params[1] * x[1];
		}
		bool hasAnalyticalJacobian() const override { return true; }
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
			matrix = Matrix<Real>(2, 2, {0, 1, -_params[0] * _params[0], -_params[1]});
		}
		bool isDissipative() const override { return _params[1] > 0; }
	};
	DampedOscillator oscillator;
	Vector<Real> derivative(2);
	oscillator.derivs(0, Vector<Real>{1, 0}, derivative);
	Matrix<Real> matrix;
	oscillator.jacobian(0, Vector<Real>{1, 0}, matrix);
	std::cout << "Custom oscillator derivative: " << derivative << "\n";
	if (derivative != Vector<Real>{0, -1} || matrix.rows() != 2 || matrix(0, 1) != 1)
		throw std::logic_error("Custom oscillator example failed");

	class LotkaVolterraSystem : public DynamicalSystemBase<2, 4> {
	public:
		LotkaVolterraSystem(Real alpha = 1, Real beta = 0.1, Real gamma = 1.5, Real delta = 0.075) {
			setParams(Vector<Real>{alpha, beta, gamma, delta});
			_stateNames = {"prey", "predator"};
			_paramNames = {"alpha", "beta", "gamma", "delta"};
		}
		void derivs(Real, const Vector<Real>& x, Vector<Real>& derivative) const override {
			derivative[0] = _params[0] * x[0] - _params[1] * x[0] * x[1];
			derivative[1] = _params[3] * x[0] * x[1] - _params[2] * x[1];
		}
		bool hasAnalyticalJacobian() const override { return true; }
		void jacobian(Real, const Vector<Real>& x, Matrix<Real>& matrix) const override {
			matrix = Matrix<Real>(2, 2, {
				_params[0] - _params[1] * x[1], -_params[1] * x[0],
				_params[3] * x[1], _params[3] * x[0] - _params[2]});
		}
	};
	LotkaVolterraSystem populations;
	auto equilibrium = FixedPointFinder::Find(populations, Vector<Real>{18, 12});
	if (equilibrium.status != AlgorithmStatus::Success)
		throw std::runtime_error(equilibrium.error_message);
	std::cout << "Equilibrium populations: " << equilibrium.location << "\n";
	std::cout << "Linearized type: " << to_string(equilibrium.type) << "\n";
	if (!IsWithinAbs(equilibrium.location[0], Real(20), Precision::CoarseTolerance) ||
		!IsWithinAbs(equilibrium.location[1], Real(10), Precision::CoarseTolerance))
		throw std::logic_error("Population equilibrium example failed");
}

void Docs_Demo_SystemsApiMaps()
{
	LogisticMap logistic(3.9);
	auto orbit = logistic.orbit(Vector<Real>{0.3}, 32);
	auto lyap = DiscreteMapLyapunov<1>::Compute(logistic, Vector<Real>{0.3}, 10000, 1000);
	std::cout << "Logistic samples: " << orbit.size() << "\n";
	std::cout << "Finite-iteration growth estimate: " << lyap.maxExponent << "\n";
	logistic.setParam(0, 4);
	std::cout << "Analytical rate at r=4: " << logistic.analyticalLyapunov() << "\n";
	if (orbit.size() != 32 || lyap.exponents.size() != 1 || !std::isfinite(lyap.maxExponent))
		throw std::logic_error("Logistic API example failed");

	HenonMap henon(1.4, 0.3);
	auto henonOrbit = henon.orbit(Vector<Real>{0, 0}, 32);
	auto henonLyap = DiscreteMapLyapunov<2>::Compute(henon, Vector<Real>{0, 0}, 10000, 1000);
	std::cout << "Henon samples: " << henonOrbit.size() << "\n";
	std::cout << "Henon finite-iteration spectrum: " << henonLyap.exponents << "\n";
	if (henonOrbit.size() != 32 || henonLyap.exponents.size() != 2)
		throw std::logic_error("Henon API example failed");

	StandardMap standard(1);
	Vector<Real> state{0.5, 1};
	for (int iteration = 0; iteration < 32; ++iteration)
		state = standard.iterate(state);
	Matrix<Real> jacobian;
	standard.jacobian(state, jacobian);
	Real determinant = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
	std::cout << "Standard-map Jacobian determinant: " << determinant << "\n";
	if (!IsWithinAbs(determinant, Real(1), Precision::LooseTolerance))
		throw std::logic_error("Standard-map determinant example failed");

	for (Real parameter : {Real(2.8), Real(3.2), Real(3.5), Real(3.9)}) {
		LogisticMap sample(parameter);
		Vector<Real> x{0.3};
		for (int iteration = 0; iteration < 256; ++iteration)
			x = sample.iterate(x);
		for (int iteration = 0; iteration < 8; ++iteration) {
			x = sample.iterate(x);
			std::cout << parameter << " " << x[0] << "\n";
			if (!std::isfinite(x[0]))
				throw std::logic_error("Discrete parameter-sampling example failed");
		}
	}
}

void Docs_Demo_SystemsApiContinuous()
{
	LorenzSystem lorenz(10, 28, 8.0 / 3.0);
	auto lyap = LyapunovAnalyzer::Compute(lorenz, Vector<Real>{1, 1, 1}, 20, 0.25, 0.01);
	std::cout << "Lorenz finite-time maximum: " << lyap.maxExponent << "\n";
	std::cout << "Positive-growth heuristic: " << lyap.isChaotic << "\n";

	RosslerSystem rossler(0.2, 0.2, 5.7);
	auto rosslerLyap = LyapunovAnalyzer::Compute(rossler, Vector<Real>{1, 1, 0}, 20, 0.25, 0.01);
	std::cout << "Rossler finite-time maximum: " << rosslerLyap.maxExponent << "\n";

	VanDerPolSystem vdp(1);
	Vector<Real> initial{2, 0};
	DynamicalSystemAnalyzer<> vdpAnalyzer(vdp);
	auto points = vdpAnalyzer.IntegrateTrajectory(initial, 10, 0.1, 0.01);
	vdp.setParam(0, 5);
	std::cout << "Van der Pol trajectory points: " << points.size() << "\n";

	DuffingSystem duffing(0.3, -1, 1, 0.5, 1.2); // delta, alpha, beta, gamma, omega
	Vector<Real> duffingDerivative(3);
	duffing.derivs(0, Vector<Real>{1, 0, 0}, duffingDerivative);
	std::cout << "Duffing phase rate: " << duffingDerivative[2] << "\n";
	std::cout << "Current autonomy metadata: " << duffing.isAutonomous() << "\n";

	ChuaCircuit chua(15.6, 28, -1.143, -0.714);
	Vector<Real> chuaDerivative(3);
	chua.derivs(0, Vector<Real>{0.1, 0.1, 0.1}, chuaDerivative);
	std::cout << "Chua derivative: " << chuaDerivative << "\n";
	if (lyap.exponents.size() != 3 || rosslerLyap.exponents.size() != 3 ||
		points.empty() || duffingDerivative[2] != duffing.getParam(4))
		throw std::logic_error("Continuous model API examples failed");
	for (Real value : chuaDerivative)
		if (!std::isfinite(value))
			throw std::logic_error("Chua API example failed");

	HenonHeilesSystem hamiltonian;
	auto state = hamiltonian.getDefaultInitialCondition();
	std::cout << hamiltonian.getInvariantName(0) << ": " << hamiltonian.computeInvariant(0, state) << "\n";
	std::cout << "Independent E parameter: " << hamiltonian.getParam(0) << "\n";
	DoublePendulumSystem pendulum;
	std::cout << "Pendulum energy: " << pendulum.computeInvariant(0, pendulum.getDefaultInitialCondition()) << "\n";
	if (!std::isfinite(hamiltonian.computeInvariant(0, state)) ||
		!std::isfinite(pendulum.computeInvariant(0, pendulum.getDefaultInitialCondition())))
		throw std::logic_error("Invariant API examples failed");
}

void Docs_Demo_SystemsApiWorkflow()
{
	LorenzSystem lorenz(10, 28, 8.0 / 3.0);
	DynamicalSystemAnalyzer<> analyzer(lorenz);
	auto fixedPoint = analyzer.FindFixedPoint(Vector<Real>{0, 0, 0});
	if (fixedPoint.status != AlgorithmStatus::Success)
		throw std::runtime_error(fixedPoint.error_message);
	std::cout << "Origin linearized type: " << to_string(fixedPoint.type) << "\n";
	auto lyap = analyzer.ComputeLyapunov(Vector<Real>{1, 1, 1}, 20, 0.25, 0.01);
	auto bifurcation = analyzer.ComputeBifurcation(1, 24, 28, 3, Vector<Real>{1, 1, 1}, 2, 2, 4);
	auto trajectory = analyzer.IntegrateTrajectory(Vector<Real>{1, 1, 1}, 2, 0.1, 0.01);
	std::cout << "Spectrum: " << lyap.exponents << "\n";
	std::cout << "Sweep samples: " << bifurcation.parameterValues.size() << "\n";
	std::cout << "Trajectory samples: " << trajectory.size() << "\n";
	if (bifurcation.parameterValues.size() != 3 || trajectory.size() != 21 || lorenz.getParam(1) != 28)
		throw std::logic_error("Bounded analysis workflow example failed");

	PoincareSection<Real> section(2, 27, +1);
	auto crossings = analyzer.ComputePoincareSection(Vector<Real>{1, 1, 20}, section, 4, 0.01, 50, 1);
	std::cout << "Bounded Poincare crossings: " << crossings.size() << "\n";
	if (crossings.size() != 4)
		throw std::logic_error("Poincare API example failed");

	VanDerPolSystem vdp(0.5);
	auto diagram = BifurcationAnalyzer::Sweep(vdp, 0, 0.5, 1.5, 3, Vector<Real>{2, 0}, 0, 5, 10, 0.02);
	for (std::size_t sample = 0; sample < diagram.parameterValues.size(); ++sample)
		for (Real maximum : diagram.attractorValues[sample])
			std::cout << diagram.parameterValues[sample] << " " << maximum << "\n";
	if (diagram.parameterValues.size() != 3 || vdp.getParam(0) != Real(0.5))
		throw std::logic_error("Continuous parameter-sweep example failed");

	RungeKutta4_StepCalculator rk4;
	ODESystemFixedStepSolver fixedSolver(lorenz, rk4);
	auto fixedSolution = fixedSolver.integrate(Vector<Real>{1, 1, 1}, 0, 1, 100);
	ODEAdaptiveIntegrator<DormandPrince5_Stepper> adaptive(lorenz);
	ODEIntegratorConfig config;
	config.initial_step_size = 0.01;
	config.output_interval = 0.1;
	auto adaptiveSolution = adaptive.integrate(Vector<Real>{1, 1, 1}, 0, 1, config);
	std::cout << "Fixed-step stored points: " << fixedSolution.size() << "\n";
	std::cout << "Adaptive stored points: " << adaptiveSolution.size() << "\n";
	if (fixedSolution.size() != 101 || adaptiveSolution.size() != 11)
		throw std::logic_error("ODE solver API examples failed");
}

void Docs_Demo_SystemsApiOverview()
{
	VanDerPolSystem vdp(1);
	DynamicalSystemAnalyzer<> analyzer(vdp);
	auto fixedPoint = analyzer.FindFixedPoint(Vector<Real>{0.1, 0.1});
	if (fixedPoint.status != AlgorithmStatus::Success)
		throw std::runtime_error(fixedPoint.error_message);
	std::cout << "Van der Pol linearized type: " << to_string(fixedPoint.type) << "\n";
	auto points = analyzer.FindFixedPointsInBox(Vector<Real>{-2, -2}, Vector<Real>{2, 2}, 3);
	std::cout << "Distinct successful box roots: " << points.size() << "\n";
	if (points.empty())
		throw std::logic_error("Overview box search example found no equilibrium");
	std::vector<Vector<Real>> guesses = {
		Vector<Real>({0.0, 0.0}),
		Vector<Real>({1.0, 1.0}),
		Vector<Real>({-1.0, -1.0}),
		Vector<Real>({1.0, -1.0})
	};
	auto fixedPoints = FixedPointFinder::FindMultiple(vdp, guesses);
	std::cout << "Found " << fixedPoints.size() << " unique fixed points:\n";
	for (const auto& fp : fixedPoints) {
		std::cout << "  Location: " << fp.location
			<< " Type: " << to_string(fp.type) << "\n";
	}
	if (fixedPoints.size() != 1)
		throw std::logic_error("Multiple-guess API example failed");

	LorenzSystem lorenz;
	DynamicalSystemAnalyzer<> lorenzAnalyzer(lorenz);
	Vector<Real> initial{1, 1, 1};
	auto lyap = lorenzAnalyzer.ComputeLyapunov(initial, 20, 0.25, 0.01);
	bool positiveGrowth = lorenzAnalyzer.IsChaotic(initial, 20);
	Real growthTime = lorenzAnalyzer.GetLyapunovTime(initial, 20);
	std::cout << "Finite-time spectrum: " << lyap.exponents << "\n";
	std::cout << "Positive-growth heuristic: " << positiveGrowth << "\n";
	std::cout << "Reciprocal detected growth rate: " << growthTime << "\n";
	if (!std::isfinite(lyap.maxExponent))
		throw std::logic_error("Overview finite-time example failed");
}

void Docs_Demo_DynamicalSystem()
{
	std::cout << "\n";
	std::cout << "##########################################################################\n";
	std::cout << "###              DynamicalSystem Documentation Demos                  ###\n";
	std::cout << "##########################################################################\n";
	std::cout << "\nDynamicalSystem.h provides comprehensive tools for analyzing\n";
	std::cout << "continuous and discrete dynamical systems, including:\n";
	std::cout << "- Built-in classic systems (Lorenz, Rössler, Van der Pol, etc.)\n";
	std::cout << "- Fixed point finding and classification\n";
	std::cout << "- Lyapunov exponent computation\n";
	std::cout << "- Bifurcation analysis\n";
	std::cout << "- Discrete maps (Logistic, Hénon, Standard, Tent)\n";
	
	Docs_Demo_DynamicalSystem_BasicSystems();
	Docs_Demo_DynamicalSystem_Trajectories();
	Docs_Demo_DynamicalSystem_FixedPoints();
	Docs_Demo_DynamicalSystem_FixedPointTypes();
	Docs_Demo_DynamicalSystem_LyapunovExponents();
	Docs_Demo_DynamicalSystem_LogisticMap();
	Docs_Demo_DynamicalSystem_HenonMap();
	Docs_Demo_DynamicalSystem_StandardMap();
	Docs_Demo_DynamicalSystem_TentMap();
	Docs_Demo_DynamicalSystem_Bifurcation();
	Docs_Demo_DynamicalSystem_Parameters();
	Docs_Demo_DynamicalSystem_Jacobian();
	Docs_Demo_DynamicalSystem_ChaosDetection();
	
	std::cout << "\n";
	std::cout << "##########################################################################\n";
	std::cout << "###            DynamicalSystem Demos Completed                        ###\n";
	std::cout << "##########################################################################\n";
}
