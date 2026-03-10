#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "systems/DynamicalSystem.h"

#include "algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#include "algorithms/ODESolvers/ODEFixedStepIntegrators.h"
#include "algorithms/ODESolvers/ODESystemStepCalculators.h"
#endif

#include <set>

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
	
	std::cout << "Initial guess: "; guess.Print(std::cout, 6, 2); std::cout << "\n";
	std::cout << "Fixed point found: "; fp.location.Print(std::cout, 10, 6); std::cout << "\n";
	std::cout << "Convergence residual: " << std::scientific << fp.convergenceResidual << "\n";
	std::cout << "Newton iterations: " << fp.iterations << "\n";
	std::cout << "Type: " << ToString(fp.type) << "\n";
	std::cout << "Stable: " << (fp.isStable ? "yes" : "no") << "\n";
	
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
	std::cout << " - " << ToString(fp_origin.type) << "\n";
	
	// C+ fixed point (around sqrt(β(ρ-1)), sqrt(β(ρ-1)), ρ-1)
	Real val = std::sqrt(8.0/3.0 * (28.0 - 1.0));
	auto fp_cplus = FixedPointFinder::Find(lorenz, Vector<Real>({val, val, 27.0}));
	std::cout << "C+ point: "; fp_cplus.location.Print(std::cout, 8, 4);
	std::cout << " - " << ToString(fp_cplus.type) << "\n";
	
	// C- fixed point
	auto fp_cminus = FixedPointFinder::Find(lorenz, Vector<Real>({-val, -val, 27.0}));
	std::cout << "C- point: "; fp_cminus.location.Print(std::cout, 8, 4);
	std::cout << " - " << ToString(fp_cminus.type) << "\n";
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
		          << ": " << std::setw(16) << ToString(fp.type);
		
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
	
	std::cout << "\nNote: For all μ > 0, the origin is an unstable focus.\n";
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
	
	std::cout << "\nLyapunov exponents quantify chaos by measuring exponential divergence.\n";
	std::cout << "- λ₁ > 0: Chaotic (nearby trajectories diverge exponentially)\n";
	std::cout << "- λ₁ = 0: Marginally stable or quasi-periodic\n";
	std::cout << "- λ₁ < 0: Stable (trajectories converge)\n";
	std::cout << "- Sum < 0: Dissipative system (phase space contracts)\n";
	
	// Lorenz system - famous chaotic attractor
	std::cout << "\n--- Lorenz System (σ=10, ρ=28, β=8/3) ---\n";
	LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
	
	std::cout << "Computing Lyapunov exponents (this may take a moment)...\n";
	auto lyap_lorenz = LyapunovAnalyzer::Compute(lorenz,
		Vector<Real>({1.0, 1.0, 1.0}),  // Initial condition
		500.0,   // Total time
		1.0,     // Orthonormalization interval
		0.01     // Time step
	);
	
	std::cout << "\nLyapunov exponents: ";
	lyap_lorenz.exponents.Print(std::cout, 10, 4);
	std::cout << "\n";
	std::cout << "Max exponent (λ₁): " << std::fixed << std::setprecision(4) 
	          << lyap_lorenz.maxExponent << "\n";
	std::cout << "Sum of exponents: " << lyap_lorenz.sum << "\n";
	std::cout << "Kaplan-Yorke dimension: " << lyap_lorenz.kaplanYorkeDimension << "\n";
	std::cout << "Chaotic: " << (lyap_lorenz.isChaotic ? "YES" : "no") << "\n";
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
	std::cout << "Chaotic: " << (lyap_rossler.isChaotic ? "YES" : "no") << "\n";
	std::cout << "Kaplan-Yorke dimension: " << lyap_rossler.kaplanYorkeDimension << "\n";
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
	auto lyap = DiscreteMapLyapunov<1>::compute(logistic4, Vector<Real>({0.3}), 10000, 1000);
	
	std::cout << "r = 4.0: λ = " << std::setprecision(6) << lyap.maxExponent << "\n";
	std::cout << "Analytical: λ = log(2) = " << std::log(2.0) << "\n";
	std::cout << "Chaotic: " << (lyap.isChaotic ? "YES" : "no") << "\n";
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
	auto lyap = DiscreteMapLyapunov<2>::compute(henon, Vector<Real>({0.0, 0.0}), 50000, 1000);
	
	std::cout << "Lyapunov exponents: λ₁ = " << std::setprecision(4) << lyap.exponents[0]
	          << ", λ₂ = " << lyap.exponents[1] << "\n";
	std::cout << "Sum: " << (lyap.exponents[0] + lyap.exponents[1]) << "\n";
	std::cout << "Chaotic: " << (lyap.isChaotic ? "YES" : "no") << "\n";
	
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
	std::cout << "K ≈ 0.97: Last KAM torus destroyed (global chaos onset)\n";
	std::cout << "K > 1: Widespread chaos with islands\n";
	
	std::vector<Real> K_values = {0.5, 0.9716, 2.0};
	
	for (Real K : K_values) {
		StandardMap stdmap(K);
		
		std::cout << "\nK = " << std::fixed << std::setprecision(4) << K << ":\n";
		std::cout << "  Area preserving: " << (stdmap.isAreaPreserving() ? "yes" : "no") << "\n";
		
		// Compute Lyapunov exponent
		auto lyap = DiscreteMapLyapunov<2>::compute(stdmap, 
			Vector<Real>({0.5, 1.0}), 10000, 500);
		
		std::cout << "  Max Lyapunov exponent: " << std::setprecision(4) << lyap.maxExponent << "\n";
		std::cout << "  Chaotic: " << (lyap.isChaotic ? "yes" : "no") << "\n";
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
		auto lyap = DiscreteMapLyapunov<1>::compute(tent, Vector<Real>({0.3}), 10000, 500);
		
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
	
	std::cout << "\nBifurcation diagrams show how attractors change with parameters.\n";
	
	// Simple demonstration with Logistic map
	std::cout << "\n--- Logistic Map Period-Doubling Route to Chaos ---\n";
	std::cout << "\nSampling attractor at different r values:\n";
	
	std::cout << std::setw(6) << "r" << std::setw(15) << "Attractor type" << "  Values\n";
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
		if (unique_vals.size() <= 1) type = "fixed point";
		else if (unique_vals.size() <= 2) type = "period-2";
		else if (unique_vals.size() <= 4) type = "period-4";
		else if (unique_vals.size() <= 8) type = "period-8";
		else type = "chaos";
		
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
		          << " (range: [" << minVal << ", " << maxVal << "])\n";
	}
	
	// Modify parameters
	std::cout << "\nModifying parameters:\n";
	lorenz.setParam(1, 24.0);  // Change ρ below chaos threshold
	std::cout << "  Set ρ = " << lorenz.getParam(1) << " (below chaos threshold ~24.74)\n";
	
	// Get/set all parameters at once
	Vector<Real> params = lorenz.getParams();
	std::cout << "\nAll parameters: "; params.Print(std::cout, 8, 2); std::cout << "\n";
	
	// Restore chaotic parameters
	lorenz.setParams(Vector<Real>({10.0, 28.0, 8.0/3.0}));
	std::cout << "Restored to chaotic regime: σ=" << lorenz.getParam(0) 
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
///                         CHAOS DETECTION WORKFLOW                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DynamicalSystem_ChaosDetection()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Complete Chaos Detection Workflow\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nA system is chaotic if:\n";
	std::cout << "1. Largest Lyapunov exponent > 0\n";
	std::cout << "2. Bounded trajectories (attractor exists)\n";
	std::cout << "3. Sensitive dependence on initial conditions\n";
	
	LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
	
	std::cout << "\n--- Testing Lorenz System ---\n";
	
	// 1. Check fixed points
	std::cout << "\n1. Fixed Point Analysis:\n";
	auto fp = FixedPointFinder::Find(lorenz, Vector<Real>({0.0, 0.0, 0.0}));
	std::cout << "   Origin: " << ToString(fp.type) << " (unstable)\n";
	
	Real val = std::sqrt(8.0/3.0 * 27.0);
	auto fp2 = FixedPointFinder::Find(lorenz, Vector<Real>({val, val, 27.0}));
	std::cout << "   C+ point: " << ToString(fp2.type) << "\n";
	
	// 2. Compute Lyapunov exponents
	std::cout << "\n2. Lyapunov Exponents:\n";
	auto lyap = LyapunovAnalyzer::Compute(lorenz,
		Vector<Real>({1.0, 1.0, 1.0}), 200.0, 1.0, 0.01);
	
	std::cout << "   λ₁ = " << std::fixed << std::setprecision(4) << lyap.exponents[0] 
	          << (lyap.exponents[0] > 0 ? " > 0 ✓" : " ≤ 0") << "\n";
	std::cout << "   λ₂ = " << lyap.exponents[1] << "\n";
	std::cout << "   λ₃ = " << lyap.exponents[2] << "\n";
	
	// 3. Check boundedness (sum < 0 for dissipative)
	std::cout << "\n3. Boundedness Check:\n";
	std::cout << "   Sum of exponents: " << lyap.sum << "\n";
	std::cout << "   Dissipative (sum < 0): " << (lyap.sum < 0 ? "yes ✓" : "no") << "\n";
	
	// 4. Kaplan-Yorke dimension
	std::cout << "\n4. Attractor Dimension:\n";
	std::cout << "   Kaplan-Yorke dimension: " << lyap.kaplanYorkeDimension << "\n";
	std::cout << "   (Fractional dimension indicates strange attractor)\n";
	
	// Verdict
	std::cout << "\n--- Verdict ---\n";
	std::cout << "CHAOTIC: " << (lyap.isChaotic ? "YES ✓" : "NO") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN ENTRY POINT                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

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
