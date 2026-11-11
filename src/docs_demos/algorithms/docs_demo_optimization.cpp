#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/Function.h>
#include <mml/algorithms/Optimization/Optimization.h>
#include <mml/algorithms/Optimization/OptimizationMultidim.h>
#include <mml/algorithms/Optimization/LinearProgramming.h>
#include <mml/algorithms/Optimization/Constraints/AugmentedLagrangian.h>
#include <mml/algorithms/Optimization/Constraints/BoxConstrainedPowell.h>
#endif

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string_view>

using namespace MML;

void Docs_Demo_Optimization_BoundsContracts()
{
	Optimization::BoundConstraints unit(1, 0, 1);
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	bool invalidRejected = false;
	try {
		const auto unexpected = unit.IsFeasible(Vector<Real>{nan});
		(void)unexpected;
	} catch (const Optimization::BoundConstraintsError&) {
		invalidRejected = true;
	}
	if (!invalidRejected)
		throw std::runtime_error("NaN feasibility query was not rejected");

	const auto unbounded = Optimization::BoundConstraints::Unbounded(1);
	if (!unbounded.IsFeasible(VectorN<Real, 1>{2}) || unbounded.MaxViolation(Vector<Real>{2}) != 0)
		throw std::runtime_error("Unbounded finite-coordinate feasibility changed");

	const Real inf = std::numeric_limits<Real>::infinity();
	for (Real invalid : {nan, inf, -inf}) {
		invalidRejected = false;
		try {
			const auto unexpected = unbounded.MaxViolation(VectorN<Real, 1>{invalid});
			(void)unexpected;
		} catch (const Optimization::BoundConstraintsError&) {
			invalidRejected = true;
		}
		if (!invalidRejected)
			throw std::runtime_error("Nonfinite unbounded violation query was not rejected");
	}
	const Optimization::BoundConstraints box(2, 0, 1);
	const auto projected = box.Project(Vector<Real>{-2, 3});
	if (projected[0] != 0 || projected[1] != 1 || !box.IsFeasible(projected))
		throw std::runtime_error("Finite bound projection changed");
	std::cout << "Bounds: nonfinite query coordinates rejected; unbounded markers and projection preserved\n";
}

void Docs_Demo_Optimization_CallbackAccounting()
{
	int scalarCalls = 0;
	auto scalar = [&](Real x) {
		++scalarCalls;
		const Real offset = x - Real(1) / 3;
		return 1 + offset * offset;
	};
	Minimization1DConfig config;
	config.max_iterations = 1;
	const auto minimum = Minimization::BrentMinimize(scalar, -1, 0, config);
	if (minimum.function_evaluations != scalarCalls || scalarCalls != 5 || minimum.iterations != 1)
		throw std::runtime_error("Brent bracketing/refinement accounting failed");
	std::cout << "Brent: " << minimum.iterations << " iteration, " << scalarCalls << " objective calls\n";

	scalarCalls = 0;
	const auto golden = Minimization::GoldenSectionSearch(scalar, -1, 0, config);
	if (golden.function_evaluations != scalarCalls || scalarCalls != 6)
		throw std::runtime_error("Golden-section accounting failed");
	int derivativeCalls = 0;
	auto derivative = [&](Real x) { ++derivativeCalls; return 2 * (x - Real(1) / 3); };
	scalarCalls = 0;
	const auto differentiable = Minimization::BrentMinimizeWithDeriv(scalar, derivative, -1, 0, config);
	if (differentiable.function_evaluations != scalarCalls || scalarCalls != 5 || derivativeCalls != 2)
		throw std::runtime_error("Derivative Brent objective/derivative accounting failed");
	std::cout << "Derivative Brent: " << scalarCalls << " objective calls, "
	          << derivativeCalls << " separate derivative calls\n";

	int failedCalls = 0;
	auto decreasing = [&](Real x) { ++failedCalls; return -std::asinh(x); };
	const auto unbracketed = Minimization::BrentMinimize(decreasing, 0, 1);
	if (unbracketed.converged || unbracketed.status != AlgorithmStatus::AlgorithmSpecificFailure
		|| unbracketed.function_evaluations != failedCalls || failedCalls == 0)
		throw std::runtime_error("Failed-bracketing accounting failed");
	std::cout << "Failed bracket: " << failedCalls << " objective calls retained\n";

	scalarCalls = 0;
	Minimization1DConfig successConfig;
	const auto successful = Minimization::BrentMinimize(scalar, -1, 0, successConfig);
	if (!successful.converged || successful.status != AlgorithmStatus::Success
		|| !successful.IsSuccess() || !static_cast<bool>(successful)
		|| successful.function_evaluations != scalarCalls)
		throw std::runtime_error("Brent success-status contract failed");
	scalarCalls = derivativeCalls = 0;
	const auto successfulDerivative = Minimization::BrentMinimizeWithDeriv(scalar, derivative, -1, 0, successConfig);
	if (!successfulDerivative.converged || successfulDerivative.status != AlgorithmStatus::Success
		|| !successfulDerivative.IsSuccess() || !static_cast<bool>(successfulDerivative)
		|| successfulDerivative.function_evaluations != scalarCalls || derivativeCalls == 0)
		throw std::runtime_error("Derivative Brent success-status contract failed");
	if (minimum.status != AlgorithmStatus::MaxIterationsExceeded || minimum.IsSuccess()
		|| static_cast<bool>(minimum) || unbracketed.IsSuccess() || static_cast<bool>(unbracketed))
		throw std::runtime_error("Brent failure-status contract failed");
	std::cout << "Brent and derivative Brent: converged results report usable Success\n";

	int objectiveCalls = 0;
	ScalarFunctionFromStdFunc<2> objective([&](const VectorN<Real, 2>& x) {
		++objectiveCalls;
		return 1 + x[0] * x[0] + 2 * x[1] * x[1];
	});
	Optimization::Powell optimizer(Precision::OptimizationTolerance, 1);
	const auto powell = optimizer.Minimize(objective, VectorN<Real, 2>{-5, 5});
	if (powell.function_evaluations != objectiveCalls || optimizer.getNumFuncEvals() != objectiveCalls
		|| powell.iterations != 1 || powell.converged)
		throw std::runtime_error("Powell nested-callback accounting failed");
	std::cout << "Powell: " << powell.iterations << " iteration, " << objectiveCalls << " objective calls\n";

	objectiveCalls = 0;
	ScalarFunctionFromStdFunc<2> failed([&](const VectorN<Real, 2>& x) {
		if (++objectiveCalls == 2) throw std::logic_error("demo Powell callback failure");
		return 1 + x[0] * x[0] + 2 * x[1] * x[1];
	});
	bool propagated = false;
	try {
		const auto unexpected = optimizer.Minimize(failed, VectorN<Real, 2>{-5, 5});
		(void)unexpected;
	} catch (const std::logic_error& error) {
		if (std::string_view(error.what()) != "demo Powell callback failure") throw;
		propagated = true;
	}
	if (!propagated || objectiveCalls != 2 || optimizer.getNumFuncEvals() != objectiveCalls)
		throw std::runtime_error("Powell throwing-callback accounting failed");
	std::cout << "Powell propagated failure: " << optimizer.getNumFuncEvals() << " callback attempts\n";

	objectiveCalls = 0;
	Optimization::BoundConstraints bounds(Vector<Real>{-2, -2}, Vector<Real>{2, 2});
	Optimization::MultidimOptimizationConfig boundedConfig;
	boundedConfig.max_iterations = 1;
	boundedConfig.tolerance = Precision::OptimizationTolerance;
	const auto bounded = Optimization::BoxConstrainedPowellMinimize(objective, VectorN<Real, 2>{1, 1}, bounds, boundedConfig);
	if (bounded.function_evaluations != objectiveCalls || bounded.iterations != 1)
		throw std::runtime_error("Bounded Powell final-check accounting failed");
	std::cout << "Bounded Powell: " << objectiveCalls << " objective calls including final projected check\n";

	objectiveCalls = 0;
	Matrix<Real> wrongDirections(3, 3, Real(0));
	for (int i = 0; i < 3; ++i) wrongDirections(i, i) = 1;
	const auto before = wrongDirections;
	bool rejected = false;
	try {
		const auto unexpected = optimizer.Minimize(objective, VectorN<Real, 2>{-5, 5}, wrongDirections);
		(void)unexpected;
	} catch (const Optimization::MultidimOptimizationInputError& error) {
		if (std::string_view(error.what()).find("expected 2 x 2, got 3 x 3") == std::string_view::npos) throw;
		rejected = true;
	}
	if (!rejected || objectiveCalls != 0 || optimizer.getNumFuncEvals() != 0 || wrongDirections != before)
		throw std::runtime_error("Powell direction-shape preflight failed");
	std::cout << "Powell: wrong direction dimensions rejected before callbacks or mutation\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    BRACKETING A MINIMUM                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Bracketing()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Bracketing a Minimum\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nBefore minimizing, we need to bracket the minimum: find (a,b,c)\n";
	std::cout << "such that f(b) < f(a) and f(b) < f(c), guaranteeing min in [a,c].\n";
	
	// Simple parabola with minimum at x=2
	RealFunctionFromStdFunc func([](Real x) { return (x - 2.0) * (x - 2.0) + 1.0; });
	
	std::cout << "\n--- f(x) = (x-2)² + 1  (minimum at x=2) ---\n";
	
	// Find bracket starting from [0, 1]
	MinimumBracket bracket = Minimization::BracketMinimum(func, 0.0, 1.0);
	
	std::cout << "\nBracket found:\n";
	std::cout << "  ax = " << bracket.ax << ", f(ax) = " << bracket.fa << std::endl;
	std::cout << "  bx = " << bracket.bx << ", f(bx) = " << bracket.fb << std::endl;
	std::cout << "  cx = " << bracket.cx << ", f(cx) = " << bracket.fc << std::endl;
	std::cout << "  Valid: " << (bracket.valid ? "yes" : "no") << std::endl;
	
	std::cout << "\nNote: bx is the lowest point, guaranteeing min in [ax, cx]\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    GOLDEN SECTION SEARCH                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_GoldenSection()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Golden Section Search\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nGolden section: simple, robust, linear convergence rate.\n";
	std::cout << "Uses golden ratio φ = (1+√5)/2 ≈ 1.618 for optimal bracketing.\n";
	
	// Function with minimum at x = 3
	RealFunctionFromStdFunc func([](Real x) { 
		return (x - 3.0) * (x - 3.0) + 2.0; 
	});
	
	std::cout << "\n--- f(x) = (x-3)² + 2  (minimum at x=3, f=2) ---\n";
	
	// Bracket first
	MinimumBracket bracket = Minimization::BracketMinimum(func, 0.0, 1.0);
	
	// Golden section search
	Minimization1DResult result = Minimization::GoldenSectionSearch(func, bracket);

	std::cout << "\nGolden Section Result:\n";
	std::cout << "  x_min = " << result.xmin << " (exact: 3.0)\n";
	std::cout << "  f_min = " << result.fmin << " (exact: 2.0)\n";
	std::cout << "  Iterations: " << result.iterations << std::endl;
	std::cout << "  Converged: " << (result.converged ? "yes" : "no") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    BRENT'S METHOD                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Brent()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Brent's Method (Parabolic Interpolation)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nBrent's method: combines golden section with parabolic interpolation.\n";
	std::cout << "Superlinear convergence for smooth functions, never worse than golden.\n";
	
	// More complex function
	RealFunctionFromStdFunc func([](Real x) { 
		return std::sin(x) + 0.1 * x * x - 2.0;  // min around x ≈ -0.45
	});
	
	std::cout << "\n--- f(x) = sin(x) + 0.1*x² - 2 ---\n";
	
	// Bracket in [-2, 2]
	MinimumBracket bracket = Minimization::BracketMinimum(func, -2.0, 2.0);
	
	// Brent's method
	Minimization1DResult result = Minimization::BrentMinimize(func, bracket);

	std::cout << "\nBrent's Method Result:\n";
	std::cout << "  x_min = " << result.xmin << std::endl;
	std::cout << "  f_min = " << result.fmin << std::endl;
	std::cout << "  Iterations: " << result.iterations << std::endl;
	
	// Verify
	std::cout << "\nVerification: f'(x_min) ≈ " 
	          << (func(result.xmin + 1e-6) - func(result.xmin - 1e-6)) / 2e-6 
	          << " (should be ~0)\n";
}

void Docs_Demo_Optimization_Brent_WithDeriv()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Brent's Method with Derivatives\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nWhen derivatives are available, convergence is faster.\n";
	
	// Function and its derivative
	RealFunctionFromStdFunc func([](Real x) { 
		return x * x * x * x - 3.0 * x * x + 2.0;  // x⁴ - 3x² + 2
	});
	RealFunctionFromStdFunc dfunc([](Real x) { 
		return 4.0 * x * x * x - 6.0 * x;  // 4x³ - 6x
	});
	
	std::cout << "\n--- f(x) = x⁴ - 3x² + 2, f'(x) = 4x³ - 6x ---\n";
	std::cout << "Minima at x = ±√(3/2) ≈ ±1.225\n";
	
	// Bracket around positive minimum
	MinimumBracket bracket = Minimization::BracketMinimum(func, 0.5, 2.0);
	
	// Brent with derivatives
	Minimization1DResult result = Minimization::BrentMinimizeWithDeriv(func, dfunc, bracket);

	std::cout << "\nBrent with Derivatives Result:\n";
	std::cout << "  x_min = " << result.xmin << " (exact: " << std::sqrt(1.5) << ")\n";
	std::cout << "  f_min = " << result.fmin << std::endl;
	std::cout << "  f'(x_min) = " << dfunc(result.xmin) << " (should be ~0)\n";
	std::cout << "  Iterations: " << result.iterations << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    1D MINIMIZATION WITH BOUNDS                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Bounded()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Bounded Minimization\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFind minimum within specified bounds [a, b].\n";
	
	// Function with global minimum outside bounds
	RealFunctionFromStdFunc func([](Real x) { 
		return (x - 5.0) * (x - 5.0);  // min at x=5
	});
	
	std::cout << "\n--- f(x) = (x-5)²  (global min at x=5) ---\n";
	std::cout << "Searching in [0, 3] - minimum is outside this range!\n";
	
	// Minimize in [0, 3]
	MinimumBracket bracket;
	bracket.ax = 0.0;
	bracket.cx = 3.0;
	bracket.bx = 1.5;
	bracket.fa = func(bracket.ax);
	bracket.fb = func(bracket.bx);
	bracket.fc = func(bracket.cx);
	bracket.valid = true;

	Minimization1DResult result = Minimization::GoldenSectionSearch(func, bracket, 1e-8);

	std::cout << "\nConstrained Result:\n";
	std::cout << "  x_min = " << result.xmin << " (at right boundary ≈ 3)\n";
	std::cout << "  f_min = " << result.fmin << " (compared to global min = 0 at x=5)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    COMPARING METHODS                                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Comparison()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Method Comparison\n";
	std::cout << "==========================================================================\n";
	
	// Rosenbrock-like 1D function (steep valley)
	RealFunctionFromStdFunc func([](Real x) { 
		return 100.0 * std::pow(1.0 - x, 2) + std::pow(x - 1.0, 4);
	});
	
	std::cout << "\n--- Rosenbrock-like: 100*(1-x)² + (x-1)⁴ ---\n";
	std::cout << "Minimum at x = 1\n";
	
	MinimumBracket bracket = Minimization::BracketMinimum(func, -5.0, 5.0);
	
	// Golden section
	Minimization1DResult golden = Minimization::GoldenSectionSearch(func, bracket);

	// Brent
	Minimization1DResult brent = Minimization::BrentMinimize(func, bracket);

	std::cout << "\nComparison:\n";
	std::cout << "                Golden Section     Brent\n";
	std::cout << "  x_min:        " << golden.xmin << "              " << brent.xmin << std::endl;
	std::cout << "  f_min:        " << golden.fmin << "           " << brent.fmin << std::endl;
	std::cout << "  Iterations:   " << golden.iterations << "                    " << brent.iterations << std::endl;
	
	std::cout << "\nBrent typically needs fewer iterations due to parabolic acceleration.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    MULTIDIMENSIONAL OPTIMIZATION (PREVIEW)                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Multidim()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Multidimensional Optimization Overview\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFor N-dimensional minimization, see Optimization/OptimizationMultidim.h:\n";
	std::cout << "  - Coordinate descent (simple, dimension-by-dimension)\n";
	std::cout << "  - Powell's method (conjugate directions)\n";
	std::cout << "  - Nelder-Mead (simplex, derivative-free)\n";
	std::cout << "  - Gradient descent (with line search)\n";
	std::cout << "  - BFGS (quasi-Newton, fast convergence)\n";
	
	std::cout << "\nExample setup:\n";
	std::cout << "  ScalarFunction<2> f([](const VectorN<Real,2>& x) {\n";
	std::cout << "      return (x[0]-1)*(x[0]-1) + (x[1]-2)*(x[1]-2);\n";
	std::cout << "  });\n";
	std::cout << "  VectorN<Real,2> x0 = {0, 0};\n";
	std::cout << "  auto result = MultidimMinimization::NelderMead(f, x0);\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    LINEAR PROGRAMMING                                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_LinearProgramming()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Dense Linear Programming\n";
	std::cout << "==========================================================================\n";

	std::cout << "\nMaximize 3x + 2y subject to x + y <= 4, 2x + y <= 5, x,y >= 0.\n";

	Optimization::LinearProgram lp(2, "ProductionPlan");
	lp.SetVariableNames({ "x", "y" });
	lp.SetObjective({ static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0) }, Optimization::LPObjective::Maximize);
	lp.AddConstraint({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0) }, Optimization::LPConstraintType::LessEqual, static_cast<MML::Real>(4.0), "capacity1");
	lp.AddConstraint({ static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0) }, Optimization::LPConstraintType::LessEqual, static_cast<MML::Real>(5.0), "capacity2");

	Optimization::LinearProgrammingResult result = Optimization::SolveLP(lp);
	std::cout << "\nStatus: " << result.statusMessage() << std::endl;
	if (result.IsOptimal()) {
		std::cout << "  x = " << result.x[0] << std::endl;
		std::cout << "  y = " << result.x[1] << std::endl;
		std::cout << "  objective = " << result.objectiveValue << std::endl;
		std::cout << "  iterations = " << result.iterations << std::endl;
	}
}

void Docs_Demo_LP_Sensitivity()
{
	using namespace Optimization;
	LinearProgram lp(1);
	lp.SetObjective({2}, LPObjective::Maximize);
	lp.AddConstraint({1}, LPConstraintType::LessEqual, 3);
	lp.AddConstraint({1}, LPConstraintType::LessEqual, 5);
	const auto result = SolveLPWithSensitivity(lp);
	if (!result.IsOptimal() || result.rhsRanges[1].first != 3 ||
		result.rhsRanges[1].second != std::numeric_limits<Real>::infinity() ||
		result.dualValues[0] != 2 || result.dualValues[1] != 0)
		throw LinearProgrammingError("Sensitivity demo: wrong original RHS range or shadow price");
	std::cout << "Second RHS range: [" << result.rhsRanges[1].first << ", "
		<< result.rhsRanges[1].second << "]; shadow prices: "
		<< result.dualValues[0] << ", " << result.dualValues[1] << '\n';

	LinearProgram scaled(2);
	scaled.SetObjective({1, 1}, LPObjective::Maximize);
	scaled.AddConstraint({1, 0}, LPConstraintType::LessEqual, 2);
	scaled.AddConstraint({0, Real(1e-6)}, LPConstraintType::LessEqual, Real(3e-6));
	const auto conditioned = SolveLPWithSensitivity(scaled);
	if (!conditioned.IsOptimal() ||
		std::abs(conditioned.basisConditionNumber / Real(1e6) - 1) > Precision::DefaultToleranceStrict)
		throw LinearProgrammingError("Sensitivity demo: original basis condition was lost during pivoting");
	std::cout << "Original diag(1, 1e-6) basis condition: " << conditioned.basisConditionNumber << '\n';
}

void Docs_Demo_LP_Duality()
{
	using namespace Optimization;
	LinearProgram primal(1);
	primal.SetObjective({-1});
	primal.AddConstraint({1}, LPConstraintType::Equal, 2);
	LinearProgram dual = primal.Dual();
	LinearProgrammingResult result = SolveLP(dual);
	if (!result.IsOptimal())
		throw LinearProgrammingError("Duality demo: equality dual did not solve");
	Real multiplier = result.x[0] - result.x[1];
	std::cout << result.objectiveValue << " " << multiplier << '\n';
	if (std::abs(result.objectiveValue + 2) > Precision::DefaultToleranceStrict ||
		std::abs(multiplier + 1) > Precision::DefaultToleranceStrict)
		throw LinearProgrammingError("Duality demo: wrong equality multiplier or objective");

	for (const auto type : {LPConstraintType::LessEqual, LPConstraintType::GreaterEqual}) {
		LinearProgram inequality(1);
		inequality.SetObjective({1});
		inequality.AddConstraint({1}, type, 1);
		const auto p = SolveLP(inequality);
		const auto d = SolveLP(inequality.Dual());
		const Real expected = type == LPConstraintType::LessEqual ? Real(0) : Real(1);
		if (!p.IsOptimal() || !d.IsOptimal() ||
			std::abs(p.objectiveValue - expected) > Precision::DefaultToleranceStrict ||
			std::abs(d.objectiveValue - expected) > Precision::DefaultToleranceStrict)
			throw LinearProgrammingError("Duality demo: inequality primal/dual mismatch");
	}
	std::cout << "LP dual row signs, free multiplier coordinates and strong duality verified.\n";
}

void Docs_Demo_LP_ModelContracts()
{
	using namespace Optimization;
	LinearProgram invalid(1);
	invalid.SetObjective({1});
	invalid.AddConstraint({1}, LPConstraintType::LessEqual, std::numeric_limits<Real>::quiet_NaN());
	std::string error;
	if (invalid.Validate(error) || error.empty())
		throw LinearProgrammingError("Model demo: NaN RHS passed validation");
	std::cout << error << '\n';
	bool rejected = false;
	try {
		(void)SolveLP(invalid);
	} catch (const LinearProgrammingError&) {
		rejected = true;
	}
	if (!rejected)
		throw LinearProgrammingError("Model demo: solver accepted NaN RHS");

	LinearProgram valid(1);
	valid.SetObjective({1});
	valid.AddConstraint({-1}, LPConstraintType::LessEqual, -1);
	if (!valid.Validate(error) || !error.empty())
		throw LinearProgrammingError("Model demo: valid negative RHS failed validation");
	const auto result = SolveLP(valid);
	if (!result.IsOptimal() || std::abs(result.x[0] - 1) > Precision::DefaultToleranceStrict ||
		std::abs(result.objectiveValue - 1) > Precision::DefaultToleranceStrict)
		throw LinearProgrammingError("Model demo: normalized negative RHS gave a wrong optimum");
	std::cout << "LP finite-model validation and negative-RHS conversion verified.\n";
}

void Docs_Demo_ALM_Contracts()
{
	using namespace Optimization;
	ALMConfig config;
	config.max_outer_iterations = 4;
	config.max_inner_iterations = 256;
	config.inner_simplex_delta = Real(3) / 10;
	AugmentedLagrangian<2> solver(config);
	solver.SetBounds(2, 3);
	auto result = solver.Optimize(
		[](const VectorN<Real, 2>& x) { return 1 + x[0] * x[0] + x[1] * x[1]; },
		VectorN<Real, 2>{Real(5) / 2, Real(12) / 5});
	std::cout << MML::to_string(result.status) << " " << result.objective_value << '\n';
	if (!result.converged || result.x[0] < 2 || result.x[0] > 3 || result.x[1] < 2 || result.x[1] > 3 ||
		std::abs(result.objective_value - 9) > Precision::LooseTolerance)
		throw MultidimOptimizationError("ALM demo: wrong box-constrained result");

	config.max_inner_iterations = 2;
	AugmentedLagrangian<1> exhausted(config);
	int calls = 0;
	auto objective = [&](const VectorN<Real, 1>& x) { ++calls; return x[0] * x[0]; };
	const auto partial = exhausted.Optimize(objective, VectorN<Real, 1>{1});
	if (partial.IsSuccess() || partial.converged || partial.outer_iterations != 1 ||
		partial.status != AlgorithmStatus::MaxIterationsExceeded || partial.function_evaluations != calls ||
		!std::isnan(partial.stationarity_measure))
		throw MultidimOptimizationError("ALM demo: inner failure or objective count was lost");
	std::cout << "ALM bounds, inner failure and objective-call accounting verified.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization()
{
	std::cout << "\n##########################################################################\n";
	std::cout << "#                     OPTIMIZATION DEMOS                                 #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_Optimization_Bracketing();
	Docs_Demo_Optimization_GoldenSection();
	Docs_Demo_Optimization_Brent();
	Docs_Demo_Optimization_Brent_WithDeriv();
	Docs_Demo_Optimization_Bounded();
	Docs_Demo_Optimization_Comparison();
	Docs_Demo_Optimization_Multidim();
	Docs_Demo_Optimization_LinearProgramming();
	Docs_Demo_LP_Duality();
	Docs_Demo_LP_ModelContracts();
	Docs_Demo_ALM_Contracts();
	
	std::cout << "\n=== All Optimization Demos Complete ===\n";
}
