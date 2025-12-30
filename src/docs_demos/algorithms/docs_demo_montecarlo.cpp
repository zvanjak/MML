#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/Function.h"
#include "interfaces/IFunction.h"
#include "algorithms/MonteCarloIntegration.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                    BASIC MONTE CARLO INTEGRATION                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MonteCarlo_Basic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Basic Monte Carlo Integration\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nMonte Carlo: estimate ∫f(x)dx by random sampling.\n";
	std::cout << "Error = O(1/√N) - independent of dimension!\n";
	
	// 1D integral: ∫₀¹ x² dx = 1/3
	std::cout << "\n--- 1D: ∫₀¹ x² dx = 1/3 ≈ 0.333... ---\n";
	
	// Create scalar function using lambda wrapper
	ScalarFunctionFromStdFunc<1> f1([](const VectorN<Real,1>& x) {
		return x[0] * x[0];
	});
	
	VectorN<Real,1> lower1 = {0.0};
	VectorN<Real,1> upper1 = {1.0};
	
	MonteCarloIntegrator<1> mc1(42);  // fixed seed for reproducibility
	
	MonteCarloConfig config;
	config.num_samples = 10000;
	
	auto result1 = mc1.integrate(f1, lower1, upper1, config);
	
	std::cout << "Estimated: " << result1.value << std::endl;
	std::cout << "Error estimate: " << result1.error_estimate << std::endl;
	std::cout << "Exact: " << 1.0/3.0 << std::endl;
	std::cout << "Samples used: " << result1.samples_used << std::endl;
}

void Docs_Demo_MonteCarlo_2D()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: 2D Monte Carlo Integration\n";
	std::cout << "==========================================================================\n";
	
	// Area of unit circle: ∫∫_{x²+y²≤1} 1 dA = π
	std::cout << "\n--- Area of unit circle: ∫∫_{x²+y²≤1} 1 dA = π ---\n";
	
	ScalarFunctionFromStdFunc<2> circle([](const VectorN<Real,2>& x) {
		return (x[0]*x[0] + x[1]*x[1] <= 1.0) ? 1.0 : 0.0;
	});
	
	VectorN<Real,2> lower = {-1.0, -1.0};
	VectorN<Real,2> upper = {1.0, 1.0};
	
	MonteCarloIntegrator<2> mc2(42);
	
	MonteCarloConfig config;
	config.num_samples = 100000;
	
	auto result = mc2.integrate(circle, lower, upper, config);
	
	std::cout << "Estimated π: " << result.value << std::endl;
	std::cout << "Actual π:    " << Constants::PI << std::endl;
	std::cout << "Error: " << std::abs(result.value - Constants::PI) << std::endl;
	std::cout << "Statistical error estimate: " << result.error_estimate << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    HIGH-DIMENSIONAL INTEGRATION                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MonteCarlo_HighDim()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: High-Dimensional Integration\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nMonte Carlo shines in high dimensions where grid methods fail.\n";
	std::cout << "Grid with 10 pts/dim in 10D needs 10^10 points!\n";
	
	// Volume of 5D unit hypersphere: V₅ = 8π²/15 ≈ 5.264
	std::cout << "\n--- Volume of 5D unit hypersphere: V₅ = 8π²/15 ≈ 5.264 ---\n";
	
	ScalarFunctionFromStdFunc<5> hypersphere([](const VectorN<Real,5>& x) {
		Real sum = 0;
		for (int i = 0; i < 5; i++) sum += x[i] * x[i];
		return (sum <= 1.0) ? 1.0 : 0.0;
	});
	
	VectorN<Real,5> lower = {-1, -1, -1, -1, -1};
	VectorN<Real,5> upper = {1, 1, 1, 1, 1};
	
	MonteCarloIntegrator<5> mc5(42);
	
	MonteCarloConfig config;
	config.num_samples = 500000;
	
	auto result = mc5.integrate(hypersphere, lower, upper, config);
	
	Real exact = 8.0 * Constants::PI * Constants::PI / 15.0;
	
	std::cout << "Estimated: " << result.value << std::endl;
	std::cout << "Exact:     " << exact << std::endl;
	std::cout << "Relative error: " << std::abs(result.value - exact)/exact * 100 << "%" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    VARIANCE REDUCTION: ANTITHETIC VARIATES                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MonteCarlo_Antithetic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Antithetic Variates\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nAntithetic variates: use x and (1-x) together to reduce variance.\n";
	std::cout << "Works when f(x) and f(1-x) are negatively correlated.\n";
	
	// Test on ∫₀¹ exp(x) dx = e - 1 ≈ 1.718
	ScalarFunctionFromStdFunc<1> f([](const VectorN<Real,1>& x) {
		return std::exp(x[0]);
	});
	
	VectorN<Real,1> lower = {0.0};
	VectorN<Real,1> upper = {1.0};
	
	MonteCarloIntegrator<1> mc(42);
	
	// Without antithetic
	MonteCarloConfig config1;
	config1.num_samples = 10000;
	config1.use_antithetic = false;
	auto result1 = mc.integrate(f, lower, upper, config1);
	
	// With antithetic
	MonteCarloConfig config2;
	config2.num_samples = 10000;
	config2.use_antithetic = true;
	auto result2 = mc.integrate(f, lower, upper, config2);
	
	Real exact = std::exp(1.0) - 1.0;
	
	std::cout << "\n--- ∫₀¹ exp(x) dx = e-1 ≈ 1.718 ---\n";
	std::cout << "\nWithout antithetic:\n";
	std::cout << "  Value: " << result1.value << ", Error: " << result1.error_estimate << std::endl;
	
	std::cout << "\nWith antithetic:\n";
	std::cout << "  Value: " << result2.value << ", Error: " << result2.error_estimate << std::endl;
	
	std::cout << "\nExact: " << exact << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    CONVERGENCE STUDY                                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MonteCarlo_Convergence()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Convergence Study\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nMonte Carlo error decreases as O(1/√N).\n";
	
	// Simple integral: ∫₀¹ x dx = 0.5
	ScalarFunctionFromStdFunc<1> f([](const VectorN<Real,1>& x) {
		return x[0];
	});
	
	VectorN<Real,1> lower = {0.0};
	VectorN<Real,1> upper = {1.0};
	
	MonteCarloIntegrator<1> mc(42);
	
	std::cout << "\n--- ∫₀¹ x dx = 0.5 ---\n";
	std::cout << "\n  N          Estimate     Error       Error*√N\n";
	
	for (int logN = 2; logN <= 6; logN++) {
		int N = static_cast<int>(std::pow(10, logN));
		MonteCarloConfig config;
		config.num_samples = N;
		
		auto result = mc.integrate(f, lower, upper, config);
		Real error = std::abs(result.value - 0.5);
		
		std::cout << "  " << N << "        " << result.value 
		          << "      " << error 
		          << "      " << error * std::sqrt(static_cast<Real>(N)) << std::endl;
	}
	
	std::cout << "\nNote: Error*√N should be roughly constant (confirming O(1/√N)).\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    PRACTICAL EXAMPLE: GAUSSIAN INTEGRAL                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MonteCarlo_Gaussian()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Gaussian Integral\n";
	std::cout << "==========================================================================\n";
	
	// ∫∫ exp(-(x²+y²)) dxdy over [-3,3]² ≈ π (since full integral is π)
	ScalarFunctionFromStdFunc<2> gaussian([](const VectorN<Real,2>& x) {
		return std::exp(-(x[0]*x[0] + x[1]*x[1]));
	});
	
	VectorN<Real,2> lower = {-3.0, -3.0};
	VectorN<Real,2> upper = {3.0, 3.0};
	
	MonteCarloIntegrator<2> mc(42);
	
	MonteCarloConfig config;
	config.num_samples = 100000;
	
	auto result = mc.integrate(gaussian, lower, upper, config);
	
	std::cout << "\n--- ∫∫_{[-3,3]²} exp(-(x²+y²)) dxdy ≈ π ---\n";
	std::cout << "Estimated: " << result.value << std::endl;
	std::cout << "Exact (full plane): " << Constants::PI << std::endl;
	std::cout << "Error estimate: " << result.error_estimate << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MonteCarlo()
{
	std::cout << "\n##########################################################################\n";
	std::cout << "#                 MONTE CARLO INTEGRATION DEMOS                          #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_MonteCarlo_Basic();
	Docs_Demo_MonteCarlo_2D();
	Docs_Demo_MonteCarlo_HighDim();
	Docs_Demo_MonteCarlo_Antithetic();
	Docs_Demo_MonteCarlo_Convergence();
	Docs_Demo_MonteCarlo_Gaussian();
	
	std::cout << "\n=== All Monte Carlo Demos Complete ===\n";
}
