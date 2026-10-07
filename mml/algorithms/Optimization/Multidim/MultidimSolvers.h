///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Multidim/MultidimSolvers.h                                           ///
///  Description: Multidimensional optimization convenience functions                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_MULTIDIM_SOLVERS_H
#define MML_OPTIMIZATION_MULTIDIM_SOLVERS_H

#include <mml/algorithms/Optimization/Multidim/NelderMead.h>
#include <mml/algorithms/Optimization/Multidim/Powell.h>
#include <mml/algorithms/Optimization/Multidim/QuasiNewton.h>

namespace MML::Optimization {


	/**
     * @brief Minimize using Powell's method (no derivatives needed)
     */
	template<int N>
	MultidimMinimizationResult PowellMinimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start, Real ftol = PrecisionValues<Real>::OptimizationTolerance) {
		Powell optimizer(ftol);
		return optimizer.Minimize(func, start);
	}

	/**
     * @brief Minimize using conjugate gradient (requires derivatives)
     */
	template<int N>
	MultidimMinimizationResult ConjugateGradientMinimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start,
														 Real ftol = PrecisionValues<Real>::OptimizationTolerance,
														 ConjugateGradient::Method method = ConjugateGradient::Method::PolakRibiere) {
		ConjugateGradient optimizer(ftol, PrecisionValues<Real>::OptimizationGradientTolerance, 200, method);
		return optimizer.Minimize(func, start);
	}

	/**
     * @brief Minimize using BFGS quasi-Newton method (requires derivatives)
     */
	template<int N>
	MultidimMinimizationResult BFGSMinimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start,
											Real ftol = PrecisionValues<Real>::OptimizationTolerance) {
		BFGS optimizer(ftol);
		return optimizer.Minimize(func, start);
	}

	/**
     * @brief Minimize using L-BFGS limited-memory quasi-Newton (requires derivatives)
     * 
     * Preferred over BFGS for large-scale problems (N > 100) due to O(m×N) memory
     * instead of O(N²). Uses only the last m correction pairs to approximate
     * the inverse Hessian.
     * 
     * @tparam N Dimension of the problem
     * @param func Differentiable scalar function to minimize
     * @param start Starting point
     * @param ftol Function value tolerance (default: 3e-8)
     * @param memorySize Number of correction pairs to store (default: 10)
     * @return Optimization result
     */
	template<int N>
	MultidimMinimizationResult LBFGSMinimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start,
											 Real ftol = PrecisionValues<Real>::OptimizationTolerance, int memorySize = 10) {
		LBFGS optimizer(ftol, PrecisionValues<Real>::OptimizationGradientTolerance, 1000, memorySize);
		return optimizer.Minimize(func, start);
	}

	///////////////////////////////////////////////////////////////////////////
	///      Config-Based Overloads for Powell/CG/BFGS (API Standardization) ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Minimize using Powell's method with config
     */
	template<int N>
	MultidimMinimizationResult PowellMinimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start,
											  const MultidimOptimizationConfig& config) {
		AlgorithmTimer timer;
		Powell optimizer(config.tolerance, config.max_iterations);
		MultidimMinimizationResult result = optimizer.Minimize(func, start);
		
		result.algorithm_name = "Powell";
		result.elapsed_time_ms = timer.elapsed_ms();
		if (!result.converged) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "Powell method did not converge within " + std::to_string(config.max_iterations) + " iterations";
		}
		return result;
	}

	/**
     * @brief Minimize using conjugate gradient with config
     */
	template<int N>
	MultidimMinimizationResult ConjugateGradientMinimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start,
														 const MultidimOptimizationConfig& config,
														 ConjugateGradient::Method method = ConjugateGradient::Method::PolakRibiere) {
		AlgorithmTimer timer;
		ConjugateGradient optimizer(config.tolerance, PrecisionValues<Real>::OptimizationGradientTolerance, config.max_iterations, method);
		MultidimMinimizationResult result = optimizer.Minimize(func, start);
		
		result.algorithm_name = "ConjugateGradient";
		result.elapsed_time_ms = timer.elapsed_ms();
		if (!result.converged) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "Conjugate gradient did not converge within " + std::to_string(config.max_iterations) + " iterations";
		}
		return result;
	}

	/**
     * @brief Minimize using BFGS with config
     */
	template<int N>
	MultidimMinimizationResult BFGSMinimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start,
											const MultidimOptimizationConfig& config) {
		AlgorithmTimer timer;
		BFGS optimizer(config.tolerance, PrecisionValues<Real>::OptimizationGradientTolerance, config.max_iterations);
		MultidimMinimizationResult result = optimizer.Minimize(func, start);
		
		result.algorithm_name = "BFGS";
		result.elapsed_time_ms = timer.elapsed_ms();
		if (!result.converged) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "BFGS did not converge within " + std::to_string(config.max_iterations) + " iterations";
		}
		return result;
	}

	/**
     * @brief Minimize using L-BFGS with config
     * 
     * For large-scale problems where BFGS memory requirements are prohibitive.
     * Memory usage: O(m×N) where m = config.lbfgs_memory_size (default 10).
     */
	template<int N>
	MultidimMinimizationResult LBFGSMinimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start,
											 const MultidimOptimizationConfig& config) {
		AlgorithmTimer timer;
		LBFGS optimizer(config.tolerance, PrecisionValues<Real>::OptimizationGradientTolerance, config.max_iterations, config.lbfgs_memory_size);
		MultidimMinimizationResult result = optimizer.Minimize(func, start);
		
		result.algorithm_name = "L-BFGS";
		result.elapsed_time_ms = timer.elapsed_ms();
		if (!result.converged) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "L-BFGS did not converge within " + std::to_string(config.max_iterations) + " iterations";
		}
		return result;
	}


} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_MULTIDIM_SOLVERS_H
