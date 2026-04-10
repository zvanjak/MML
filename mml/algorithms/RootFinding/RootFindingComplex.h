///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFindingComplex.h                                                ///
///  Description: Root-finding algorithms for complex functions f:C->C                ///
///               Newton's method in C, Muller's method                               ///
///                                                                                   ///
///  REFERENCES:                                                                      ///
///    [NR3]    Press et al., Numerical Recipes 3rd ed., Ch. 9                       ///
///    [M56]    Muller, D.E. (1956). "A Method for Solving Algebraic Equations        ///
///             Using an Automatic Computer." Mathematical Tables and Other Aids       ///
///             to Computation, 10(56), pp. 208-215                                   ///
///    [AH74]   Ahlfors, L.V. (1979). Complex Analysis, 3rd ed. McGraw-Hill          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_COMPLEX_H
#define MML_ROOTFINDING_COMPLEX_H

#include "mml/MMLBase.h"
#include "mml/core/AlgorithmTypes.h"
#include "mml/core/NumericValidation.h"
#include "mml/interfaces/IComplexFunction.h"
#include "mml/core/Derivation/DerivationComplex.h"

namespace MML {
	namespace RootFinding {

		/*********************************************************************/
		/*****           Complex Root Result Type                        *****/
		/*********************************************************************/

		/// Result of a complex root-finding operation.
		///
		/// Mirrors RootFindingResult but with Complex root and function value.
		struct ComplexRootFindingResult {
			/// The computed root value
			Complex root = Complex(0.0, 0.0);

			/// Function value at root: f(root), should be near zero if converged
			Complex function_value = Complex(0.0, 0.0);

			/// Number of iterations actually used
			int iterations_used = 0;

			/// True if algorithm converged within tolerance and iteration limits
			bool converged = false;

			/// Actual achieved tolerance |f(root)|
			Real achieved_tolerance = 0.0;

			/// Algorithm termination status
			AlgorithmStatus status = AlgorithmStatus::Success;

			/// Error message if not converged (empty string if successful)
			std::string error_message;

			/// Name of the algorithm used
			std::string algorithm_name;

			/// Elapsed wall-clock time in milliseconds
			double elapsed_time_ms = 0.0;
		};

		/*********************************************************************/
		/*****           Complex Newton-Raphson Method                   *****/
		/*********************************************************************/

		/// Find root of a complex function using Newton-Raphson iteration.
		///
		/// ALGORITHM:
		/// - z_{n+1} = z_n - f(z_n) / f'(z_n)
		/// - Derivative computed numerically via 4-point central difference
		/// - Converges when |f(z)| < tolerance
		///
		/// CONVERGENCE:
		/// - Quadratic convergence near simple roots
		/// - May diverge near multiple roots or saddle points
		/// - No bracket constraint (complex plane has no natural ordering)
		///
		/// @param func        Complex function f:C→C
		/// @param z0          Initial guess
		/// @param tolerance   Convergence tolerance on |f(z)|  
		/// @param max_iterations Maximum iterations (default: 100)
		/// @return ComplexRootFindingResult with root and diagnostics
		static ComplexRootFindingResult FindRootNewtonComplex(
			const IComplexFunction& func,
			Complex z0,
			Real tolerance = Precision::ComplexRootFindingTolerance,
			int max_iterations = 100)
		{
			ComplexRootFindingResult result;
			result.algorithm_name = "NewtonComplex";
			AlgorithmTimer timer;

			Complex z = z0;

			for (int j = 0; j < max_iterations; j++)
			{
				Complex fz = func(z);
				Real abs_fz = std::abs(fz);

				// Check convergence
				if (abs_fz < tolerance)
				{
					result.root = z;
					result.function_value = fz;
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = abs_fz;
					result.status = AlgorithmStatus::Success;
					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				}

				// Compute derivative numerically
				Complex dfz = Derivation::NDer4Complex(func, z);

				// Check for near-zero derivative
				Real abs_df = std::abs(dfz);
				if (abs_df < std::sqrt(Constants::Eps))
				{
					result.root = z;
					result.function_value = fz;
					result.iterations_used = j + 1;
					result.converged = false;
					result.achieved_tolerance = abs_fz;
					result.status = AlgorithmStatus::NumericalInstability;
					result.error_message = "Near-zero derivative in NewtonComplex (possible saddle point or multiple root)";
					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				}

				// Newton step
				Complex dz = fz / dfz;
				z -= dz;

				// Check for non-finite result
				if (!std::isfinite(z.real()) || !std::isfinite(z.imag()))
				{
					result.root = z;
					result.function_value = fz;
					result.iterations_used = j + 1;
					result.converged = false;
					result.status = AlgorithmStatus::NumericalInstability;
					result.error_message = "Non-finite iterate in NewtonComplex";
					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				}
			}

			// Maximum iterations exceeded
			Complex fz = func(z);
			result.root = z;
			result.function_value = fz;
			result.iterations_used = max_iterations;
			result.converged = false;
			result.achieved_tolerance = std::abs(fz);
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "NewtonComplex: max iterations (" + std::to_string(max_iterations)
			                     + ") exceeded, |f(z)| = " + std::to_string(std::abs(fz));
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		/*********************************************************************/
		/*****           Muller's Method                                 *****/
		/*********************************************************************/

		/// Find root of a complex function using Muller's method.
		///
		/// Muller's method is a generalization of the secant method that uses
		/// quadratic (rather than linear) interpolation through three points.
		/// Its key advantage is that it naturally finds **complex roots** even
		/// from real starting points, because quadratic interpolation can produce
		/// complex intermediate values.
		///
		/// ALGORITHM:
		/// Given three points (z0, f0), (z1, f1), (z2, f2):
		/// 1. Fit a quadratic through these three points
		/// 2. Find the root of the quadratic closest to z2
		/// 3. Replace the oldest point and repeat
		///
		/// The quadratic is constructed via divided differences:
		///   q = (z2 - z1) / (z1 - z0)
		///   A = q·f2 - q·(1+q)·f1 + q²·f0
		///   B = (2q+1)·f2 - (1+q)²·f1 + q²·f0
		///   C = (1+q)·f2
		///
		/// CONVERGENCE:
		/// - Superlinear: order ≈ 1.84 (faster than secant's 1.618)
		/// - Robust: works well even with poor initial guesses
		/// - Naturally handles complex arithmetic
		///
		/// @param func           Complex function f:C→C
		/// @param z0             First initial point
		/// @param z1             Second initial point
		/// @param z2             Third initial point (closest to expected root)
		/// @param tolerance      Convergence tolerance on |f(z)|
		/// @param max_iterations Maximum iterations (default: 100)
		/// @return ComplexRootFindingResult with root and diagnostics
		static ComplexRootFindingResult FindRootMuller(
			const IComplexFunction& func,
			Complex z0,
			Complex z1,
			Complex z2,
			Real tolerance = Precision::ComplexRootFindingTolerance,
			int max_iterations = 100)
		{
			ComplexRootFindingResult result;
			result.algorithm_name = "Muller";
			AlgorithmTimer timer;

			Complex f0 = func(z0);
			Complex f1 = func(z1);
			Complex f2 = func(z2);

			for (int j = 0; j < max_iterations; j++)
			{
				// Check convergence on the most recent point
				Real abs_f2 = std::abs(f2);
				if (abs_f2 < tolerance)
				{
					result.root = z2;
					result.function_value = f2;
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = abs_f2;
					result.status = AlgorithmStatus::Success;
					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				}

				// Divided differences for quadratic interpolation
				Complex q = (z2 - z1) / (z1 - z0);
				Complex A = q * f2 - q * (REAL(1.0) + q) * f1 + q * q * f0;
				Complex B = (REAL(2.0) * q + REAL(1.0)) * f2
				          - (REAL(1.0) + q) * (REAL(1.0) + q) * f1
				          + q * q * f0;
				Complex C = (REAL(1.0) + q) * f2;

				// Discriminant
				Complex disc = B * B - REAL(4.0) * A * C;
				Complex sqrt_disc = std::sqrt(disc);

				// Choose denominator with larger magnitude for numerical stability
				Complex denom_plus = B + sqrt_disc;
				Complex denom_minus = B - sqrt_disc;
				Complex denom = (std::abs(denom_plus) >= std::abs(denom_minus))
				              ? denom_plus : denom_minus;

				// Check for degenerate denominator
				if (std::abs(denom) < Constants::Eps)
				{
					result.root = z2;
					result.function_value = f2;
					result.iterations_used = j + 1;
					result.converged = false;
					result.achieved_tolerance = abs_f2;
					result.status = AlgorithmStatus::Stalled;
					result.error_message = "Muller: degenerate quadratic (zero denominator)";
					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				}

				// New estimate
				Complex z3 = z2 - (z2 - z1) * (REAL(2.0) * C / denom);

				// Check for non-finite result
				if (!std::isfinite(z3.real()) || !std::isfinite(z3.imag()))
				{
					result.root = z2;
					result.function_value = f2;
					result.iterations_used = j + 1;
					result.converged = false;
					result.status = AlgorithmStatus::NumericalInstability;
					result.error_message = "Non-finite iterate in Muller";
					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				}

				// Shift points: drop oldest, add newest
				z0 = z1; f0 = f1;
				z1 = z2; f1 = f2;
				z2 = z3; f2 = func(z2);
			}

			// Maximum iterations exceeded
			result.root = z2;
			result.function_value = f2;
			result.iterations_used = max_iterations;
			result.converged = false;
			result.achieved_tolerance = std::abs(f2);
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "Muller: max iterations (" + std::to_string(max_iterations)
			                     + ") exceeded, |f(z)| = " + std::to_string(std::abs(f2));
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		/// @brief Convenience overload: Muller's method with evenly-spaced starting points.
		/// @details Creates three starting points z0-h, z0, z0+h around initial guess z0.
		///          The step h defaults to 0.5 to provide reasonable initial spacing.
		static ComplexRootFindingResult FindRootMuller(
			const IComplexFunction& func,
			Complex z0,
			Real tolerance = Precision::ComplexRootFindingTolerance,
			int max_iterations = 100)
		{
			Real h = 0.5;
			return FindRootMuller(func, z0 - h, z0, z0 + h, tolerance, max_iterations);
		}

	} // namespace RootFinding
} // namespace MML

#endif // MML_ROOTFINDING_COMPLEX_H
