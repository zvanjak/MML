///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationComplexStep.h                                             ///
///  Description: Complex-step derivative method                                      ///
///               f'(x) = Im[f(x + ih)] / h                                          ///
///               Works with any callable accepting Complex arguments                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_COMPLEX_STEP_H
#define MML_DERIVATION_COMPLEX_STEP_H

#include "MMLBase.h"

#include <complex>
#include <type_traits>

namespace MML
{
	namespace Derivation
	{
		///////////////////////////////////////////////////////////////////////////////////
		///                     Complex-step derivative method                           ///
		///                                                                              ///
		///  The complex-step method computes derivatives using:                         ///
		///      f'(x) = Im[f(x + ih)] / h                                              ///
		///                                                                              ///
		///  Advantages over finite differences:                                         ///
		///  - No subtractive cancellation (no f(x+h) - f(x))                           ///
		///  - Step size h can be extremely small (e.g. 1e-100)                          ///
		///  - Accuracy limited only by function's analytic continuation                 ///
		///  - Essentially exact for analytic functions                                  ///
		///                                                                              ///
		///  Limitations:                                                                ///
		///  - Requires function to accept std::complex<Real> arguments                  ///
		///  - Function must be analytic (satisfy Cauchy-Riemann equations)               ///
		///  - Cannot use with IRealFunction (Real-only interface)                       ///
		///  - Operations inside f must handle complex arithmetic (use std::sin,         ///
		///    std::cos, std::exp, etc. — NOT C math.h functions)                        ///
		///                                                                              ///
		///  Usage:                                                                      ///
		///      auto f = [](Complex z) { return z*z*z + Complex(2.0)*z; };              ///
		///      Real deriv = Derivation::ComplexStep(f, 1.0);  // f'(1) = 5            ///
		///                                                                              ///
		///  Reference: Squire & Trapp (1998), "Using Complex Variables to Estimate      ///
		///             Derivatives of Real Functions"                                   ///
		///////////////////////////////////////////////////////////////////////////////////

		/// @brief Compute first derivative using the complex-step method.
		/// @tparam F Callable type accepting Complex and returning Complex
		/// @param f Function to differentiate (must accept Complex arguments)
		/// @param x Point at which to evaluate the derivative
		/// @param h Step size (can be extremely small for double)
		/// @return Approximation of f'(x)
		template<typename F>
		static Real ComplexStep(F&& f, Real x, Real h = Precision::ComplexStepH)
		{
			Complex z(x, h);
			Complex fz = f(z);
			return fz.imag() / h;
		}

		/// @brief Compute second derivative using the complex-step method.
		///
		/// Uses the formula:
		///   f''(x) = 2 * [f(x) - Re(f(x + ih))] / h^2
		///
		/// @tparam F Callable type accepting Complex and returning Complex
		/// @param f Function to differentiate (must accept Complex arguments)
		/// @param x Point at which to evaluate the second derivative
		/// @param h Step size (optimal is larger than for first derivative)
		/// @return Approximation of f''(x)
		template<typename F>
		static Real ComplexStep2(F&& f, Real x, Real h = Precision::ComplexStepH2)
		{
			Real fx = std::real(f(Complex(x, 0.0)));
			Complex fz = f(Complex(x, h));
			return 2.0 * (fx - fz.real()) / (h * h);
		}

	} // namespace Derivation
} // namespace MML

#endif // MML_DERIVATION_COMPLEX_STEP_H
