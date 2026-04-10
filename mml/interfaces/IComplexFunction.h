///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IComplexFunction.h                                                  ///
///  Description: Complex function interface hierarchy                                ///
///               Abstract base classes for complex-valued functions                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IComplexFunction.h
 * @brief Function interfaces for complex-valued functions.
 * 
 * **Complex Functions (C → C):**
 * - IComplexFunction - Basic complex-to-complex function
 * 
 * **Real-to-Complex Functions (R → C):**
 * - IRealToComplexFunction - Maps real parameter to complex value
 *   Used for contour parameterizations, signals, wavefunctions
 */

#if !defined MML_ICOMPLEX_FUNCTION_H
#define MML_ICOMPLEX_FUNCTION_H

#include "MMLBase.h"
#include "interfaces/IFunction.h"

namespace MML
{
	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for complex-valued functions of a complex variable.
	 * 
	 * Represents functions f: C → C. Used in complex analysis algorithms:
	 * contour integration, residue computation, complex root finding.
	 * 
	 * @note Implementations should ensure operator() is const and thread-safe.
	 */
	class IComplexFunction : public IFunction<Complex, Complex>
	{
	public:
		/**
		 * @brief Evaluate the function at complex point z.
		 * @param z The complex input value
		 * @return f(z) as a complex value
		 */
		virtual Complex operator()(Complex z) const = 0;

		virtual ~IComplexFunction() {}
	};

	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for complex-valued functions of a real variable.
	 * 
	 * Represents functions f: R → C. Used for:
	 * - Contour parameterizations γ(t) in complex analysis
	 * - Complex signals x(t) = A(t)·exp(iφ(t))
	 * - Quantum wavefunctions ψ(x,t) along a line
	 */
	class IRealToComplexFunction : public IFunction<Complex, Real>
	{
	public:
		/**
		 * @brief Evaluate the function at real parameter t.
		 * @param t The real input parameter
		 * @return f(t) as a complex value
		 */
		virtual Complex operator()(Real t) const = 0;

		virtual ~IRealToComplexFunction() {}
	};

} // namespace MML

#endif // MML_ICOMPLEX_FUNCTION_H
