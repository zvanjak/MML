///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ComplexFunction.h                                                   ///
///  Description: Concrete wrappers for complex-valued functions                      ///
///               Function pointer and lambda adapters for C→C and R→C               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMPLEX_FUNCTION_H
#define MML_COMPLEX_FUNCTION_H

#include "MMLBase.h"
#include "interfaces/IComplexFunction.h"

#include <functional>

namespace MML
{
	///////////////////////////     COMPLEX FUNCTION  C → C     ////////////////////////////

	/// @brief Wrapper for a complex function C → C using function pointer.
	/// @details Wraps a C-style function pointer for use with complex analysis algorithms.
	class ComplexFunction : public IComplexFunction
	{
		Complex(*_func)(Complex);   ///< Function pointer to evaluate
	public:
		/// @brief Constructs wrapper from function pointer.
		/// @param inFunc Function pointer Complex → Complex
		ComplexFunction(Complex(*inFunc)(Complex)) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param z Complex input value
		/// @return f(z)
		Complex operator()(Complex z) const override { return _func(z); }
	};

	/// @brief Wrapper for a complex function C → C using std::function.
	/// @details Allows using lambdas, functors, and bound functions.
	/// @warning **Lambda capture lifetime hazard:** The std::function is copied into this
	///          wrapper, but if the source lambda captures local variables **by reference**,
	///          those references will dangle if this object outlives the captured variables.
	///          Safe:   `ComplexFunctionFromStdFunc f([](Complex z) { return z*z; });`
	///          Safe:   `ComplexFunctionFromStdFunc f([a](Complex z) { return a*z; });`  (capture by value)
	///          UNSAFE: `ComplexFunctionFromStdFunc f([&a](Complex z) { return a*z; });` (if a goes out of scope)
	class ComplexFunctionFromStdFunc : public IComplexFunction
	{
		std::function<Complex(Complex)> _func;   ///< Callable object to evaluate
	public:
		/// @brief Constructs wrapper from std::function.
		/// @param inFunc Callable object Complex → Complex (lambda, functor, etc.)
		ComplexFunctionFromStdFunc(const std::function<Complex(Complex)>& inFunc) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param z Complex input value
		/// @return f(z)
		Complex operator()(Complex z) const override { return _func(z); }
	};

	///////////////////////////     REAL TO COMPLEX  R → C     ////////////////////////////

	/// @brief Wrapper for a real-to-complex function R → C using function pointer.
	/// @details Used for contour parameterizations, time-dependent complex signals, etc.
	class RealToComplexFunction : public IRealToComplexFunction
	{
		Complex(*_func)(Real);   ///< Function pointer to evaluate
	public:
		/// @brief Constructs wrapper from function pointer.
		/// @param inFunc Function pointer Real → Complex
		RealToComplexFunction(Complex(*inFunc)(Real)) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param t Real input parameter
		/// @return f(t)
		Complex operator()(Real t) const override { return _func(t); }
	};

	/// @brief Wrapper for a real-to-complex function R → C using std::function.
	/// @details Allows using lambdas, functors, and bound functions.
	/// @warning **Lambda capture lifetime hazard:** The std::function is copied into this
	///          wrapper, but if the source lambda captures local variables **by reference**,
	///          those references will dangle if this object outlives the captured variables.
	///          Safe:   `RealToComplexFunctionFromStdFunc f([](Real t) { return Complex(cos(t), sin(t)); });`
	///          UNSAFE: `RealToComplexFunctionFromStdFunc f([&r](Real t) { return Complex(r, t); });`
	class RealToComplexFunctionFromStdFunc : public IRealToComplexFunction
	{
		std::function<Complex(Real)> _func;   ///< Callable object to evaluate
	public:
		/// @brief Constructs wrapper from std::function.
		/// @param inFunc Callable object Real → Complex (lambda, functor, etc.)
		RealToComplexFunctionFromStdFunc(const std::function<Complex(Real)>& inFunc) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param t Real input parameter
		/// @return f(t)
		Complex operator()(Real t) const override { return _func(t); }
	};

} // namespace MML

#endif // MML_COMPLEX_FUNCTION_H
