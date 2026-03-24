///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionHelpers.h                                                   ///
///  Description: Function utilities (composition, scaling, parameter binding)        ///
///               Adapters for converting between function types                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_FUNCTION_HELPERS_H
#define MML_FUNCTION_HELPERS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector/VectorN.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "base/Polynom.h"

#include "core/Derivation.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////////////////
	///                              TAYLOR SERIES EXPANSION                                ///
	///////////////////////////////////////////////////////////////////////////////////////////
	// Compute polynomial approximation of f(x) around point a using Taylor expansion:
	// f(x) ≈ f(a) + f'(a)(x-a) + f''(a)(x-a)²/2! + f'''(a)(x-a)³/3! + ...
	// 
	// Note on derivative accuracy: Higher derivatives use lower-order methods because
	// they require more function evaluations. This balances accuracy with efficiency.
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Compute 2nd order Taylor series approximation: f(x) ≈ f(a) + f'(a)(x-a) + f''(a)(x-a)²/2!
	/// @param f Function to approximate
	/// @param a Expansion point
	/// @return Polynomial approximation of degree 2 in standard form
	/// @warning Derivative uses NDer6 (6th order) and NSecDer4 (4th order) for accuracy
	inline PolynomRealFunc TaylorSeries2(const IRealFunction& f, Real a)
	{
		PolynomRealFunc ret(2);

		Real val = f(a);
		Real coef1 = Derivation::NDer6(f, a);           // f'(a), 6th order accuracy
		Real coef2 = Derivation::NSecDer4(f, a) / 2.0;  // f''(a)/2!, 4th order accuracy

		// Convert from Taylor form to standard polynomial: p(x) = c0 + c1*x + c2*x²
		ret[0] = val - coef1 * a + coef2 * POW2(a);
		ret[1] = coef1 - 2.0 * coef2 * a;
		ret[2] = coef2;

		return ret;
	}

	/// @brief Compute 3rd order Taylor series approximation: f(x) ≈ f(a) + f'(a)(x-a) + f''(a)(x-a)²/2! + f'''(a)(x-a)³/3!
	/// @param f Function to approximate
	/// @param a Expansion point
	/// @return Polynomial approximation of degree 3 in standard form
	/// @warning Uses NDer6 (1st), NSecDer4 (2nd), and NThirdDer2 (3rd) for balanced accuracy vs efficiency
	inline PolynomRealFunc TaylorSeries3(const IRealFunction& f, Real a)
	{
		PolynomRealFunc ret(3);

		Real val = f(a);
		Real coef1 = Derivation::NDer6(f, a);            // f'(a), 6th order accuracy
		Real coef2 = Derivation::NSecDer4(f, a) / 2.0;   // f''(a)/2!, 4th order accuracy
		Real coef3 = Derivation::NThirdDer2(f, a) / 6.0; // f'''(a)/3!, 2nd order accuracy

		// Convert from Taylor form to standard polynomial: p(x) = c0 + c1*x + c2*x² + c3*x³
		ret[0] = val - coef1 * a + coef2 * POW2(a) - coef3 * pow(a, 3);
		ret[1] = coef1 - 2.0 * coef2 * a + 3.0 * coef3 * POW2(a);
		ret[2] = coef2 - 3.0 * coef3 * a;
		ret[3] = coef3;

		return ret;
	}

	// Note: TaylorSeries4 is not implemented because Derivation doesn't have NFourthDer
	// If fourth derivative support is added to Derivation.h, TaylorSeries4 can be added here.

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         FIRST DERIVATIVE WRAPPER CLASSES                            ///
	///////////////////////////////////////////////////////////////////////////////////////////
	// These classes wrap a function and return its numerical derivative.
	// The number suffix indicates the order of accuracy of the numerical method.
	// 
	// WARNING: These store a reference to the original function. Ensure the original
	// function outlives these wrapper objects to avoid dangling references.
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief First derivative wrapper with 1st order accuracy (forward/backward difference)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncDerived1 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		/// @brief Constructor with automatic step size selection
		RealFuncDerived1(const IRealFunction& f) : _f(f), _step(0.0) {}
		/// @brief Constructor with custom step size
		RealFuncDerived1(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer1(_f, x, _step) : Derivation::NDer1(_f, x);
		}
	};

	/// @brief First derivative wrapper with 2nd order accuracy (central difference)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncDerived2 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		/// @brief Constructor with automatic step size selection
		RealFuncDerived2(const IRealFunction& f) : _f(f), _step(0.0) {}
		/// @brief Constructor with custom step size
		RealFuncDerived2(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer2(_f, x, _step) : Derivation::NDer2(_f, x);
		}
	};

	/// @brief First derivative wrapper with 4th order accuracy (5-point stencil)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncDerived4 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		/// @brief Constructor with automatic step size selection
		RealFuncDerived4(const IRealFunction& f) : _f(f), _step(0.0) {}
		/// @brief Constructor with custom step size
		RealFuncDerived4(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer4(_f, x, _step) : Derivation::NDer4(_f, x);
		}
	};

	/// @brief First derivative wrapper with 6th order accuracy (7-point stencil)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncDerived6 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		/// @brief Constructor with automatic step size selection
		RealFuncDerived6(const IRealFunction& f) : _f(f), _step(0.0) {}
		/// @brief Constructor with custom step size
		RealFuncDerived6(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer6(_f, x, _step) : Derivation::NDer6(_f, x);
		}
	};

	/// @brief First derivative wrapper with 8th order accuracy (9-point stencil)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncDerived8 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		/// @brief Constructor with automatic step size selection
		RealFuncDerived8(const IRealFunction& f) : _f(f), _step(0.0) {}
		/// @brief Constructor with custom step size
		RealFuncDerived8(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer8(_f, x, _step) : Derivation::NDer8(_f, x);
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                        SECOND DERIVATIVE WRAPPER CLASSES                            ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Second derivative wrapper with 2nd order accuracy (3-point stencil)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncSecondDerived2 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		/// @brief Constructor with automatic step size selection
		RealFuncSecondDerived2(const IRealFunction& f) : _f(f), _step(0.0) {}
		/// @brief Constructor with custom step size
		RealFuncSecondDerived2(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NSecDer2(_f, x, _step) : Derivation::NSecDer2(_f, x);
		}
	};

	/// @brief Second derivative wrapper with 4th order accuracy (5-point stencil)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncSecondDerived4 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		/// @brief Constructor with automatic step size selection
		RealFuncSecondDerived4(const IRealFunction& f) : _f(f), _step(0.0) {}
		/// @brief Constructor with custom step size
		RealFuncSecondDerived4(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NSecDer4(_f, x, _step) : Derivation::NSecDer4(_f, x);
		}
	};

	/// @brief Second derivative wrapper with 4th order accuracy (5-point stencil)
	/// @note Named "6" for future NSecDer6 implementation; currently uses NSecDer4
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncSecondDerived6 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		/// @brief Constructor with automatic step size selection
		RealFuncSecondDerived6(const IRealFunction& f) : _f(f), _step(0.0) {}
		/// @brief Constructor with custom step size
		RealFuncSecondDerived6(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NSecDer4(_f, x, _step) : Derivation::NSecDer4(_f, x);
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                      FUNCTION ARITHMETIC HELPER CLASSES                             ///
	///////////////////////////////////////////////////////////////////////////////////////////
	// Classes for performing arithmetic operations on functions.
	// All classes store references - ensure original functions outlive these wrappers.
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Sum of two functions: h(x) = f(x) + g(x)
	/// @warning Stores references - ensure f1, f2 outlive this object
	class RealFuncSum : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		/// @brief Constructor
		/// @param f1 First function
		/// @param f2 Second function
		RealFuncSum(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) + _f2(x); }
	};

	/// @brief Difference of two functions: h(x) = f(x) - g(x)
	/// @warning Stores references - ensure f1, f2 outlive this object
	class RealFuncDiff : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		/// @brief Constructor
		/// @param f1 First function (minuend)
		/// @param f2 Second function (subtrahend)
		RealFuncDiff(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) - _f2(x); }
	};

	/// @brief Product of two functions: h(x) = f(x) * g(x)
	/// @warning Stores references - ensure f1, f2 outlive this object
	class RealFuncProduct : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		/// @brief Constructor
		/// @param f1 First function
		/// @param f2 Second function
		RealFuncProduct(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) * _f2(x); }
	};

	/// @brief Quotient of two functions: h(x) = f(x) / g(x)
	/// @warning No division-by-zero check - ensure f2(x) ≠ 0 in your domain
	/// @warning Stores references - ensure f1, f2 outlive this object
	class RealFuncQuotient : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		/// @brief Constructor
		/// @param f1 Numerator function
		/// @param f2 Denominator function (must not be zero)
		RealFuncQuotient(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) / _f2(x); }
	};

	/// @brief Function composition: h(x) = f(g(x))
	/// @warning Stores references - ensure outer, inner outlive this object
	class RealFuncCompose : public IRealFunction
	{
		const IRealFunction& _outer;  // f
		const IRealFunction& _inner;  // g
	public:
		/// @brief Constructor
		/// @param outer Outer function f
		/// @param inner Inner function g
		RealFuncCompose(const IRealFunction& outer, const IRealFunction& inner) 
			: _outer(outer), _inner(inner) {}
		Real operator()(Real x) const { return _outer(_inner(x)); }
	};

	/// @brief Scalar multiplication: h(x) = c * f(x)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncScale : public IRealFunction
	{
		const IRealFunction& _f;
		Real _scale;
	public:
		/// @brief Constructor
		/// @param f Function to scale
		/// @param scale Scaling factor c
		RealFuncScale(const IRealFunction& f, Real scale) : _f(f), _scale(scale) {}
		Real operator()(Real x) const { return _scale * _f(x); }
	};

	/// @brief Scalar addition (vertical shift): h(x) = f(x) + c
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncShift : public IRealFunction
	{
		const IRealFunction& _f;
		Real _shift;
	public:
		/// @brief Constructor
		/// @param f Function to shift
		/// @param shift Vertical shift amount c
		RealFuncShift(const IRealFunction& f, Real shift) : _f(f), _shift(shift) {}
		Real operator()(Real x) const { return _f(x) + _shift; }
	};

	/// @brief Function negation: h(x) = -f(x)
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncNegate : public IRealFunction
	{
		const IRealFunction& _f;
	public:
		/// @brief Constructor
		/// @param f Function to negate
		RealFuncNegate(const IRealFunction& f) : _f(f) {}
		Real operator()(Real x) const { return -_f(x); }
	};

	/// @brief Absolute value: h(x) = |f(x)|
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncAbs : public IRealFunction
	{
		const IRealFunction& _f;
	public:
		/// @brief Constructor
		/// @param f Function to take absolute value of
		RealFuncAbs(const IRealFunction& f) : _f(f) {}
		Real operator()(Real x) const { return std::abs(_f(x)); }
	};

	/// @brief Power function: h(x) = [f(x)]^n
	/// @warning Stores reference to f - ensure f outlives this object
	class RealFuncPow : public IRealFunction
	{
		const IRealFunction& _f;
		Real _power;
	public:
		/// @brief Constructor
		/// @param f Base function
		/// @param power Exponent n
		RealFuncPow(const IRealFunction& f, Real power) : _f(f), _power(power) {}
		Real operator()(Real x) const { return std::pow(_f(x), _power); }
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                        COMPARISON / DISTANCE HELPER CLASSES                         ///
	///////////////////////////////////////////////////////////////////////////////////////////
	// Useful for computing norms and comparing functions.
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Absolute difference: h(x) = |f(x) - g(x)|
	/// @note Useful for L¹ norm and maximum absolute error computation
	/// @warning Stores references - ensure f1, f2 outlive this object
	class RealFuncAbsDiff : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		/// @brief Constructor
		/// @param f1 First function
		/// @param f2 Second function
		RealFuncAbsDiff(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return std::abs(_f1(x) - _f2(x)); }
	};

	/// @brief Squared difference: h(x) = [f(x) - g(x)]²
	/// @note Useful for L² norm computation and least-squares fitting
	/// @warning Stores references - ensure f1, f2 outlive this object
	class RealFuncDiffSqr : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		/// @brief Constructor
		/// @param f1 First function
		/// @param f2 Second function
		RealFuncDiffSqr(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return POW2(_f1(x) - _f2(x)); }
	};

} // end namespace

#endif