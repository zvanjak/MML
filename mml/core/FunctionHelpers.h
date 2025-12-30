///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionHelpers.h                                                   ///
///  Description: Function utilities (composition, scaling, parameter binding)        ///
///               Adapters for converting between function types                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_FUNCTION_HELPERS_H
#define MML_FUNCTION_HELPERS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/VectorN.h"
#include "base/Vector.h"
#include "base/Matrix.h"
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

	class RealFuncDerived1 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		RealFuncDerived1(const IRealFunction& f) : _f(f), _step(0.0) {}
		RealFuncDerived1(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer1(_f, x, _step) : Derivation::NDer1(_f, x);
		}
	};

	class RealFuncDerived2 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		RealFuncDerived2(const IRealFunction& f) : _f(f), _step(0.0) {}
		RealFuncDerived2(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer2(_f, x, _step) : Derivation::NDer2(_f, x);
		}
	};

	class RealFuncDerived4 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		RealFuncDerived4(const IRealFunction& f) : _f(f), _step(0.0) {}
		RealFuncDerived4(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer4(_f, x, _step) : Derivation::NDer4(_f, x);
		}
	};

	class RealFuncDerived6 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		RealFuncDerived6(const IRealFunction& f) : _f(f), _step(0.0) {}
		RealFuncDerived6(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer6(_f, x, _step) : Derivation::NDer6(_f, x);
		}
	};

	class RealFuncDerived8 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		RealFuncDerived8(const IRealFunction& f) : _f(f), _step(0.0) {}
		RealFuncDerived8(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NDer8(_f, x, _step) : Derivation::NDer8(_f, x);
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                        SECOND DERIVATIVE WRAPPER CLASSES                            ///
	///////////////////////////////////////////////////////////////////////////////////////////

	class RealFuncSecondDerived2 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		RealFuncSecondDerived2(const IRealFunction& f) : _f(f), _step(0.0) {}
		RealFuncSecondDerived2(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NSecDer2(_f, x, _step) : Derivation::NSecDer2(_f, x);
		}
	};

	class RealFuncSecondDerived4 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		RealFuncSecondDerived4(const IRealFunction& f) : _f(f), _step(0.0) {}
		RealFuncSecondDerived4(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NSecDer4(_f, x, _step) : Derivation::NSecDer4(_f, x);
		}
	};

	class RealFuncSecondDerived6 : public IRealFunction
	{
		const IRealFunction& _f;
		Real _step;
	public:
		RealFuncSecondDerived6(const IRealFunction& f) : _f(f), _step(0.0) {}
		RealFuncSecondDerived6(const IRealFunction& f, Real step) : _f(f), _step(step) {}

		Real operator()(Real x) const {
			return _step != 0.0 ? Derivation::NSecDer4(_f, x, _step, nullptr) : Derivation::NSecDer4(_f, x, nullptr);
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                      FUNCTION ARITHMETIC HELPER CLASSES                             ///
	///////////////////////////////////////////////////////////////////////////////////////////
	// Classes for performing arithmetic operations on functions.
	// All classes store references - ensure original functions outlive these wrappers.
	///////////////////////////////////////////////////////////////////////////////////////////

	// Sum: h(x) = f(x) + g(x)
	class RealFuncSum : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		RealFuncSum(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) + _f2(x); }
	};

	// Difference: h(x) = f(x) - g(x)
	class RealFuncDiff : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		RealFuncDiff(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) - _f2(x); }
	};
	// Backward compatibility alias
	using RealFuncDiffHelper = RealFuncDiff;

	// Product: h(x) = f(x) * g(x)
	class RealFuncProduct : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		RealFuncProduct(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) * _f2(x); }
	};

	// Quotient: h(x) = f(x) / g(x)
	// WARNING: No division-by-zero check for performance. Ensure g(x) ≠ 0 in your domain.
	class RealFuncQuotient : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		RealFuncQuotient(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) / _f2(x); }
	};

	// Composition: h(x) = f(g(x))
	class RealFuncCompose : public IRealFunction
	{
		const IRealFunction& _outer;  // f
		const IRealFunction& _inner;  // g
	public:
		RealFuncCompose(const IRealFunction& outer, const IRealFunction& inner) 
			: _outer(outer), _inner(inner) {}
		Real operator()(Real x) const { return _outer(_inner(x)); }
	};

	// Scalar multiplication: h(x) = c * f(x)
	class RealFuncScale : public IRealFunction
	{
		const IRealFunction& _f;
		Real _scale;
	public:
		RealFuncScale(const IRealFunction& f, Real scale) : _f(f), _scale(scale) {}
		Real operator()(Real x) const { return _scale * _f(x); }
	};

	// Scalar addition (shift): h(x) = f(x) + c
	class RealFuncShift : public IRealFunction
	{
		const IRealFunction& _f;
		Real _shift;
	public:
		RealFuncShift(const IRealFunction& f, Real shift) : _f(f), _shift(shift) {}
		Real operator()(Real x) const { return _f(x) + _shift; }
	};

	// Negation: h(x) = -f(x)
	class RealFuncNegate : public IRealFunction
	{
		const IRealFunction& _f;
	public:
		RealFuncNegate(const IRealFunction& f) : _f(f) {}
		Real operator()(Real x) const { return -_f(x); }
	};

	// Absolute value: h(x) = |f(x)|
	class RealFuncAbs : public IRealFunction
	{
		const IRealFunction& _f;
	public:
		RealFuncAbs(const IRealFunction& f) : _f(f) {}
		Real operator()(Real x) const { return std::abs(_f(x)); }
	};

	// Power: h(x) = f(x)^n
	class RealFuncPow : public IRealFunction
	{
		const IRealFunction& _f;
		Real _power;
	public:
		RealFuncPow(const IRealFunction& f, Real power) : _f(f), _power(power) {}
		Real operator()(Real x) const { return std::pow(_f(x), _power); }
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                        COMPARISON / DISTANCE HELPER CLASSES                         ///
	///////////////////////////////////////////////////////////////////////////////////////////
	// Useful for computing norms and comparing functions.
	///////////////////////////////////////////////////////////////////////////////////////////

	// Absolute difference: h(x) = |f(x) - g(x)|
	class RealFuncAbsDiff : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		RealFuncAbsDiff(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return std::abs(_f1(x) - _f2(x)); }
	};
	// Backward compatibility alias
	using RealFuncAbsDiffHelper = RealFuncAbsDiff;

	// Squared difference: h(x) = (f(x) - g(x))²
	// Useful for L² norm computation
	class RealFuncDiffSqr : public IRealFunction
	{
		const IRealFunction& _f1;
		const IRealFunction& _f2;
	public:
		RealFuncDiffSqr(const IRealFunction& f1, const IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return POW2(_f1(x) - _f2(x)); }
	};
	// Backward compatibility alias
	using RealFuncDiffSqrHelper = RealFuncDiffSqr;

} // end namespace

#endif