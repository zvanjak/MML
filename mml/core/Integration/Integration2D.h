///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Integration2D.h                                                     ///
///  Description: 2D numerical integration (double integrals)                         ///
///               Variable limits: ∬_D f(x,y) dx dy where D = [x1,x2] × [y1(x),y2(x)] ///
///                                                                                   ///
///  Usage:       // Integrate over triangular region                                 ///
///               auto y_lo = [](Real x) { return 0.0; };                             ///
///               auto y_hi = [](Real x) { return 1.0 - x; };                         ///
///               auto r = Integrate2DCore(f, SIMPSON, 0, 1, y_lo, y_hi);                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_INTEGRATION_2D_H
#define MML_INTEGRATION_2D_H

#include <mml/MMLBase.h>

#include <mml/interfaces/IFunction.h>
#include "Integration1D.h"

namespace MML
{
	/// @brief Helper: integrates f(x,y) over y for fixed x
	struct Integral2DInner : public IRealFunction
	{
		mutable Real				_currX;
		const IScalarFunction<2>& _funcToIntegrate;
		mutable int _evaluations = 0;

		Integral2DInner(const IScalarFunction<2>& func) : _currX(0.0), _funcToIntegrate(func) {}

		Real operator()(const Real y) const
		{
			VectorN<Real, 2> v{ _currX, y };
			if (_evaluations == std::numeric_limits<int>::max())
				throw ResourceLimitError("Nested integration evaluation counter exhausted");
			++_evaluations;
			return _funcToIntegrate(v);
		}
	};

	/// @brief Helper: outer integrator that calls inner for each x value
	struct Integral2DOuter : public IRealFunction
	{
		mutable Integral2DInner _fInner;

		IntegrationMethod			_integrMethod;
		const IScalarFunction<2>&		_funcToIntegrate;

		Real(*_yRangeLow)(Real);
		Real(*_yRangeUpp)(Real);
		IntegrationConfig _config;
		mutable IntegrationDetail::NestedErrors _errors;

		Integral2DOuter(const IScalarFunction<2>& func, IntegrationMethod inMethod, 
										Real yy1(Real), Real yy2(Real), const IntegrationConfig& config = {})
				: _fInner(func), _integrMethod(inMethod), _funcToIntegrate(func),
					_yRangeLow(yy1), _yRangeUpp(yy2), _config(config) {
			if (!yy1 || !yy2) throw NotImplementedError("Nested integration requires both bound callbacks");
		}

		// for given x, will return (ie. integrate function over Y-range)
		Real operator()(const Real x) const
		{
			_fInner._currX = x;
			const auto inner = IntegrationDetail::NestedAxis(_fInner, _yRangeLow(x), _yRangeUpp(x), _integrMethod, _config);
			_errors.Observe(inner);
			return inner.value;
		}
	};
	static void RunIntegration2D(IntegrationResult& result, const IScalarFunction<2>& function, IntegrationMethod method,
		Real a, Real b, Real lower(Real), Real upper(Real), const IntegrationConfig& config = {}) {
		Integral2DOuter outer(function, method, lower, upper, config);
		IntegrationDetail::EvaluationCountGuard count{result, outer._fInner._evaluations};
		auto integrated = IntegrationDetail::NestedAxis(outer, a, b, method, config, &result);
		outer._errors.Combine(integrated, a, b, method);
		IntegrationDetail::Adopt(result, integrated, "Integrate2D");
	}
	/// @brief Double integral ∬_D f(x,y) dx dy with variable y-limits
	/// @param func Scalar function f: ℝ² → ℝ to integrate
	/// @param method Integration method (TRAP, SIMPSON, ROMBERG, GAUSS10)
	/// @param x1,x2 Fixed x-integration bounds
	/// @param y1,y2 Functions defining y-bounds: y ∈ [y1(x), y2(x)]
	/// @return IntegrationResult with value, error_estimate, iterations, converged
	static IntegrationResult Integrate2DCore(const IScalarFunction<2>& func, IntegrationMethod method, 
																			 const Real x1, const Real x2, 
																			 Real y1(Real), Real y2(Real))
	{
		IntegrationResult result;
		RunIntegration2D(result, func, method, x1, x2, y1, y2);
		return result;
	}

	static inline Real Integrate2D(const IScalarFunction<2>& func, IntegrationMethod method,
																			 const Real x1, const Real x2, 
																			 Real y1(Real), Real y2(Real))
	{
		return Integrate2DCore(func, method, x1, x2, y1, y2).ValueOrThrow("Integrate2D");
	}

	/******************************************************************************/
	/*****              Detailed API - 2D Integration                         *****/
	/******************************************************************************/

	/// 2D integration with full diagnostics
	static inline IntegrationDetailedResult Integrate2DDetailed(
		const IScalarFunction<2>& func, IntegrationMethod method,
		const Real x1, const Real x2,
		Real y1(Real), Real y2(Real),
		const IntegrationConfig& config = {})
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"Integrate2D", config,
			[&](IntegrationDetailedResult& result) {
				RunIntegration2D(result, func, method, x1, x2, y1, y2, config);
			});
	}
}

#endif // MML_INTEGRATION_2D_H