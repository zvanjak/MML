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
///               auto r = Integrate2D(f, SIMPSON, 0, 1, y_lo, y_hi);                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_INTEGRATION_2D_H
#define MML_INTEGRATION_2D_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"


namespace MML
{
	/// @brief Helper: integrates f(x,y) over y for fixed x
	struct Integral2DInner : public IRealFunction
	{
		mutable Real				_currX;
		const IScalarFunction<2>& _funcToIntegrate;

		Integral2DInner(const IScalarFunction<2>& func) : _currX(0.0), _funcToIntegrate(func) {}

		Real operator()(const Real y) const
		{
			VectorN<Real, 2> v{ _currX, y };
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

		Integral2DOuter(const IScalarFunction<2>& func, IntegrationMethod inMethod, 
										Real yy1(Real), Real yy2(Real)) 
				: _yRangeLow(yy1), _yRangeUpp(yy2), 
					_fInner(func), _funcToIntegrate(func), _integrMethod(inMethod)	{	}

		// for given x, will return (ie. integrate function over Y-range)
		Real operator()(const Real x) const
		{
			_fInner._currX = x;
			switch (_integrMethod)
			{
				case SIMPSON:
					return IntegrateSimpson(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
				case ROMBERG:
					return IntegrateRomberg(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
				case GAUSS10:
					return IntegrateGauss10(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
				case GAUSS10KRONROD21:
					return IntegrateGK21(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
				default:
					return IntegrateTrap(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
			}
		}
	};
	/// @brief Double integral ∬_D f(x,y) dx dy with variable y-limits
	/// @param func Scalar function f: ℝ² → ℝ to integrate
	/// @param method Integration method (TRAP, SIMPSON, ROMBERG, GAUSS10)
	/// @param x1,x2 Fixed x-integration bounds
	/// @param y1,y2 Functions defining y-bounds: y ∈ [y1(x), y2(x)]
	/// @return IntegrationResult with value, error_estimate, iterations, converged
	static IntegrationResult Integrate2D(const IScalarFunction<2>& func, IntegrationMethod method, 
																			 const Real x1, const Real x2, 
																			 Real y1(Real), Real y2(Real))
	{
		Integral2DOuter f1(func, method, y1, y2);

		switch (method)
		{
			case SIMPSON:
				return IntegrateSimpson(f1, x1, x2);
			case ROMBERG:
				return IntegrateRomberg(f1, x1, x2);
			case GAUSS10:
				return IntegrateGauss10(f1, x1, x2);  // Now returns IntegrationResult directly
			case GAUSS10KRONROD21:
				return IntegrateGK21(f1, x1, x2);
			default:
				return IntegrateTrap(f1, x1, x2);
		}
	}

	/******************************************************************************/
	/*****              Detailed API - 2D Integration                         *****/
	/******************************************************************************/

	/// 2D integration with full diagnostics
	static IntegrationDetailedResult Integrate2DDetailed(
		const IScalarFunction<2>& func, IntegrationMethod method,
		const Real x1, const Real x2,
		Real y1(Real), Real y2(Real),
		const IntegrationConfig& config = {})
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"Integrate2D", config,
			[&](IntegrationDetailedResult& result) {
				auto r = Integrate2D(func, method, x1, x2, y1, y2);
				IntegrationDetail::PopulateFromSimple(result, r, "Integrate2D");
			});
	}
}

#endif // MML_INTEGRATION_2D_H