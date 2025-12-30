///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Integration2D.h                                                     ///
///  Description: 2D numerical integration (double integrals)                         ///
///               Rectangular and triangular domain integration                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_INTEGRATION_2D_H
#define MML_INTEGRATION_2D_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"


namespace MML
{
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
					return IntegrateSimpson(_fInner, _yRangeLow(x), _yRangeUpp(x));
				case ROMBERG:
					return IntegrateRomberg(_fInner, _yRangeLow(x), _yRangeUpp(x));
				case GAUSS10:
					return IntegrateGauss10(_fInner, _yRangeLow(x), _yRangeUpp(x));
				default:
					return IntegrateTrap(_fInner, _yRangeLow(x), _yRangeUpp(x));
			}
		}
	};

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
			default:
				return IntegrateTrap(f1, x1, x2);
		}
	}
}

#endif // MML_INTEGRATION_2D_H