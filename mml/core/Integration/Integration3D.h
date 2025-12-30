///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Integration3D.h                                                     ///
///  Description: 3D numerical integration (triple integrals)                         ///
///               Box, spherical, and cylindrical domain integration                  ///
///               Supports TRAP, SIMPSON, ROMBERG, and GAUSS10 methods                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_INTEGRATION_3D_H
#define MML_INTEGRATION_3D_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"


namespace MML
{
	struct Integral3DInnermost : public IRealFunction
	{
		mutable Real					_currX, _currY;
		const IScalarFunction<3>&		_funcToIntegrate;
		IntegrationMethod				_integrMethod;

		Integral3DInnermost(const IScalarFunction<3>& func, IntegrationMethod method = GAUSS10) 
			: _funcToIntegrate(func), _currX{ 0 }, _currY{ 0 }, _integrMethod(method) {}
		
		Real operator()(const Real z) const
		{
			VectorN<Real, 3> v{ _currX, _currY, z };

			return _funcToIntegrate(v);
		}
	};

	struct Integral3DInner : public IRealFunction
	{
		mutable Integral3DInnermost _fInnermost;

		const IScalarFunction<3>& _funcToIntegrate;
		IntegrationMethod			_integrMethod;

		Real(*_zRangeLow)(Real, Real);
		Real(*_zRangeUpp)(Real, Real);

		Integral3DInner(const IScalarFunction<3>& func, IntegrationMethod method,
										Real zz1(Real, Real), Real zz2(Real, Real)) 
			: _zRangeLow(zz1), _zRangeUpp(zz2), _funcToIntegrate(func), 
			  _fInnermost(func, method), _integrMethod(method) {}

		Real operator()(const Real y) const
		{
			_fInnermost._currY = y;

			switch (_integrMethod)
			{
				case SIMPSON:
					return IntegrateSimpson(_fInnermost, 
											_zRangeLow(_fInnermost._currX, y), 
											_zRangeUpp(_fInnermost._currX, y));
				case ROMBERG:
					return IntegrateRomberg(_fInnermost, 
											_zRangeLow(_fInnermost._currX, y), 
											_zRangeUpp(_fInnermost._currX, y));
				case GAUSS10:
					return IntegrateGauss10(_fInnermost, 
											_zRangeLow(_fInnermost._currX, y), 
											_zRangeUpp(_fInnermost._currX, y));
				default:  // TRAP
					return IntegrateTrap(_fInnermost, 
										 _zRangeLow(_fInnermost._currX, y), 
										 _zRangeUpp(_fInnermost._currX, y));
			}
		}
	};

	struct Integral3DOuter : public IRealFunction
	{
		mutable Integral3DInner _fInner;

		const IScalarFunction<3>& _funcToIntegrate;
		IntegrationMethod			_integrMethod;
		Real(*_yRangeLow)(Real);
		Real(*_yRangeUpp)(Real);

		Integral3DOuter(const IScalarFunction<3>& func, IntegrationMethod method,
						Real yy1(Real), Real yy2(Real), 
						Real z1(Real, Real), Real z2(Real, Real)) 
			: _yRangeLow(yy1), _yRangeUpp(yy2), _funcToIntegrate(func), 
			  _fInner(func, method, z1, z2), _integrMethod(method)
		{
		}

		Real operator()(const Real x) const
		{
			_fInner._fInnermost._currX = x;

			switch (_integrMethod)
			{
				case SIMPSON:
					return IntegrateSimpson(_fInner, _yRangeLow(x), _yRangeUpp(x));
				case ROMBERG:
					return IntegrateRomberg(_fInner, _yRangeLow(x), _yRangeUpp(x));
				case GAUSS10:
					return IntegrateGauss10(_fInner, _yRangeLow(x), _yRangeUpp(x));
				default:  // TRAP
					return IntegrateTrap(_fInner, _yRangeLow(x), _yRangeUpp(x));
			}
		}
	};

	/// 3D integration with selectable method (TRAP, SIMPSON, ROMBERG, GAUSS10)
	static IntegrationResult Integrate3D(const IScalarFunction<3>& func, 
										 IntegrationMethod method,
										 const Real x1, const Real x2, 
										 Real y1(Real), Real y2(Real),
										 Real z1(Real, Real), Real z2(Real, Real))
	{
		Integral3DOuter f1(func, method, y1, y2, z1, z2);

		switch (method)
		{
			case SIMPSON:
				return IntegrateSimpson(f1, x1, x2);
			case ROMBERG:
				return IntegrateRomberg(f1, x1, x2);
			case GAUSS10:
				return IntegrateGauss10(f1, x1, x2);
			default:  // TRAP
				return IntegrateTrap(f1, x1, x2);
		}
	}

	/// 3D integration with default GAUSS10 method (backward compatible)
	static IntegrationResult Integrate3D(const IScalarFunction<3>& func, 
										 const Real x1, const Real x2, 
										 Real y1(Real), Real y2(Real),
										 Real z1(Real, Real), Real z2(Real, Real))
	{
		return Integrate3D(func, GAUSS10, x1, x2, y1, y2, z1, z2);
	}	
}

#endif // MML_INTEGRATION_3D_H
