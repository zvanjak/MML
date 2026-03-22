///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Integration3D.h                                                     ///
///  Description: 3D numerical integration (triple integrals)                         ///
///               Variable limits: ∭_V f(x,y,z) dx dy dz                              ///
///               where V = [x1,x2] × [y1(x),y2(x)] × [z1(x,y),z2(x,y)]               ///
///                                                                                   ///
///  Usage:       auto r = Integrate3D(f, GAUSS10, x1, x2, y1, y2, z1, z2);           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_INTEGRATION_3D_H
#define MML_INTEGRATION_3D_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"


namespace MML
{
	/// @brief Helper: integrates f(x,y,z) over z for fixed x,y
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

	/// @brief Helper: integrates over y, calling innermost for each y
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
											_zRangeUpp(_fInnermost._currX, y)).value;
				case ROMBERG:
					return IntegrateRomberg(_fInnermost, 
											_zRangeLow(_fInnermost._currX, y), 
											_zRangeUpp(_fInnermost._currX, y)).value;
				case GAUSS10:
					return IntegrateGauss10(_fInnermost, 
											_zRangeLow(_fInnermost._currX, y), 
											_zRangeUpp(_fInnermost._currX, y)).value;
				case GAUSS10KRONROD21:
					return IntegrateGK21(_fInnermost, 
											_zRangeLow(_fInnermost._currX, y), 
											_zRangeUpp(_fInnermost._currX, y)).value;
				default:  // TRAP
					return IntegrateTrap(_fInnermost, 
										 _zRangeLow(_fInnermost._currX, y), 
										 _zRangeUpp(_fInnermost._currX, y)).value;
			}
		}
	};

	/// @brief Helper: outer integrator that drives the nested integration
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
					return IntegrateSimpson(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
				case ROMBERG:
					return IntegrateRomberg(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
				case GAUSS10:
					return IntegrateGauss10(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
				case GAUSS10KRONROD21:
					return IntegrateGK21(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
				default:  // TRAP
					return IntegrateTrap(_fInner, _yRangeLow(x), _yRangeUpp(x)).value;
			}
		}
	};

	/// @brief Triple integral ∭_V f(x,y,z) dx dy dz with variable limits
	/// @param func Scalar function f: ℝ³ → ℝ to integrate
	/// @param method Integration method (TRAP, SIMPSON, ROMBERG, GAUSS10)
	/// @param x1,x2 Fixed x-integration bounds
	/// @param y1,y2 Functions: y ∈ [y1(x), y2(x)]
	/// @param z1,z2 Functions: z ∈ [z1(x,y), z2(x,y)]
	/// @return IntegrationResult with value, error_estimate, iterations, converged
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
			case GAUSS10KRONROD21:
				return IntegrateGK21(f1, x1, x2);
			default:  // TRAP
				return IntegrateTrap(f1, x1, x2);
		}
	}

	/// @brief Triple integral with default GAUSS10 method
	/// @note Convenience overload - uses GAUSS10 for speed
	static IntegrationResult Integrate3D(const IScalarFunction<3>& func, 
										 const Real x1, const Real x2, 
										 Real y1(Real), Real y2(Real),
										 Real z1(Real, Real), Real z2(Real, Real))
	{
		return Integrate3D(func, GAUSS10, x1, x2, y1, y2, z1, z2);
	}

	/******************************************************************************/
	/*****              Detailed API - 3D Integration                         *****/
	/******************************************************************************/

	/// 3D integration with method selection and full diagnostics
	static IntegrationDetailedResult Integrate3DDetailed(
		const IScalarFunction<3>& func, IntegrationMethod method,
		const Real x1, const Real x2,
		Real y1(Real), Real y2(Real),
		Real z1(Real, Real), Real z2(Real, Real),
		const IntegrationConfig& config = {})
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"Integrate3D", config,
			[&](IntegrationDetailedResult& result) {
				auto r = Integrate3D(func, method, x1, x2, y1, y2, z1, z2);
				IntegrationDetail::PopulateFromSimple(result, r, "Integrate3D");
			});
	}

	/// 3D integration with default GAUSS10 method and full diagnostics
	static IntegrationDetailedResult Integrate3DDetailed(
		const IScalarFunction<3>& func,
		const Real x1, const Real x2,
		Real y1(Real), Real y2(Real),
		Real z1(Real, Real), Real z2(Real, Real),
		const IntegrationConfig& config = {})
	{
		return Integrate3DDetailed(func, GAUSS10, x1, x2, y1, y2, z1, z2, config);
	}

} // end namespace MML

#endif // MML_INTEGRATION_3D_H
