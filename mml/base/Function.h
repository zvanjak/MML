///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Function.h                                                          ///
///  Description: Function wrappers for scalar and vector functions                   ///
///               Lambda adapters, composition, and evaluation utilities              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_FUNCTION_H
#define MML_FUNCTION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/VectorN.h"
#include "base/Vector.h"
#include "base/Matrix.h"

#include <functional>

namespace MML
{
	///////////////////////////     REAL FUNCTION      ////////////////////////////////////
	class RealFunction : public IRealFunction
	{
		Real(*_func)(const Real);
	public:
		RealFunction(Real(*inFunc)(const Real)) : _func(inFunc) {}

		Real operator()(const Real x) const { return _func(x); }
	};
	class RealFunctionFromStdFunc : public IRealFunction
	{
		std::function<Real(const Real)> _func;
	public:
		RealFunctionFromStdFunc(const std::function<Real(const Real)> &inFunc) : _func(inFunc) {}

		Real operator()(const Real x) const { return _func(x); }
	};

	///////////////////////////     SCALAR FUNCTION       //////////////////////////////////
	template<int N>
	class ScalarFunction : public IScalarFunction<N>
	{
		Real(*_func)(const VectorN<Real, N>&);
	public:
		ScalarFunction(Real(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	template<int N>
	class ScalarFunctionFromStdFunc : public IScalarFunction<N>
	{
		std::function<Real(const VectorN<Real, N>&)> _func;
	public:
		ScalarFunctionFromStdFunc(const std::function<Real(const VectorN<Real, N>&)> inFunc) : _func(inFunc) {}

		Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/////////////////////////    VECTOR FUNCTION N -> N      ///////////////////////////////////
	template<int N>
	class VectorFunction : public IVectorFunction<N>
	{
		VectorN<Real, N>(*_func)(const VectorN<Real, N>&);
	public:
		VectorFunction(VectorN<Real, N>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	template<int N>
	class VectorFunctionFromStdFunc : public IVectorFunction<N>
	{
		std::function<VectorN<Real, N>(const VectorN<Real, N>&)> _func;
	public:
		VectorFunctionFromStdFunc(const std::function<VectorN<Real, N>(const VectorN<Real, N>&)>& inFunc) 
			: _func(inFunc) {}

		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/////////////////////////    VECTOR FUNCTION N -> M      ///////////////////////////////////
	template<int N, int M>
	class VectorFunctionNM : public IVectorFunctionNM<N, M>
	{
		VectorN<Real, M>(*_func)(const VectorN<Real, N>&);
	public:
		VectorFunctionNM(VectorN<Real, M>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	template<int N, int M>
	class VectorFunctionNMFromStdFunc : public IVectorFunctionNM<N, M>
	{
		std::function<VectorN<Real, M>(const VectorN<Real, N>&)> _func;
	public:
		VectorFunctionNMFromStdFunc(const std::function<VectorN<Real, M>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	//////////////////////     PARAMETRIC CURVE             ///////////////////////////////////
	template<int N>
	class   ParametricCurve : public IParametricCurve<N>
	{
		Real _minT;
		Real _maxT;
		VectorN<Real, N>(*_func)(Real);
	public:
		ParametricCurve(VectorN<Real, N>(*inFunc)(Real)) 
			: _func(inFunc), _minT(Constants::NegInf), _maxT(Constants::PosInf) {}
		ParametricCurve(Real minT, Real maxT, VectorN<Real, N>(*inFunc)(Real)) 
			: _func(inFunc), _minT(minT), _maxT(maxT) {}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		virtual VectorN<Real, N> operator()(Real x) const { return _func(x); }
	};

	template<int N>
	class ParametricCurveFromStdFunc : public IParametricCurve<N>
	{
		Real _minT;
		Real _maxT;
		std::function<VectorN<Real, N>(Real)> _func;
	public:
		ParametricCurveFromStdFunc(const std::function<VectorN<Real, N>(Real)>& inFunc) 
			: _func(inFunc), _minT(Constants::NegInf), _maxT(Constants::PosInf) {}
		ParametricCurveFromStdFunc(Real minT, Real maxT, const std::function<VectorN<Real, N>(Real)>& inFunc) 
			: _func(inFunc), _minT(minT), _maxT(maxT) {}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		VectorN<Real, N> operator()(Real x) const { return _func(x); }
	};

	/////////////////////       PARAMETRIC SURFACE         //////////////////////////////////
	template<int N>
	class ParametricSurface : public IParametricSurface<N>
	{
		static_assert(N >= 3, "ParametricSurface requires N >= 3");

		Real _minX, _maxX;
		Real(*_y1)(Real);			// lower y boundary as function of x
		Real(*_y2)(Real);			// upper y boundary as function of x
		VectorN<Real, N>(*_func)(Real x, Real y);

	public:
		ParametricSurface(VectorN<Real, N>(*inFunc)(Real x, Real y),
											Real minX, Real maxX,
											Real(*y1)(Real), Real(*y2)(Real))
			: _func(inFunc), _minX(minX), _maxX(maxX), _y1(y1), _y2(y2)
		{	}

		VectorN<Real, N> operator()(Real x, Real y) const override
		{
			if (x < _minX || x > _maxX)
				throw std::domain_error("ParametricSurface: x out of domain");
			Real y_min = _y1(x), y_max = _y2(x);
			if (y < y_min || y > y_max)
				throw std::domain_error("ParametricSurface: y out of domain for given x");
			return _func(x, y);
		}

		virtual Real getMinU() const override { return _minX; }
		virtual Real getMaxU() const override { return _maxX; }
		virtual Real getMinW(Real x) const override { return _y1(x); }
		virtual Real getMaxW(Real x) const override { return _y2(x); }
	};
	
	template<int N>
	class ParametricSurfaceRect : public IParametricSurfaceRect<N>
	{
		static_assert(N >= 3, "ParametricSurfaceRect requires N >= 3");

		Real _minX;
		Real _maxX;
		Real _minY;
		Real _maxY;
		VectorN<Real, N>(*_func)(Real u, Real w);

	public:
		ParametricSurfaceRect(VectorN<Real, N>(*inFunc)(Real u, Real w)) 
				: _func(inFunc), _minX(Constants::NegInf), _maxX(Constants::PosInf), 
					_minY(Constants::NegInf), _maxY(Constants::PosInf) {}
		ParametricSurfaceRect(VectorN<Real, N>(*inFunc)(Real u, Real w), 
													Real minX, Real maxX, Real minY, Real maxY) 
				: _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

		VectorN<Real, N> operator()(Real u, Real w) const override { return _func(u, w); }

		virtual Real getMinU() const override { return _minX; }
		virtual Real getMaxU() const override { return _maxX; }
		virtual Real getMinW() const override { return _minY; }
		virtual Real getMaxW() const override { return _maxY; }
	};

	// imati cemo i surface Discrete, kreiran od triangles?

	template<int N>
	class ParametricSurfaceFromStdFunc : public IParametricSurfaceRect<N>
	{
		static_assert(N >= 3, "ParametricSurfaceFromStdFunc requires N >= 3");

		Real _minX;
		Real _maxX;
		Real _minY;
		Real _maxY;
		std::function<VectorN<Real, N>(Real u, Real w)> _func;
	public:
		ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc) : _func(inFunc), _minX(Constants::NegInf), _maxX(Constants::PosInf), _minY(Constants::NegInf), _maxY(Constants::PosInf) {}
		ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc, Real minX, Real maxX, Real minY, Real maxY) : _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

		VectorN<Real, N> operator()(Real u, Real w) const { return _func(u, w); }

		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW() const { return _minY; }
		virtual Real getMaxW() const { return _maxY; }
	};
} // end namespace

#endif