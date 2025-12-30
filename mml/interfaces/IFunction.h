///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IFunction.h                                                         ///
///  Description: Function interface hierarchy (IRealFunction, IVectorFunction, etc.) ///
///               Abstract base classes for all function types in MML                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_IFUNCTION_H
#define MML_IFUNCTION_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorN.h"

namespace MML
{
	//////////////////////////////////////////////////////////////////////
	template<typename _RetType, typename _ArgType>
	class IFunction
	{
	public:
		virtual _RetType operator()(_ArgType) const = 0;

		virtual ~IFunction() {}
	};

	//////////////////////////////////////////////////////////////////////
	class IRealFunction : public IFunction<Real, Real>
	{
	public:
		virtual Real operator()(Real) const = 0;
		
		void GetValues(Real x1, Real x2, int numPnt, Vector<Real>& outX, Vector<Real>& outY) const
		{
			outX.Resize(numPnt);
			outY.Resize(numPnt);

			for (int i = 0; i < numPnt; i++) {
				outX[i] = x1 + i * (x2 - x1) / (numPnt - 1);
				outY[i] = (*this)(outX[i]);
			}
		}

		virtual ~IRealFunction() {}
	};

	class IRealFunctionParametrized : public IRealFunction
	{
	public:
		virtual int		getNumParam() const = 0;
		virtual Real	getParam(int i) const = 0;
		virtual void	setParam(int i, Real val) = 0;

		virtual Vector<Real>	getParams() const = 0;
		virtual void					setParams(const Vector<Real>&) = 0;
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IScalarFunction : public IFunction<Real, const VectorN<Real, N>&>
	{
	public:
		virtual Real operator()(const VectorN<Real, N>& x) const = 0;

		virtual ~IScalarFunction() {}
	};

	template<int N>
	class IScalarFunctionParametrized : public IScalarFunction<N>
	{
	public:
		virtual int		getNumParam() const = 0;
		virtual Real	getParam(int i) const = 0;
		virtual void	setParam(int i, Real val) = 0;

		virtual Vector<Real>	getParams() const = 0;
		virtual void					setParams(const Vector<Real>&) = 0;
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IRealToVectorFunction : public IFunction<VectorN<Real, N>, Real>
	{
	public:
		virtual VectorN<Real, N> operator()(Real x) const = 0;

		virtual ~IRealToVectorFunction() {}
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IVectorFunction : public IFunction<VectorN<Real, N>, const VectorN<Real, N>&>
	{
	public:
		virtual VectorN<Real, N> operator()(const VectorN<Real, N>& x) const = 0;

		virtual ~IVectorFunction() {}
	};

	template<int N>
	class IVectorFunctionParametrized : public IVectorFunction<N>
	{
	public:
		virtual int		getNumParam() const = 0;
		virtual Real	getParam(int i) const = 0;
		virtual void	setParam(int i, Real val) = 0;

		virtual Vector<Real>	getParams() const = 0;
		virtual void					setParams(const Vector<Real>&) = 0;
	};

	//////////////////////////////////////////////////////////////////////
	template<int N, int M>
	class IVectorFunctionNM : public IFunction<VectorN<Real, M>, const VectorN<Real, N>&>
	{
	public:
		virtual VectorN<Real, M> operator()(const VectorN<Real, N>& x) const = 0;
		
		virtual Real operator()(const VectorN<Real, N>& x, int component) const
		{
			VectorN<Real, M> val = (*this)(x);
			return val[component];
		}
		
		virtual ~IVectorFunctionNM() {}
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IParametricCurve : public IRealToVectorFunction<N>
	{
	public:
		virtual Real getMinT() const = 0;
		virtual Real getMaxT() const = 0;

		std::vector<VectorN<Real, N>> GetTrace(double t1, double t2, int numPoints) const
		{
			std::vector<VectorN<Real, N>> ret;
                        Real deltaT = (t2 - t1) / (numPoints - 1);
		}

		virtual ~IParametricCurve() {}
	};

	template<int N>
	class IParametricCurveParametrized : public IParametricCurve<N>
	{
	public:
			virtual int		getNumParam() const = 0;
			virtual Real	getParam(int i) const = 0;
			virtual void	setParam(int i, Real val) = 0;

			virtual VectorN<Real, N> operator()(Real t, const Vector<Real> &params) = 0;
	};

	//////////////////////////////////////////////////////////////////////
	// complex surface, with fixed u limits, but variable w limits (dependent on u)
	template<int N>
	class IParametricSurface : public IVectorFunctionNM<2, N>
	{
	public:
		virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

		virtual Real getMinU() const = 0;
		virtual Real getMaxU() const = 0;
		virtual Real getMinW(Real u) const = 0;
		virtual Real getMaxW(Real u) const = 0;

		virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
		{
			return operator()(coord[0], coord[1]);
		}

		virtual ~IParametricSurface() {}
	};

	// simple regular surface, defined on rectangular coordinate patch
	template<int N>
	class IParametricSurfaceRect : public IVectorFunctionNM<2, N>
	{
	public:
		virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

		virtual Real getMinU() const = 0;
		virtual Real getMaxU() const = 0;
		virtual Real getMinW() const = 0;
		virtual Real getMaxW() const = 0;

		virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
		{
			return operator()(coord[0], coord[1]);
		}

		virtual ~IParametricSurfaceRect() {}
	};
}
#endif