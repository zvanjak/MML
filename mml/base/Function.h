///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Function.h                                                          ///
///  Description: Function wrappers for scalar and vector functions                   ///
///               Lambda adapters, composition, and evaluation utilities              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_FUNCTION_H
#define MML_FUNCTION_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "interfaces/IFunction.h"

#include "base/Vector/VectorN.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

#include <functional>

namespace MML
{
	///////////////////////////     REAL FUNCTION      ////////////////////////////////////

	/// @brief Wrapper for a real function R → R using function pointer.
	/// @details Wraps a C-style function pointer for use with numerical algorithms.
	class RealFunction : public IRealFunction
	{
		Real(*_func)(const Real);   ///< Function pointer to evaluate
	public:
		/// @brief Constructs wrapper from function pointer.
		/// @param inFunc Function pointer Real → Real
		RealFunction(Real(*inFunc)(const Real)) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param x Input value
		/// @return f(x)
		Real operator()(const Real x) const { return _func(x); }
	};

	/// @brief Wrapper for a real function R → R using std::function.
	/// @details Allows using lambdas, functors, and bound functions.
	/// @warning **Lambda capture lifetime hazard:** The std::function is copied into this
	///          wrapper, but if the source lambda captures local variables **by reference**,
	///          those references will dangle if this object outlives the captured variables.
	///          Safe:   `RealFunctionFromStdFunc f([](Real x) { return x*x; });`
	///          Safe:   `RealFunctionFromStdFunc f([a](Real x) { return a*x; });`  (capture by value)
	///          UNSAFE: `RealFunctionFromStdFunc f([&a](Real x) { return a*x; });` (if a goes out of scope)
	class RealFunctionFromStdFunc : public IRealFunction
	{
		std::function<Real(const Real)> _func;   ///< Callable object to evaluate
	public:
		/// @brief Constructs wrapper from std::function.
		/// @param inFunc Callable object Real → Real (lambda, functor, etc.)
		RealFunctionFromStdFunc(const std::function<Real(const Real)> &inFunc) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param x Input value
		/// @return f(x)
		Real operator()(const Real x) const { return _func(x); }
	};

	///////////////////////////     SCALAR FUNCTION       //////////////////////////////////

	/// @brief Wrapper for a scalar function R^N → R using function pointer.
	/// @tparam N Dimension of input space
	template<int N>
	class ScalarFunction : public IScalarFunction<N>
	{
		Real(*_func)(const VectorN<Real, N>&);   ///< Function pointer to evaluate
	public:
		/// @brief Constructs wrapper from function pointer.
		/// @param inFunc Function pointer VectorN<N> → Real
		ScalarFunction(Real(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param x Input vector
		/// @return f(x)
		Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/// @brief Wrapper for a scalar function R^N → R using std::function.
	/// @tparam N Dimension of input space
	/// @warning **Lambda capture lifetime hazard:** If constructed from a lambda capturing
	///          local variables by reference, ensure those variables outlive this object.
	template<int N>
	class ScalarFunctionFromStdFunc : public IScalarFunction<N>
	{
		std::function<Real(const VectorN<Real, N>&)> _func;   ///< Callable object
	public:
		/// @brief Constructs wrapper from std::function.
		/// @param inFunc Callable object VectorN<N> → Real
		ScalarFunctionFromStdFunc(const std::function<Real(const VectorN<Real, N>&)> inFunc) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param x Input vector
		/// @return f(x)
		Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/////////////////////////    VECTOR FUNCTION N -> N      ///////////////////////////////////

	/// @brief Wrapper for a vector function R^N → R^N using function pointer.
	/// @tparam N Dimension of input and output space
	template<int N>
	class VectorFunction : public IVectorFunction<N>
	{
		VectorN<Real, N>(*_func)(const VectorN<Real, N>&);   ///< Function pointer
	public:
		/// @brief Constructs wrapper from function pointer.
		/// @param inFunc Function pointer VectorN<N> → VectorN<N>
		VectorFunction(VectorN<Real, N>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param x Input vector
		/// @return f(x)
		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/// @brief Wrapper for a vector function R^N → R^N using std::function.
	/// @tparam N Dimension of input and output space
	/// @warning **Lambda capture lifetime hazard:** If constructed from a lambda capturing
	///          local variables by reference, ensure those variables outlive this object.
	template<int N>
	class VectorFunctionFromStdFunc : public IVectorFunction<N>
	{
		std::function<VectorN<Real, N>(const VectorN<Real, N>&)> _func;   ///< Callable object
	public:
		/// @brief Constructs wrapper from std::function.
		/// @param inFunc Callable object VectorN<N> → VectorN<N>
		VectorFunctionFromStdFunc(const std::function<VectorN<Real, N>(const VectorN<Real, N>&)>& inFunc) 
			: _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param x Input vector
		/// @return f(x)
		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/////////////////////////    VECTOR FUNCTION N -> M      ///////////////////////////////////

	/// @brief Wrapper for a vector function R^N → R^M using function pointer.
	/// @tparam N Dimension of input space
	/// @tparam M Dimension of output space
	template<int N, int M>
	class VectorFunctionNM : public IVectorFunctionNM<N, M>
	{
		VectorN<Real, M>(*_func)(const VectorN<Real, N>&);   ///< Function pointer
	public:
		/// @brief Constructs wrapper from function pointer.
		/// @param inFunc Function pointer VectorN<N> → VectorN<M>
		VectorFunctionNM(VectorN<Real, M>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param x Input vector
		/// @return f(x)
		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/// @brief Wrapper for a vector function R^N → R^M using std::function.
	/// @tparam N Dimension of input space
	/// @tparam M Dimension of output space
	/// @warning **Lambda capture lifetime hazard:** If constructed from a lambda capturing
	///          local variables by reference, ensure those variables outlive this object.
	template<int N, int M>
	class VectorFunctionNMFromStdFunc : public IVectorFunctionNM<N, M>
	{
		std::function<VectorN<Real, M>(const VectorN<Real, N>&)> _func;   ///< Callable object
	public:
		/// @brief Constructs wrapper from std::function.
		/// @param inFunc Callable object VectorN<N> → VectorN<M>
		VectorFunctionNMFromStdFunc(const std::function<VectorN<Real, M>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

		/// @brief Evaluates the wrapped function.
		/// @param x Input vector
		/// @return f(x)
		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	//////////////////////     PARAMETRIC CURVE             ///////////////////////////////////

	/// @brief Parametric curve γ: R → R^N using function pointer.
	/// @details Represents a path in N-dimensional space parameterized by t.
	/// @tparam N Dimension of output space (typically 2 or 3)
	template<int N>
	class   ParametricCurve : public IParametricCurve<N>
	{
		Real _minT;   ///< Minimum parameter value
		Real _maxT;   ///< Maximum parameter value
		VectorN<Real, N>(*_func)(Real);   ///< Function pointer t → point
	public:
		/// @brief Constructs unbounded parametric curve.
		/// @param inFunc Function pointer Real → VectorN<N>
		ParametricCurve(VectorN<Real, N>(*inFunc)(Real)) 
			: _func(inFunc), _minT(Constants::NegInf), _maxT(Constants::PosInf) {}

		/// @brief Constructs bounded parametric curve.
		/// @param minT Minimum parameter value
		/// @param maxT Maximum parameter value
		/// @param inFunc Function pointer Real → VectorN<N>
		ParametricCurve(Real minT, Real maxT, VectorN<Real, N>(*inFunc)(Real)) 
			: _func(inFunc), _minT(minT), _maxT(maxT) {}

		/// @brief Gets minimum parameter value.
		/// @return Minimum t
		Real getMinT() const { return _minT; }

		/// @brief Gets maximum parameter value.
		/// @return Maximum t
		Real getMaxT() const { return _maxT; }

		/// @brief Evaluates curve at parameter t.
		/// @param x Parameter value
		/// @return Point on curve γ(t)
		virtual VectorN<Real, N> operator()(Real x) const { return _func(x); }
	};

	/// @brief Parametric curve γ: R → R^N using std::function.
	/// @tparam N Dimension of output space
	/// @warning **Lambda capture lifetime hazard:** If constructed from a lambda capturing
	///          local variables by reference, ensure those variables outlive this object.
	template<int N>
	class ParametricCurveFromStdFunc : public IParametricCurve<N>
	{
		Real _minT;   ///< Minimum parameter value
		Real _maxT;   ///< Maximum parameter value
		std::function<VectorN<Real, N>(Real)> _func;   ///< Callable object
	public:
		/// @brief Constructs unbounded parametric curve.
		/// @param inFunc Callable object Real → VectorN<N>
		ParametricCurveFromStdFunc(const std::function<VectorN<Real, N>(Real)>& inFunc) 
			: _func(inFunc), _minT(Constants::NegInf), _maxT(Constants::PosInf) {}

		/// @brief Constructs bounded parametric curve.
		/// @param minT Minimum parameter value
		/// @param maxT Maximum parameter value
		/// @param inFunc Callable object Real → VectorN<N>
		ParametricCurveFromStdFunc(Real minT, Real maxT, const std::function<VectorN<Real, N>(Real)>& inFunc) 
			: _func(inFunc), _minT(minT), _maxT(maxT) {}

		/// @brief Gets minimum parameter value.
		Real getMinT() const { return _minT; }

		/// @brief Gets maximum parameter value.
		Real getMaxT() const { return _maxT; }

		/// @brief Evaluates curve at parameter t.
		/// @param x Parameter value
		/// @return Point on curve γ(t)
		VectorN<Real, N> operator()(Real x) const { return _func(x); }
	};

	/////////////////////       PARAMETRIC SURFACE         //////////////////////////////////

	/// @brief Parametric surface with general (non-rectangular) domain.
	/// @details Domain is defined by x ∈ [minX, maxX] and y ∈ [y1(x), y2(x)].
	/// @tparam N Dimension of output space (must be >= 3)
	template<int N>
	class ParametricSurface : public IParametricSurface<N>
	{
		static_assert(N >= 3, "ParametricSurface requires N >= 3");

		Real _minX, _maxX;                      ///< X-parameter bounds
		Real(*_y1)(Real);                       ///< Lower y boundary y1(x)
		Real(*_y2)(Real);                       ///< Upper y boundary y2(x)
		VectorN<Real, N>(*_func)(Real x, Real y);   ///< Surface function

	public:
		/// @brief Constructs parametric surface with general domain.
		/// @param inFunc Surface function (x,y) → point
		/// @param minX Minimum x parameter
		/// @param maxX Maximum x parameter
		/// @param y1 Lower y boundary as function of x
		/// @param y2 Upper y boundary as function of x
		ParametricSurface(VectorN<Real, N>(*inFunc)(Real x, Real y),
											Real minX, Real maxX,
											Real(*y1)(Real), Real(*y2)(Real))
			: _func(inFunc), _minX(minX), _maxX(maxX), _y1(y1), _y2(y2)
		{	}

		/// @brief Evaluates surface at (x, y).
		/// @param x First parameter
		/// @param y Second parameter
		/// @return Point on surface σ(x, y)
		/// @throws DomainError if parameters are out of domain
		VectorN<Real, N> operator()(Real x, Real y) const override
		{
			if (x < _minX || x > _maxX)
				throw DomainError("ParametricSurface: x out of domain");
			Real y_min = _y1(x), y_max = _y2(x);
			if (y < y_min || y > y_max)
				throw DomainError("ParametricSurface: y out of domain for given x");
			return _func(x, y);
		}

		/// @brief Gets minimum x parameter.
		virtual Real getMinU() const override { return _minX; }
		/// @brief Gets maximum x parameter.
		virtual Real getMaxU() const override { return _maxX; }
		/// @brief Gets minimum y parameter for given x.
		virtual Real getMinW(Real x) const override { return _y1(x); }
		/// @brief Gets maximum y parameter for given x.
		virtual Real getMaxW(Real x) const override { return _y2(x); }
	};
	
	/// @brief Parametric surface with rectangular domain.
	/// @details Domain is x ∈ [minX, maxX], y ∈ [minY, maxY].
	/// @tparam N Dimension of output space (must be >= 3)
	template<int N>
	class ParametricSurfaceRect : public IParametricSurfaceRect<N>
	{
		static_assert(N >= 3, "ParametricSurfaceRect requires N >= 3");

		Real _minX;   ///< Minimum x parameter
		Real _maxX;   ///< Maximum x parameter
		Real _minY;   ///< Minimum y parameter
		Real _maxY;   ///< Maximum y parameter
		VectorN<Real, N>(*_func)(Real u, Real w);   ///< Surface function

	public:
		/// @brief Constructs unbounded rectangular parametric surface.
		/// @param inFunc Surface function (u,w) → point
		ParametricSurfaceRect(VectorN<Real, N>(*inFunc)(Real u, Real w)) 
				: _func(inFunc), _minX(Constants::NegInf), _maxX(Constants::PosInf), 
					_minY(Constants::NegInf), _maxY(Constants::PosInf) {}

		/// @brief Constructs bounded rectangular parametric surface.
		/// @param inFunc Surface function (u,w) → point
		/// @param minX Minimum u parameter
		/// @param maxX Maximum u parameter
		/// @param minY Minimum w parameter
		/// @param maxY Maximum w parameter
		ParametricSurfaceRect(VectorN<Real, N>(*inFunc)(Real u, Real w), 
													Real minX, Real maxX, Real minY, Real maxY) 
				: _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

		/// @brief Evaluates surface at (u, w).
		/// @param u First parameter
		/// @param w Second parameter
		/// @return Point on surface σ(u, w)
		/// @brief Evaluates surface at (u, w).
		/// @param u First parameter
		/// @param w Second parameter
		/// @return Point on surface σ(u, w)
		VectorN<Real, N> operator()(Real u, Real w) const override { return _func(u, w); }

		/// @brief Gets minimum u parameter.
		virtual Real getMinU() const override { return _minX; }
		/// @brief Gets maximum u parameter.
		virtual Real getMaxU() const override { return _maxX; }
		/// @brief Gets minimum w parameter.
		virtual Real getMinW() const override { return _minY; }
		/// @brief Gets maximum w parameter.
		virtual Real getMaxW() const override { return _maxY; }
	};

	/// @brief Parametric surface with rectangular domain using std::function.
	/// @tparam N Dimension of output space (must be >= 3)
	/// @warning **Lambda capture lifetime hazard:** If constructed from a lambda capturing
	///          local variables by reference, ensure those variables outlive this object.
	template<int N>
	class ParametricSurfaceFromStdFunc : public IParametricSurfaceRect<N>
	{
		static_assert(N >= 3, "ParametricSurfaceFromStdFunc requires N >= 3");

		Real _minX;   ///< Minimum u parameter
		Real _maxX;   ///< Maximum u parameter
		Real _minY;   ///< Minimum w parameter
		Real _maxY;   ///< Maximum w parameter
		std::function<VectorN<Real, N>(Real u, Real w)> _func;   ///< Callable surface function
	public:
		/// @brief Constructs unbounded surface from std::function.
		/// @param inFunc Callable object (u,w) → VectorN<N>
		ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc) : _func(inFunc), _minX(Constants::NegInf), _maxX(Constants::PosInf), _minY(Constants::NegInf), _maxY(Constants::PosInf) {}

		/// @brief Constructs bounded surface from std::function.
		/// @param inFunc Callable object (u,w) → VectorN<N>
		/// @param minX Minimum u parameter
		/// @param maxX Maximum u parameter
		/// @param minY Minimum w parameter
		/// @param maxY Maximum w parameter
		ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc, Real minX, Real maxX, Real minY, Real maxY) : _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

		/// @brief Evaluates surface at (u, w).
		/// @param u First parameter
		/// @param w Second parameter
		/// @return Point on surface σ(u, w)
		VectorN<Real, N> operator()(Real u, Real w) const { return _func(u, w); }

		/// @brief Gets minimum u parameter.
		virtual Real getMinU() const { return _minX; }
		/// @brief Gets maximum u parameter.
		virtual Real getMaxU() const { return _maxX; }
		/// @brief Gets minimum w parameter.
		virtual Real getMinW() const { return _minY; }
		/// @brief Gets maximum w parameter.
		virtual Real getMaxW() const { return _maxY; }
	};
} // end namespace

#endif