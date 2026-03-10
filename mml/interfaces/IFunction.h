///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IFunction.h                                                         ///
///  Description: Function interface hierarchy (IRealFunction, IVectorFunction, etc.) ///
///               Abstract base classes for all function types in MML                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IFunction.h
 * @brief Core function interface hierarchy for mathematical functions.
 * 
 * This file defines the complete hierarchy of function interfaces used throughout MML:
 * 
 * **Base Template:**
 * - IFunction<RetType, ArgType> - Generic function interface
 * 
 * **Real Functions (R → R):**
 * - IRealFunction - Basic real-to-real function
 * - IRealFunctionParametrized - Real function with adjustable parameters
 * 
 * **Scalar Functions (R^N → R):**
 * - IScalarFunction<N> - Multi-variable scalar function
 * - IScalarFunctionParametrized<N> - Parametrized scalar function
 * 
 * **Vector Functions:**
 * - IRealToVectorFunction<N> - R → R^N mapping
 * - IVectorFunction<N> - R^N → R^N mapping
 * - IVectorFunctionParametrized<N> - Parametrized vector function
 * - IVectorFunctionNM<N,M> - R^N → R^M mapping
 * 
 * **Parametric Geometry:**
 * - IParametricCurve<N> - Curve in N-dimensional space
 * - IParametricCurveParametrized<N> - Parametrized curve
 * - IParametricSurface<N> - Surface with variable domain
 * - IParametricSurfaceRect<N> - Surface on rectangular domain
 */

#if !defined  MML_IFUNCTION_H
#define MML_IFUNCTION_H

#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "interfaces/IParametrized.h"

#include <vector>

namespace MML
{
	/**
	 * @brief Generic function interface template.
	 * 
	 * Base template for all function types in MML. Defines a callable
	 * interface mapping from _ArgType to _RetType.
	 * 
	 * @tparam _RetType Return type of the function
	 * @tparam _ArgType Argument type of the function
	 */
	template<typename _RetType, typename _ArgType>
	class IFunction
	{
	public:
		/**
		 * @brief Evaluate the function at the given argument.
		 * @param arg The input argument
		 * @return The function value at the argument
		 */
		virtual _RetType operator()(_ArgType) const = 0;

		virtual ~IFunction() {}
	};

	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for real-valued functions of a single real variable.
	 * 
	 * Represents functions f: R → R. This is the most fundamental function
	 * type, used extensively in root finding, integration, and differentiation.
	 * 
	 * @note Implementations should ensure operator() is const and thread-safe
	 *       for parallel evaluation.
	 */
	class IRealFunction : public IFunction<Real, Real>
	{
	public:
		/**
		 * @brief Evaluate the function at point x.
		 * @param x The input value
		 * @return f(x)
		 */
		virtual Real operator()(Real) const = 0;
		
		/**
		 * @brief Sample the function over a range of values.
		 * 
		 * Generates numPnt equally-spaced sample points and their function values,
		 * useful for plotting or numerical analysis.
		 * 
		 * @param x1 Start of sampling range
		 * @param x2 End of sampling range
		 * @param numPnt Number of sample points
		 * @param[out] outX Vector of x coordinates
		 * @param[out] outY Vector of corresponding y = f(x) values
		 */
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

	/**
	 * @brief Real function with adjustable parameters.
	 * 
	 * Extends IRealFunction to support functions with internal parameters
	 * that can be get/set at runtime. Useful for curve fitting, optimization,
	 * and parametric studies.
	 * 
	 * Example: f(x; a, b) = a*sin(b*x) where a, b are parameters.
	 */
	class IRealFunctionParametrized : public IRealFunction, public IParametrized
	{
	public:
		virtual ~IRealFunctionParametrized() = default;
	};

	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for scalar functions of N real variables.
	 * 
	 * Represents functions f: R^N → R. Used in optimization, gradient
	 * computation, and multi-variable calculus.
	 * 
	 * @tparam N Dimension of the input space
	 */
	template<int N>
	class IScalarFunction : public IFunction<Real, const VectorN<Real, N>&>
	{
	public:
		/**
		 * @brief Evaluate the function at point x.
		 * @param x N-dimensional input vector
		 * @return f(x)
		 */
		virtual Real operator()(const VectorN<Real, N>& x) const = 0;

		virtual ~IScalarFunction() {}
	};

	/**
	 * @brief Scalar function with adjustable parameters.
	 * 
	 * @tparam N Dimension of the input space
	 */
	template<int N>
	class IScalarFunctionParametrized : public IScalarFunction<N>, public IParametrized
	{
	public:
		virtual ~IScalarFunctionParametrized() = default;
	};

	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for vector-valued functions of a single real variable.
	 * 
	 * Represents functions f: R → R^N. Commonly used for parametric curves
	 * and time-dependent vector quantities.
	 * 
	 * @tparam N Dimension of the output vector
	 */
	template<int N>
	class IRealToVectorFunction : public IFunction<VectorN<Real, N>, Real>
	{
	public:
		/**
		 * @brief Evaluate the function at parameter x.
		 * @param x The input parameter
		 * @return N-dimensional output vector
		 */
		virtual VectorN<Real, N> operator()(Real x) const = 0;

		virtual ~IRealToVectorFunction() {}
	};

	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for vector-valued functions of N real variables.
	 * 
	 * Represents functions f: R^N → R^N (same input/output dimension).
	 * Used for vector fields, coordinate transformations, and systems
	 * of equations.
	 * 
	 * @tparam N Dimension of both input and output spaces
	 */
	template<int N>
	class IVectorFunction : public IFunction<VectorN<Real, N>, const VectorN<Real, N>&>
	{
	public:
		/**
		 * @brief Evaluate the function at point x.
		 * @param x N-dimensional input vector
		 * @return N-dimensional output vector
		 */
		virtual VectorN<Real, N> operator()(const VectorN<Real, N>& x) const = 0;

		/**
		 * @brief Evaluate a single component of the output.
		 * @param x N-dimensional input vector
		 * @param component Index of the output component (0 to N-1)
		 * @return The specified component of f(x)
		 * 
		 * @note This method delegates to evaluateComponent() which can be
		 *       overridden for efficient per-component evaluation. The default
		 *       implementation computes the full vector, which may be O(N).
		 */
		Real operator()(const VectorN<Real, N>& x, int component) const
		{
			return evaluateComponent(x, component);
		}

		virtual ~IVectorFunction() {}

	protected:
		/**
		 * @brief Override point for efficient single-component evaluation.
		 * 
		 * Derived classes can override this method to provide O(1) per-component
		 * evaluation instead of the default O(N) full-vector computation.
		 * This is particularly beneficial when:
		 * - Computing Jacobians (avoids O(N²) becoming O(N))
		 * - Optimization algorithms need individual gradient components
		 * - Large dimensions (N >> 1)
		 * 
		 * @param x N-dimensional input vector
		 * @param component Index of the output component (0 to N-1)
		 * @return The specified component of f(x)
		 */
		virtual Real evaluateComponent(const VectorN<Real, N>& x, int component) const
		{
			VectorN<Real, N> val = (*this)(x);
			return val[component];
		}
	};

	/**
	 * @brief Vector function with adjustable parameters.
	 * 
	 * @tparam N Dimension of input and output spaces
	 */
	template<int N>
	class IVectorFunctionParametrized : public IVectorFunction<N>, public IParametrized
	{
	public:
		virtual ~IVectorFunctionParametrized() = default;
	};

	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for vector-valued functions with different input/output dimensions.
	 * 
	 * Represents functions f: R^N → R^M where N ≠ M in general.
	 * Useful for projections, embeddings, and general mappings.
	 * 
	 * @tparam N Dimension of the input space
	 * @tparam M Dimension of the output space
	 */
	template<int N, int M>
	class IVectorFunctionNM : public IFunction<VectorN<Real, M>, const VectorN<Real, N>&>
	{
	public:
		/**
		 * @brief Evaluate the function at point x.
		 * @param x N-dimensional input vector
		 * @return M-dimensional output vector
		 */
		virtual VectorN<Real, M> operator()(const VectorN<Real, N>& x) const = 0;
		
		/**
		 * @brief Evaluate a single component of the output.
		 * @param x N-dimensional input vector
		 * @param component Index of the output component (0 to M-1)
		 * @return The specified component of f(x)
		 * 
		 * @note This method delegates to evaluateComponent() which can be
		 *       overridden for efficient per-component evaluation. The default
		 *       implementation computes the full vector, which may be O(M).
		 */
		Real operator()(const VectorN<Real, N>& x, int component) const
		{
			return evaluateComponent(x, component);
		}
		
		virtual ~IVectorFunctionNM() {}

	protected:
		/**
		 * @brief Override point for efficient single-component evaluation.
		 * 
		 * Derived classes can override this method to provide O(1) per-component
		 * evaluation instead of the default O(M) full-vector computation.
		 * This is particularly beneficial when:
		 * - Computing Jacobians (avoids O(M²) becoming O(M))
		 * - Optimization algorithms need individual gradient components
		 * - Large output dimensions (M >> 1)
		 * 
		 * @param x N-dimensional input vector
		 * @param component Index of the output component (0 to M-1)
		 * @return The specified component of f(x)
		 */
		virtual Real evaluateComponent(const VectorN<Real, N>& x, int component) const
		{
			VectorN<Real, M> val = (*this)(x);
			return val[component];
		}
	};

	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for parametric curves in N-dimensional space.
	 * 
	 * Represents a curve γ: [tMin, tMax] → R^N. The parameter t typically
	 * represents arc length or time. Used for path integrals, trajectory
	 * analysis, and geometric modeling.
	 * 
	 * @tparam N Dimension of the ambient space
	 */
	template<int N>
	class IParametricCurve : public IRealToVectorFunction<N>
	{
	public:
		/** @brief Get the minimum parameter value. */
		virtual Real getMinT() const = 0;
		
		/** @brief Get the maximum parameter value. */
		virtual Real getMaxT() const = 0;

		/**
		 * @brief Generate points along the curve for visualization.
		 * @param t1 Start parameter
		 * @param t2 End parameter
		 * @param numPoints Number of sample points
		 * @return Vector of points on the curve
		 */
		std::vector<VectorN<Real, N>> GetTrace(double t1, double t2, int numPoints) const
		{
			std::vector<VectorN<Real, N>> ret;
			if (numPoints <= 0)
				return ret;
			ret.reserve(static_cast<size_t>(numPoints));

			const Real tStart = static_cast<Real>(t1);
			if (numPoints == 1)
			{
				ret.push_back((*this)(tStart));
				return ret;
			}

			const Real tEnd = static_cast<Real>(t2);
			const Real deltaT = (tEnd - tStart) / (numPoints - 1);
			for (int i = 0; i < numPoints; i++)
			{
				const Real t = tStart + static_cast<Real>(i) * deltaT;
				ret.push_back((*this)(t));
			}
			return ret;
		}

		virtual ~IParametricCurve() {}
	};

	/**
	 * @brief Parametric curve with adjustable parameters.
	 * 
	 * @tparam N Dimension of the ambient space
	 */
	template<int N>
	class IParametricCurveParametrized : public IParametricCurve<N>, public IParametrized
	{
	public:
			/**
			 * @brief Evaluate the curve with explicit parameters.
			 * @param t Curve parameter
			 * @param params External parameter vector
			 * @return Point on the curve
			 */
			virtual VectorN<Real, N> operator()(Real t, const Vector<Real> &params) const = 0;
	};

	//////////////////////////////////////////////////////////////////////
	/**
	 * @brief Interface for parametric surfaces with variable domain.
	 * 
	 * Represents a surface S: D → R^N where D is a domain in R^2 with
	 * u ∈ [uMin, uMax] and w ∈ [wMin(u), wMax(u)]. The w-limits can
	 * depend on u, allowing non-rectangular domains (e.g., triangular,
	 * circular patches).
	 * 
	 * @tparam N Dimension of the ambient space (typically 3)
	 */
	template<int N>
	class IParametricSurface : public IVectorFunctionNM<2, N>
	{
	public:
		/**
		 * @brief Evaluate the surface at parameter (u, w).
		 * @param u First surface parameter
		 * @param w Second surface parameter
		 * @return Point on the surface
		 */
		virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

		/** @brief Get minimum u parameter. */
		virtual Real getMinU() const = 0;
		
		/** @brief Get maximum u parameter. */
		virtual Real getMaxU() const = 0;
		
		/**
		 * @brief Get minimum w parameter at given u.
		 * @param u The u-coordinate
		 * @return Minimum w value at this u
		 */
		virtual Real getMinW(Real u) const = 0;
		
		/**
		 * @brief Get maximum w parameter at given u.
		 * @param u The u-coordinate
		 * @return Maximum w value at this u
		 */
		virtual Real getMaxW(Real u) const = 0;

		/**
		 * @brief Evaluate using a 2D coordinate vector.
		 * @param coord Vector containing (u, w)
		 * @return Point on the surface
		 */
		virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
		{
			return operator()(coord[0], coord[1]);
		}

		virtual ~IParametricSurface() {}
	};

	/**
	 * @brief Interface for parametric surfaces on rectangular domains.
	 * 
	 * Simplified surface interface where the domain is a rectangle:
	 * u ∈ [uMin, uMax], w ∈ [wMin, wMax] with constant bounds.
	 * 
	 * @tparam N Dimension of the ambient space (typically 3)
	 */
	template<int N>
	class IParametricSurfaceRect : public IVectorFunctionNM<2, N>
	{
	public:
		/**
		 * @brief Evaluate the surface at parameter (u, w).
		 * @param u First surface parameter
		 * @param w Second surface parameter
		 * @return Point on the surface
		 */
		virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

		/** @brief Get minimum u parameter. */
		virtual Real getMinU() const = 0;
		
		/** @brief Get maximum u parameter. */
		virtual Real getMaxU() const = 0;
		
		/** @brief Get minimum w parameter (constant). */
		virtual Real getMinW() const = 0;
		
		/** @brief Get maximum w parameter (constant). */
		virtual Real getMaxW() const = 0;

		/**
		 * @brief Evaluate using a 2D coordinate vector.
		 * @param coord Vector containing (u, w)
		 * @return Point on the surface
		 */
		virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
		{
			return operator()(coord[0], coord[1]);
		}

		virtual ~IParametricSurfaceRect() {}
	};
}
#endif