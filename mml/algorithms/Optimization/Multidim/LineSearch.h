///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Multidim/LineSearch.h                                                 ///
///  Description: Line-function wrappers and line minimization helper                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_LINE_SEARCH_H
#define MML_OPTIMIZATION_LINE_SEARCH_H

#include <mml/algorithms/Optimization/Multidim/MultidimTypes.h>
#include <mml/algorithms/Optimization/Optimization.h>

namespace MML::Optimization {

	///                     LINE SEARCH METHODS                       ///
	/////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	///            Interface for differentiable scalar functions            ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Interface for N-dimensional scalar function with gradient
     * 
     * Functions that provide gradient information can use more efficient
     * optimization methods like conjugate gradient or BFGS.
     */
	template<int N>
	class IDifferentiableScalarFunction : public IScalarFunction<N> {
	public:
		/**
         * @brief Compute the gradient at point x
         * @param x Point at which to evaluate gradient
         * @param grad Output gradient vector (filled by this method)
         */
		virtual void Gradient(const VectorN<Real, N>& x, VectorN<Real, N>& grad) const = 0;

		virtual ~IDifferentiableScalarFunction() = default;
	};

	///////////////////////////////////////////////////////////////////////////
	///                    1D function along a line                         ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Wrapper to convert N-dim function to 1D function along a line
     * 
     * Given f(x) in N dimensions, creates g(t) = f(p + t*xi) for line search.
     * This is used by all line-search based optimizers.
     */
	template<int N>
	class LineFunction : public IRealFunction {
	private:
		const IScalarFunction<N>& _func;
		const VectorN<Real, N>& _p;	 // Base point
		const VectorN<Real, N>& _xi; // Direction

	public:
		LineFunction(const IScalarFunction<N>& func, const VectorN<Real, N>& p, const VectorN<Real, N>& xi)
			: _func(func)
			, _p(p)
			, _xi(xi) {}

		Real operator()(Real t) const override {
			VectorN<Real, N> xt;
			for (int j = 0; j < N; ++j)
				xt[j] = _p[j] + t * _xi[j];
			return _func(xt);
		}
	};

	/**
     * @brief 1D function with derivative along a line (for gradient-based methods)
     */
	template<int N>
	class DLineFunction : public IRealFunction {
	private:
		const IDifferentiableScalarFunction<N>& _func;
		VectorN<Real, N> _p;			// Base point (mutable for evaluation)
		VectorN<Real, N> _xi;			// Direction
		mutable VectorN<Real, N> _xt;	// Current point
		mutable VectorN<Real, N> _grad; // Gradient at current point

	public:
		DLineFunction(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& p, const VectorN<Real, N>& xi)
			: _func(func)
			, _p(p)
			, _xi(xi) {}

		void updateBasePoint(const VectorN<Real, N>& p) { _p = p; }
		void updateDirection(const VectorN<Real, N>& xi) { _xi = xi; }

		Real operator()(Real t) const override {
			for (int j = 0; j < N; ++j)
				_xt[j] = _p[j] + t * _xi[j];
			return _func(_xt);
		}

		/**
         * @brief Compute directional derivative df/dt = grad(f) · xi
         */
		Real derivative(Real t) const {
			for (int j = 0; j < N; ++j)
				_xt[j] = _p[j] + t * _xi[j];
			_func.Gradient(_xt, _grad);

			Real df = 0.0;
			for (int j = 0; j < N; ++j)
				df += _grad[j] * _xi[j];
			return df;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                      Line Minimization                              ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Line minimization helper for N-dimensional optimization
     * 
     * Given a point p and direction xi, finds the minimum along the line
     * p + t*xi using Brent's method, then updates p and xi.
     */
	class LineMinimizer {
	public:
		/**
         * @brief Minimize function along line p + t*xi
         * @tparam N Problem dimension
         * @param func Function to minimize
         * @param p Current point (updated to minimum along line)
         * @param xi Search direction (updated to displacement vector)
         * @return Function value at minimum
         */
		template<int N>
		static Real Minimize(const IScalarFunction<N>& func, VectorN<Real, N>& p, VectorN<Real, N>& xi, Real tol = PrecisionValues<Real>::OptimizationTolerance) {
			LineFunction<N> f1dim(func, p, xi);

			// Bracket the minimum
			auto bracket = Minimization::BracketMinimum(f1dim, 0.0, 1.0);
			if (!bracket.valid) {
				// Try with different initial interval
				bracket = Minimization::BracketMinimum(f1dim, 0.0, 0.1);
			}

			// Find minimum using Brent's method
			auto result = Minimization::BrentMinimize(f1dim, bracket, tol);
			Real xmin = result.xmin;

			// Update p and xi
			for (int j = 0; j < N; ++j) {
				xi[j] *= xmin;
				p[j] += xi[j];
			}

			return result.fmin;
		}
	};


} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_LINE_SEARCH_H
