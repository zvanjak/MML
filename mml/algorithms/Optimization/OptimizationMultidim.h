///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        OptimizationMultidim.h                                              ///
///  Description: Multidimensional optimization (Nelder-Mead, Powell, Conjugate Grad) ///
///               BFGS, Levenberg-Marquardt for non-linear least squares             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_MULTIDIM_H
#define MML_OPTIMIZATION_MULTIDIM_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "interfaces/IFunction.h"
#include "base/Function.h"
#include "base/Vector.h"
#include "base/Matrix.h"

#include "algorithms/Optimization.h"

namespace MML
{
    /////////////////////////////////////////////////////////////////////
    ///                 MULTIDIMENSIONAL OPTIMIZATION                 ///
    /////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    ///                     MultidimOptimizationError                       ///
    ///////////////////////////////////////////////////////////////////////////
    class MultidimOptimizationError : public std::runtime_error
    {
    public:
        explicit MultidimOptimizationError(const std::string& message)
            : std::runtime_error("MultidimOptimizationError: " + message) {}
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                  MultidimMinimizationResult                        ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Result structure for multidimensional minimization
     */
    struct MultidimMinimizationResult
    {
        Vector<Real> xmin;      ///< Location of minimum
        Real         fmin;       ///< Function value at minimum
        int          iterations; ///< Number of iterations (or function evaluations)
        bool         converged;  ///< True if converged within tolerance

        MultidimMinimizationResult()
            : fmin(0.0), iterations(0), converged(false) {}

        MultidimMinimizationResult(const Vector<Real>& x, Real f, int iter, bool conv)
            : xmin(x), fmin(f), iterations(iter), converged(conv) {}
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                 Simplex (Nelder-Mead) Minimization                 ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Downhill simplex (Nelder-Mead) method for multidimensional minimization
     * 
     * The Nelder-Mead simplex method is a direct search method that does not require 
     * derivatives. It works by maintaining a simplex (a geometric figure with N+1 vertices 
     * in N dimensions) and iteratively replacing the worst vertex.
     * 
     * Reference: Numerical Recipes Chapter 10.5, Nelder & Mead (1965)
     * 
     * The algorithm uses four basic operations:
     * - Reflection: Reflect the worst point through the centroid
     * - Expansion: If reflection is good, try going further
     * - Contraction: If reflection is bad, try a point between worst and centroid
     * - Shrink: If all else fails, shrink the simplex toward the best point
     * 
     * Standard coefficients:
     * - alpha = 1.0 (reflection coefficient)
     * - gamma = 2.0 (expansion coefficient)
     * - rho = 0.5 (contraction coefficient)  
     * - sigma = 0.5 (shrink coefficient)
     */
    class NelderMead
    {
    public:
        // Simplex coefficients
        static constexpr Real ALPHA = 1.0;   // Reflection coefficient
        static constexpr Real GAMMA = 2.0;   // Expansion coefficient
        static constexpr Real RHO   = 0.5;   // Contraction coefficient
        static constexpr Real SIGMA = 0.5;   // Shrink coefficient

    private:
        Real _ftol;         // Fractional tolerance for convergence
        int  _maxIter;      // Maximum iterations
        int  _nfunc;        // Number of function evaluations
        int  _ndim;         // Problem dimension
        int  _mpts;         // Number of simplex points (ndim + 1)

        Matrix<Real> _p;    // Simplex vertices [mpts x ndim]
        Vector<Real> _y;    // Function values at vertices [mpts]

    public:
        /**
         * @brief Construct Nelder-Mead optimizer
         * @param ftol Fractional convergence tolerance (default 1e-8)
         * @param maxIter Maximum iterations (default 5000)
         */
        NelderMead(Real ftol = 1.0e-8, int maxIter = 5000)
            : _ftol(ftol), _maxIter(maxIter), _nfunc(0), _ndim(0), _mpts(0) {}

        Real getFtol() const { return _ftol; }
        void setFtol(Real ftol) { _ftol = ftol; }

        int getMaxIter() const { return _maxIter; }
        void setMaxIter(int maxIter) { _maxIter = maxIter; }

        int getNumFuncEvals() const { return _nfunc; }

        /**
         * @brief Get the final simplex (for analysis/debugging)
         */
        const Matrix<Real>& getSimplex() const { return _p; }

        /**
         * @brief Get function values at simplex vertices
         */
        const Vector<Real>& getSimplexValues() const { return _y; }

        ///////////////////////////////////////////////////////////////////////////
        ///                      Minimize with uniform delta                    ///
        ///////////////////////////////////////////////////////////////////////////
        /**
         * @brief Minimize using a starting point with uniform perturbation
         * @tparam N Dimension of the problem
         * @param func Scalar function to minimize
         * @param start Starting point
         * @param delta Uniform perturbation for creating initial simplex (default 1.0)
         * @return Minimization result
         * 
         * Creates initial simplex by adding delta to each coordinate of start
         */
        template<int N>
        MultidimMinimizationResult Minimize(const IScalarFunction<N>& func,
                                            const VectorN<Real, N>& start,
                                            Real delta = 1.0)
        {
            Vector<Real> deltas(N, delta);
            return Minimize(func, start, deltas);
        }

        ///////////////////////////////////////////////////////////////////////////
        ///                    Minimize with per-dimension deltas               ///
        ///////////////////////////////////////////////////////////////////////////
        /**
         * @brief Minimize using a starting point with per-dimension perturbations
         * @tparam N Dimension of the problem
         * @param func Scalar function to minimize
         * @param start Starting point
         * @param deltas Per-dimension perturbations for initial simplex
         * @return Minimization result
         * 
         * Creates initial simplex by adding deltas[i] to coordinate i of start
         */
        template<int N>
        MultidimMinimizationResult Minimize(const IScalarFunction<N>& func,
                                            const VectorN<Real, N>& start,
                                            const Vector<Real>& deltas)
        {
            if (deltas.size() != N)
                throw MultidimOptimizationError("Deltas vector dimension mismatch");

            // Create initial simplex: N+1 vertices
            _ndim = N;
            _mpts = N + 1;
            _p = Matrix<Real>(_mpts, _ndim);
            _y = Vector<Real>(_mpts);

            // First vertex is the starting point
            for (int j = 0; j < _ndim; ++j)
                _p(0, j) = start[j];

            // Other vertices: perturb one coordinate at a time
            for (int i = 1; i <= _ndim; ++i)
            {
                for (int j = 0; j < _ndim; ++j)
                    _p(i, j) = start[j];
                _p(i, i - 1) += deltas[i - 1];
            }

            return MinimizeFromSimplex(func);
        }

        ///////////////////////////////////////////////////////////////////////////
        ///                    Minimize from explicit simplex                   ///
        ///////////////////////////////////////////////////////////////////////////
        /**
         * @brief Minimize from an explicitly specified initial simplex
         * @tparam N Dimension of the problem
         * @param func Scalar function to minimize
         * @param simplex Initial simplex as (N+1) x N matrix (rows are vertices)
         * @return Minimization result
         */
        template<int N>
        MultidimMinimizationResult Minimize(const IScalarFunction<N>& func,
                                            const Matrix<Real>& simplex)
        {
            if (simplex.RowNum() != N + 1 || simplex.ColNum() != N)
                throw MultidimOptimizationError("Simplex dimensions invalid: expected " +
                    std::to_string(N + 1) + " x " + std::to_string(N));

            _ndim = N;
            _mpts = N + 1;
            _p = simplex;
            _y = Vector<Real>(_mpts);

            return MinimizeFromSimplex(func);
        }

    private:
        ///////////////////////////////////////////////////////////////////////////
        ///                   Core Nelder-Mead Algorithm                        ///
        ///////////////////////////////////////////////////////////////////////////
        template<int N>
        MultidimMinimizationResult MinimizeFromSimplex(const IScalarFunction<N>& func)
        {
            const Real TINY = 1.0e-10;
            int ihi, ilo, inhi;
            VectorN<Real, N> x;
            Vector<Real> psum(_ndim);

            // Evaluate function at all simplex vertices
            for (int i = 0; i < _mpts; ++i)
            {
                for (int j = 0; j < _ndim; ++j)
                    x[j] = _p(i, j);
                _y[i] = func(x);
            }
            _nfunc = _mpts;

            // Compute initial psum (sum of vertex coordinates)
            GetPsum(psum);

            // Main iteration loop
            for (int iter = 0; iter < _maxIter; ++iter)
            {
                // Find lowest (best), highest (worst), and next-highest
                ilo = 0;
                if (_y[0] > _y[1])
                {
                    ihi = 0;
                    inhi = 1;
                }
                else
                {
                    ihi = 1;
                    inhi = 0;
                }

                for (int i = 0; i < _mpts; ++i)
                {
                    if (_y[i] <= _y[ilo]) ilo = i;
                    if (_y[i] > _y[ihi])
                    {
                        inhi = ihi;
                        ihi = i;
                    }
                    else if (_y[i] > _y[inhi] && i != ihi)
                    {
                        inhi = i;
                    }
                }

                // Check convergence
                Real rtol = 2.0 * std::abs(_y[ihi] - _y[ilo]) /
                            (std::abs(_y[ihi]) + std::abs(_y[ilo]) + TINY);

                if (rtol < _ftol)
                {
                    // Converged - swap best to position 0
                    std::swap(_y[0], _y[ilo]);
                    for (int j = 0; j < _ndim; ++j)
                    {
                        std::swap(_p(0, j), _p(ilo, j));
                        x[j] = _p(0, j);
                    }

                    Vector<Real> result(_ndim);
                    for (int j = 0; j < _ndim; ++j)
                        result[j] = x[j];

                    return MultidimMinimizationResult(result, _y[0], _nfunc, true);
                }

                // Check iteration limit
                if (_nfunc >= _maxIter)
                {
                    Vector<Real> result(_ndim);
                    for (int j = 0; j < _ndim; ++j)
                        result[j] = _p(ilo, j);
                    return MultidimMinimizationResult(result, _y[ilo], _nfunc, false);
                }

                // Try reflection
                Real ytry = Amotry<N>(psum, ihi, -ALPHA, func);
                _nfunc++;

                if (ytry <= _y[ilo])
                {
                    // Reflection is better than best - try expansion
                    ytry = Amotry<N>(psum, ihi, GAMMA, func);
                    _nfunc++;
                }
                else if (ytry >= _y[inhi])
                {
                    // Reflection is worse than second-worst
                    Real ysave = _y[ihi];
                    ytry = Amotry<N>(psum, ihi, RHO, func);
                    _nfunc++;

                    if (ytry >= ysave)
                    {
                        // Contraction failed - do shrink
                        for (int i = 0; i < _mpts; ++i)
                        {
                            if (i != ilo)
                            {
                                for (int j = 0; j < _ndim; ++j)
                                {
                                    _p(i, j) = SIGMA * (_p(i, j) + _p(ilo, j));
                                    x[j] = _p(i, j);
                                }
                                _y[i] = func(x);
                            }
                        }
                        _nfunc += _ndim;
                        GetPsum(psum);
                    }
                }
            }

            // Max iterations reached without convergence
            int ilo_final = 0;
            for (int i = 1; i < _mpts; ++i)
                if (_y[i] < _y[ilo_final]) ilo_final = i;

            Vector<Real> result(_ndim);
            for (int j = 0; j < _ndim; ++j)
                result[j] = _p(ilo_final, j);

            return MultidimMinimizationResult(result, _y[ilo_final], _nfunc, false);
        }

        ///////////////////////////////////////////////////////////////////////////
        ///                           Helper Functions                          ///
        ///////////////////////////////////////////////////////////////////////////
        
        /**
         * @brief Compute sum of vertex coordinates (for centroid calculation)
         */
        void GetPsum(Vector<Real>& psum)
        {
            for (int j = 0; j < _ndim; ++j)
            {
                Real sum = 0.0;
                for (int i = 0; i < _mpts; ++i)
                    sum += _p(i, j);
                psum[j] = sum;
            }
        }

        /**
         * @brief Extrapolate by a factor through the face opposite the high point
         * @tparam N Problem dimension
         * @param psum Sum of vertex coordinates
         * @param ihi Index of highest (worst) point
         * @param fac Extrapolation factor
         * @param func Function to evaluate
         * @return Function value at trial point
         * 
         * This implements the core simplex transformation. The trial point is:
         * ptry = psum * (1-fac)/ndim - p[ihi] * ((1-fac)/ndim - fac)
         *      = centroid_of_face - fac * (centroid_of_face - p[ihi])
         * 
         * For fac = -1: reflection through centroid
         * For fac = 2: expansion beyond reflection point
         * For fac = 0.5: contraction toward centroid
         */
        template<int N>
        Real Amotry(Vector<Real>& psum, int ihi, Real fac, const IScalarFunction<N>& func)
        {
            VectorN<Real, N> ptry;
            Real fac1 = (1.0 - fac) / _ndim;
            Real fac2 = fac1 - fac;

            for (int j = 0; j < _ndim; ++j)
                ptry[j] = psum[j] * fac1 - _p(ihi, j) * fac2;

            Real ytry = func(ptry);

            if (ytry < _y[ihi])
            {
                // Accept the new point
                _y[ihi] = ytry;
                for (int j = 0; j < _ndim; ++j)
                {
                    psum[j] += ptry[j] - _p(ihi, j);
                    _p(ihi, j) = ptry[j];
                }
            }
            return ytry;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                   Convenience Wrapper Functions                     ///
    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Minimize a scalar function using Nelder-Mead
     * @tparam N Dimension of the problem
     * @param func Function to minimize
     * @param start Starting point
     * @param delta Initial simplex size (default 1.0)
     * @param ftol Convergence tolerance (default 1e-8)
     * @return Minimization result
     */
    template<int N>
    MultidimMinimizationResult NelderMeadMinimize(
        const IScalarFunction<N>& func,
        const VectorN<Real, N>& start,
        Real delta = 1.0,
        Real ftol = 1.0e-8)
    {
        NelderMead optimizer(ftol);
        return optimizer.Minimize(func, start, delta);
    }

    /**
     * @brief Minimize a scalar function using Nelder-Mead with custom deltas
     * @tparam N Dimension of the problem
     * @param func Function to minimize
     * @param start Starting point
     * @param deltas Per-dimension simplex sizes
     * @param ftol Convergence tolerance (default 1e-8)
     * @return Minimization result
     */
    template<int N>
    MultidimMinimizationResult NelderMeadMinimize(
        const IScalarFunction<N>& func,
        const VectorN<Real, N>& start,
        const Vector<Real>& deltas,
        Real ftol = 1.0e-8)
    {
        NelderMead optimizer(ftol);
        return optimizer.Minimize(func, start, deltas);
    }

    ///////////////////////////////////////////////////////////////////////////
    ///                     Maximization Wrapper                            ///
    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Helper class to negate a function for maximization
     */
    template<int N>
    class NegatedScalarFunction : public IScalarFunction<N>
    {
    private:
        const IScalarFunction<N>& _func;

    public:
        explicit NegatedScalarFunction(const IScalarFunction<N>& func) : _func(func) {}

        Real operator()(const VectorN<Real, N>& x) const override
        {
            return -_func(x);
        }
    };

    /**
     * @brief Maximize a scalar function using Nelder-Mead
     * @tparam N Dimension of the problem
     * @param func Function to maximize
     * @param start Starting point
     * @param delta Initial simplex size (default 1.0)
     * @param ftol Convergence tolerance (default 1e-8)
     * @return Maximization result (fmin is negated to give actual maximum value)
     */
    template<int N>
    MultidimMinimizationResult NelderMeadMaximize(
        const IScalarFunction<N>& func,
        const VectorN<Real, N>& start,
        Real delta = 1.0,
        Real ftol = 1.0e-8)
    {
        NegatedScalarFunction<N> negFunc(func);
        NelderMead optimizer(ftol);
        auto result = optimizer.Minimize(negFunc, start, delta);
        result.fmin = -result.fmin;  // Convert back to maximum value
        return result;
    }

    /////////////////////////////////////////////////////////////////////
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
    class IDifferentiableScalarFunction : public IScalarFunction<N>
    {
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
    class LineFunction : public IRealFunction
    {
    private:
        const IScalarFunction<N>& _func;
        const VectorN<Real, N>& _p;       // Base point
        const VectorN<Real, N>& _xi;      // Direction

    public:
        LineFunction(const IScalarFunction<N>& func,
                     const VectorN<Real, N>& p,
                     const VectorN<Real, N>& xi)
            : _func(func), _p(p), _xi(xi) {}

        Real operator()(Real t) const override
        {
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
    class DLineFunction : public IRealFunction
    {
    private:
        const IDifferentiableScalarFunction<N>& _func;
        VectorN<Real, N> _p;       // Base point (mutable for evaluation)
        VectorN<Real, N> _xi;      // Direction
        mutable VectorN<Real, N> _xt;      // Current point
        mutable VectorN<Real, N> _grad;    // Gradient at current point

    public:
        DLineFunction(const IDifferentiableScalarFunction<N>& func,
                      const VectorN<Real, N>& p,
                      const VectorN<Real, N>& xi)
            : _func(func), _p(p), _xi(xi) {}

        void updateBasePoint(const VectorN<Real, N>& p) { _p = p; }
        void updateDirection(const VectorN<Real, N>& xi) { _xi = xi; }

        Real operator()(Real t) const override
        {
            for (int j = 0; j < N; ++j)
                _xt[j] = _p[j] + t * _xi[j];
            return _func(_xt);
        }

        /**
         * @brief Compute directional derivative df/dt = grad(f) · xi
         */
        Real derivative(Real t) const
        {
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
    class LineMinimizer
    {
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
        static Real Minimize(const IScalarFunction<N>& func,
                             VectorN<Real, N>& p,
                             VectorN<Real, N>& xi,
                             Real tol = 3.0e-8)
        {
            LineFunction<N> f1dim(func, p, xi);

            // Bracket the minimum
            auto bracket = Minimization::BracketMinimum(f1dim, 0.0, 1.0);
            if (!bracket.valid)
            {
                // Try with different initial interval
                bracket = Minimization::BracketMinimum(f1dim, 0.0, 0.1);
            }

            // Find minimum using Brent's method
            auto result = Minimization::BrentMinimize(f1dim, bracket, tol);
            Real xmin = result.xmin;

            // Update p and xi
            for (int j = 0; j < N; ++j)
            {
                xi[j] *= xmin;
                p[j] += xi[j];
            }

            return result.fmin;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                      Powell's Method                                ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Powell's direction set method for multidimensional minimization
     * 
     * Powell's method is a derivative-free optimization method that performs
     * successive line minimizations along a set of directions that are updated
     * to become mutually conjugate.
     * 
     * Reference: Numerical Recipes Chapter 10.7
     * 
     * The algorithm:
     * 1. Start with N unit vectors as search directions
     * 2. Minimize along each direction in sequence
     * 3. Construct new direction from total displacement
     * 4. Replace direction with largest decrease with new direction
     * 5. Repeat until convergence
     * 
     * This method is particularly effective when:
     * - Derivatives are not available
     * - Function is smooth and well-behaved
     * - Problem dimension is moderate (N < 20)
     */
    class Powell
    {
    private:
        Real _ftol;         // Fractional tolerance
        int  _maxIter;      // Maximum iterations
        int  _iter;         // Iteration counter
        Real _fret;         // Current function value

    public:
        Powell(Real ftol = 3.0e-8, int maxIter = 200)
            : _ftol(ftol), _maxIter(maxIter), _iter(0), _fret(0.0) {}

        Real getFtol() const { return _ftol; }
        void setFtol(Real ftol) { _ftol = ftol; }

        int getMaxIter() const { return _maxIter; }
        void setMaxIter(int maxIter) { _maxIter = maxIter; }

        int getIterations() const { return _iter; }
        Real getCurrentFValue() const { return _fret; }

        /**
         * @brief Minimize using starting point with identity direction matrix
         */
        template<int N>
        MultidimMinimizationResult Minimize(const IScalarFunction<N>& func,
                                            const VectorN<Real, N>& start)
        {
            // Initialize direction matrix to identity
            Matrix<Real> ximat(N, N, 0.0);
            for (int i = 0; i < N; ++i)
                ximat(i, i) = 1.0;

            return Minimize(func, start, ximat);
        }

        /**
         * @brief Minimize with custom initial direction matrix
         */
        template<int N>
        MultidimMinimizationResult Minimize(const IScalarFunction<N>& func,
                                            const VectorN<Real, N>& start,
                                            Matrix<Real>& ximat)
        {
            const Real TINY = 1.0e-25;

            VectorN<Real, N> p = start;
            VectorN<Real, N> pt, ptt, xi;

            _fret = func(p);

            // Save initial point
            for (int j = 0; j < N; ++j)
                pt[j] = p[j];

            for (_iter = 0; _iter < _maxIter; ++_iter)
            {
                Real fp = _fret;
                int ibig = 0;
                Real del = 0.0;

                // Minimize along each direction
                for (int i = 0; i < N; ++i)
                {
                    // Extract i-th direction from matrix columns
                    for (int j = 0; j < N; ++j)
                        xi[j] = ximat(j, i);

                    Real fptt = _fret;
                    _fret = LineMinimizer::Minimize(func, p, xi, _ftol);

                    // Track direction with largest decrease
                    if (fptt - _fret > del)
                    {
                        del = fptt - _fret;
                        ibig = i + 1;
                    }
                }

                // Check convergence
                if (2.0 * (fp - _fret) <= _ftol * (std::abs(fp) + std::abs(_fret)) + TINY)
                {
                    Vector<Real> result(N);
                    for (int j = 0; j < N; ++j)
                        result[j] = p[j];
                    return MultidimMinimizationResult(result, _fret, _iter + 1, true);
                }

                // Construct extrapolated point and new direction
                for (int j = 0; j < N; ++j)
                {
                    ptt[j] = 2.0 * p[j] - pt[j];  // Extrapolated point
                    xi[j] = p[j] - pt[j];          // New direction
                    pt[j] = p[j];                  // Save current point
                }

                Real fptt = func(ptt);

                if (fptt < fp)
                {
                    Real t = 2.0 * (fp - 2.0 * _fret + fptt) *
                             (fp - _fret - del) * (fp - _fret - del) -
                             del * (fp - fptt) * (fp - fptt);

                    if (t < 0.0)
                    {
                        // Replace direction ibig with new direction
                        _fret = LineMinimizer::Minimize(func, p, xi, _ftol);

                        for (int j = 0; j < N; ++j)
                        {
                            ximat(j, ibig - 1) = ximat(j, N - 1);
                            ximat(j, N - 1) = xi[j];
                        }
                    }
                }
            }

            // Max iterations reached
            Vector<Real> result(N);
            for (int j = 0; j < N; ++j)
                result[j] = p[j];
            return MultidimMinimizationResult(result, _fret, _iter, false);
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///            Conjugate Gradient Method (Fletcher-Reeves-Polak-Ribiere)///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Conjugate gradient method for minimization with derivatives
     * 
     * The conjugate gradient method is one of the most effective methods for
     * minimizing smooth functions when gradients are available. It generates
     * search directions that are conjugate with respect to the Hessian.
     * 
     * Reference: Numerical Recipes Chapter 10.8
     * 
     * Two variants are implemented:
     * - Fletcher-Reeves: gamma = (g_new · g_new) / (g_old · g_old)
     * - Polak-Ribiere:   gamma = (g_new · (g_new - g_old)) / (g_old · g_old)
     * 
     * Polak-Ribiere is generally preferred as it automatically resets to
     * steepest descent when progress stalls.
     * 
     * The algorithm:
     * 1. Compute gradient at starting point
     * 2. Set initial search direction to negative gradient (steepest descent)
     * 3. Minimize along search direction
     * 4. Compute new gradient
     * 5. Update search direction using conjugate gradient formula
     * 6. Repeat until convergence
     */
    class ConjugateGradient
    {
    public:
        enum class Method
        {
            FletcherReeves,  // Original CG formula
            PolakRibiere     // Generally preferred variant
        };

    private:
        Real   _ftol;       // Function tolerance
        Real   _gtol;       // Gradient tolerance
        int    _maxIter;    // Maximum iterations
        int    _iter;       // Iteration counter
        Real   _fret;       // Current function value
        Method _method;     // CG variant

    public:
        ConjugateGradient(Real ftol = 3.0e-8, Real gtol = 1.0e-8, 
                          int maxIter = 200, Method method = Method::PolakRibiere)
            : _ftol(ftol), _gtol(gtol), _maxIter(maxIter), 
              _iter(0), _fret(0.0), _method(method) {}

        Real getFtol() const { return _ftol; }
        void setFtol(Real ftol) { _ftol = ftol; }

        Real getGtol() const { return _gtol; }
        void setGtol(Real gtol) { _gtol = gtol; }

        int getMaxIter() const { return _maxIter; }
        void setMaxIter(int maxIter) { _maxIter = maxIter; }

        Method getMethod() const { return _method; }
        void setMethod(Method method) { _method = method; }

        int getIterations() const { return _iter; }
        Real getCurrentFValue() const { return _fret; }

        /**
         * @brief Minimize a differentiable function using conjugate gradient
         */
        template<int N>
        MultidimMinimizationResult Minimize(const IDifferentiableScalarFunction<N>& func,
                                            const VectorN<Real, N>& start)
        {
            const Real EPS = 1.0e-18;

            VectorN<Real, N> p = start;
            VectorN<Real, N> g, h, xi;

            // Evaluate function and gradient at starting point
            Real fp = func(p);
            func.Gradient(p, xi);

            // Initialize: g = -gradient, h = xi = g (steepest descent)
            for (int j = 0; j < N; ++j)
            {
                g[j] = -xi[j];
                xi[j] = h[j] = g[j];
            }

            for (_iter = 0; _iter < _maxIter; ++_iter)
            {
                // Line minimization along xi
                _fret = LineMinimizer::Minimize(func, p, xi, _ftol);

                // Check function convergence
                if (2.0 * std::abs(_fret - fp) <= _ftol * (std::abs(_fret) + std::abs(fp) + EPS))
                {
                    Vector<Real> result(N);
                    for (int j = 0; j < N; ++j)
                        result[j] = p[j];
                    return MultidimMinimizationResult(result, _fret, _iter + 1, true);
                }

                fp = _fret;

                // Compute new gradient
                func.Gradient(p, xi);

                // Check gradient convergence
                Real test = 0.0;
                Real den = std::max(std::abs(fp), REAL(1.0));
                for (int j = 0; j < N; ++j)
                {
                    Real temp = std::abs(xi[j]) * std::max(std::abs(p[j]), REAL(1.0)) / den;
                    if (temp > test) test = temp;
                }
                if (test < _gtol)
                {
                    Vector<Real> result(N);
                    for (int j = 0; j < N; ++j)
                        result[j] = p[j];
                    return MultidimMinimizationResult(result, _fret, _iter + 1, true);
                }

                // Compute gamma (CG update coefficient)
                Real gg = 0.0, dgg = 0.0;
                for (int j = 0; j < N; ++j)
                {
                    gg += g[j] * g[j];

                    if (_method == Method::FletcherReeves)
                    {
                        dgg += xi[j] * xi[j];  // Fletcher-Reeves
                    }
                    else
                    {
                        dgg += (xi[j] + g[j]) * xi[j];  // Polak-Ribiere
                    }
                }

                if (gg == 0.0)
                {
                    // Gradient is zero - we're done
                    Vector<Real> result(N);
                    for (int j = 0; j < N; ++j)
                        result[j] = p[j];
                    return MultidimMinimizationResult(result, _fret, _iter + 1, true);
                }

                Real gam = dgg / gg;

                // Update search direction
                for (int j = 0; j < N; ++j)
                {
                    g[j] = -xi[j];
                    xi[j] = h[j] = g[j] + gam * h[j];
                }
            }

            // Max iterations reached
            Vector<Real> result(N);
            for (int j = 0; j < N; ++j)
                result[j] = p[j];
            return MultidimMinimizationResult(result, _fret, _iter, false);
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///               BFGS Quasi-Newton Method                              ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief BFGS (Broyden-Fletcher-Goldfarb-Shanno) quasi-Newton method
     * 
     * BFGS is one of the most effective methods for smooth unconstrained
     * optimization. It builds an approximation to the inverse Hessian matrix
     * using only gradient information.
     * 
     * The algorithm:
     * 1. Start with identity approximation to inverse Hessian
     * 2. Compute search direction p = -H * gradient
     * 3. Line search to find step size
     * 4. Update inverse Hessian approximation using BFGS formula
     * 5. Repeat until convergence
     * 
     * Advantages:
     * - Superlinear convergence rate
     * - Only requires gradient (not Hessian)
     * - Self-correcting: recovers from poor initial H estimate
     * 
     * Disadvantages:
     * - O(N²) storage for inverse Hessian approximation
     * - May be expensive for very large N
     */
    class BFGS
    {
    private:
        Real _ftol;         // Function tolerance
        Real _gtol;         // Gradient tolerance
        int  _maxIter;      // Maximum iterations
        int  _iter;         // Iteration counter
        Real _fret;         // Current function value

    public:
        BFGS(Real ftol = 3.0e-8, Real gtol = 1.0e-8, int maxIter = 200)
            : _ftol(ftol), _gtol(gtol), _maxIter(maxIter), _iter(0), _fret(0.0) {}

        Real getFtol() const { return _ftol; }
        void setFtol(Real ftol) { _ftol = ftol; }

        Real getGtol() const { return _gtol; }
        void setGtol(Real gtol) { _gtol = gtol; }

        int getMaxIter() const { return _maxIter; }
        void setMaxIter(int maxIter) { _maxIter = maxIter; }

        int getIterations() const { return _iter; }
        Real getCurrentFValue() const { return _fret; }

        /**
         * @brief Minimize using BFGS quasi-Newton method
         */
        template<int N>
        MultidimMinimizationResult Minimize(const IDifferentiableScalarFunction<N>& func,
                                            const VectorN<Real, N>& start)
        {
            const Real EPS = 1.0e-18;
            const Real STPMX = 100.0;  // Maximum step size

            VectorN<Real, N> p = start;
            VectorN<Real, N> g, dg, hdg, pnew, xi;
            Matrix<Real> hessin(N, N, 0.0);

            // Initialize inverse Hessian to identity
            for (int i = 0; i < N; ++i)
                hessin(i, i) = 1.0;

            // Evaluate function and gradient
            _fret = func(p);
            func.Gradient(p, g);

            // Initial search direction is steepest descent
            for (int j = 0; j < N; ++j)
                xi[j] = -g[j];

            // Compute maximum step size
            Real sum = 0.0;
            for (int j = 0; j < N; ++j)
                sum += p[j] * p[j];
            Real stpmax = STPMX * std::max(std::sqrt(sum), static_cast<Real>(N));

            for (_iter = 0; _iter < _maxIter; ++_iter)
            {
                // Line search along xi
                Real fp = _fret;
                VectorN<Real, N> pold = p;

                // Perform line minimization
                _fret = LineSearchBacktrack(func, p, g, xi, stpmax);

                // Check for convergence on function value
                if (2.0 * std::abs(_fret - fp) <= _ftol * (std::abs(_fret) + std::abs(fp) + EPS))
                {
                    Vector<Real> result(N);
                    for (int j = 0; j < N; ++j)
                        result[j] = p[j];
                    return MultidimMinimizationResult(result, _fret, _iter + 1, true);
                }

                // Update xi to be the actual step taken
                for (int j = 0; j < N; ++j)
                    xi[j] = p[j] - pold[j];

                // Save old gradient and compute new gradient
                for (int j = 0; j < N; ++j)
                    dg[j] = g[j];
                func.Gradient(p, g);

                // Check gradient convergence
                Real test = 0.0;
                Real den = std::max(std::abs(_fret), REAL(1.0));
                for (int j = 0; j < N; ++j)
                {
                    Real temp = std::abs(g[j]) * std::max(std::abs(p[j]), REAL(1.0)) / den;
                    if (temp > test) test = temp;
                }
                if (test < _gtol)
                {
                    Vector<Real> result(N);
                    for (int j = 0; j < N; ++j)
                        result[j] = p[j];
                    return MultidimMinimizationResult(result, _fret, _iter + 1, true);
                }

                // Compute gradient difference
                for (int j = 0; j < N; ++j)
                    dg[j] = g[j] - dg[j];

                // Compute H * dg
                for (int i = 0; i < N; ++i)
                {
                    hdg[i] = 0.0;
                    for (int j = 0; j < N; ++j)
                        hdg[i] += hessin(i, j) * dg[j];
                }

                // BFGS update of inverse Hessian
                Real fac = 0.0, fae = 0.0, sumdg = 0.0, sumxi = 0.0;
                for (int j = 0; j < N; ++j)
                {
                    fac += dg[j] * xi[j];
                    fae += dg[j] * hdg[j];
                    sumdg += dg[j] * dg[j];
                    sumxi += xi[j] * xi[j];
                }

                if (fac > std::sqrt(EPS * sumdg * sumxi))
                {
                    fac = 1.0 / fac;
                    Real fad = 1.0 / fae;

                    // Vector that makes BFGS different from DFP
                    for (int j = 0; j < N; ++j)
                        dg[j] = fac * xi[j] - fad * hdg[j];

                    // BFGS formula for updating inverse Hessian
                    for (int i = 0; i < N; ++i)
                    {
                        for (int j = i; j < N; ++j)
                        {
                            hessin(i, j) += fac * xi[i] * xi[j] - 
                                            fad * hdg[i] * hdg[j] +
                                            fae * dg[i] * dg[j];
                            hessin(j, i) = hessin(i, j);
                        }
                    }
                }

                // Compute new search direction
                for (int i = 0; i < N; ++i)
                {
                    xi[i] = 0.0;
                    for (int j = 0; j < N; ++j)
                        xi[i] -= hessin(i, j) * g[j];
                }
            }

            // Max iterations reached
            Vector<Real> result(N);
            for (int j = 0; j < N; ++j)
                result[j] = p[j];
            return MultidimMinimizationResult(result, _fret, _iter, false);
        }

    private:
        /**
         * @brief Backtracking line search with Armijo condition
         */
        template<int N>
        Real LineSearchBacktrack(const IDifferentiableScalarFunction<N>& func,
                                 VectorN<Real, N>& p,
                                 const VectorN<Real, N>& g,
                                 VectorN<Real, N>& xi,
                                 Real stpmax)
        {
            const Real ALF = 1.0e-4;   // Ensures sufficient decrease
            const Real TOLX = 1.0e-12; // Convergence criterion on x

            // Scale if step is too big
            Real sum = 0.0;
            for (int j = 0; j < N; ++j)
                sum += xi[j] * xi[j];
            sum = std::sqrt(sum);

            if (sum > stpmax)
            {
                Real scale = stpmax / sum;
                for (int j = 0; j < N; ++j)
                    xi[j] *= scale;
            }

            // Compute slope
            Real slope = 0.0;
            for (int j = 0; j < N; ++j)
                slope += g[j] * xi[j];

            if (slope >= 0.0)
            {
                // Reset to steepest descent if slope is not negative
                for (int j = 0; j < N; ++j)
                    xi[j] = -g[j];
                slope = 0.0;
                for (int j = 0; j < N; ++j)
                    slope += g[j] * xi[j];
            }

            // Compute lambda_min
            Real test = 0.0;
            for (int j = 0; j < N; ++j)
            {
                Real temp = std::abs(xi[j]) / std::max(std::abs(p[j]), REAL(1.0));
                if (temp > test) test = temp;
            }
            Real alamin = TOLX / test;

            Real alam = 1.0;  // Always try full Newton step first
            Real f = func(p);
            Real alam2 = 0.0, f2 = 0.0;

            for (;;)
            {
                VectorN<Real, N> pnew;
                for (int j = 0; j < N; ++j)
                    pnew[j] = p[j] + alam * xi[j];

                Real fnew = func(pnew);

                if (alam < alamin)
                {
                    // Convergence on delta x
                    for (int j = 0; j < N; ++j)
                        p[j] = pnew[j];
                    return fnew;
                }
                else if (fnew <= f + ALF * alam * slope)
                {
                    // Sufficient decrease - accept step
                    for (int j = 0; j < N; ++j)
                        p[j] = pnew[j];
                    return fnew;
                }
                else
                {
                    // Backtrack
                    Real tmplam;
                    if (alam == 1.0)
                    {
                        // First backtrack: quadratic
                        tmplam = -slope / (2.0 * (fnew - f - slope));
                    }
                    else
                    {
                        // Subsequent backtracks: cubic
                        Real rhs1 = fnew - f - alam * slope;
                        Real rhs2 = f2 - f - alam2 * slope;
                        Real a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) /
                                 (alam - alam2);
                        Real b = (-alam2 * rhs1 / (alam * alam) +
                                  alam * rhs2 / (alam2 * alam2)) / (alam - alam2);

                        if (a == 0.0)
                        {
                            tmplam = -slope / (2.0 * b);
                        }
                        else
                        {
                            Real disc = b * b - 3.0 * a * slope;
                            if (disc < 0.0)
                                tmplam = 0.5 * alam;
                            else if (b <= 0.0)
                                tmplam = (-b + std::sqrt(disc)) / (3.0 * a);
                            else
                                tmplam = -slope / (b + std::sqrt(disc));
                        }

                        // Limit step size reduction
                        if (tmplam > 0.5 * alam)
                            tmplam = 0.5 * alam;
                    }

                    alam2 = alam;
                    f2 = fnew;
                    alam = std::max(tmplam, REAL(0.1) * alam);
                }
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                   Convenience Wrapper Functions                     ///
    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Minimize using Powell's method (no derivatives needed)
     */
    template<int N>
    MultidimMinimizationResult PowellMinimize(
        const IScalarFunction<N>& func,
        const VectorN<Real, N>& start,
        Real ftol = 3.0e-8)
    {
        Powell optimizer(ftol);
        return optimizer.Minimize(func, start);
    }

    /**
     * @brief Minimize using conjugate gradient (requires derivatives)
     */
    template<int N>
    MultidimMinimizationResult ConjugateGradientMinimize(
        const IDifferentiableScalarFunction<N>& func,
        const VectorN<Real, N>& start,
        Real ftol = 3.0e-8,
        ConjugateGradient::Method method = ConjugateGradient::Method::PolakRibiere)
    {
        ConjugateGradient optimizer(ftol, 1.0e-8, 200, method);
        return optimizer.Minimize(func, start);
    }

    /**
     * @brief Minimize using BFGS quasi-Newton method (requires derivatives)
     */
    template<int N>
    MultidimMinimizationResult BFGSMinimize(
        const IDifferentiableScalarFunction<N>& func,
        const VectorN<Real, N>& start,
        Real ftol = 3.0e-8)
    {
        BFGS optimizer(ftol);
        return optimizer.Minimize(func, start);
    }

} // namespace MML

#endif // MML_OPTIMIZATION_MULTIDIM_H