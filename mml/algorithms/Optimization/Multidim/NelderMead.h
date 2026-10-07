///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Multidim/NelderMead.h                                                ///
///  Description: Nelder-Mead simplex minimization and related convenience wrappers        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_NELDER_MEAD_H
#define MML_OPTIMIZATION_NELDER_MEAD_H

#include <mml/algorithms/Optimization/Multidim/MultidimTypes.h>

#include <algorithm>
#include <cmath>

namespace MML::Optimization {

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
	class NelderMead {
	public:
		// Simplex coefficients
		static constexpr Real ALPHA = 1.0; // Reflection coefficient
		static constexpr Real GAMMA = 2.0; // Expansion coefficient
		static constexpr Real RHO = 0.5;   // Contraction coefficient
		static constexpr Real SIGMA = 0.5; // Shrink coefficient

	private:
		Real _ftol;	  // Fractional tolerance for convergence
		int _maxIter; // Maximum number of function evaluations (also caps outer iterations)
		int _nfunc;	  // Number of function evaluations
		int _ndim;	  // Problem dimension
		int _mpts;	  // Number of simplex points (ndim + 1)

		Matrix<Real> _p; // Simplex vertices [mpts x ndim]
		Vector<Real> _y; // Function values at vertices [mpts]

	public:
		/**
         * @brief Construct Nelder-Mead optimizer
         * @param ftol Fractional convergence tolerance (default 1e-8)
         * @param maxIter Maximum number of function evaluations (default 5000).
         *                Note: this budget counts function evaluations, not simplex
         *                iterations; each iteration performs 1-2+ evaluations.
         */
		NelderMead(Real ftol = PrecisionValues<Real>::OptimizationGradientTolerance, int maxIter = 5000)
			: _ftol(ftol)
			, _maxIter(maxIter)
			, _nfunc(0)
			, _ndim(0)
			, _mpts(0) {}

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
         * @throws MultidimOptimizationInputError if inputs are invalid
         * 
         * Creates initial simplex by adding delta to each coordinate of start
         */
		template<int N>
		MultidimMinimizationResult Minimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start, Real delta = 1.0) {
			// Validate inputs
			ValidateVectorFinite<N>(start, "NelderMead::Minimize starting point");
			if (!std::isfinite(delta) || delta == 0.0)
				throw MultidimOptimizationInputError("NelderMead::Minimize: delta must be finite and non-zero");
			
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
         * @throws MultidimOptimizationInputError if inputs are invalid
         * 
         * Creates initial simplex by adding deltas[i] to coordinate i of start
         */
		template<int N>
		MultidimMinimizationResult Minimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start, const Vector<Real>& deltas) {
			if (deltas.size() != N)
				throw MultidimOptimizationError("Deltas vector dimension mismatch");

			// Validate inputs
			ValidateVectorFinite<N>(start, "NelderMead::Minimize starting point");
			for (int i = 0; i < N; ++i) {
				if (!std::isfinite(deltas[i]) || deltas[i] == 0.0)
					throw MultidimOptimizationInputError("NelderMead::Minimize: deltas[" + std::to_string(i) + 
						"] must be finite and non-zero");
			}

			// Create initial simplex: N+1 vertices
			_ndim = N;
			_mpts = N + 1;
			_p = Matrix<Real>(_mpts, _ndim);
			_y = Vector<Real>(_mpts);

			// First vertex is the starting point
			for (int j = 0; j < _ndim; ++j)
				_p(0, j) = start[j];

			// Other vertices: perturb one coordinate at a time
			for (int i = 1; i <= _ndim; ++i) {
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
         * @throws MultidimOptimizationInputError if simplex contains NaN/Inf values
         */
		template<int N>
		MultidimMinimizationResult Minimize(const IScalarFunction<N>& func, const Matrix<Real>& simplex) {
			if (simplex.rows() != N + 1 || simplex.cols() != N)
				throw MultidimOptimizationError("Simplex dimensions invalid: expected " + std::to_string(N + 1) + " x " +
												std::to_string(N));

			// Validate all simplex entries are finite
			for (int i = 0; i < N + 1; ++i) {
				for (int j = 0; j < N; ++j) {
					if (!std::isfinite(simplex(i, j)))
						throw MultidimOptimizationInputError("NelderMead::Minimize: simplex(" + std::to_string(i) + 
							"," + std::to_string(j) + ") is " + (std::isnan(simplex(i, j)) ? "NaN" : "Inf"));
				}
			}

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
		MultidimMinimizationResult MinimizeFromSimplex(const IScalarFunction<N>& func) {
			// Validate tolerance
			ValidateMultidimTolerance(_ftol, "NelderMead");
			if (_maxIter <= 0)
				throw MultidimOptimizationInputError("NelderMead: maxIter must be positive");
			
			const Real tiny = std::max(Real(1.0e-10), std::numeric_limits<Real>::epsilon());
			int ihi, ilo, inhi;
			VectorN<Real, N> x;
			Vector<Real> psum(_ndim);

			// Evaluate function at all simplex vertices
			for (int i = 0; i < _mpts; ++i) {
				for (int j = 0; j < _ndim; ++j)
					x[j] = _p(i, j);
				_y[i] = func(x);
				ValidateMultidimFunctionValue(_y[i], "NelderMead initial evaluation");
			}
			_nfunc = _mpts;

			// Compute initial psum (sum of vertex coordinates)
			GetPsum(psum);

			// Main iteration loop
			for (int iter = 0; iter < _maxIter; ++iter) {
				// Find lowest (best), highest (worst), and next-highest
				ilo = 0;
				if (_y[0] > _y[1]) {
					ihi = 0;
					inhi = 1;
				} else {
					ihi = 1;
					inhi = 0;
				}

				for (int i = 0; i < _mpts; ++i) {
					if (_y[i] <= _y[ilo])
						ilo = i;
					if (_y[i] > _y[ihi]) {
						inhi = ihi;
						ihi = i;
					} else if (_y[i] > _y[inhi] && i != ihi) {
						inhi = i;
					}
				}

				// Check convergence
				Real rtol = 2.0 * std::abs(_y[ihi] - _y[ilo]) / (std::abs(_y[ihi]) + std::abs(_y[ilo]) + tiny);

				if (rtol < _ftol) {
					// Converged - swap best to position 0
					std::swap(_y[0], _y[ilo]);
					for (int j = 0; j < _ndim; ++j) {
						std::swap(_p(0, j), _p(ilo, j));
						x[j] = _p(0, j);
					}

					Vector<Real> result(_ndim);
					for (int j = 0; j < _ndim; ++j)
						result[j] = x[j];

					return MultidimMinimizationResult(result, _y[0], _nfunc, true);
				}

				// Check iteration limit
				if (_nfunc >= _maxIter) {
					Vector<Real> result(_ndim);
					for (int j = 0; j < _ndim; ++j)
						result[j] = _p(ilo, j);
					return MultidimMinimizationResult(result, _y[ilo], _nfunc, false);
				}

				// Try reflection
				Real ytry = Amotry<N>(psum, ihi, -ALPHA, func);
				_nfunc++;

				if (ytry <= _y[ilo]) {
					// Reflection is better than best - try expansion
					ytry = Amotry<N>(psum, ihi, GAMMA, func);
					_nfunc++;
				} else if (ytry >= _y[inhi]) {
					// Reflection is worse than second-worst
					Real ysave = _y[ihi];
					ytry = Amotry<N>(psum, ihi, RHO, func);
					_nfunc++;

					if (ytry >= ysave) {
						// Contraction failed - do shrink
						for (int i = 0; i < _mpts; ++i) {
							if (i != ilo) {
								for (int j = 0; j < _ndim; ++j) {
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
				if (_y[i] < _y[ilo_final])
					ilo_final = i;

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
		void GetPsum(Vector<Real>& psum) {
			for (int j = 0; j < _ndim; ++j) {
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
		Real Amotry(Vector<Real>& psum, int ihi, Real fac, const IScalarFunction<N>& func) {
			VectorN<Real, N> ptry;
			Real fac1 = (1.0 - fac) / _ndim;
			Real fac2 = fac1 - fac;

			for (int j = 0; j < _ndim; ++j)
				ptry[j] = psum[j] * fac1 - _p(ihi, j) * fac2;

			Real ytry = func(ptry);

			if (ytry < _y[ihi]) {
				// Accept the new point
				_y[ihi] = ytry;
				for (int j = 0; j < _ndim; ++j) {
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
	MultidimMinimizationResult NelderMeadMinimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start, Real delta = 1.0,
												  Real ftol = PrecisionValues<Real>::OptimizationGradientTolerance) {
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
	MultidimMinimizationResult NelderMeadMinimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start, const Vector<Real>& deltas,
												  Real ftol = PrecisionValues<Real>::OptimizationGradientTolerance) {
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
	class NegatedScalarFunction : public IScalarFunction<N> {
	private:
		const IScalarFunction<N>& _func;

	public:
		explicit NegatedScalarFunction(const IScalarFunction<N>& func)
			: _func(func) {}

		Real operator()(const VectorN<Real, N>& x) const override { return -_func(x); }
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
	MultidimMinimizationResult NelderMeadMaximize(const IScalarFunction<N>& func, const VectorN<Real, N>& start, Real delta = 1.0,
												  Real ftol = PrecisionValues<Real>::OptimizationGradientTolerance) {
		NegatedScalarFunction<N> negFunc(func);
		NelderMead optimizer(ftol);
		auto result = optimizer.Minimize(negFunc, start, delta);
		result.fmin = -result.fmin; // Convert back to maximum value
		return result;
	}

	///////////////////////////////////////////////////////////////////////////
	///              Config-Based Overloads (API Standardization)           ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Minimize a scalar function using Nelder-Mead with config
     * @tparam N Dimension of the problem
     * @param func Function to minimize
     * @param start Starting point
     * @param config Algorithm configuration
     * @return Minimization result with enhanced diagnostics
     */
	template<int N>
	MultidimMinimizationResult NelderMeadMinimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start,
												  const MultidimOptimizationConfig& config) {
		AlgorithmTimer timer;
		NelderMead optimizer(config.tolerance, config.max_iterations);
		MultidimMinimizationResult result = optimizer.Minimize(func, start, config.initial_delta);
		
		// Fill enhanced diagnostic fields
		result.algorithm_name = "NelderMead";
		result.elapsed_time_ms = timer.elapsed_ms();
		result.function_evaluations = optimizer.getNumFuncEvals();
		if (!result.converged) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "Nelder-Mead did not converge within " + std::to_string(config.max_iterations) + " iterations";
		}
		return result;
	}

	/**
     * @brief Maximize a scalar function using Nelder-Mead with config
     * @tparam N Dimension of the problem
     * @param func Function to maximize
     * @param start Starting point
     * @param config Algorithm configuration
     * @return Maximization result with enhanced diagnostics (fmin is the maximum value)
     */
	template<int N>
	MultidimMinimizationResult NelderMeadMaximize(const IScalarFunction<N>& func, const VectorN<Real, N>& start,
												  const MultidimOptimizationConfig& config) {
		AlgorithmTimer timer;
		NegatedScalarFunction<N> negFunc(func);
		NelderMead optimizer(config.tolerance, config.max_iterations);
		MultidimMinimizationResult result = optimizer.Minimize(negFunc, start, config.initial_delta);
		
		result.fmin = -result.fmin; // Convert back to maximum value
		result.algorithm_name = "NelderMead";
		result.elapsed_time_ms = timer.elapsed_ms();
		result.function_evaluations = optimizer.getNumFuncEvals();
		if (!result.converged) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "Nelder-Mead did not converge within " + std::to_string(config.max_iterations) + " iterations";
		}
		return result;
	}


} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_NELDER_MEAD_H
