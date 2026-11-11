///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Multidim/Powell.h                                                     ///
///  Description: Powell direction-set minimization                                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_POWELL_H
#define MML_OPTIMIZATION_POWELL_H

#include <mml/algorithms/Optimization/Multidim/LineSearch.h>

#include <mml/base/Function.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace MML::Optimization {

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
	class Powell {
	private:
		Real _ftol;	  // Fractional tolerance
		int _maxIter; // Maximum iterations
		int _iter;	  // Iteration counter
		Real _fret;	  // Current function value
		int _nfunc;	  // Objective callback attempts in the latest Minimize call

	public:
		Powell(Real ftol = PrecisionValues<Real>::OptimizationTolerance, int maxIter = 200)
			: _ftol(ftol)
			, _maxIter(maxIter)
			, _iter(0)
			, _fret(0.0)
			, _nfunc(0) {}

		Real getFtol() const { return _ftol; }
		void setFtol(Real ftol) { _ftol = ftol; }

		int getMaxIter() const { return _maxIter; }
		void setMaxIter(int maxIter) { _maxIter = maxIter; }

		int getIterations() const { return _iter; }
		Real getCurrentFValue() const { return _fret; }
		/// Includes nested line-search and throwing objective calls; reset before validation.
		int getNumFuncEvals() const { return _nfunc; }

		/**
         * @brief Minimize using starting point with identity direction matrix
         * @throws MultidimOptimizationInputError if inputs are invalid
         */
		template<int N>
		MultidimOptimizationResult Minimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start) {
			_nfunc = 0;
			// Validate inputs
			ValidateVectorFinite<N>(start, "Powell::Minimize starting point");
			ValidateMultidimTolerance(_ftol, "Powell");
			if (_maxIter <= 0)
				throw MultidimOptimizationInputError("Powell: maxIter must be positive");
			
			// Initialize direction matrix to identity
			Matrix<Real> ximat(N, N, 0.0);
			for (int i = 0; i < N; ++i)
				ximat(i, i) = 1.0;

			return Minimize(func, start, ximat);
		}

		/**
         * @brief Minimize with custom initial direction matrix
         * @param ximat Finite N-by-N initial directions (columns); updated during minimization
         * @throws MultidimOptimizationInputError if inputs are invalid
         */
		template<int N>
		MultidimOptimizationResult Minimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start, Matrix<Real>& ximat) {
			_nfunc = 0;
			// Validate inputs
			ValidateVectorFinite<N>(start, "Powell::Minimize starting point");
			ValidateMultidimTolerance(_ftol, "Powell");
			if (_maxIter <= 0)
				throw MultidimOptimizationInputError("Powell: maxIter must be positive");
			// Validate direction matrix
			if (ximat.rows() != N || ximat.cols() != N)
				throw MultidimOptimizationInputError("Powell::Minimize: direction matrix dimensions invalid: expected " +
					std::to_string(N) + " x " + std::to_string(N) + ", got " +
					std::to_string(ximat.rows()) + " x " + std::to_string(ximat.cols()));
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					if (!std::isfinite(ximat(i, j)))
						throw MultidimOptimizationInputError("Powell::Minimize: direction matrix contains NaN/Inf");
				}
			}
			
			const Real TINY = Precision::DivisionSafetyThreshold;
			ScalarFunctionFromStdFunc<N> countedObjective([&](const VectorN<Real, N>& point) {
				if (_nfunc == std::numeric_limits<int>::max())
					throw MultidimOptimizationError("Powell: function evaluation counter overflow");
				++_nfunc;
				return func(point);
			});

			VectorN<Real, N> p = start;
			VectorN<Real, N> pt, ptt, xi;

			_fret = countedObjective(p);
			ValidateMultidimFunctionValue(_fret, "Powell initial evaluation");

			// Save initial point
			for (int j = 0; j < N; ++j)
				pt[j] = p[j];

			for (_iter = 0; _iter < _maxIter; ++_iter) {
				Real fp = _fret;
				int ibig = 0;
				Real del = 0.0;

				// Minimize along each direction
				for (int i = 0; i < N; ++i) {
					// Extract i-th direction from matrix columns
					for (int j = 0; j < N; ++j)
						xi[j] = ximat(j, i);

					Real fptt = _fret;
					_fret = LineMinimizer::Minimize(countedObjective, p, xi, _ftol);

					// Track direction with largest decrease
					if (fptt - _fret > del) {
						del = fptt - _fret;
						ibig = i + 1;
					}
				}

				// Check convergence
				if (2.0 * (fp - _fret) <= _ftol * (std::abs(fp) + std::abs(_fret)) + TINY) {
					Vector<Real> result(N);
					for (int j = 0; j < N; ++j)
						result[j] = p[j];
					MultidimOptimizationResult output(result, _fret, _iter + 1, true);
					output.function_evaluations = _nfunc;
					return output;
				}

				// Construct extrapolated point and new direction
				for (int j = 0; j < N; ++j) {
					ptt[j] = MML::Real(2.0) * p[j] - pt[j]; // Extrapolated point
					xi[j] = p[j] - pt[j];		 // New direction
					pt[j] = p[j];				 // Save current point
				}

				Real fptt = countedObjective(ptt);

				if (fptt < fp) {
					Real t = MML::Real(2.0) * (fp - MML::Real(2.0) * _fret + fptt) * (fp - _fret - del) * (fp - _fret - del) - del * (fp - fptt) * (fp - fptt);

					if (t < 0.0) {
						// Replace direction ibig with new direction
						_fret = LineMinimizer::Minimize(countedObjective, p, xi, _ftol);

						for (int j = 0; j < N; ++j) {
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
			MultidimOptimizationResult output(result, _fret, _iter, false);
			output.function_evaluations = _nfunc;
			return output;
		}
	};


} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_POWELL_H
