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

#include <algorithm>
#include <cmath>

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

	public:
		Powell(Real ftol = PrecisionValues<Real>::OptimizationTolerance, int maxIter = 200)
			: _ftol(ftol)
			, _maxIter(maxIter)
			, _iter(0)
			, _fret(0.0) {}

		Real getFtol() const { return _ftol; }
		void setFtol(Real ftol) { _ftol = ftol; }

		int getMaxIter() const { return _maxIter; }
		void setMaxIter(int maxIter) { _maxIter = maxIter; }

		int getIterations() const { return _iter; }
		Real getCurrentFValue() const { return _fret; }

		/**
         * @brief Minimize using starting point with identity direction matrix
         * @throws MultidimOptimizationInputError if inputs are invalid
         */
		template<int N>
		MultidimMinimizationResult Minimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start) {
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
         * @throws MultidimOptimizationInputError if inputs are invalid
         */
		template<int N>
		MultidimMinimizationResult Minimize(const IScalarFunction<N>& func, const VectorN<Real, N>& start, Matrix<Real>& ximat) {
			// Validate inputs
			ValidateVectorFinite<N>(start, "Powell::Minimize starting point");
			ValidateMultidimTolerance(_ftol, "Powell");
			if (_maxIter <= 0)
				throw MultidimOptimizationInputError("Powell: maxIter must be positive");
			// Validate direction matrix
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					if (!std::isfinite(ximat(i, j)))
						throw MultidimOptimizationInputError("Powell::Minimize: direction matrix contains NaN/Inf");
				}
			}
			
			const Real TINY = 1.0e-25;

			VectorN<Real, N> p = start;
			VectorN<Real, N> pt, ptt, xi;

			_fret = func(p);
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
					_fret = LineMinimizer::Minimize(func, p, xi, _ftol);

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
					return MultidimMinimizationResult(result, _fret, _iter + 1, true);
				}

				// Construct extrapolated point and new direction
				for (int j = 0; j < N; ++j) {
					ptt[j] = 2.0 * p[j] - pt[j]; // Extrapolated point
					xi[j] = p[j] - pt[j];		 // New direction
					pt[j] = p[j];				 // Save current point
				}

				Real fptt = func(ptt);

				if (fptt < fp) {
					Real t = 2.0 * (fp - 2.0 * _fret + fptt) * (fp - _fret - del) * (fp - _fret - del) - del * (fp - fptt) * (fp - fptt);

					if (t < 0.0) {
						// Replace direction ibig with new direction
						_fret = LineMinimizer::Minimize(func, p, xi, _ftol);

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
			return MultidimMinimizationResult(result, _fret, _iter, false);
		}
	};


} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_POWELL_H
