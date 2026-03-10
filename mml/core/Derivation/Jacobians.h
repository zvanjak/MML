///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Jacobians.h                                                         ///
///  Description: Jacobian matrix computations for vector functions                   ///
///               Supports both static (VectorN) and dynamic (Vector) sizes           ///
///               Numerical differentiation using 4th-order central differences       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_JACOBIANS_H
#define MML_DERIVATION_JACOBIANS_H

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "base/Vector/VectorN.h"
#include "base/Matrix/MatrixNM.h"

#include "core/Derivation/DerivationVectorFunction.h"

#include <functional>
#include <cmath>
#include <limits>

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                          STATIC-SIZE JACOBIANS (VectorN, MatrixNM)                                 ********/
		/********************************************************************************************************************/

		/// @brief Calculate Jacobian matrix for a square vector function f: R^N -> R^N
		/// @tparam N Dimension of input and output
		/// @param func The vector function to differentiate
		/// @param pos Point at which to evaluate the Jacobian
		/// @return N×N Jacobian matrix where J(i,j) = ∂f_i/∂x_j
		/// Complexity: O(N²) partial derivatives, each using 4th-order stencil (4 f-evals).
		///            Total: 4×N² scalar function evaluations.
		template<int N>
		static MatrixNM<Real, N, N> calcJacobian(const IVectorFunction<N>& func, const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, N, N> jac;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = NDer4Partial(func, i, j, pos);
				}

			return jac;
		}

		/// @brief Calculate Jacobian matrix for a general vector function f: R^N -> R^M
		/// @tparam N Dimension of input
		/// @tparam M Dimension of output
		/// @param func The vector function to differentiate
		/// @param pos Point at which to evaluate the Jacobian
		/// @return M×N Jacobian matrix where J(i,j) = ∂f_i/∂x_j
		template<int N, int M>
		static MatrixNM<Real, M, N> calcJacobian(const IVectorFunctionNM<N, M>& func, const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, M, N> jac;

			for (int i = 0; i < M; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(func, i, j, pos);
				}

			return jac;
		}

		/********************************************************************************************************************/
		/********                          DYNAMIC-SIZE JACOBIANS (Vector, Matrix)                                   ********/
		/********************************************************************************************************************/

		/// @brief Calculate numerical partial derivative for a dynamic vector function
		/// @param func Function taking Vector<Real> and returning Vector<Real>
		/// @param funcIndex Output component index (which f_i to differentiate)
		/// @param derivIndex Input variable index (differentiate with respect to x_j)
		/// @param pos Point at which to evaluate the derivative
		/// @param h Step size for finite differences (optional, defaults to optimal)
		/// @return ∂f_funcIndex/∂x_derivIndex at pos
		/// 
		/// Uses 4th-order central difference formula:
		/// f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
		static Real NDer4PartialDyn(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			int funcIndex, int derivIndex,
			const Vector<Real>& pos,
			Real h = 0.0)
		{
			// Optimal step size for 4th order method: h ≈ ε^(1/5) where ε is machine epsilon
			if (h == 0.0)
				h = std::pow(std::numeric_limits<Real>::epsilon(), 0.2);

			Vector<Real> x = pos;
			Real x_orig = x[derivIndex];

			// 4th-order central difference coefficients: (-1, +8, -8, +1) / 12h
			x[derivIndex] = x_orig + 2 * h;
			Real f_plus2h = func(x)[funcIndex];

			x[derivIndex] = x_orig + h;
			Real f_plus_h = func(x)[funcIndex];

			x[derivIndex] = x_orig - h;
			Real f_minus_h = func(x)[funcIndex];

			x[derivIndex] = x_orig - 2 * h;
			Real f_minus2h = func(x)[funcIndex];

			return (-f_plus2h + 8 * f_plus_h - 8 * f_minus_h + f_minus2h) / (12 * h);
		}

		/// @brief Calculate numerical partial derivative using 2nd-order central difference
		/// @param func Function taking Vector<Real> and returning Vector<Real>
		/// @param funcIndex Output component index
		/// @param derivIndex Input variable index
		/// @param pos Point at which to evaluate
		/// @param h Step size (optional)
		/// @return ∂f_funcIndex/∂x_derivIndex at pos
		/// 
		/// Uses central difference: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
		static Real NDer2PartialDyn(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			int funcIndex, int derivIndex,
			const Vector<Real>& pos,
			Real h = 0.0)
		{
			// Optimal step size for 2nd order: h ≈ ε^(1/3)
			if (h == 0.0)
				h = std::pow(std::numeric_limits<Real>::epsilon(), 1.0 / 3.0);

			Vector<Real> x = pos;
			Real x_orig = x[derivIndex];

			x[derivIndex] = x_orig + h;
			Real f_plus = func(x)[funcIndex];

			x[derivIndex] = x_orig - h;
			Real f_minus = func(x)[funcIndex];

			return (f_plus - f_minus) / (2 * h);
		}

		/// @brief Calculate Jacobian matrix for a dynamic-size vector function (square case)
		/// @param func Function f: R^n -> R^n taking and returning Vector<Real>
		/// @param pos Point at which to evaluate the Jacobian
		/// @param h Step size for finite differences (0 = auto-select optimal)
		/// @return n×n Jacobian matrix where J(i,j) = ∂f_i/∂x_j
		/// Complexity: 4×n vector function evaluations (4th-order central differences).
		///            Each vector f-eval costs O(n), so total is O(4n × n) = O(n² × f-eval cost).
		/// 
		/// This is the main entry point for dynamic-size Jacobian computation.
		/// Uses 4th-order accurate central differences for high precision.
		static Matrix<Real> calcJacobianDyn(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			Real h = 0.0)
		{
			int n = static_cast<int>(pos.size());
			Matrix<Real> jac(n, n);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
				{
					jac(i, j) = NDer4PartialDyn(func, i, j, pos, h);
				}

			return jac;
		}

		/// @brief Calculate Jacobian matrix for a non-square dynamic vector function
		/// @param func Function f: R^n -> R^m
		/// @param pos Point at which to evaluate (determines n)
		/// @param outputDim Number of output components m
		/// @param h Step size for finite differences (0 = auto-select)
		/// @return m×n Jacobian matrix where J(i,j) = ∂f_i/∂x_j
		static Matrix<Real> calcJacobianDyn(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			int outputDim,
			Real h = 0.0)
		{
			int n = static_cast<int>(pos.size());
			int m = outputDim;
			Matrix<Real> jac(m, n);

			for (int i = 0; i < m; ++i)
				for (int j = 0; j < n; ++j)
				{
					jac(i, j) = NDer4PartialDyn(func, i, j, pos, h);
				}

			return jac;
		}

		/// @brief Calculate Jacobian in-place for maximum efficiency
		/// @param func Function f: R^n -> R^n
		/// @param pos Point at which to evaluate
		/// @param[out] jac Pre-allocated Jacobian matrix (must be n×n)
		/// @param h Step size (0 = auto)
		/// 
		/// This version avoids allocation by writing directly to a pre-sized matrix.
		/// Useful in iterative algorithms (Newton-Raphson, etc.) where Jacobian
		/// is recomputed many times.
		static void calcJacobianDynInPlace(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			Matrix<Real>& jac,
			Real h = 0.0)
		{
			int n = static_cast<int>(pos.size());
			jac.Resize(n, n);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
				{
					jac(i, j) = NDer4PartialDyn(func, i, j, pos, h);
				}
		}

		/// @brief Calculate Jacobian using 2nd-order differences (faster, less accurate)
		/// @param func Function f: R^n -> R^n
		/// @param pos Point at which to evaluate
		/// @param h Step size (0 = auto)
		/// @return n×n Jacobian matrix
		/// Complexity: 2×n vector function evaluations (central differences).
		///            Half the cost of 4th-order method.
		/// 
		/// Use this when speed is more important than precision, or when the
		/// function itself has limited accuracy.
		static Matrix<Real> calcJacobianDyn2(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			Real h = 0.0)
		{
			int n = static_cast<int>(pos.size());
			Matrix<Real> jac(n, n);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
				{
					jac(i, j) = NDer2PartialDyn(func, i, j, pos, h);
				}

			return jac;
		}

		/********************************************************************************************************************/
		/********                       ODE SYSTEM JACOBIAN (derivs-style interface)                                 ********/
		/********************************************************************************************************************/

		/// @brief Calculate Jacobian for an ODE system with derivs(t, x, dxdt) interface
		/// @param derivs The ODE right-hand side function
		/// @param t Current time
		/// @param x Current state vector
		/// @param[out] jac Output Jacobian matrix (will be resized to n×n)
		/// @param h Step size (0 = auto-select based on x values)
		/// Complexity: 2×n derivs() evaluations (central differences).
		///            For stiff solvers, this is called each step — dominates cost for large n.
		/// 
		/// This is specifically designed for ODE systems where the derivative
		/// function has the signature: void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		/// 
		/// The Jacobian is computed as J(i,j) = ∂(dxdt_i)/∂(x_j), holding t constant.
		static void calcJacobianODE(
			const std::function<void(Real, const Vector<Real>&, Vector<Real>&)>& derivs,
			Real t,
			const Vector<Real>& x,
			Matrix<Real>& jac,
			Real h = 0.0)
		{
			int n = static_cast<int>(x.size());
			jac.Resize(n, n);

			Vector<Real> f_plus(n), f_minus(n);
			Vector<Real> xp = x;

			// Auto-select step size if not specified
			Real eps = (h > 0.0) ? h : std::sqrt(std::numeric_limits<Real>::epsilon());

			for (int j = 0; j < n; ++j) {
				// Adaptive step size based on x[j] magnitude
				Real hj = eps * std::max(std::abs(x[j]), 1.0);

				// Central difference
				xp[j] = x[j] + hj;
				derivs(t, xp, f_plus);

				xp[j] = x[j] - hj;
				derivs(t, xp, f_minus);

				xp[j] = x[j];  // Restore

				for (int i = 0; i < n; ++i) {
					jac(i, j) = (f_plus[i] - f_minus[i]) / (2.0 * hj);
				}
			}
		}

		/// @brief Calculate Jacobian for an ODE system, returning the matrix
		/// @param derivs The ODE right-hand side function
		/// @param t Current time
		/// @param x Current state vector
		/// @param h Step size (0 = auto)
		/// @return n×n Jacobian matrix
		static Matrix<Real> calcJacobianODE(
			const std::function<void(Real, const Vector<Real>&, Vector<Real>&)>& derivs,
			Real t,
			const Vector<Real>& x,
			Real h = 0.0)
		{
			Matrix<Real> jac;
			calcJacobianODE(derivs, t, x, jac, h);
			return jac;
		}

		/********************************************************************************************************************/
		/********                          DISCRETE MAP JACOBIAN (x_next = f(x))                                     ********/
		/********************************************************************************************************************/

		/// @brief Calculate Jacobian for a discrete map x_{n+1} = f(x_n)
		/// @param mapFunc The map function
		/// @param x Current state
		/// @param[out] jac Output Jacobian (will be resized)
		/// @param h Step size (0 = auto)
		/// 
		/// For discrete dynamical systems and iterated maps.
		static void calcJacobianMap(
			const std::function<Vector<Real>(const Vector<Real>&)>& mapFunc,
			const Vector<Real>& x,
			Matrix<Real>& jac,
			Real h = 0.0)
		{
			int n = static_cast<int>(x.size());
			jac.Resize(n, n);

			Vector<Real> xp = x;
			Real eps = (h > 0.0) ? h : std::sqrt(std::numeric_limits<Real>::epsilon());

			for (int j = 0; j < n; ++j) {
				Real hj = eps * std::max(std::abs(x[j]), 1.0);

				xp[j] = x[j] + hj;
				Vector<Real> f_plus = mapFunc(xp);

				xp[j] = x[j] - hj;
				Vector<Real> f_minus = mapFunc(xp);

				xp[j] = x[j];

				for (int i = 0; i < n; ++i) {
					jac(i, j) = (f_plus[i] - f_minus[i]) / (2.0 * hj);
				}
			}
		}

	}  // namespace Derivation
}  // namespace MML

#endif