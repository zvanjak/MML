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

#include <mml/MMLBase.h>
#include <algorithm>

#include "DerivationBase.h"

#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/detail/NumericalJacobian.h>

#include <mml/core/Derivation/DerivationVectorFunction.h>

#include <functional>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace MML
{
	namespace Derivation
	{
		namespace Detail
		{
			using NumericalJacobianDetail::ValidateJacobianStep;
			using NumericalJacobianDetail::ValidateJacobianPoint;
			using NumericalJacobianDetail::ValidateJacobianStencil;
			using NumericalJacobianDetail::ValidateJacobianSample;
			using NumericalJacobianDetail::FinishJacobianQuotient;
			using NumericalJacobianDetail::JacobianDifference2;
			using NumericalJacobianDetail::JacobianCentralHalfSpan;

			static Real JacobianDifference4(Real plus2, Real plus1, Real minus1, Real minus2, Real h)
			{
				// Cancel common offsets before binary scaling; scale samples only if a difference overflows.
				Real inner = plus1 - minus1, outer = minus2 - plus2;
				int exponent;
				if (std::isfinite(inner) && std::isfinite(outer)) {
					std::frexp(std::max(std::abs(inner), std::abs(outer)), &exponent);
					inner = std::scalbn(inner, -exponent);
					outer = std::scalbn(outer, -exponent);
				}
				else {
					std::frexp(std::max({std::abs(plus2), std::abs(plus1), std::abs(minus1), std::abs(minus2)}), &exponent);
					inner = std::scalbn(plus1, -exponent) - std::scalbn(minus1, -exponent);
					outer = std::scalbn(minus2, -exponent) - std::scalbn(plus2, -exponent);
				}
				return FinishJacobianQuotient(Real(8) * inner + outer, exponent, h, Real(12));
			}

			template<int Order>
			static Real EvaluateJacobianPartialDyn(
				const std::function<Vector<Real>(const Vector<Real>&)>& func,
				int funcIndex, int derivIndex, const Vector<Real>& pos, Real h)
			{
				ValidateJacobianStep(h);
				if (!func || funcIndex < 0 || derivIndex < 0 || derivIndex >= pos.size())
					throw ArgumentError("Jacobian partial: invalid function or index");
				ValidateJacobianPoint(pos);
				if (h == 0)
					h = std::pow(std::numeric_limits<Real>::epsilon(), Real(1) / (Order + 1));
				ValidateJacobianStencil(pos[derivIndex], h, Order == 4 ? 2 : 1);
				Vector<Real> point = pos;
				int outputSize = -1;
				auto evaluate = [&](int offset) {
					point[derivIndex] = pos[derivIndex] + offset * h;
					const Vector<Real> value = func(point);
					if (funcIndex >= value.size())
						throw VectorDimensionError("Jacobian partial: missing output component", funcIndex + 1, value.size());
					if (outputSize < 0) outputSize = value.size();
					ValidateJacobianSample(value, outputSize);
					return value[funcIndex];
				};
				if constexpr (Order == 4) {
					const Real plus2 = evaluate(2), plus1 = evaluate(1), minus1 = evaluate(-1), minus2 = evaluate(-2);
					return JacobianDifference4(plus2, plus1, minus1, minus2, h);
				}
				else {
					const Real plus = evaluate(1), minus = evaluate(-1);
					return JacobianDifference2(plus, minus, h);
				}
			}
		}

		/********************************************************************************************************************/
		/********                          STATIC-SIZE JACOBIANS (VectorN, MatrixNM)                                 ********/
		/********************************************************************************************************************/

		/// @brief Calculate Jacobian matrix for a square vector function f: R^N -> R^N
		/// @tparam N Dimension of input and output
		/// @param func The vector function to differentiate
		/// @param pos Point at which to evaluate the Jacobian
		/// @return N×N Jacobian matrix where J(i,j) = ∂f_i/∂x_j
		/// Complexity: 4×N full-vector evaluations using one stencil per input column.
		template<int N>
		static MatrixNM<Real, N, N> CalcJacobian(const IVectorFunction<N>& func,
			const VectorN<Real, N>& pos, Real h = 0.0, int* functionEvaluations = nullptr)
		{
			MatrixNM<Real, N, N> jac;
			Detail::ValidateJacobianStep(h);
			Detail::ValidateJacobianPoint(pos);
			VectorN<Real, N> perturbed = pos;
			for (int column = 0; column < N; ++column) {
				const Real original = pos[column];
				const Real step = h > 0.0 ? h
					: std::pow(std::numeric_limits<Real>::epsilon(), Real(0.2))
						* std::max(Real(1.0), std::abs(original));
				Detail::ValidateJacobianStencil(original, step, 2);
				auto evaluate = [&](Real value) {
					perturbed[column] = value;
					if (functionEvaluations != nullptr) ++(*functionEvaluations);
					VectorN<Real, N> result = func(perturbed);
					Detail::ValidateJacobianSample(result, N);
					return result;
				};
				VectorN<Real, N> plus2 = evaluate(original + Real(2.0) * step);
				VectorN<Real, N> plus1 = evaluate(original + step);
				VectorN<Real, N> minus1 = evaluate(original - step);
				VectorN<Real, N> minus2 = evaluate(original - Real(2.0) * step);
				perturbed[column] = original;
				for (int row = 0; row < N; ++row)
					jac(row, column) = Detail::JacobianDifference4(plus2[row], plus1[row], minus1[row], minus2[row], step);
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
		static MatrixNM<Real, M, N> CalcJacobian(const IVectorFunctionNM<N, M>& func,
			const VectorN<Real, N>& pos, Real h = 0.0, int* functionEvaluations = nullptr)
		{
			MatrixNM<Real, M, N> jac;
			Detail::ValidateJacobianStep(h);
			Detail::ValidateJacobianPoint(pos);
			VectorN<Real, N> perturbed = pos;
			for (int column = 0; column < N; ++column) {
				const Real original = pos[column];
				const Real step = h > 0.0 ? h
					: std::pow(std::numeric_limits<Real>::epsilon(), Real(0.2))
						* std::max(Real(1.0), std::abs(original));
				Detail::ValidateJacobianStencil(original, step, 2);
				auto evaluate = [&](Real value) {
					perturbed[column] = value;
					if (functionEvaluations != nullptr) ++(*functionEvaluations);
					VectorN<Real, M> result = func(perturbed);
					Detail::ValidateJacobianSample(result, M);
					return result;
				};
				VectorN<Real, M> plus2 = evaluate(original + Real(2.0) * step);
				VectorN<Real, M> plus1 = evaluate(original + step);
				VectorN<Real, M> minus1 = evaluate(original - step);
				VectorN<Real, M> minus2 = evaluate(original - Real(2.0) * step);
				perturbed[column] = original;
				for (int row = 0; row < M; ++row)
					jac(row, column) = Detail::JacobianDifference4(plus2[row], plus1[row], minus1[row], minus2[row], step);
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
		/// @param h Absolute step; zero selects eps^(1/5) without coordinate scaling
		/// @return ∂f_funcIndex/∂x_derivIndex at pos
		/// 
		/// Uses 4th-order central difference formula:
		/// f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
		static inline Real NDer4PartialDyn(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			int funcIndex, int derivIndex,
			const Vector<Real>& pos,
			Real h = 0.0)
		{
			return Detail::EvaluateJacobianPartialDyn<4>(func, funcIndex, derivIndex, pos, h);
		}

		/// @brief Calculate numerical partial derivative using 2nd-order central difference
		/// @param func Function taking Vector<Real> and returning Vector<Real>
		/// @param funcIndex Output component index
		/// @param derivIndex Input variable index
		/// @param pos Point at which to evaluate
		/// @param h Absolute step; zero selects eps^(1/3) without coordinate scaling
		/// @return ∂f_funcIndex/∂x_derivIndex at pos
		/// 
		/// Uses central difference: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
		static Real NDer2PartialDyn(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			int funcIndex, int derivIndex,
			const Vector<Real>& pos,
			Real h = 0.0)
		{
			return Detail::EvaluateJacobianPartialDyn<2>(func, funcIndex, derivIndex, pos, h);
		}

		static void CalcJacobianColumns(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos, Matrix<Real>& jac, int outputDim,
			Real h, int* functionEvaluations);

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
		static Matrix<Real> CalcJacobian(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			Real h = 0.0,
			int* functionEvaluations = nullptr)
		{
			Matrix<Real> jac;
			CalcJacobianColumns(func, pos, jac, static_cast<int>(pos.size()), h,
				functionEvaluations);
			return jac;
		}

		/// @brief Calculate Jacobian matrix for a non-square dynamic vector function
		/// @param func Function f: R^n -> R^m
		/// @param pos Point at which to evaluate (determines n)
		/// @param outputDim Number of output components m
		/// @param h Step size for finite differences (0 = auto-select)
		/// @return m×n Jacobian matrix where J(i,j) = ∂f_i/∂x_j
		static Matrix<Real> CalcJacobian(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			int outputDim,
			Real h = 0.0,
			int* functionEvaluations = nullptr)
		{
			Matrix<Real> jac;
			CalcJacobianColumns(func, pos, jac, outputDim, h, functionEvaluations);
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
		static void CalcJacobianInPlace(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			Matrix<Real>& jac,
			Real h = 0.0,
			int* functionEvaluations = nullptr)
		{
			CalcJacobianColumns(func, pos, jac, static_cast<int>(pos.size()), h,
				functionEvaluations);
		}

		/// @brief Internal dynamic m-by-n Jacobian kernel using four vector evaluations per input column.
		/// @param func Function f: R^n -> R^m
		/// @param pos Point at which to evaluate
		/// @param[out] jac Output Jacobian, resized to m-by-n
		/// @param outputDim Expected output dimension m
		/// @param h Explicit absolute step; zero selects eps^(1/5)*max(1,|x_j|) per column
		/// @param functionEvaluations Optional accumulated count of attempted vector evaluations
		static void CalcJacobianColumns(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			Matrix<Real>& jac,
			int outputDim,
			Real h,
			int* functionEvaluations)
		{
			const int inputDim = static_cast<int>(pos.size());
			if (!func || inputDim <= 0 || outputDim <= 0 || !std::isfinite(h) || h < 0.0)
				throw ArgumentError("CalcJacobianInPlace: invalid function, dimensions, or step");
			Detail::ValidateJacobianPoint(pos);
			jac.resize(outputDim, inputDim);
			Vector<Real> perturbed = pos;
			for (int column = 0; column < inputDim; ++column) {
				const Real original = pos[column];
				const Real step = h > 0.0 ? h
					: std::pow(std::numeric_limits<Real>::epsilon(), Real(0.2))
						* std::max(Real(1.0), std::abs(original));
				Detail::ValidateJacobianStencil(original, step, 2);

				auto evaluate = [&](Real value) {
					perturbed[column] = value;
					if (functionEvaluations != nullptr) ++(*functionEvaluations);
					Vector<Real> result = func(perturbed);
					Detail::ValidateJacobianSample(result, outputDim);
					return result;
				};

				Vector<Real> plus2 = evaluate(original + Real(2.0) * step);
				Vector<Real> plus1 = evaluate(original + step);
				Vector<Real> minus1 = evaluate(original - step);
				Vector<Real> minus2 = evaluate(original - Real(2.0) * step);
				perturbed[column] = original;
				for (int row = 0; row < outputDim; ++row)
					jac(row, column) = Detail::JacobianDifference4(plus2[row], plus1[row], minus1[row], minus2[row], step);
			}
		}

		/// @brief Calculate Jacobian using per-entry 2nd-order differences
		/// @param func Function f: R^n -> R^n
		/// @param pos Point at which to evaluate
		/// @param h Step size (0 = auto)
		/// @return n×n Jacobian matrix
		/// Complexity: 2×n×n vector function evaluations (two per entry).
		/// 
		/// Useful when a lower-order formula is preferred for limited-accuracy callbacks.
		static Matrix<Real> CalcJacobianOrder2(
			const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos,
			Real h = 0.0)
		{
			int n = static_cast<int>(pos.size());
			Detail::ValidateJacobianStep(h);
			if (!func || n <= 0)
				throw ArgumentError("CalcJacobianOrder2: invalid function or dimensions");
			Detail::ValidateJacobianPoint(pos);
			Matrix<Real> jac(n, n);
			const std::function<Vector<Real>(const Vector<Real>&)> checkedFunction = [&](const Vector<Real>& x) {
				Vector<Real> value = func(x);
				Detail::ValidateJacobianSample(value, n);
				return value;
			};

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
				{
					jac(i, j) = NDer2PartialDyn(checkedFunction, i, j, pos, h);
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
		/// @param h Relative step factor; hj = h*max(1,|x_j|), zero selects cbrt(eps)
		/// Complexity: 2×n derivs() evaluations (central differences).
		///            For stiff solvers, this is called each step — dominates cost for large n.
		/// 
		/// This is specifically designed for ODE systems where the derivative
		/// function has the signature: void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		/// 
		/// The Jacobian is computed as J(i,j) = ∂(dxdt_i)/∂(x_j), holding t constant.
		static void CalcJacobianODE(
			const std::function<void(Real, const Vector<Real>&, Vector<Real>&)>& derivs,
			Real t,
			const Vector<Real>& x,
			Matrix<Real>& jac,
			Real h = 0.0)
		{
			if (!derivs)
				throw ArgumentError("CalcJacobianODE: derivative callback is empty");
			NumericalJacobianDetail::JacobianODE(derivs, t, x, jac, h);
		}

		/// @brief Calculate Jacobian for an ODE system, returning the matrix
		/// @param derivs The ODE right-hand side function
		/// @param t Current time
		/// @param x Current state vector
		/// @param h Relative step factor; zero selects cbrt(eps), then scales by max(1,|x_j|)
		/// @return n×n Jacobian matrix
		static Matrix<Real> CalcJacobianODE(
			const std::function<void(Real, const Vector<Real>&, Vector<Real>&)>& derivs,
			Real t,
			const Vector<Real>& x,
			Real h = 0.0)
		{
			Matrix<Real> jac;
			CalcJacobianODE(derivs, t, x, jac, h);
			return jac;
		}

		/********************************************************************************************************************/
		/********                          DISCRETE MAP JACOBIAN (x_next = f(x))                                     ********/
		/********************************************************************************************************************/

		/// @brief Calculate Jacobian for a discrete map x_{n+1} = f(x_n)
		/// @param mapFunc The map function
		/// @param x Current state
		/// @param[out] jac Output Jacobian (will be resized)
		/// @param h Relative step factor; zero selects cbrt(eps), then scales by max(1,|x_j|)
		/// 
		/// For discrete dynamical systems and iterated maps.
		static void CalcJacobianMap(
			const std::function<Vector<Real>(const Vector<Real>&)>& mapFunc,
			const Vector<Real>& x,
			Matrix<Real>& jac,
			Real h = 0.0)
		{
			int n = static_cast<int>(x.size());
			Detail::ValidateJacobianStep(h);
			if (!mapFunc || n <= 0)
				throw ArgumentError("CalcJacobianMap: invalid function or dimensions");
			Detail::ValidateJacobianPoint(x);
			jac.resize(n, n);

			Vector<Real> xp = x;
			Real eps = (h > 0.0) ? h : Precision::CbrtEps;

			for (int j = 0; j < n; ++j) {
				Real hj = eps * std::max<Real>(std::abs(x[j]), Real(1.0));
				Detail::ValidateJacobianStencil(x[j], hj, 1);
				const Real halfSpan = Detail::JacobianCentralHalfSpan(x[j], hj);

				xp[j] = x[j] + hj;
				Vector<Real> f_plus = mapFunc(xp);
				Detail::ValidateJacobianSample(f_plus, n);

				xp[j] = x[j] - hj;
				Vector<Real> f_minus = mapFunc(xp);
				Detail::ValidateJacobianSample(f_minus, n);

				xp[j] = x[j];

				for (int i = 0; i < n; ++i) {
					jac(i, j) = Detail::JacobianDifference2(f_plus[i], f_minus[i], halfSpan);
				}
			}
		}

		/********************************************************************************************************************/
		/********              Deprecated 2.0 spellings (CONTRACT_POLICY section 11; one minor release)              ********/
		/********************************************************************************************************************/

		template<int N>
		[[deprecated("use CalcJacobian()")]]
		static MatrixNM<Real, N, N> calcJacobian(const IVectorFunction<N>& func,
			const VectorN<Real, N>& pos, Real h = 0.0, int* functionEvaluations = nullptr)
		{ return CalcJacobian<N>(func, pos, h, functionEvaluations); }

		template<int N, int M>
		[[deprecated("use CalcJacobian()")]]
		static MatrixNM<Real, M, N> calcJacobian(const IVectorFunctionNM<N, M>& func,
			const VectorN<Real, N>& pos, Real h = 0.0, int* functionEvaluations = nullptr)
		{ return CalcJacobian<N, M>(func, pos, h, functionEvaluations); }

		[[deprecated("use CalcJacobian() - dynamic size is selected by the argument types")]]
		static inline Matrix<Real> calcJacobianDyn(const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos, Real h = 0.0, int* functionEvaluations = nullptr)
		{ return CalcJacobian(func, pos, h, functionEvaluations); }

		[[deprecated("use CalcJacobian() - dynamic size is selected by the argument types")]]
		static inline Matrix<Real> calcJacobianDyn(const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos, int outputDim, Real h = 0.0, int* functionEvaluations = nullptr)
		{ return CalcJacobian(func, pos, outputDim, h, functionEvaluations); }

		[[deprecated("use CalcJacobianInPlace()")]]
		static inline void calcJacobianDynInPlace(const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos, Matrix<Real>& jac, Real h = 0.0, int* functionEvaluations = nullptr)
		{ CalcJacobianInPlace(func, pos, jac, h, functionEvaluations); }

		[[deprecated("use CalcJacobianColumns()")]]
		static inline void calcJacobianDynColumns(const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos, Matrix<Real>& jac, int outputDim, Real h, int* functionEvaluations)
		{ CalcJacobianColumns(func, pos, jac, outputDim, h, functionEvaluations); }

		[[deprecated("use CalcJacobianOrder2()")]]
		static inline Matrix<Real> calcJacobianDyn2(const std::function<Vector<Real>(const Vector<Real>&)>& func,
			const Vector<Real>& pos, Real h = 0.0)
		{ return CalcJacobianOrder2(func, pos, h); }

		[[deprecated("use CalcJacobianODE()")]]
		static inline void calcJacobianODE(const std::function<void(Real, const Vector<Real>&, Vector<Real>&)>& derivs,
			Real t, const Vector<Real>& x, Matrix<Real>& jac, Real h = 0.0)
		{ CalcJacobianODE(derivs, t, x, jac, h); }

		[[deprecated("use CalcJacobianODE()")]]
		static inline Matrix<Real> calcJacobianODE(const std::function<void(Real, const Vector<Real>&, Vector<Real>&)>& derivs,
			Real t, const Vector<Real>& x, Real h = 0.0)
		{ return CalcJacobianODE(derivs, t, x, h); }

		[[deprecated("use CalcJacobianMap()")]]
		static inline void calcJacobianMap(const std::function<Vector<Real>(const Vector<Real>&)>& mapFunc,
			const Vector<Real>& x, Matrix<Real>& jac, Real h = 0.0)
		{ CalcJacobianMap(mapFunc, x, jac, h); }

	}  // namespace Derivation
}  // namespace MML

#endif