///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Hessians.h                                                          ///
///  Description: Hessian matrix computations for scalar functions                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_HESSIANS_H
#define MML_DERIVATION_HESSIANS_H

#include <mml/MMLBase.h>

#include <mml/base/Matrix/MatrixNM.h>
#include <mml/core/Derivation/DerivationScalarFunction.h>

#include <cmath>

namespace MML
{
	namespace Derivation
	{
		/// @brief Calculate Hessian matrix in-place for a scalar function f: R^N -> R.
		/// @tparam N Dimension of the input space
		/// @param func The scalar function to differentiate
		/// @param point Point at which to evaluate the Hessian
		/// @param[out] hessian Output N x N Hessian matrix
		/// @param h Step size for finite differences (0 = existing per-entry automatic step)
		template<int N>
		static void calcHessian(const IScalarFunction<N>& func,
			const VectorN<Real, N>& point,
			MatrixNM<Real, N, N>& hessian,
			Real h = REAL(0.0))
		{
			if (!std::isfinite(h) || h < REAL(0.0))
				throw ArgumentError("calcHessian: invalid step");

			for (int row = 0; row < N; ++row)
			{
				for (int col = row; col < N; ++col)
				{
					const Real value = h > REAL(0.0)
						? NSecDer4Partial(func, row, col, point, h)
						: NSecDer4Partial(func, row, col, point);
					hessian(row, col) = value;
					if (row != col)
						hessian(col, row) = value;
				}
			}
		}

		/// @brief Calculate Hessian matrix for a scalar function f: R^N -> R.
		/// @tparam N Dimension of the input space
		/// @param func The scalar function to differentiate
		/// @param point Point at which to evaluate the Hessian
		/// @param h Step size for finite differences (0 = existing per-entry automatic step)
		/// @return N x N Hessian matrix where H(i,j) = d^2 f / dx_i dx_j
		template<int N>
		static MatrixNM<Real, N, N> calcHessian(const IScalarFunction<N>& func,
			const VectorN<Real, N>& point, Real h = REAL(0.0))
		{
			MatrixNM<Real, N, N> hessian;
			calcHessian(func, point, hessian, h);
			return hessian;
		}
	}
}

#endif // MML_DERIVATION_HESSIANS_H