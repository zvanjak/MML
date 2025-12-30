///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Jacobians.h                                                         ///
///  Description: Jacobian matrix computations for coordinate transformations         ///
///               Numerical and analytical Jacobian evaluation                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_JACOBIANS_H
#define MML_DERIVATION_JACOBIANS_H

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

#include "core/Derivation/DerivationVectorFunction.h"

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                                      CALCULATING JACOBIANS                                         ********/
		/********************************************************************************************************************/
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
  };
}

#endif