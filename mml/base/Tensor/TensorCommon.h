///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        TensorCommon.h                                                     ///
///  Description: Shared implementation helpers for concrete tensors                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_TENSOR_COMMON_H
#define MML_TENSOR_COMMON_H

#include <array>

#include <mml/MMLBase.h>
#include <mml/interfaces/ITensor.h>
#include <mml/base/Matrix/MatrixNM.h>

namespace MML
{
	namespace detail
	{
		template<int Rank>
		void ValidateTensorVarianceCounts(int covariant, int contravariant)
		{
			if (covariant < 0 || covariant > Rank || contravariant < 0 || contravariant > Rank ||
				contravariant != Rank - covariant)
				throw TensorCovarContravarNumError("Tensor variance counts must be nonnegative and sum to the rank",
					contravariant, covariant);
		}

		template<int Rank>
		inline void ValidateTensorContractionIndices(int first, int second)
		{
			if (first < 0 || first >= Rank || second < 0 || second >= Rank)
				throw TensorIndexError("Tensor contraction - index position out of range");
			if (first == second)
				throw TensorIndexError("Tensor contraction - index positions must differ");
		}

		template<int N>
		Real DeterminantForLeviCivitaTensor(const MatrixNM<Real, N, N>& matrix)
		{
			static_assert(N == 3 || N == 4, "DeterminantForLeviCivitaTensor is implemented only for N=3 and N=4");
			if constexpr (N == 3)
			{
				return matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1))
					 - matrix(0, 1) * (matrix(1, 0) * matrix(2, 2) - matrix(1, 2) * matrix(2, 0))
					 + matrix(0, 2) * (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0));
			}
			else if constexpr (N == 4)
			{
				Real det = static_cast<MML::Real>(0.0);
				for (int col = 0; col < 4; col++)
				{
					MatrixNM<Real, 3, 3> minor;
					for (int i = 1; i < 4; i++)
					{
						int minorCol = 0;
						for (int j = 0; j < 4; j++)
						{
							if (j == col)
								continue;
							minor(i - 1, minorCol++) = matrix(i, j);
						}
					}
					Real cofactorSign = (col % 2 == 0) ? static_cast<MML::Real>(1.0) : -static_cast<MML::Real>(1.0);
					det += cofactorSign * matrix(0, col) * DeterminantForLeviCivitaTensor<3>(minor);
				}
				return det;
			}
		}

		template<int Rank>
		void ValidateTensorPermutation(const std::array<int, Rank>& permutation)
		{
			bool seen[Rank] = { false };
			for (int i = 0; i < Rank; i++)
			{
				int sourceIndex = permutation[i];
				if (sourceIndex < 0 || sourceIndex >= Rank)
					throw TensorIndexError("Tensor index permutation contains an out-of-range index");
				if (seen[sourceIndex])
					throw TensorIndexError("Tensor index permutation contains a duplicate index");
				seen[sourceIndex] = true;
			}
		}
	}
}
#endif
