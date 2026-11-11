///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FieldMatrixAlgorithms.h                                             ///
///  Description: Exact fixed-size matrix algorithms over fields                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FIELD_MATRIX_ALGORITHMS_H
#define MML_FIELD_MATRIX_ALGORITHMS_H

#include <mml/MMLExceptions.h>
#include <mml/base/Matrix/MatrixNM.h>

#include <utility>

namespace MML
{
	namespace Algebra
	{
		template<class Field, int N>
		Field FieldMatrixDeterminant(MatrixNM<Field, N, N> matrix)
		{
			Field determinant(1);
			for (int column = 0; column < N; column++) {
				int pivot = column;
				while (pivot < N && matrix(pivot, column) == Field(0))
					pivot++;
				if (pivot == N)
					return Field(0);

				if (pivot != column) {
					for (int index = 0; index < N; index++)
						std::swap(matrix(column, index), matrix(pivot, index));
					determinant = -determinant;
				}

				const Field pivotValue = matrix(column, column);
				determinant *= pivotValue;
				for (int row = column + 1; row < N; row++) {
					const Field factor = matrix(row, column) / pivotValue;
					for (int index = column; index < N; index++)
						matrix(row, index) -= factor * matrix(column, index);
				}
			}
			return determinant;
		}

		template<class Field, int N>
		MatrixNM<Field, N, N> FieldMatrixInverse(MatrixNM<Field, N, N> matrix)
		{
			MatrixNM<Field, N, N> inverse;
			for (int index = 0; index < N; index++)
				inverse(index, index) = Field(1);

			for (int column = 0; column < N; column++) {
				int pivot = column;
				while (pivot < N && matrix(pivot, column) == Field(0))
					pivot++;
				if (pivot == N)
					throw DomainError("Field matrix is singular");

				if (pivot != column)
					for (int index = 0; index < N; index++) {
						std::swap(matrix(column, index), matrix(pivot, index));
						std::swap(inverse(column, index), inverse(pivot, index));
					}

				const Field pivotInverse = Field(1) / matrix(column, column);
				for (int index = 0; index < N; index++) {
					matrix(column, index) *= pivotInverse;
					inverse(column, index) *= pivotInverse;
				}

				for (int row = 0; row < N; row++) {
					if (row == column)
						continue;
					const Field factor = matrix(row, column);
					for (int index = 0; index < N; index++) {
						matrix(row, index) -= factor * matrix(column, index);
						inverse(row, index) -= factor * inverse(column, index);
					}
				}
			}
			return inverse;
		}
	}
}

#endif // MML_FIELD_MATRIX_ALGORITHMS_H