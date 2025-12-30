///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixUtils.h                                                       ///
///  Description: Matrix utility functions (special matrix generation)                ///
///               Vandermonde, companion, Hilbert, rotation matrices                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIXUTILS_H
#define MML_MATRIXUTILS_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Matrix.h"
#include "base/MatrixNM.h"
#include "base/Polynom.h"

#include "core/LinAlgEqSolvers.h"

namespace MML
{
	namespace Utils
	{
		// function returning whether Matrix is nilpotent
		template<class Type>
		bool IsNilpotent(const Matrix<Type>& A)
		{
			int N = A.RowNum();
			if (N != A.ColNum())
				throw MatrixDimensionError("IsNilpotent - must be square matrix", N, A.ColNum(), -1, -1);
			
			Matrix<Type> M = A;
			for (int k = 1; k < N; k++)
			{
				M = A * M;
				if (M.IsZero())
					return true;
			}
			return false;
		}

		template<class Type>
		bool IsUnipotent(const Matrix<Type>& A)
		{
			int N = A.RowNum();
			if (N != A.ColNum())
				throw MatrixDimensionError("IsUnipotent - must be square matrix", N, A.ColNum(), -1, -1);

			Matrix<Type> B = A - Matrix<Type>::GetUnitMatrix(N);
			Matrix<Type> M = B;
			for (int k = 1; k < N; k++)
			{
				if (M.IsZero())
					return true;
				M = B * M;
			}
			return M.IsZero();
		}

		// Faddeev-Leverrier algorithm for characteristic polynomial, determinant and inverse matrix
		void FaddeevAlg(const Matrix<Real>& A, PolynomRealFunc& outCharPoly, Real& outDet, Matrix<Real>& outInv)
		{
			int N = A.RowNum();

			Matrix<Real> nullMat(N, N);		// null matrix
			auto identMat = Matrix<Real>::GetUnitMatrix(N);

			Vector<Matrix<Real>> M(N + 1, nullMat);

			outCharPoly.SetDegree(N);

			M[0] = nullMat;
			outCharPoly[N] = 1.0;
			for (int k = 1; k <= A.RowNum(); k++)
			{
				M[k] = A * M[k - 1] + outCharPoly[N - k + 1] * identMat;
				outCharPoly[N - k] = -1.0 / k * (A * M[k]).Trace();
			}

			outDet = std::pow(-1.0, N) * outCharPoly[0];
			outInv = M[N] / outDet; 
		}
	
		// Using LU decomposition to compute the determinant
		template<class Type>
		Type Det(const Matrix<Type>& A)
		{
			try
			{
				LUSolver<Type> lu(A);

				return lu.det();
			}
			catch(const SingularMatrixError&)
			{
				return 0.0;
			}
		}

		// Calculating the rank of a matrix using Gaussian elimination
		template<class Type>
		int Rank(const Matrix<Type>& A, Real EPS = Defaults::RankAlgEPS)
		{
			int rows = A.RowNum();
			int cols = A.ColNum();
			Matrix<Type> mat = A;			// Make a copy to preserve the original

			int rank = 0;
			std::vector<bool> row_selected(rows, false);

			for (int col = 0; col < cols; ++col)
			{
				int pivot_row = -1;
				for (int row = 0; row < rows; ++row)
				{
					if (!row_selected[row] && std::abs(mat(row, col)) > EPS)
					{
						pivot_row = row;
						break;
					}
				}
				if (pivot_row == -1)
					continue;

				++rank;
				row_selected[pivot_row] = true;

				// Eliminate below
				for (int row = 0; row < rows; ++row)
				{
					if (row != pivot_row)
					{
						Type factor = mat(row, col) / mat(pivot_row, col);
						for (int k = col; k < cols; ++k)
							mat(row, k) -= factor * mat(pivot_row, k);
					}
				}
			}
			return rank;
		}
	}
} // namespace MML
#endif


