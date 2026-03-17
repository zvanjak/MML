///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgQR.h                                                          ///
///  Description: QR decomposition solver using Householder reflections               ///
///               Supports least squares, determinant, and matrix inverse             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_LINEAR_ALG_QR_H
#define MML_LINEAR_ALG_QR_H

#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "base/BaseUtils.h"

namespace MML
{
	/// @brief QR decomposition solver using Householder reflections
	/// @tparam Type Numeric type (Real, Complex, etc.)
	/// @note Decomposes A=QR where Q is orthogonal and R is upper triangular
	/// @note Works for square (m=n) and overdetermined (m>n) systems
	/// @note Complexity: O(2mn²-2n³/3) decomposition, O(mn) per solve
	template<class Type>
	class QRSolver
	{
	private:
		int _m;              // number of rows
		int _n;              // number of columns
		Matrix<Type> _QR;    // Combined QR storage: R in upper triangle, Householder vectors in lower
		Vector<Type> _c;     // Diagonal of R (stored separately)
		Vector<Type> _d;     // Householder scaling factors
		bool _sing;          // Singularity flag
		int _num_reflections; // Number of Householder reflections performed (for determinant sign)

	public:
		/// @brief Get number of rows of the decomposed matrix
		int rows() const { return _m; }
		/// @brief Get number of columns of the decomposed matrix
		int cols() const { return _n; }
		/// @brief Check if the decomposed matrix is singular
		bool isSingular() const { return _sing; }

		/// @brief Constructor - performs QR decomposition using Householder reflections
		/// @param a Input matrix (m×n with m≥n)
		/// @throws MatrixDimensionError if m<n
		/// @note Stores R in upper triangle of QR, Householder vectors in lower triangle
		QRSolver(const Matrix<Type>& a) : _m(a.rows()), _n(a.cols()), _QR(a), _c(_n), _d(_n), _sing(false), _num_reflections(0)
		{
			if (_m < _n)
				throw MatrixDimensionError("QRSolver: Matrix must have m >= n (rows >= columns)", _m, _n, _m, _n);

			int i, j, k;
			Type scale, sigma, sum, tau;
			Vector<Type> vec(_m);

			// Perform Householder reduction
			// For square matrices: process n-1 columns (last has no elements below diagonal)
			// For overdetermined: process all n columns (last column still has elements below diagonal)
			int ncols = (_m == _n) ? _n - 1 : _n;
			for (k = 0; k < ncols; k++)
			{
				// Compute the norm of the k-th column below the diagonal
				scale = 0.0;
				for (i = k; i < _m; i++)
					scale = std::max(scale, Abs(_QR[i][k]));

				if (scale < std::numeric_limits<Real>::epsilon())
				{
					// Singular case: column below diagonal is effectively zero
					_sing = true;
					_c[k] = _d[k] = 0.0;
				}
				else
				{
					// Form the Householder vector
					for (i = k; i < _m; i++)
						_QR[i][k] /= scale;

					sum = 0.0;
					for (i = k; i < _m; i++)
						sum += _QR[i][k] * _QR[i][k];

					// Choose sign of sigma to avoid cancellation errors
					sigma = (_QR[k][k] >= 0.0 ? std::sqrt(sum) : -std::sqrt(sum));
					_QR[k][k] += sigma;
					_c[k] = sigma * _QR[k][k];
					_d[k] = -scale * sigma;

					// Count this reflection for determinant
					_num_reflections++;

					// Apply the transformation to remaining columns
					for (j = k + 1; j < _n; j++)
					{
						sum = 0.0;
						for (i = k; i < _m; i++)
							sum += _QR[i][k] * _QR[i][j];

						tau = sum / _c[k];

						for (i = k; i < _m; i++)
							_QR[i][j] -= tau * _QR[i][k];
					}
				}
			}

			// For square matrices, handle the last diagonal element separately
			if (_m == _n)
			{
				_d[_n - 1] = _QR[_n - 1][_n - 1];
				// Use norm-scaled threshold instead of exact zero
				Real norm_a = 0.0;
				for (int ii = 0; ii < _m; ii++)
					for (int jj = 0; jj < _n; jj++)
						if (Abs(a(ii, jj)) > norm_a)
							norm_a = Abs(a(ii, jj));
				Real threshold = std::numeric_limits<Real>::epsilon() * norm_a * _n;
				if (Abs(_d[_n - 1]) < threshold)
					_sing = true;
			}
		}

		// Solves the set of n linear equations A*x = b for a square system (m == n)
		// Uses back-substitution on R after applying Q^T to b
		void Solve(const Vector<Type>& b, Vector<Type>& x)
		{
			if (_m != _n)
				throw MatrixDimensionError("QRSolver::Solve - Use LeastSquaresSolve for overdetermined systems (m > n)", _m, _n, _m, _m);
			if (b.size() != static_cast<size_t>(_m))
				throw VectorDimensionError("QRSolver::Solve - Vector size mismatch", _m, b.size());
			if (_sing)
				throw SingularMatrixError("QRSolver::Solve - Singular matrix");

			// Apply Q^T to b
			QtMultiply(b, x);

			// Back-substitution on R
			RSolve(x, x);
		}

		Vector<Type> Solve(const Vector<Type>& b)
		{
			Vector<Type> x(_n);
			Solve(b, x);
			return x;
		}

		// Least-squares solution for overdetermined systems (m > n)
		// Minimizes ||Ax - b||₂ by solving R*x = Q^T*b
		void LeastSquaresSolve(const Vector<Type>& b, Vector<Type>& x)
		{
			if (b.size() != static_cast<size_t>(_m))
				throw VectorDimensionError("QRSolver::LeastSquaresSolve - Vector size mismatch", _m, b.size());
			if (_sing)
				throw SingularMatrixError("QRSolver::LeastSquaresSolve - Singular matrix");

			// Apply Q^T to b
			Vector<Type> qtb(_m);
			QtMultiply(b, qtb);

			// Back-substitution on R (using only first n elements of qtb)
			x.Resize(_n);
			for (int i = 0; i < _n; i++)
				x[i] = qtb[i];

			RSolve(x, x);
		}

		Vector<Type> LeastSquaresSolve(const Vector<Type>& b)
		{
			Vector<Type> x(_n);
			LeastSquaresSolve(b, x);
			return x;
		}

		// Solves R*x = b where R is upper triangular (stored in QR)
		// b and x can be the same vector (in-place operation)
		void RSolve(const Vector<Type>& b, Vector<Type>& x)
		{
			if (b.size() != static_cast<size_t>(_n))
				throw VectorDimensionError("QRSolver::RSolve - Vector size mismatch", _n, b.size());
			if (_sing)
				throw SingularMatrixError("QRSolver::RSolve - Singular matrix");

			x.Resize(_n);
			for (int i = 0; i < _n; i++)
				x[i] = b[i];

			// Back-substitution
			for (int i = _n - 1; i >= 0; i--)
			{
				Type sum = x[i];
				for (int j = i + 1; j < _n; j++)
					sum -= _QR[i][j] * x[j];
				x[i] = sum / _d[i];
			}
		}

		// Multiplies Q^T * b and stores result in qtb
		// Uses the Householder vectors stored in QR
		void QtMultiply(const Vector<Type>& b, Vector<Type>& qtb)
		{
			if (b.size() != static_cast<size_t>(_m))
				throw VectorDimensionError("QRSolver::QtMultiply - Vector size mismatch", _m, b.size());

			qtb.Resize(_m);
			for (int i = 0; i < _m; i++)
				qtb[i] = b[i];

			// Apply Householder transformations (n-1 for square, n for overdetermined)
			int ncols = (_m == _n) ? _n - 1 : _n;
			for (int k = 0; k < ncols; k++)
			{
				if (_c[k] != 0.0)
				{
					Type sum = 0.0;
					for (int i = k; i < _m; i++)
						sum += _QR[i][k] * qtb[i];

					Type tau = sum / _c[k];

					for (int i = k; i < _m; i++)
						qtb[i] -= tau * _QR[i][k];
				}
			}
		}

		// Multiplies Q * b and stores result in qb
		// Useful for computing residuals and verifying decomposition
		void QMultiply(const Vector<Type>& b, Vector<Type>& qb)
		{
			if (b.size() != static_cast<size_t>(_m))
				throw VectorDimensionError("QRSolver::QMultiply - Vector size mismatch", _m, b.size());

			qb.Resize(_m);
			for (int i = 0; i < _m; i++)
				qb[i] = b[i];

			// Apply Householder transformations in reverse order (n-1 for square, n for overdetermined)
			int ncols = (_m == _n) ? _n - 1 : _n;
			for (int k = ncols - 1; k >= 0; k--)
			{
				if (_c[k] != 0.0)
				{
					Type sum = 0.0;
					for (int i = k; i < _m; i++)
						sum += _QR[i][k] * qb[i];

					Type tau = sum / _c[k];

					for (int i = k; i < _m; i++)
						qb[i] -= tau * _QR[i][k];
				}
			}
		}

		// Extract the upper triangular matrix R
		Matrix<Type> GetR() const
		{
			Matrix<Type> R(_n, _n);
			for (int i = 0; i < _n; i++)
			{
				for (int j = 0; j < _n; j++)
				{
					if (j >= i)
						R[i][j] = (i == j) ? _d[i] : _QR[i][j];
					else
						R[i][j] = 0.0;
				}
			}
			return R;
		}

		// Extract the orthogonal matrix Q (expensive - use sparingly)
		Matrix<Type> GetQ() const
		{
			Matrix<Type> Q(_m, _n);

			// Initialize Q as I (first n columns)
			for (int i = 0; i < _m; i++)
				for (int j = 0; j < _n; j++)
					Q[i][j] = (i == j) ? 1.0 : 0.0;

			// Apply Householder transformations in reverse order
			for (int k = _n - 1; k >= 0; k--)
			{
				if (_c[k] != 0.0)
				{
					for (int j = 0; j < _n; j++)
					{
						Type sum = 0.0;
						for (int i = k; i < _m; i++)
							sum += _QR[i][k] * Q[i][j];

						Type tau = sum / _c[k];

						for (int i = k; i < _m; i++)
							Q[i][j] -= tau * _QR[i][k];
					}
				}
			}
			return Q;
		}

		// Compute determinant (only for square matrices)
		// det(A) = det(Q) * det(R), where det(Q) = (-1)^num_reflections and det(R) = product of diagonal
		Type det() const
		{
			if (_m != _n)
				throw MatrixDimensionError("QRSolver::det - Determinant only defined for square matrices", _m, _n, _m, _m);
			if (_sing)
				return 0.0;

			Type dd = 1.0;
			for (int i = 0; i < _n; i++)
				dd *= _d[i];

			// Each Householder reflection has determinant -1
			// So det(Q) = (-1)^num_reflections
			if (_num_reflections % 2 == 1)
				dd = -dd;

			return dd;
		}

		// Check if matrix is singular
		bool IsSingular() const { return _sing; }

		// Matrix inversion using QR decomposition (only for square matrices)
		void inverse(Matrix<Type>& ainv)
		{
			if (_m != _n)
				throw MatrixDimensionError("QRSolver::inverse - Inverse only defined for square matrices", _m, _n, _m, _m);
			if (_sing)
				throw SingularMatrixError("QRSolver::inverse - Cannot invert singular matrix");

			ainv.Resize(_n, _n);
			Vector<Type> ei(_n), col(_n);

			// Solve A * ainv[,j] = e_j for each column of the inverse
			for (int j = 0; j < _n; j++)
			{
				// Set up unit vector e_j
				for (int i = 0; i < _n; i++)
					ei[i] = (i == j) ? 1.0 : 0.0;

				// Solve A * col = e_j
				Solve(ei, col);

				// Store result in j-th column of ainv
				for (int i = 0; i < _n; i++)
					ainv[i][j] = col[i];
			}
		}
	};

} // namespace MML

#endif // MML_LINEAR_ALG_QR_H
