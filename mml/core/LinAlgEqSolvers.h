///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgEqSolvers.h                                                   ///
///  Description: Direct linear system solvers (LU, Cholesky, QR, SVD)                ///
///               Matrix decompositions and equation solving algorithms               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_LINEAR_ALG_EQ_SOLVERS_H
#define MML_LINEAR_ALG_EQ_SOLVERS_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixBandDiag.h"
#include "base/BaseUtils.h"

namespace MML
{
	///////////////////////   GAUSS-JORDAN SOLVER    /////////////////////////////
	template<class Type>
	class GaussJordanSolver
	{
	public:
		// Solving with Matrix RHS (i.e. solving simultaneously for multiple RHS)
		// Throws SingularMatrixError or MatrixNumericalError on failure
		static void SolveInPlace(Matrix<Type>& a, Matrix<Type>& b)
		{
			int i, icol = 0, irow = 0, j, k, l, ll;
			Real big;
			Type dum, pivinv;

			int n = a.RowNum();
			int m = b.ColNum();
			std::vector<int> indxc(n), indxr(n), ipiv(n);
			for (j = 0; j < n; j++) ipiv[j] = 0;
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (Abs(a[j][k]) >= big) {
									big = Abs(a[j][k]);
									irow = j;
									icol = k;
								}
							}
						}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l = 0; l < n; l++) std::swap(a[irow][l], a[icol][l]);
					for (l = 0; l < m; l++) std::swap(b[irow][l], b[icol][l]);
				}
				indxr[i] = irow;
				indxc[i] = icol;

				if (a[icol][icol] == Real{ 0.0 })
					throw SingularMatrixError("GaussJordanSolver::SolveInPlace - Singular Matrix");

				pivinv = Real{ 1.0 } / a[icol][icol];
				
				// Check for non-finite pivot inverse (overflow from small pivot, Real types only)
				if constexpr (std::is_same_v<Type, Real>) {
					if (std::isnan(pivinv) || std::isinf(pivinv))
						throw MatrixNumericalError("Non-finite pivot inverse in Gauss-Jordan elimination");
				}
				
				a[icol][icol] = 1.0;
				for (l = 0; l < n; l++) a[icol][l] *= pivinv;
				for (l = 0; l < m; l++) b[icol][l] *= pivinv;
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a[ll][icol];
						a[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
						for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
					}
			}
			for (l = n - 1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						std::swap(a[k][indxr[l]], a[k][indxc[l]]);
			}
		}
		
		// Solving for a given RHS vector
		// Throws SingularMatrixError or MatrixNumericalError on failure
		static void SolveInPlace(Matrix<Type>& a, Vector<Type>& b)
		{
			auto bmat = Utils::ColumnMatrixFromVector(b);
			SolveInPlace(a, bmat);
			b = bmat.VectorFromColumn(0);
		}

		// Solving for a given RHS vector, returns result
		// Throws SingularMatrixError or MatrixNumericalError on failure
		static Vector<Type> Solve(Matrix<Type>& a, const Vector<Type>& b)
		{
			Matrix<Type> bmat = Utils::ColumnMatrixFromVector(b);
			SolveInPlace(a, bmat);
			return bmat.VectorFromColumn(0);
		}

		// solving for a given RHS vector, but with return value
		// both matrix and vector are const and are not changed
		// (in case of singular matrix 'a', exception is thrown)
		static Vector<Type> SolveConst(const Matrix<Type>& a, const Vector<Type>& b)
		{
			if (a.RowNum() != a.ColNum())
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix must be square", a.RowNum(), a.ColNum(), -1, -1);
			if (a.RowNum() != b.size())
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix and vector must have same size", a.RowNum(), a.ColNum(), b.size(), -1);
			if (b.size() == 0)
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - vector must be non-empty", a.RowNum(), a.ColNum(), b.size(), -1);
			if (a.RowNum() == 0 || a.ColNum() == 0 )
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix must be non-empty", a.RowNum(), a.ColNum(), -1, -1);

			Matrix<Real> mat(a);
			Matrix<Type> bmat = Utils::ColumnMatrixFromVector(b);
			SolveInPlace(mat, bmat);
			return bmat.VectorFromColumn(0);
		}
	};

	//////////////////////   LU DECOMPOSITION SOLVER    ///////////////////////////
	template<class Type>
	class LUSolver
	{
	private:
		int _n;
		// Note: _refOrig removed - was unused and created dangling reference risk
		// The LU decomposition is stored in _lu (copy of input matrix)

		Matrix<Type> _lu;
		std::vector<int> _indx;
		Real _d;

	public:
		// constructor - performs the LU decomposition
		// The input matrix is copied, so temporaries are safe to pass
		// Algorithm uses partial pivoting with row interchanges
		LUSolver(const Matrix<Type>& inMatRef) 
			: _n(inMatRef.RowNum()), _lu(inMatRef), _indx(_n)
		{
			int i, imax, j, k;
			Real big, temp;
			Type temp2;
			Vector<Type> vv(_n);
			
			_d = 1.0;
			// finding biggest element in each row, and saving its invers in vv[]
			for (i = 0; i < _n; i++) 
			{
				big = 0.0;
				for (j = 0; j < _n; j++)
					if ((temp = Abs(_lu[i][j])) > big) 
						big = temp;
				
				// Check for non-finite values in matrix (Real types only)
				if constexpr (std::is_same_v<Type, Real>) {
					if (std::isnan(big) || std::isinf(big))
						throw MatrixNumericalError("Non-finite values detected in LUSolver input matrix");
				}
				
				if (big == 0.0)
					throw SingularMatrixError("LUSolver::ctor - Singular Matrix");

				vv[i] = 1.0 / big;
			}

			// main loop
			for (k = 0; k < _n; k++) 
			{
				big = 0.0;
				imax = k;
				for (i = k; i < _n; i++) 
				{
					temp = Abs(vv[i] * _lu[i][k]);
					if (temp > big) {
						big = temp;
						imax = i;
					}
				}

				if (k != imax) 
				{
					for (j = 0; j < _n; j++) 
					{
						temp2 = _lu[imax][j];
						_lu[imax][j] = _lu[k][j];
						_lu[k][j] = temp2;
					}
					_d = -_d;
					vv[imax] = vv[k];
				}

				_indx[k] = imax;
				if (_lu[k][k] == Real{ 0.0 }) 
					throw SingularMatrixError("LUSolver::ctor - Singular Matrix");
				
				for (i = k + 1; i < _n; i++) 
				{
					temp2 = _lu[i][k] /= _lu[k][k];
					
					// Check for non-finite multiplier (can happen from overflow/underflow, Real types only)
					if constexpr (std::is_same_v<Type, Real>) {
						if (std::isnan(temp2) || std::isinf(temp2))
							throw MatrixNumericalError("Non-finite multiplier in LU decomposition");
					}
					
					for (j = k + 1; j < _n; j++)
						_lu[i][j] -= temp2 * _lu[k][j];
				}
			}
		}

		// solving with Matrix RHS (ie. solving simultaneously for multiple RHS)
		// while using the LU decomposition already performed in the constructor
		void Solve(Matrix<Type>& matRHS, Matrix<Type> &outSol)
		{
			int i, j;

			int n = matRHS.RowNum();
			int m = matRHS.ColNum();

			if (matRHS.RowNum() != _n )
				throw("LUSolver::solve bad sizes");

			Vector<Type> xx(n);

			for (j = 0; j < m; j++) 
			{
				for (i = 0; i < n; i++)
					xx[i] = matRHS[i][j];

				Solve(xx, xx);

				for (i = 0; i < n; i++)
					outSol[i][j] = xx[i];
			}
		}

		// solving for a given RHS vector
		// while using the LU decomposition already performed in the constructor
		bool Solve(const Vector<Type>& b, Vector<Type>& x)
		{
			// _lu, _n, and _indx are not modified by this routine
			// and can be left in place for successive calls with different right-hand sides b. This routine takes
			// into account the possibility that b will begin with many zero elements, so it is efficient for use
			// in Matrix<Real> inversion
			int i, ii = 0, ip, j;
			Type sum;

			if (b.size() != _n || x.size() != _n)
				return false;

			for (i = 0; i < _n; i++)
				x[i] = b[i];

			for (i = 0; i < _n; i++) {
				ip = _indx[i];
				sum = x[ip];
				x[ip] = x[i];
				if (ii != 0)
					for (j = ii - 1; j < i; j++) 
						sum -= _lu[i][j] * x[j];
				else if (sum != Real{ 0.0 })
					ii = i + 1;
				x[i] = sum;
			}

			for (i = _n - 1; i >= 0; i--) {
				sum = x[i];
				for (j = i + 1; j < _n; j++) 
					sum -= _lu[i][j] * x[j];
				x[i] = sum / _lu[i][i];
			}
			return true;
		}

		Vector<Type> Solve(const Vector<Type>& b)
		{
			Vector<Type> x(b.size());
			if (Solve(b, x) == true)
				return x;
			else
				throw VectorDimensionError("LUDecompositionSolver::Solve - bad sizes", b.size(), _n);
		}

		// Using the stored LU decomposition, return in ainv the matrix inverse 
		void inverse(Matrix<Type>& ainv)
		{
			int i, j;
			ainv.Resize(_n, _n);
			for (i = 0; i < _n; i++) {
				for (j = 0; j < _n; j++) ainv[i][j] = 0.;
				ainv[i][i] = 1.;
			}
			Solve(ainv, ainv);
		}

		Type det()
		{
			Type dd = _d;
			for (int i = 0; i < _n; i++)
				dd = dd * _lu[i][i];
			return dd;
		}
	};


	template<class Type>
	class LUSolverInPlace
	{
	private:
		int _n;
		Matrix<Type> &_lu;
		std::vector<int> _indx;
		Real _d;

	public:
		// constructor - performs the LU decomposition
		// Algorithm uses partial pivoting with row interchanges
		LUSolverInPlace(Matrix<Type>& inMatRef)
			: _n(inMatRef.RowNum()), _lu(inMatRef), _indx(_n)
		{
			const Real TINY = 1.0e-40;
			int i, imax, j, k;
			Real big, temp;
			Type temp2;
			Vector<Type> vv(_n);
			_d = 1.0;
			for (i = 0; i < _n; i++) {
				big = 0.0;
				for (j = 0; j < _n; j++)
					if ((temp = Abs(_lu[i][j])) > big) big = temp;
				if (big == 0.0)
					throw SingularMatrixError("LUDecompositionSolver::ctor - Singular Matrix");

				vv[i] = 1.0 / big;
			}
			for (k = 0; k < _n; k++) {
				big = 0.0;
				imax = k;
				for (i = k; i < _n; i++) {
					temp = Abs(vv[i] * _lu[i][k]);
					if (temp > big) {
						big = temp;
						imax = i;
					}
				}
				if (k != imax) {
					for (j = 0; j < _n; j++) {
						temp2 = _lu[imax][j];
						_lu[imax][j] = _lu[k][j];
						_lu[k][j] = temp2;
					}
					_d = -_d;
					vv[imax] = vv[k];
				}
				_indx[k] = imax;
				if (_lu[k][k] == Real{ 0.0 })
					_lu[k][k] = TINY;

				for (i = k + 1; i < _n; i++) {
					temp2 = _lu[i][k] /= _lu[k][k];
					for (j = k + 1; j < _n; j++)
						_lu[i][j] -= temp2 * _lu[k][j];
				}
			}
		}
		// solving with Matrix RHS (ie. solving simultaneously for multiple RHS)
		// while using the LU decomposition already performed in the constructor

		Matrix<Type> Solve(Matrix<Type>& matRHS)
		{
			int i, j;

			int n = matRHS.RowNum();
			int m = matRHS.ColNum();

			if (matRHS.RowNum() != n)
				throw("LUSolver::solve bad sizes");

			Matrix<Type> outSol(n, m);
			Vector<Type> xx(n);

			for (j = 0; j < m; j++)
			{
				for (i = 0; i < n; i++)
					xx[i] = matRHS[i][j];

				Solve(xx, xx);

				for (i = 0; i < n; i++)
					outSol[i][j] = xx[i];
			}

			return outSol;
		}

		// solving for a given RHS vector
		// while using the LU decomposition already performed in the constructor
		bool Solve(const Vector<Type>& b, Vector<Type>& x)
		{
			// _lu, _n, and _indx are not modified by this routine
			// and can be left in place for successive calls with different right-hand sides b. This routine takes
			// into account the possibility that b will begin with many zero elements, so it is efficient for use
			// in Matrix<Real> inversion
			int i, ii = 0, ip, j;
			Type sum;

			if (b.size() != _n || x.size() != _n)
				return false;

			for (i = 0; i < _n; i++)
				x[i] = b[i];

			for (i = 0; i < _n; i++) {
				ip = _indx[i];
				sum = x[ip];
				x[ip] = x[i];
				if (ii != 0)
					for (j = ii - 1; j < i; j++)
						sum -= _lu[i][j] * x[j];
				else if (sum != Real{ 0.0 })
					ii = i + 1;
				x[i] = sum;
			}

			for (i = _n - 1; i >= 0; i--) {
				sum = x[i];
				for (j = i + 1; j < _n; j++)
					sum -= _lu[i][j] * x[j];
				x[i] = sum / _lu[i][i];
			}
			return true;
		}

		Vector<Type> Solve(const Vector<Type>& b)
		{
			Vector<Type> x(b.size());
			if (Solve(b, x) == true)
				return x;
			else
				throw VectorDimensionError("LUDecompositionSolver::Solve - bad sizes", b.size(), _n);
		}
	};


	///////////////////////   BAND DIAGONAL SOLVER    /////////////////////////////
	/**
	 * @brief Solver for band diagonal linear systems using LU decomposition
	 * 
	 * This solver is optimized for band diagonal matrices, exploiting their
	 * sparse structure for O(n*m^2) complexity instead of O(n^3) for general matrices,
	 * where n is the matrix dimension and m is the bandwidth.
	 * 
	 * Based on the algorithm from Numerical Recipes (bandec/banbks).
	 * 
	 * Usage:
	 * @code
	 *   BandDiagonalMatrix A(n, m1, m2, data);  // n×n with lower bandwidth m1, upper bandwidth m2
	 *   Vector<Real> b = ...;                    // Right-hand side
	 *   BandDiagonalSolver solver(A);           // Performs LU decomposition
	 *   Vector<Real> x = solver.Solve(b);       // Solve Ax = b
	 * @endcode
	 */
	class BandDiagonalSolver
	{
	private:
		int _n;           // Matrix dimension
		int _m1;          // Lower bandwidth (number of subdiagonals)
		int _m2;          // Upper bandwidth (number of superdiagonals)
		Matrix<Real> _au; // Upper triangular factor with fill-in
		Matrix<Real> _al; // Lower triangular factor (stores multipliers)
		std::vector<int> _indx;  // Row permutation indices
		Real _d;          // +1 or -1 depending on row exchanges

	public:
		/**
		 * @brief Constructor - performs LU decomposition of band diagonal matrix
		 * @param a Band diagonal matrix to decompose
		 * @throws SingularMatrixError if matrix is singular
		 * 
		 * The decomposition stores L and U factors for subsequent solve operations.
		 * The algorithm uses partial pivoting within the band structure.
		 */
		BandDiagonalSolver(const BandDiagonalMatrix& a)
			: _n(a.GetDimension()), 
			  _m1(a.GetLowerBandwidth()), 
			  _m2(a.GetUpperBandwidth()),
			  _au(_n, _m1 + _m2 + 1),
			  _al(_n, _m1),
			  _indx(_n)
		{
			const Real TINY = 1.0e-40;
			int mm = _m1 + _m2 + 1;
			
			// Copy band diagonal matrix to working storage
			// _au stores the upper triangle with potential fill-in from pivoting
			for (int i = 0; i < _n; i++)
			{
				for (int j = 0; j < mm; j++)
				{
					int col = j - _m1 + i;
					if (col >= 0 && col < _n)
						_au[i][j] = a(i, col);
					else
						_au[i][j] = 0.0;
				}
			}
			
			_d = 1.0;
			int l = _m1;
			
			// Rearrange the storage for the initial rows
			for (int i = 0; i < _m1; i++)
			{
				for (int j = _m1 - i; j < mm; j++)
					_au[i][j - l] = _au[i][j];
				l--;
				for (int j = mm - l - 1; j < mm; j++)
					_au[i][j] = 0.0;
			}
			
			l = _m1;
			
			// LU decomposition with partial pivoting
			for (int k = 0; k < _n; k++)
			{
				Real dum = _au[k][0];
				int i = k;
				
				// Find pivot
				if (l < _n) l++;
				for (int j = k + 1; j < l; j++)
				{
					if (Abs(_au[j][0]) > Abs(dum))
					{
						dum = _au[j][0];
						i = j;
					}
				}
				
				_indx[k] = i + 1;  // Store 1-based index for compatibility
				
				if (dum == 0.0)
					_au[k][0] = TINY;  // Matrix is singular but proceed with tiny pivot
				
				// Interchange rows if necessary
				if (i != k)
				{
					_d = -_d;
					for (int j = 0; j < mm; j++)
						std::swap(_au[k][j], _au[i][j]);
				}
				
				// Compute multipliers and eliminate
				for (int i2 = k + 1; i2 < l; i2++)
				{
					dum = _au[i2][0] / _au[k][0];
					_al[k][i2 - k - 1] = dum;
					for (int j = 1; j < mm; j++)
						_au[i2][j - 1] = _au[i2][j] - dum * _au[k][j];
					_au[i2][mm - 1] = 0.0;
				}
			}
		}

		/**
		 * @brief Solve the system Ax = b using the stored LU decomposition
		 * @param b Right-hand side vector
		 * @param x Solution vector (output)
		 * @return true if solve succeeded, false if dimension mismatch
		 */
		bool Solve(const Vector<Real>& b, Vector<Real>& x) const
		{
			if (b.size() != static_cast<size_t>(_n) || x.size() != static_cast<size_t>(_n))
				return false;

			int mm = _m1 + _m2 + 1;
			int l = _m1;
			
			// Copy b to x
			for (int i = 0; i < _n; i++)
				x[i] = b[i];
			
			// Forward substitution
			for (int k = 0; k < _n; k++)
			{
				int j = _indx[k] - 1;  // Convert back to 0-based
				if (j != k)
					std::swap(x[k], x[j]);
				
				if (l < _n) l++;
				for (int j2 = k + 1; j2 < l; j2++)
					x[j2] -= _al[k][j2 - k - 1] * x[k];
			}
			
			l = 1;
			
			// Back substitution
			for (int i = _n - 1; i >= 0; i--)
			{
				Real dum = x[i];
				for (int k = 1; k < l; k++)
					dum -= _au[i][k] * x[k + i];
				x[i] = dum / _au[i][0];
				if (l < mm) l++;
			}
			
			return true;
		}

		/**
		 * @brief Solve the system Ax = b, returning the solution
		 * @param b Right-hand side vector
		 * @return Solution vector x
		 * @throws VectorDimensionError if dimension mismatch
		 */
		Vector<Real> Solve(const Vector<Real>& b) const
		{
			Vector<Real> x(_n);
			if (Solve(b, x))
				return x;
			else
				throw VectorDimensionError("BandDiagonalSolver::Solve - dimension mismatch", b.size(), _n);
		}

		/**
		 * @brief Solve multiple right-hand sides simultaneously
		 * @param B Matrix where each column is a right-hand side
		 * @return Solution matrix where each column is the corresponding solution
		 */
		Matrix<Real> Solve(const Matrix<Real>& B) const
		{
			if (B.RowNum() != _n)
				throw MatrixDimensionError("BandDiagonalSolver::Solve - row count mismatch", 
				                          B.RowNum(), B.ColNum(), _n, -1);
			
			Matrix<Real> X(_n, B.ColNum());
			Vector<Real> b(_n), x(_n);
			
			for (int j = 0; j < B.ColNum(); j++)
			{
				for (int i = 0; i < _n; i++)
					b[i] = B[i][j];
				
				Solve(b, x);
				
				for (int i = 0; i < _n; i++)
					X[i][j] = x[i];
			}
			
			return X;
		}

		/**
		 * @brief Get the determinant of the matrix
		 * @return Determinant value
		 */
		Real Det() const
		{
			Real dd = _d;
			for (int i = 0; i < _n; i++)
				dd *= _au[i][0];
			return dd;
		}

		/**
		 * @brief Get matrix dimension
		 */
		int GetDimension() const { return _n; }

		/**
		 * @brief Get lower bandwidth
		 */
		int GetLowerBandwidth() const { return _m1; }

		/**
		 * @brief Get upper bandwidth
		 */
		int GetUpperBandwidth() const { return _m2; }
	};


	///////////////////////   CHOLESKY DECOMPOSITION    /////////////////////////////
	template<class Type>
	class CholeskySolver
	{
	public:
		int n;
		Matrix<Type> el;

	public:
    // Given a positive-definite symmetric Matrix<Type> a[1..n][1..n], this routine constructs its Cholesky
    // decomposition A = L * L^T. On input, only the upper triangle of a need be given; it is not
    // modified. The Cholesky factor L is returned in the lower triangle of el, including diagonal elements.
		CholeskySolver(const Matrix<Type>& a) : n(a.RowNum()), el(a)
		{
			// Perform Cholesky decomposition
			for (int i = 0; i < n; i++)
			{
				for (int j = i; j < n; j++)
				{
					Type sum = el[i][j];
					
					// Subtract contributions from previous columns
					for (int k = 0; k < i; k++)
						sum -= el[i][k] * el[j][k];
					
					if (i == j)
					{
						// Diagonal element
						if (sum <= 0.0)
							throw SingularMatrixError("CholeskySolver: Matrix is not positive definite");
						
						el[i][i] = std::sqrt(sum);
					}
					else
					{
						// Off-diagonal element
						el[j][i] = sum / el[i][i];
					}
				}
			}
			
			// Zero out upper triangle (keeping only lower triangular L)
			for (int i = 0; i < n; i++)
				for (int j = i + 1; j < n; j++)
					el[i][j] = 0.0;
		}
		
		void Solve(const Vector<Type>& b, Vector<Type>& x)
		{
			// Solves the set of n linear equations A * x = b, where A is a positive-definite symmetric Matrix<Type>.
			// b[1..n] is input as the right-hand side Vector<Type>. 
			// The solution Vector<Type> is returned in x[1..n].
			// Uses forward and backward substitution: A*x = b => L*L^T*x = b
			// First solve L*y = b (forward substitution), then solve L^T*x = y (backward substitution)
			
			if (b.size() != n)
				throw VectorDimensionError("CholeskySolver::Solve - Vector size mismatch", n, b.size());
			
			x.Resize(n);
			
			// Forward substitution: solve L*y = b for y
			for (int i = 0; i < n; i++)
			{
				Type sum = b[i];
				for (int k = 0; k < i; k++)
					sum -= el[i][k] * x[k];
				x[i] = sum / el[i][i];
			}
			
			// Backward substitution: solve L^T*x = y for x
			for (int i = n - 1; i >= 0; i--)
			{
				Type sum = x[i];
				for (int k = i + 1; k < n; k++)
					sum -= el[k][i] * x[k];
				x[i] = sum / el[i][i];
			}
		}
		
		Vector<Type> Solve(const Vector<Type>& b)
		{
			Vector<Type> x(b.size());
			Solve(b, x);
			return x;
		}
		
		void inverse(Matrix<Type>& ainv)
		{
			// Computes the inverse of the original matrix A using the Cholesky decomposition
			// A^-1 = (L*L^T)^-1 = (L^T)^-1 * L^-1
			
			ainv.Resize(n, n);
			Vector<Type> ei(n), col(n);
			
			// Solve A * ainv[,j] = e_j for each column of the inverse
			for (int j = 0; j < n; j++)
			{
				// Set up unit vector e_j
				for (int i = 0; i < n; i++)
					ei[i] = (i == j) ? 1.0 : 0.0;
				
				// Solve A * col = e_j
				Solve(ei, col);
				
				// Store result in j-th column of ainv
				for (int i = 0; i < n; i++)
					ainv[i][j] = col[i];
			}
		}
		
		Type logdet()
		{
			// Computes log(det(A)) where A = L*L^T
			// det(A) = det(L*L^T) = det(L)^2 = (product of diagonal elements of L)^2
			// log(det(A)) = 2 * log(det(L)) = 2 * sum(log(L[i][i]))
			Type sum = 0.0;
			for (int i = 0; i < n; i++)
				sum += std::log(el[i][i]);
			return 2.0 * sum;
		}
	};

	///////////////////////   QR DECOMPOSITION    /////////////////////////////
	// QR decomposition using Householder reflections
	// Decomposes A = Q*R where Q is orthogonal and R is upper triangular
	// Works for both square and overdetermined (m >= n) systems
	template<class Type>
	class QRSolver
	{
	public:
		int m;              // number of rows
		int n;              // number of columns
		Matrix<Type> QR;    // Combined QR storage: R in upper triangle, Householder vectors in lower
		Vector<Type> c;     // Diagonal of R (stored separately)
		Vector<Type> d;     // Householder scaling factors
		bool sing;          // Singularity flag
		int num_reflections; // Number of Householder reflections performed (for determinant sign)

	public:
		// Constructs the QR decomposition of a[m][n] (m >= n)
		// Uses Householder reflections to compute Q and R such that A = Q*R
		// On output, the upper triangular matrix R is stored in the upper triangle and diagonal of QR
		// The Householder vectors are stored in the lower triangle of QR
		QRSolver(const Matrix<Type>& a) : m(a.RowNum()), n(a.ColNum()), QR(a), c(n), d(n), sing(false), num_reflections(0)
		{
			if (m < n)
				throw MatrixDimensionError("QRSolver: Matrix must have m >= n (rows >= columns)", m, n, m, n);

			int i, j, k;
			Type scale, sigma, sum, tau;
			Vector<Type> vec(m);

			// Perform Householder reduction
			// For square matrices: process n-1 columns (last has no elements below diagonal)
			// For overdetermined: process all n columns (last column still has elements below diagonal)
			int ncols = (m == n) ? n - 1 : n;
			for (k = 0; k < ncols; k++)
			{
				// Compute the norm of the k-th column below the diagonal
				scale = 0.0;
				for (i = k; i < m; i++)
					scale = std::max(scale, Abs(QR[i][k]));

				if (scale == 0.0)
				{
					// Singular case
					sing = true;
					c[k] = d[k] = 0.0;
				}
				else
				{
					// Form the Householder vector
					for (i = k; i < m; i++)
						QR[i][k] /= scale;

					sum = 0.0;
					for (i = k; i < m; i++)
						sum += QR[i][k] * QR[i][k];

					// Choose sign of sigma to avoid cancellation errors
					sigma = (QR[k][k] >= 0.0 ? std::sqrt(sum) : -std::sqrt(sum));
					QR[k][k] += sigma;
					c[k] = sigma * QR[k][k];
					d[k] = -scale * sigma;

					// Count this reflection for determinant
					num_reflections++;

					// Apply the transformation to remaining columns
					for (j = k + 1; j < n; j++)
					{
						sum = 0.0;
						for (i = k; i < m; i++)
							sum += QR[i][k] * QR[i][j];

						tau = sum / c[k];

						for (i = k; i < m; i++)
							QR[i][j] -= tau * QR[i][k];
					}
				}
			}

			// For square matrices, handle the last diagonal element separately
			if (m == n)
			{
				d[n - 1] = QR[n - 1][n - 1];
				c[n - 1] = 0.0;  // No Householder reflection for last column in square matrix
				if (d[n - 1] == 0.0)
					sing = true;
			}
		}

		// Solves the set of n linear equations A*x = b for a square system (m == n)
		// Uses back-substitution on R after applying Q^T to b
		void Solve(const Vector<Type>& b, Vector<Type>& x)
		{
			if (m != n)
				throw MatrixDimensionError("QRSolver::Solve - Use LeastSquaresSolve for overdetermined systems (m > n)", m, n, m, m);
			if (b.size() != m)
				throw VectorDimensionError("QRSolver::Solve - Vector size mismatch", m, b.size());
			if (sing)
				throw SingularMatrixError("QRSolver::Solve - Singular matrix");

			// Apply Q^T to b
			QtMultiply(b, x);

			// Back-substitution on R
			RSolve(x, x);
		}

		Vector<Type> Solve(const Vector<Type>& b)
		{
			Vector<Type> x(n);
			Solve(b, x);
			return x;
		}

		// Least-squares solution for overdetermined systems (m > n)
		// Minimizes ||Ax - b||₂ by solving R*x = Q^T*b
		void LeastSquaresSolve(const Vector<Type>& b, Vector<Type>& x)
		{
			if (b.size() != m)
				throw VectorDimensionError("QRSolver::LeastSquaresSolve - Vector size mismatch", m, b.size());
			if (sing)
				throw SingularMatrixError("QRSolver::LeastSquaresSolve - Singular matrix");

			// Apply Q^T to b
			Vector<Type> qtb(m);
			QtMultiply(b, qtb);

			// Back-substitution on R (using only first n elements of qtb)
			x.Resize(n);
			for (int i = 0; i < n; i++)
				x[i] = qtb[i];
			
			RSolve(x, x);
		}

		Vector<Type> LeastSquaresSolve(const Vector<Type>& b)
		{
			Vector<Type> x(n);
			LeastSquaresSolve(b, x);
			return x;
		}

		// Solves R*x = b where R is upper triangular (stored in QR)
		// b and x can be the same vector (in-place operation)
		void RSolve(const Vector<Type>& b, Vector<Type>& x)
		{
			if (b.size() != n)
				throw VectorDimensionError("QRSolver::RSolve - Vector size mismatch", n, b.size());
			if (sing)
				throw SingularMatrixError("QRSolver::RSolve - Singular matrix");

			x.Resize(n);
			for (int i = 0; i < n; i++)
				x[i] = b[i];

			// Back-substitution
			for (int i = n - 1; i >= 0; i--)
			{
				Type sum = x[i];
				for (int j = i + 1; j < n; j++)
					sum -= QR[i][j] * x[j];
				x[i] = sum / d[i];
			}
		}

		// Multiplies Q^T * b and stores result in qtb
		// Uses the Householder vectors stored in QR
		void QtMultiply(const Vector<Type>& b, Vector<Type>& qtb)
		{
			if (b.size() != m)
				throw VectorDimensionError("QRSolver::QtMultiply - Vector size mismatch", m, b.size());

			qtb.Resize(m);
			for (int i = 0; i < m; i++)
				qtb[i] = b[i];

			// Apply Householder transformations (n-1 for square, n for overdetermined)
			int ncols = (m == n) ? n - 1 : n;
			for (int k = 0; k < ncols; k++)
			{
				if (c[k] != 0.0)
				{
					Type sum = 0.0;
					for (int i = k; i < m; i++)
						sum += QR[i][k] * qtb[i];

					Type tau = sum / c[k];

					for (int i = k; i < m; i++)
						qtb[i] -= tau * QR[i][k];
				}
			}
		}

		// Multiplies Q * b and stores result in qb
		// Useful for computing residuals and verifying decomposition
		void QMultiply(const Vector<Type>& b, Vector<Type>& qb)
		{
			if (b.size() != m)
				throw VectorDimensionError("QRSolver::QMultiply - Vector size mismatch", m, b.size());

			qb.Resize(m);
			for (int i = 0; i < m; i++)
				qb[i] = b[i];

			// Apply Householder transformations in reverse order (n-1 for square, n for overdetermined)
			int ncols = (m == n) ? n - 1 : n;
			for (int k = ncols - 1; k >= 0; k--)
			{
				if (c[k] != 0.0)
				{
					Type sum = 0.0;
					for (int i = k; i < m; i++)
						sum += QR[i][k] * qb[i];

					Type tau = sum / c[k];

					for (int i = k; i < m; i++)
						qb[i] -= tau * QR[i][k];
				}
			}
		}

		// Extract the upper triangular matrix R
		Matrix<Type> GetR() const
		{
			Matrix<Type> R(n, n);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					if (j >= i)
						R[i][j] = (i == j) ? d[i] : QR[i][j];
					else
						R[i][j] = 0.0;
				}
			}
			return R;
		}

		// Extract the orthogonal matrix Q (expensive - use sparingly)
		Matrix<Type> GetQ() const
		{
			Matrix<Type> Q(m, n);
			
			// Initialize Q as I (first n columns)
			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
					Q[i][j] = (i == j) ? 1.0 : 0.0;

			// Apply Householder transformations in reverse order
			// For square matrices: only n-1 reflections were performed
			// For overdetermined: all n reflections were performed
			int ncols = (m == n) ? n - 1 : n;
			for (int k = ncols - 1; k >= 0; k--)
			{
				if (c[k] != 0.0)
				{
					for (int j = 0; j < n; j++)
					{
						Type sum = 0.0;
						for (int i = k; i < m; i++)
							sum += QR[i][k] * Q[i][j];

						Type tau = sum / c[k];

						for (int i = k; i < m; i++)
							Q[i][j] -= tau * QR[i][k];
					}
				}
			}
			return Q;
		}

		// Compute determinant (only for square matrices)
		// det(A) = det(Q) * det(R), where det(Q) = (-1)^num_reflections and det(R) = product of diagonal
		Type det() const
		{
			if (m != n)
				throw MatrixDimensionError("QRSolver::det - Determinant only defined for square matrices", m, n, m, m);
			if (sing)
				return 0.0;

			Type dd = 1.0;
			for (int i = 0; i < n; i++)
				dd *= d[i];
			
			// Each Householder reflection has determinant -1
			// So det(Q) = (-1)^num_reflections
			if (num_reflections % 2 == 1)
				dd = -dd;
				
			return dd;
		}

		// Matrix inversion using QR decomposition (only for square matrices)
		void inverse(Matrix<Type>& ainv)
		{
			if (m != n)
				throw MatrixDimensionError("QRSolver::inverse - Inverse only defined for square matrices", m, n, m, m);
			if (sing)
				throw SingularMatrixError("QRSolver::inverse - Cannot invert singular matrix");

			ainv.Resize(n, n);
			Vector<Type> ei(n), col(n);

			// Solve A * ainv[,j] = e_j for each column of the inverse
			for (int j = 0; j < n; j++)
			{
				// Set up unit vector e_j
				for (int i = 0; i < n; i++)
					ei[i] = (i == j) ? 1.0 : 0.0;

				// Solve A * col = e_j
				Solve(ei, col);

				// Store result in j-th column of ainv
				for (int i = 0; i < n; i++)
					ainv[i][j] = col[i];
			}
		}
  };

	/////////////////////////////////   SVD DECOMPOSITION      /////////////////////////////  
	class SVDecompositionSolver
	{
	private:
		int m, n;
		Matrix<Real> u, v;
		Vector<Real> w;
		Real eps, tsh;

		// Helper function: numerically stable sqrt(a^2 + b^2)
		Real pythag(const Real a, const Real b) const {
			Real absa = std::abs(a), absb = std::abs(b);
			return (absa > absb ? absa * std::sqrt(1.0 + (absb / absa) * (absb / absa)) :
				(absb == 0.0 ? 0.0 : absb * std::sqrt(1.0 + (absa / absb) * (absa / absb))));
		}

		// Core SVD decomposition using Householder reduction and QR iteration
		void decompose();

		// Reorder singular values in descending order
		void reorder();

	public:
		Vector<Real> getW() const { return w; }
		Matrix<Real> getU() const { return u; }
		Matrix<Real> getV() const { return v; }

	public:
		SVDecompositionSolver(const Matrix<Real>& a) : m(a.RowNum()), n(a.ColNum()), u(a), v(n, n), w(n)
		{
			// Given a Matrix<Real> a[m][n], this routine computes its singular value decomposition, A = U·W·V^T
			// The Matrix U replaces a on output (stored in u)
			// The diagonal matrix of singular values W is output as a vector w[n]
			// The matrix V (not the transpose V^T) is output as v[n][n]
			
			eps = std::numeric_limits<Real>::epsilon();
			decompose();
			reorder();
			tsh = 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps;
		}

		Real inv_condition() {
			return (w[0] <= 0. || w[n - 1] <= 0.) ? 0. : w[n - 1] / w[0];
		}

		void Solve(const Vector<Real>& b, Vector<Real>& x, Real thresh = -1.)
		{
			// Solve A·x = b for a vector x using the pseudoinverse of A as obtained by SVD. If positive,
			// thresh is the threshold value below which singular values are considered as zero. If thresh is
			// negative, a default based on expected roundoff error is used.
			
			Real tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			
			Vector<Real> tmp(n);
			// Calculate U^T · b
			for (int j = 0; j < n; j++)
			{
				Real s = 0.0;
				if (w[j] > tsh)  // Only include non-zero singular values
				{
					for (int i = 0; i < m; i++)
						s += u[i][j] * b[i];
					s /= w[j];  // Multiply by inverse singular value
				}
				tmp[j] = s;
			}
			
			// Calculate x = V · tmp
			for (int j = 0; j < n; j++)
			{
				Real s = 0.0;
				for (int jj = 0; jj < n; jj++)
					s += v[j][jj] * tmp[jj];
				x[j] = s;
			}
		}

		Vector<Real> Solve(const Vector<Real>& b, Real thresh = -1.)
		{
			Vector<Real> x(n);
			Solve(b, x, thresh);
			return x;
		}

		// Solves m sets of n equations A·X = B using the pseudoinverse of A. The right-hand sides are
		// input as b[m][p], while x[n][p] returns the solutions. thresh as above.
		void Solve(const Matrix<Real>& b, Matrix<Real>& x, Real thresh = -1.)
		{
			int p = b.ColNum();
			if (b.RowNum() != m || x.RowNum() != n || x.ColNum() != p)
				throw MatrixDimensionError("SVD::Solve - bad dimensions", m, n, b.RowNum(), x.RowNum());
			
			Vector<Real> bcol(m), xcol(n);
			for (int j = 0; j < p; j++)
			{
				// Extract column j from b
				for (int i = 0; i < m; i++)
					bcol[i] = b[i][j];
				
				// Solve for column j
				Solve(bcol, xcol, thresh);
				
				// Store result in column j of x
				for (int i = 0; i < n; i++)
					x[i][j] = xcol[i];
			}
		}

		// Return the rank of A, after zeroing any singular values smaller than thresh. If thresh is
		// negative, a default value based on estimated roundoff is used.        
		int Rank(Real thresh = -1.) {
			Real tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			int rank = 0;
			for (int j = 0; j < n; j++)
				if (w[j] > tsh) rank++;
			return rank;
		}

		// Return the nullity of A, after zeroing any singular values smaller than thresh. Default value as above.
		int Nullity(Real thresh = -1.) {
			Real tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			int nullity = 0;
			for (int j = 0; j < n; j++)
				if (w[j] <= tsh) nullity++;
			return nullity;
		}

		// Gives an orthonormal basis for the range of A as the columns of a returned matrix. thresh as above.
		Matrix<Real> Range(Real thresh = -1.) {
			Real tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			int rank = Rank(tsh);
			
			Matrix<Real> range(m, rank);
			int col = 0;
			for (int j = 0; j < n; j++)
			{
				if (w[j] > tsh)
				{
					for (int i = 0; i < m; i++)
						range[i][col] = u[i][j];
					col++;
				}
			}
			return range;
		}

		// Gives an orthonormal basis for the nullspace of A as the columns of a returned matrix. thresh as above
		Matrix<Real> Nullspace(Real thresh = -1.) {
			Real tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			int nullity = Nullity(tsh);
			
			Matrix<Real> nullspace(n, nullity);
			int col = 0;
			for (int j = 0; j < n; j++)
			{
				if (w[j] <= tsh)
				{
					for (int i = 0; i < n; i++)
						nullspace[i][col] = v[i][j];
					col++;
				}
			}
			return nullspace;
		}
	};

	///////////////////// SVDecompositionSolver Implementation /////////////////////

	inline void SVDecompositionSolver::decompose()
	{
		bool flag;
		int i, its, j, jj, k, l, nm;
		Real anorm, c, f, g, h, s, scale, x, y, z;
		Vector<Real> rv1(n);

		g = scale = anorm = 0.0;

		// Householder reduction to bidiagonal form
		for (i = 0; i < n; i++)
		{
			l = i + 2;
			rv1[i] = scale * g;
			g = s = scale = 0.0;

			if (i < m)
			{
				for (k = i; k < m; k++) scale += std::abs(u[k][i]);
				if (scale != 0.0)
				{
					for (k = i; k < m; k++)
					{
						u[k][i] /= scale;
						s += u[k][i] * u[k][i];
					}
					f = u[i][i];
					g = -std::copysign(std::sqrt(s), f);
					h = f * g - s;
					u[i][i] = f - g;
					for (j = l - 1; j < n; j++)
					{
						for (s = 0.0, k = i; k < m; k++) s += u[k][i] * u[k][j];
						f = s / h;
						for (k = i; k < m; k++) u[k][j] += f * u[k][i];
					}
					for (k = i; k < m; k++) u[k][i] *= scale;
				}
			}

			w[i] = scale * g;
			g = s = scale = 0.0;

			if (i + 1 <= m && i + 1 != n)
			{
				for (k = l - 1; k < n; k++) scale += std::abs(u[i][k]);
				if (scale != 0.0)
				{
					for (k = l - 1; k < n; k++)
					{
						u[i][k] /= scale;
						s += u[i][k] * u[i][k];
					}
					f = u[i][l - 1];
					g = -std::copysign(std::sqrt(s), f);
					h = f * g - s;
					u[i][l - 1] = f - g;
					for (k = l - 1; k < n; k++) rv1[k] = u[i][k] / h;
					for (j = l - 1; j < m; j++)
					{
						for (s = 0.0, k = l - 1; k < n; k++) s += u[j][k] * u[i][k];
						for (k = l - 1; k < n; k++) u[j][k] += s * rv1[k];
					}
					for (k = l - 1; k < n; k++) u[i][k] *= scale;
				}
			}
			anorm = std::max(anorm, (std::abs(w[i]) + std::abs(rv1[i])));
		}

		// Accumulation of right-hand transformations
		for (i = n - 1; i >= 0; i--)
		{
			if (i < n - 1)
			{
				if (g != 0.0)
				{
					for (j = l; j < n; j++)
						v[j][i] = (u[i][j] / u[i][l]) / g;
					for (j = l; j < n; j++)
					{
						for (s = 0.0, k = l; k < n; k++) s += u[i][k] * v[k][j];
						for (k = l; k < n; k++) v[k][j] += s * v[k][i];
					}
				}
				for (j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
			}
			v[i][i] = 1.0;
			g = rv1[i];
			l = i;
		}

		// Accumulation of left-hand transformations
		for (i = std::min(m, n) - 1; i >= 0; i--)
		{
			l = i + 1;
			g = w[i];
			for (j = l; j < n; j++) u[i][j] = 0.0;
			if (g != 0.0)
			{
				g = 1.0 / g;
				for (j = l; j < n; j++)
				{
					for (s = 0.0, k = l; k < m; k++) s += u[k][i] * u[k][j];
					f = (s / u[i][i]) * g;
					for (k = i; k < m; k++) u[k][j] += f * u[k][i];
				}
				for (j = i; j < m; j++) u[j][i] *= g;
			}
			else
				for (j = i; j < m; j++) u[j][i] = 0.0;
			++u[i][i];
		}

		// Diagonalization of the bidiagonal form via QR iteration
		for (k = n - 1; k >= 0; k--)
		{
			for (its = 0; its < 30; its++)
			{
				flag = true;
				for (l = k; l >= 0; l--)
				{
					nm = l - 1;
					if (l == 0 || std::abs(rv1[l]) <= eps * anorm)
					{
						flag = false;
						break;
					}
					if (std::abs(w[nm]) <= eps * anorm) break;
				}

				if (flag)
				{
					c = 0.0;
					s = 1.0;
					for (i = l; i < k + 1; i++)
					{
						f = s * rv1[i];
						rv1[i] = c * rv1[i];
						if (std::abs(f) <= eps * anorm) break;
						g = w[i];
						h = pythag(f, g);
						w[i] = h;
						h = 1.0 / h;
						c = g * h;
						s = -f * h;
						for (j = 0; j < m; j++)
						{
							y = u[j][nm];
							z = u[j][i];
							u[j][nm] = y * c + z * s;
							u[j][i] = z * c - y * s;
						}
					}
				}

				z = w[k];
				if (l == k)
				{
					if (z < 0.0)
					{
						w[k] = -z;
						for (j = 0; j < n; j++) v[j][k] = -v[j][k];
					}
					break;
				}

				if (its == 29)
						throw ConvergenceError("SVD: no convergence in 30 iterations", 30);
				x = w[l];
				nm = k - 1;
				y = w[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = pythag(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + std::copysign(g, f))) - h)) / x;
				c = s = 1.0;

				for (j = l; j <= nm; j++)
				{
					i = j + 1;
					g = rv1[i];
					y = w[i];
					h = s * g;
					g = c * g;
					z = pythag(f, h);
					rv1[j] = z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y *= c;
					for (jj = 0; jj < n; jj++)
					{
						x = v[jj][j];
						z = v[jj][i];
						v[jj][j] = x * c + z * s;
						v[jj][i] = z * c - x * s;
					}
					z = pythag(f, h);
					w[j] = z;
					if (z)
					{
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = c * g + s * y;
					x = c * y - s * g;
					for (jj = 0; jj < m; jj++)
					{
						y = u[jj][j];
						z = u[jj][i];
						u[jj][j] = y * c + z * s;
						u[jj][i] = z * c - y * s;
					}
				}
				rv1[l] = 0.0;
				rv1[k] = f;
				w[k] = x;
			}
		}
	}

	inline void SVDecompositionSolver::reorder()
	{
		int i, j, k, s, inc = 1;
		Real sw;
		Vector<Real> su(m), sv(n);

		// Shell sort to order singular values in descending order
		do { inc *= 3; inc++; } while (inc <= n);
		do {
			inc /= 3;
			for (i = inc; i < n; i++)
			{
				sw = w[i];
				for (k = 0; k < m; k++) su[k] = u[k][i];
				for (k = 0; k < n; k++) sv[k] = v[k][i];
				j = i;
				while (w[j - inc] < sw)
				{
					w[j] = w[j - inc];
					for (k = 0; k < m; k++) u[k][j] = u[k][j - inc];
					for (k = 0; k < n; k++) v[k][j] = v[k][j - inc];
					j -= inc;
					if (j < inc) break;
				}
				w[j] = sw;
				for (k = 0; k < m; k++) u[k][j] = su[k];
				for (k = 0; k < n; k++) v[k][j] = sv[k];
			}
		} while (inc > 1);

		// Flip signs for consistent sign convention
		for (k = 0; k < n; k++)
		{
			s = 0;
			for (i = 0; i < m; i++) if (u[i][k] < 0.) s++;
			for (j = 0; j < n; j++) if (v[j][k] < 0.) s++;
			if (s > (m + n) / 2)
			{
				for (i = 0; i < m; i++) u[i][k] = -u[i][k];
				for (j = 0; j < n; j++) v[j][k] = -v[j][k];
			}
		}
	}
}

#endif