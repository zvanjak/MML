///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgDirect.h                                                      ///
///  Description: Direct linear system solvers (Gauss-Jordan, LU, Cholesky)           ///
///               LU and band diagonal decompositions for equation solving            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_LINEAR_ALG_DIRECT_H
#define MML_LINEAR_ALG_DIRECT_H

#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixBandDiag.h"
#include "base/BaseUtils.h"
#include "core/AlgorithmTypes.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// LinearSolverConfig - Configuration for direct linear solver detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	struct LinearSolverConfig : public EvaluationConfigBase {
		// Inherits: estimate_error, check_finite, exception_policy
		// When estimate_error is true, computes residual norm ||Ax - b||
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// LinearSolverResult - Result type for direct linear solver detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	template<typename Type>
	struct LinearSolverResult : public EvaluationResultBase {
		/// The computed solution vector x
		Vector<Type> solution{};

		/// Residual norm ||Ax - b|| (populated when config.estimate_error is true)
		Real residual_norm = 0.0;
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// LinearSolverDetail - Internal helpers for Detailed API execution
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace LinearSolverDetail
	{
		/// Execute a linear solver Detailed operation with timing and exception handling.
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteLinearSolverDetailed(const char* algorithm_name,
		                                      const LinearSolverConfig& config,
		                                      ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;

				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);

				compute(result);

				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const SingularMatrixError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::SingularMatrix, ex.what(), algorithm_name);
			}
			catch (const MatrixNumericalError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::NumericalInstability, ex.what(), algorithm_name);
			}
			catch (const MatrixDimensionError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const VectorDimensionError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const std::exception& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
			}
		}
	} // namespace LinearSolverDetail

	/// @brief Gauss-Jordan elimination solver for linear systems Ax=b
	/// @tparam Type Numeric type (Real, Complex, etc.)
	/// @note Uses full pivoting for numerical stability, O(n³) complexity
	/// @warning Modifies input matrices in-place for SolveInPlace methods
	template<class Type>
	class GaussJordanSolver
	{
	public:
		/// @brief Solve Ax=B with multiple right-hand sides (in-place)
		/// @param a Coefficient matrix (modified in-place)
		/// @param b Right-hand side matrix (solution on return)
		/// @throws MatrixDimensionError if dimensions incompatible
		/// @throws SingularMatrixError if matrix is singular
		/// @throws MatrixNumericalError if pivot inverse is non-finite
		static void SolveInPlace(Matrix<Type>& a, Matrix<Type>& b)
		{
			int i, icol, irow, j, k, l, ll;
			Real big;
			Type dum, pivinv;

			int n = a.rows();
			int m = b.cols();
			
			// Dimension validation
			if (a.rows() != a.cols())
				throw MatrixDimensionError("GaussJordanSolver::SolveInPlace - matrix A must be square", a.rows(), a.cols(), -1, -1);
			if (a.rows() == 0 || a.cols() == 0)
				throw MatrixDimensionError("GaussJordanSolver::SolveInPlace - matrix A must be non-empty", a.rows(), a.cols(), -1, -1);
			if (b.rows() == 0 || b.cols() == 0)
				throw MatrixDimensionError("GaussJordanSolver::SolveInPlace - matrix B must be non-empty", a.rows(), a.cols(), b.rows(), b.cols());
			if (a.rows() != b.rows())
				throw MatrixDimensionError("GaussJordanSolver::SolveInPlace - A rows must match B rows", a.rows(), a.cols(), b.rows(), b.cols());
			std::vector<int> indxc(n), indxr(n), ipiv(n);

			// Compute infinity norm for norm-scaled singularity threshold
			Real norm_a = 0.0;
			for (int ii = 0; ii < n; ii++)
				for (int jj = 0; jj < n; jj++)
					if (Abs(a[ii][jj]) > norm_a)
						norm_a = Abs(a[ii][jj]);
			Real singularity_threshold = std::numeric_limits<Real>::epsilon() * norm_a * n;

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

				if (Abs(a[icol][icol]) < singularity_threshold)
					throw SingularMatrixError("GaussJordanSolver::SolveInPlace - Singular Matrix", Abs(a[icol][icol]));

				pivinv = Real{ 1.0 } / a[icol][icol];
				
				// Check for non-finite pivot inverse (overflow from small pivot)
				if constexpr (std::is_same_v<Type, Real>) {
					if (std::isnan(pivinv) || std::isinf(pivinv))
						throw MatrixNumericalError("Non-finite pivot inverse in Gauss-Jordan elimination");
				} else {
					Real mag = std::abs(pivinv);
					if (std::isnan(mag) || std::isinf(mag))
						throw MatrixNumericalError("Non-finite pivot inverse in Gauss-Jordan elimination (Complex)");
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
		
		/// @brief Solve Ax=b for single right-hand side vector (in-place)
		/// @param a Coefficient matrix (modified in-place)
		/// @param b Right-hand side vector (solution on return)
		/// @throws MatrixDimensionError if dimensions incompatible
		/// @throws SingularMatrixError or MatrixNumericalError on failure
		static void SolveInPlace(Matrix<Type>& a, Vector<Type>& b)
		{
			// Dimension validation
			if (a.rows() != a.cols())
				throw MatrixDimensionError("GaussJordanSolver::SolveInPlace - matrix must be square", a.rows(), a.cols(), -1, -1);
			if (a.rows() != b.size())
				throw MatrixDimensionError("GaussJordanSolver::SolveInPlace - matrix rows must match vector size", a.rows(), a.cols(), b.size(), -1);
			if (b.size() == 0)
				throw MatrixDimensionError("GaussJordanSolver::SolveInPlace - vector must be non-empty", a.rows(), a.cols(), b.size(), -1);
			
			auto bmat = Utils::ColumnMatrixFromVector(b);
			SolveInPlace(a, bmat);
			b = bmat.VectorFromColumn(0);
		}

		/// @brief Solve Ax=b and return solution vector
		/// @param a Coefficient matrix (modified in-place)
		/// @param b Right-hand side vector (not modified)
		/// @return Solution vector x
		/// @throws MatrixDimensionError if dimensions incompatible
		/// @throws SingularMatrixError or MatrixNumericalError on failure
		static Vector<Type> Solve(Matrix<Type>& a, const Vector<Type>& b)
		{
			// Dimension validation
			if (a.rows() != a.cols())
				throw MatrixDimensionError("GaussJordanSolver::Solve - matrix must be square", a.rows(), a.cols(), -1, -1);
			if (a.rows() != b.size())
				throw MatrixDimensionError("GaussJordanSolver::Solve - matrix rows must match vector size", a.rows(), a.cols(), b.size(), -1);
			if (b.size() == 0)
				throw MatrixDimensionError("GaussJordanSolver::Solve - vector must be non-empty", a.rows(), a.cols(), b.size(), -1);
			
			Matrix<Type> bmat = Utils::ColumnMatrixFromVector(b);
			SolveInPlace(a, bmat);
			return bmat.VectorFromColumn(0);
		}

		/// @brief Solve Ax=b without modifying inputs (const-correct version)
		/// @param a Coefficient matrix (not modified - internal copy made)
		/// @param b Right-hand side vector (not modified)
		/// @return Solution vector x
		/// @throws MatrixDimensionError if dimensions incompatible
		/// @throws SingularMatrixError if matrix is singular
		static Vector<Type> SolveConst(const Matrix<Type>& a, const Vector<Type>& b)
		{
			if (a.rows() != a.cols())
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix must be square", a.rows(), a.cols(), -1, -1);
			if (a.rows() != b.size())
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix and vector must have same size", a.rows(), a.cols(), b.size(), -1);
			if (b.size() == 0)
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - vector must be non-empty", a.rows(), a.cols(), b.size(), -1);
			if (a.rows() == 0 || a.cols() == 0 )
				throw MatrixDimensionError("GaussJordanSolver::SolveConst - matrix must be non-empty", a.rows(), a.cols(), -1, -1);

			Matrix<Type> mat(a);
			Matrix<Type> bmat = Utils::ColumnMatrixFromVector(b);
			SolveInPlace(mat, bmat);
			return bmat.VectorFromColumn(0);
		}
	};

	/// @brief LU decomposition solver with partial pivoting for linear systems Ax=b
	/// @tparam Type Numeric type (Real, Complex, etc.)
	/// @note Decomposes A=PLU where P is permutation, L is lower triangular, U is upper triangular
	/// @note Complexity: O(n³) decomposition, O(n²) per solve
	/// @warning Constructor performs decomposition - reuse solver for multiple RHS
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
		/// @brief Constructor - performs LU decomposition with partial pivoting
		/// @param inMatRef Coefficient matrix (copied internally, temporaries safe)
		/// @throws MatrixDimensionError if matrix is not square or is empty
		/// @throws SingularMatrixError if matrix is singular
		/// @throws MatrixNumericalError if non-finite values detected
		/// @note Uses partial pivoting with row interchanges for stability
		LUSolver(const Matrix<Type>& inMatRef) 
			: _n(inMatRef.rows()), _lu(inMatRef), _indx(_n)
		{
			// Dimension validation
			if (inMatRef.rows() != inMatRef.cols())
				throw MatrixDimensionError("LUSolver::ctor - matrix must be square", inMatRef.rows(), inMatRef.cols(), -1, -1);
			if (inMatRef.rows() == 0)
				throw MatrixDimensionError("LUSolver::ctor - matrix must be non-empty", 0, 0, -1, -1);
			int i, imax, j, k;
			Real big, temp;
			Type temp2;
			Vector<Type> vv(_n);
			
			_d = 1.0;
			// finding biggest element in each row, and saving its inverse in vv[]
			Real norm_a = 0.0;
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
					throw SingularMatrixError("LUSolver::ctor - Singular Matrix (zero row)");

				if (big > norm_a) norm_a = big;
				vv[i] = 1.0 / big;
			}
			Real singularity_threshold = std::numeric_limits<Real>::epsilon() * norm_a * _n;

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
				if (Abs(_lu[k][k]) < singularity_threshold) 
					throw SingularMatrixError("LUSolver::ctor - Singular Matrix", Abs(_lu[k][k]));
				
				for (i = k + 1; i < _n; i++) 
				{
					temp2 = _lu[i][k] /= _lu[k][k];
					
					// Check for non-finite multiplier (can happen from overflow/underflow)
					if constexpr (std::is_same_v<Type, Real>) {
						if (std::isnan(temp2) || std::isinf(temp2))
							throw MatrixNumericalError("Non-finite multiplier in LU decomposition");
					} else {
						Real mag = std::abs(temp2);
						if (std::isnan(mag) || std::isinf(mag))
							throw MatrixNumericalError("Non-finite multiplier in LU decomposition (Complex)");
					}
					
					for (j = k + 1; j < _n; j++)
						_lu[i][j] -= temp2 * _lu[k][j];
				}
			}
		}

		/// @brief Solve Ax=B with multiple right-hand sides using stored LU decomposition
		/// @param matRHS Right-hand side matrix
		/// @param outSol Solution matrix (output)
		/// @throws Exception if dimensions incompatible
		void Solve(Matrix<Type>& matRHS, Matrix<Type> &outSol)
		{
			int i, j;

			int n = matRHS.rows();
			int m = matRHS.cols();

			if (matRHS.rows() != _n)
				throw MatrixDimensionError("LUSolver::Solve - RHS row dimension mismatch", _n, m, matRHS.rows(), m);
			if (outSol.rows() != _n || outSol.cols() != m)
				throw MatrixDimensionError("LUSolver::Solve - output dimension mismatch", _n, m, outSol.rows(), outSol.cols());

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
				throw VectorDimensionError("LUSolverInPlace::Solve - bad sizes", b.size(), _n);
		}

		/// @brief Compute matrix inverse using stored LU decomposition
		/// @param ainv Matrix inverse (output, resized to n×n)
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


	/// @brief In-place LU decomposition solver - modifies input matrix directly
	/// @tparam Type Numeric type (Real, Complex, etc.)
	/// @warning Stores reference to input matrix and modifies it in-place during decomposition
	/// @note Use when you don't need to preserve original matrix (saves memory copy)
	template<class Type>
	class LUSolverInPlace
	{
	private:
		int _n;
		Matrix<Type> &_lu;
		std::vector<int> _indx;
		Real _d;

	public:
		/// @brief Constructor - performs LU decomposition in-place on input matrix
		/// @param inMatRef Coefficient matrix (modified in-place!)
		/// @throws MatrixDimensionError if matrix is not square or is empty
		/// @throws SingularMatrixError if matrix is singular
		/// @warning Matrix is modified during decomposition - original data destroyed
		LUSolverInPlace(Matrix<Type>& inMatRef)
			: _n(inMatRef.rows()), _lu(inMatRef), _indx(_n)
		{
			// Dimension validation
			if (inMatRef.rows() != inMatRef.cols())
				throw MatrixDimensionError("LUSolverInPlace::ctor - matrix must be square", inMatRef.rows(), inMatRef.cols(), -1, -1);
			if (inMatRef.rows() == 0)
				throw MatrixDimensionError("LUSolverInPlace::ctor - matrix must be non-empty", 0, 0, -1, -1);
			
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
					throw SingularMatrixError("LUSolverInPlace::ctor - Singular Matrix");

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
					throw SingularMatrixError("LUSolverInPlace::ctor - zero pivot after partial pivoting");

				for (i = k + 1; i < _n; i++) {
					temp2 = _lu[i][k] /= _lu[k][k];
					for (j = k + 1; j < _n; j++)
						_lu[i][j] -= temp2 * _lu[k][j];
				}
			}
		}
		// solving with Matrix RHS (ie. solving simultaneously for multiple RHS)
		// while using the LU decomposition already performed in the constructor

		/// @brief Solve Ax=B with multiple right-hand sides using in-place LU decomposition
		/// @param matRHS Right-hand side matrix
		/// @return Solution matrix
		Matrix<Type> Solve(Matrix<Type>& matRHS)
		{
			int i, j;

			int n = matRHS.rows();
			int m = matRHS.cols();

			if (matRHS.rows() != _n)
				throw MatrixDimensionError("LUSolver::Solve - RHS row dimension mismatch", _n, m, matRHS.rows(), m);

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

		/// @brief Solve Ax=b using in-place LU decomposition
		/// @param b Right-hand side vector
		/// @param x Solution vector (output)
		/// @return true if successful, false if dimensions incompatible
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

		/// @brief Solve Ax=b and return solution vector using in-place LU
		/// @param b Right-hand side vector
		/// @return Solution vector x
		/// @throws VectorDimensionError if dimensions incompatible
		Vector<Type> Solve(const Vector<Type>& b)
		{
			Vector<Type> x(b.size());
			if (Solve(b, x) == true)
				return x;
			else
				throw VectorDimensionError("LUSolverInPlace::Solve - bad sizes", b.size(), _n);
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
					throw SingularMatrixError("BandDiagonalSolver::ctor - Singular Matrix (zero pivot)");
				
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
			if (B.rows() != _n)
				throw MatrixDimensionError("BandDiagonalSolver::Solve - row count mismatch", 
				                          B.rows(), B.cols(), _n, -1);
			
			Matrix<Real> X(_n, B.cols());
			Vector<Real> b(_n), x(_n);
			
			for (int j = 0; j < B.cols(); j++)
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


	/// @brief Cholesky decomposition solver for symmetric positive-definite matrices
	/// @tparam Type Numeric type (Real, Complex, etc.)
	/// @note Decomposes A=LLᵀ where L is lower triangular (Cholesky factor)
	/// @note Complexity: O(n³/3) decomposition (half of LU), O(n²) per solve
	/// @warning Only works for symmetric positive-definite matrices
	template<class Type>
	class CholeskySolver
	{
	public:
		int n;
		Matrix<Type> el;

	public:
		/// @brief Constructor - performs Cholesky decomposition A=LLᵀ
		/// @param a Positive-definite symmetric matrix (only upper triangle needed)
		/// @throws MatrixDimensionError if matrix is not square or is empty
		/// @throws SingularMatrixError if matrix is not positive definite
		/// @note Cholesky factor L stored in lower triangle of el
		CholeskySolver(const Matrix<Type>& a) : n(a.rows()), el(a)
		{
			// Dimension validation
			if (a.rows() != a.cols())
				throw MatrixDimensionError("CholeskySolver::ctor - matrix must be square", a.rows(), a.cols(), -1, -1);
			if (a.rows() == 0)
				throw MatrixDimensionError("CholeskySolver::ctor - matrix must be non-empty", 0, 0, -1, -1);
			if (!a.isSymmetric())
				throw MatrixDimensionError("CholeskySolver::ctor - matrix must be symmetric", a.rows(), a.cols(), -1, -1);
			
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
		
		/// @brief Solve Ax=b using Cholesky decomposition (forward/backward substitution)
		/// @param b Right-hand side vector
		/// @param x Solution vector (output, resized to n)
		/// @throws VectorDimensionError if vector size incompatible
		/// @note Solves via Ly=b (forward), then Lᵀx=y (backward)
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
		
		/// @brief Solve Ax=b and return solution vector
		/// @param b Right-hand side vector
		/// @return Solution vector x
		Vector<Type> Solve(const Vector<Type>& b)
		{
			Vector<Type> x(b.size());
			Solve(b, x);
			return x;
		}
		
		/// @brief Compute matrix inverse A⁻¹ = (Lᵀ)⁻¹L⁻¹ using Cholesky decomposition
		/// @param ainv Matrix inverse (output, resized to n×n)
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
		
		/// @brief Compute log(det(A)) where A=LLᵀ
		/// @return log(determinant) = 2∑log(Lᵢᵢ) (numerically stable for large determinants)
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

	///////////////////////////////////////////////////////////////////////////////////////////
	// Free SolveDetailed functions - full workflow (decompose + solve) with structured results
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Gauss-Jordan solve with full diagnostics
	/// @return LinearSolverResult with solution vector, status, timing, and optional residual
	template<typename Type>
	LinearSolverResult<Type> GaussJordanSolveDetailed(
		const Matrix<Type>& A, const Vector<Type>& b,
		const LinearSolverConfig& config = {})
	{
		using ResultType = LinearSolverResult<Type>;
		return LinearSolverDetail::ExecuteLinearSolverDetailed<ResultType>(
			"GaussJordanSolver", config,
			[&](ResultType& result) {
				result.solution = GaussJordanSolver<Type>::SolveConst(A, b);
				result.function_evaluations = A.rows() * A.rows() * A.rows(); // O(n³)

				if (config.estimate_error) {
					// Compute residual norm ||Ax - b||
					Vector<Type> residual = A * result.solution - b;
					result.residual_norm = residual.NormL2();
				}
			});
	}

	/// @brief LU decomposition solve with full diagnostics
	/// @return LinearSolverResult with solution vector, status, timing, and optional residual
	template<typename Type>
	LinearSolverResult<Type> LUSolveDetailed(
		const Matrix<Type>& A, const Vector<Type>& b,
		const LinearSolverConfig& config = {})
	{
		using ResultType = LinearSolverResult<Type>;
		return LinearSolverDetail::ExecuteLinearSolverDetailed<ResultType>(
			"LUSolver", config,
			[&](ResultType& result) {
				LUSolver<Type> solver(A);
				result.solution = solver.Solve(b);
				result.function_evaluations = A.rows() * A.rows() * A.rows(); // O(n³) decomp + O(n²) solve

				if (config.estimate_error) {
					Vector<Type> residual = A * result.solution - b;
					result.residual_norm = residual.NormL2();
				}
			});
	}

	/// @brief Cholesky decomposition solve with full diagnostics (symmetric positive-definite)
	/// @return LinearSolverResult with solution vector, status, timing, and optional residual
	template<typename Type>
	LinearSolverResult<Type> CholeskySolveDetailed(
		const Matrix<Type>& A, const Vector<Type>& b,
		const LinearSolverConfig& config = {})
	{
		using ResultType = LinearSolverResult<Type>;
		return LinearSolverDetail::ExecuteLinearSolverDetailed<ResultType>(
			"CholeskySolver", config,
			[&](ResultType& result) {
				CholeskySolver<Type> solver(A);
				result.solution = solver.Solve(b);
				result.function_evaluations = A.rows() * A.rows() * A.rows() / 3; // O(n³/3)

				if (config.estimate_error) {
					Vector<Type> residual = A * result.solution - b;
					result.residual_norm = residual.NormL2();
				}
			});
	}

	/// @brief Band diagonal solve with full diagnostics
	/// @return LinearSolverResult with solution vector, status, timing, and optional residual
	inline LinearSolverResult<Real> BandDiagonalSolveDetailed(
		const BandDiagonalMatrix& A, const Vector<Real>& b,
		const LinearSolverConfig& config = {})
	{
		using ResultType = LinearSolverResult<Real>;
		return LinearSolverDetail::ExecuteLinearSolverDetailed<ResultType>(
			"BandDiagonalSolver", config,
			[&](ResultType& result) {
				BandDiagonalSolver solver(A);
				result.solution = solver.Solve(b);
				int n = A.GetDimension();
				int m = A.GetLowerBandwidth() + A.GetUpperBandwidth() + 1;
				result.function_evaluations = n * m * m; // O(n*m²)

				if (config.estimate_error) {
					Vector<Real> residual = A * result.solution - b;
					result.residual_norm = residual.NormL2();
				}
			});
	}

} // namespace MML

#endif // MML_LINEAR_ALG_DIRECT_H
