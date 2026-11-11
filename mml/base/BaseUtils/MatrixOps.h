///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixOps.h                                                         ///
///  Description: Matrix operation utilities - creation, functions, properties        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_MATRIX_OPS_H
#define MML_MATRIX_OPS_H

#include <algorithm>
#include <vector>
#include <complex>
#include <limits>

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>
#include <mml/MMLNumericValidation.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/BaseUtils/VectorOps.h>

namespace MML
{
	namespace Utils
	{
		namespace Detail
		{
			template<class Type>
			int RectangularVectorExtent(const std::vector<Vector<Type>>& vectors)
			{
				if (vectors.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
					throw MatrixAllocationError("Matrix factory - vector count exceeds integer dimensions", -1, -1);
				if (vectors.empty())
					return 0;
				const int extent = vectors.front().size();
				for (const auto& vector : vectors)
					if (vector.size() != extent)
						throw MatrixDimensionError("Matrix factory - vectors must have equal lengths",
							static_cast<int>(vectors.size()), extent, vector.size(), -1);
				return extent;
			}
		}

		// ============================================================================
		// Creating Matrix from Vector
		// ============================================================================

		/// @brief Creates 1×n row matrix from vector.
		/// @tparam Type Element type
		/// @param b Source vector
		/// @return Matrix with single row containing vector elements
		template<class Type> 
		static Matrix<Type> RowMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret(1, (int)b.size());
			for (int i = 0; i < b.size(); i++)
				ret[0][i] = b[i];

			return ret;
		}

		/// @brief Creates n×1 column matrix from vector.
		/// @tparam Type Element type
		/// @param b Source vector
		/// @return Matrix with single column containing vector elements
		template<class Type> 
		static Matrix<Type> ColumnMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret((int)b.size(), 1);
			for (int i = 0; i < b.size(); i++)
				ret[i][0] = b[i];

			return ret;
		}

		/// @brief Creates n×n diagonal matrix from vector.
		/// @tparam Type Element type
		/// @param b Source vector (becomes diagonal elements)
		/// @return Square diagonal matrix with b[i] on diagonal
		template<class Type> 
		static Matrix<Type> DiagonalMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret((int)b.size(), (int)b.size());
			for (int i = 0; i < b.size(); i++)
				ret[i][i] = b[i];

			return ret;
		}

		/// @brief Creates matrix with vectors as rows.
		/// @tparam Type Element type
		/// @param b List of vectors (each becomes a row)
		/// @return Matrix with b.size() rows
		template<class Type> 
		static Matrix<Type> MatrixFromVectorsInRows(const std::vector<Vector<Type>>& b)
		{
			const int extent = Detail::RectangularVectorExtent(b);
			Matrix<Type> ret(static_cast<int>(b.size()), extent);
			for (int i = 0; i < ret.rows(); i++)
				for (int j = 0; j < extent; j++)
					ret[i][j] = b[i][j];
			return ret;
		}

		/// @brief Creates matrix with vectors as columns.
		/// @tparam Type Element type
		/// @param b List of vectors (each becomes a column)
		/// @return Matrix with b.size() columns
		template<class Type> 
		static Matrix<Type> MatrixFromVectorsInColumns(const std::vector<Vector<Type>>& b)
		{
			const int extent = Detail::RectangularVectorExtent(b);
			Matrix<Type> ret(extent, static_cast<int>(b.size()));
			for (int i = 0; i < ret.cols(); i++)
				for (int j = 0; j < extent; j++)
					ret[j][i] = b[i][j];
			return ret;
		}

		/// @brief Creates 1×N row matrix from fixed-size vector.
		/// @tparam Type Element type
		/// @tparam N Vector dimension
		/// @param b Source vector
		/// @return MatrixNM<1,N> row matrix
		template<class Type, int N> 
		MatrixNM<Type, 1, N> RowMatrixFromVector(const VectorN<Type, N>& b)
		{
			MatrixNM<Type, 1, N>  ret;
			for (int j = 0; j < N; j++)
				ret(0, j) = b[j];

			return ret;
		}

		/// @brief Creates N×1 column matrix from fixed-size vector.
		/// @tparam Type Element type
		/// @tparam N Vector dimension
		/// @param b Source vector
		/// @return MatrixNM<N,1> column matrix
		template<class Type, int N> 
		MatrixNM<Type, N, 1> ColumnMatrixFromVector(const VectorN<Type, N>& b)
		{
			MatrixNM<Type, N, 1>  ret;
			for (int i = 0; i < N; i++)
				ret(i, 0) = b[i];
			return ret;
		}

		/// @brief Creates N×N diagonal matrix from fixed-size vector.
		/// @tparam Type Element type
		/// @tparam N Vector dimension
		/// @param b Source vector (becomes diagonal elements)
		/// @return MatrixNM<N,N> diagonal matrix
		template<class Type, int N> 
		MatrixNM<Type, N, N> DiagonalMatrixFromVector(const VectorN<Type, N>& b)
		{
			MatrixNM<Type, N, N> ret;
			for (int i = 0; i < N; i++)
				ret[i][i] = b[i];
			return ret;
		}

		// ============================================================================
		// Matrix Algebraic Operations
		// ============================================================================

		/// @brief Computes matrix commutator [A, B] = AB - BA.
		/// @tparam Type Element type
		/// @param a First matrix
		/// @param b Second matrix
		/// @return Commutator matrix AB - BA
		template<class Type> 
		static Matrix<Type> Commutator(const Matrix<Type>& a, const Matrix<Type>& b)
		{
			return a * b - b * a;
		}

		/// @brief Computes matrix anti-commutator {A, B} = AB + BA.
		/// @tparam Type Element type
		/// @param a First matrix
		/// @param b Second matrix
		/// @return Anti-commutator matrix AB + BA
		template<class Type> 
		static Matrix<Type> AntiCommutator(const Matrix<Type>& a, const Matrix<Type>& b)
		{
			return a * b + b * a;
		}

		/// @brief Decomposes matrix into symmetric and antisymmetric parts.
		/// @details For any matrix A: A = (A + Aᵀ)/2 + (A - Aᵀ)/2 = Sym + Antisym
		/// @tparam Type Element type
		/// @param orig Original square matrix
		/// @param[out] outSym Symmetric part (A + Aᵀ)/2
		/// @param[out] outAntiSym Antisymmetric part (A - Aᵀ)/2
		/// @throws MatrixDimensionError if matrix is not square
		template<class Type> 
		static void MatrixDecomposeToSymAntisym(const Matrix<Type>& orig, 
		                                        Matrix<Type>& outSym, 
		                                        Matrix<Type>& outAntiSym)
		{
			if (orig.rows() != orig.cols())
				throw MatrixDimensionError("MatrixDecompose - matrix must be square", orig.rows(), orig.cols(), -1, -1);

			auto transp = orig.Transposed();
			outSym = (orig + transp) * 0.5;
			outAntiSym = (orig - transp) * 0.5;
		}

		/// @brief Computes the maximum absolute entry-wise difference over the shared extent.
		template<class Type>
		Real MaxAbsDiff(const Matrix<Type>& left, const Matrix<Type>& right)
		{
			Real maximum = static_cast<MML::Real>(0.0);
			const int rows = std::min(left.rows(), right.rows());
			const int cols = std::min(left.cols(), right.cols());
			for (int row = 0; row < rows; ++row)
				for (int col = 0; col < cols; ++col) {
					if (!IsFiniteValue(left(row, col)) || !IsFiniteValue(right(row, col)))
						throw NumericInputError("MaxAbsDiff - entries must be finite");
					const Real difference = static_cast<Real>(std::abs(left(row, col) - right(row, col)));
					if (!std::isfinite(difference))
						throw NumericalMethodError("MaxAbsDiff - difference is not representable");
					maximum = std::max(maximum, difference);
				}
			return maximum;
		}

		/// @brief Computes Q^T A Q for a real orthogonal similarity transform.
		template<class Type>
		Matrix<Type> SimilarityTransform(const Matrix<Type>& transformation, const Matrix<Type>& matrix)
		{
			if (matrix.rows() != matrix.cols() || transformation.rows() != matrix.rows() ||
				transformation.cols() != matrix.cols())
				throw MatrixDimensionError("SimilarityTransform - equal square dimensions required",
					transformation.rows(), transformation.cols(), matrix.rows(), matrix.cols());
			const int size = matrix.rows();
			Matrix<Type> temporary(size, size);
			Matrix<Type> result(size, size);
			for (int row = 0; row < size; ++row)
				for (int col = 0; col < size; ++col)
					for (int index = 0; index < size; ++index)
						temporary(row, col) += matrix(row, index) * transformation(index, col);
			for (int row = 0; row < size; ++row)
				for (int col = 0; col < size; ++col)
					for (int index = 0; index < size; ++index)
						result(row, col) += transformation(index, row) * temporary(index, col);
			return result;
		}

		// ============================================================================
		// Matrix Functions (Power Series)
		// ============================================================================

		/// @brief Computes matrix exponential exp(A) via power series.
		/// @details exp(A) = I + A + A²/2! + A³/3! + ... (truncated at n terms)
		/// @tparam Type Element type
		/// @param a Square matrix
		/// @param n Number of terms in series expansion (default 10)
		/// @return Matrix exponential exp(A)
		template<class Type> 
		static Matrix<Type> Exp(const Matrix<Type>& a, int n = 10)
		{
			Matrix<Type> ret(a.rows(), a.cols());
			Matrix<Type> a_pow_n(a);

			Type fact = Type(1);
			ret.SetIdentity();

			for (int i = 1; i <= n; i++)
			{
				fact *= ScalarDetail::ScalarFromInteger<Type>(i);
				ret += a_pow_n / fact;

				a_pow_n = a_pow_n * a;
			}

			return ret;
		}

		/// @brief Computes matrix sine sin(A) via power series.
		/// @details sin(A) = A - A³/3! + A⁵/5! - ... (odd powers, alternating signs)
		/// @tparam Type Element type
		/// @param a Square matrix
		/// @param n Number of terms in series expansion (default 10)
		/// @return Matrix sine sin(A)
		template<class Type> 
		static Matrix<Type> Sin(const Matrix<Type>& a, int n = 10)
		{
			Matrix<Type> ret(a.rows(), a.cols());		// initialized to zero matrix!
			Matrix<Type> a_pow_n(a);
			
			Type fact = Type(1);
			
			for (int i = 1; i <= n; i += 2)
			{
				ret += a_pow_n / fact;
				
				a_pow_n = Type(-1) * a_pow_n * a * a;
				fact *= ScalarDetail::ScalarFromInteger<Type>(i + 1) * ScalarDetail::ScalarFromInteger<Type>(i + 2);
			}
			return ret;
		}

		/// @brief Computes matrix cosine cos(A) via power series.
		/// @details cos(A) = I - A²/2! + A⁴/4! - ... (even powers, alternating signs)
		/// @tparam Type Element type
		/// @param a Square matrix
		/// @param n Number of terms in series expansion (default 10)
		/// @return Matrix cosine cos(A)
		template<class Type> 
		static Matrix<Type> Cos(const Matrix<Type>& a, int n = 10)
		{
			Matrix<Type> ret(a.rows(), a.cols());
			Matrix<Type> a_pow_n(a.rows(), a.cols());

			Type fact = Type(1);
			ret.SetIdentity();
			a_pow_n.SetIdentity();

			for (int i = 1; i <= n; i+=2)
			{
				a_pow_n = Type(-1) * a_pow_n * a * a;
				fact *= ScalarDetail::ScalarFromInteger<Type>(i) * ScalarDetail::ScalarFromInteger<Type>(i + 1);

				ret += a_pow_n / fact;
			}
			return ret;
		}

		/// @brief Matrix hyperbolic sine: sinh(A) = (exp(A) - exp(-A)) / 2
		/// Power series: sinh(A) = A + A³/3! + A⁵/5! + ...
		/// @param a Square matrix
		/// @param n Number of terms in series expansion
		/// @return Matrix sinh(A)
		template<class Type> 
		static Matrix<Type> Sinh(const Matrix<Type>& a, int n = 10)
		{
			Matrix<Type> ret(a.rows(), a.cols());		// initialized to zero matrix
			Matrix<Type> a_pow_n(a);
			
			Type fact = Type(1);
			
			// sinh(x) = x + x³/3! + x⁵/5! + ... (odd powers, all positive)
			for (int i = 1; i <= n; i += 2)
			{
				ret += a_pow_n / fact;
				
				a_pow_n = a_pow_n * a * a;  // Next odd power
				fact *= ScalarDetail::ScalarFromInteger<Type>(i + 1) * ScalarDetail::ScalarFromInteger<Type>(i + 2);
			}
			return ret;
		}

		/// @brief Matrix hyperbolic cosine: cosh(A) = (exp(A) + exp(-A)) / 2
		/// Power series: cosh(A) = I + A²/2! + A⁴/4! + ...
		/// @param a Square matrix
		/// @param n Number of terms in series expansion
		/// @return Matrix cosh(A)
		template<class Type> 
		static Matrix<Type> Cosh(const Matrix<Type>& a, int n = 10)
		{
			Matrix<Type> ret(a.rows(), a.cols());
			Matrix<Type> a_pow_n(a.rows(), a.cols());

			Type fact = Type(1);
			ret.SetIdentity();
			a_pow_n.SetIdentity();

			// cosh(x) = 1 + x²/2! + x⁴/4! + ... (even powers, all positive)
			for (int i = 2; i <= n; i += 2)
			{
				a_pow_n = a_pow_n * a * a;  // Next even power
				fact *= ScalarDetail::ScalarFromInteger<Type>(i - 1) * ScalarDetail::ScalarFromInteger<Type>(i);

				ret += a_pow_n / fact;
			}
			return ret;
		}

		// ============================================================================
		// Complex Matrix Operations
		// ============================================================================

		/// @brief Extracts real part of complex matrix.
		/// @param a Complex matrix
		/// @return Real matrix with a[i][j].real()
		static inline Matrix<Real> GetRealPart(const Matrix<Complex>& a)
		{
			Matrix<Real> ret(a.rows(), a.cols());

			for (int i = 0; i < ret.rows(); i++)
				for (int j = 0; j < ret.cols(); j++) {
					ret[i][j] = a[i][j].real();
				}

			return ret;
		}

		/// @brief Extracts imaginary part of complex matrix.
		/// @param a Complex matrix
		/// @return Real matrix with a[i][j].imag()
		static inline Matrix<Real> GetImagPart(const Matrix<Complex>& a)
		{
			Matrix<Real> ret(a.rows(), a.cols());

			for (int i = 0; i < ret.rows(); i++)
				for (int j = 0; j < ret.cols(); j++) {
					ret[i][j] = a[i][j].imag();
				}

			return ret;
		}

		/// @brief Computes conjugate transpose (Hermitian adjoint) of complex matrix.
		/// @details Result[j][i] = conj(mat[i][j])
		/// @param mat Complex matrix
		/// @return Conjugate transpose matrix
		static inline Matrix<Complex> GetConjugateTranspose(const Matrix<Complex>& mat)
		{
			Matrix<Complex> ret(mat.cols(), mat.rows());

			for (int i = 0; i < mat.rows(); i++)
				for (int j = 0; j < mat.cols(); j++)
					ret[j][i] = std::conj(mat[i][j]);

			return ret;
		}

		/// @brief Creates complex matrix from real matrix (zero imaginary parts).
		/// @param mat Real matrix
		/// @return Complex matrix with mat[i][j] + 0i
		static inline Matrix<Complex> CmplxMatFromRealMat(const Matrix<Real>& mat)
		{
			Matrix<Complex> mat_cmplx(mat.rows(), mat.cols());

			for (int i = 0; i < mat.rows(); i++)
				for (int j = 0; j < mat.cols(); j++)
					mat_cmplx[i][j] = Complex(mat(i, j), 0.0);

			return mat_cmplx;
		}

		/// @brief Tests if complex matrix has all real entries (zero imaginary parts).
		/// @param mat Complex matrix
		/// @return True if all imaginary parts are exactly zero
		static inline bool IsComplexMatReal(const Matrix<Complex>& mat)
		{
			for (int i = 0; i < mat.rows(); i++)
				for (int j = 0; j < mat.cols(); j++)
					if (mat[i][j].imag() != 0.0)
						return false;

			return true;
		}

		namespace Detail
		{
			template<class TVector, class BasisComponent>
			bool OrthonormalizeColumn(TVector& vector, int basisCount, BasisComponent basisComponent, Real tolerance)
			{
				Real scale = 0;
				for (int row = 0; row < vector.size(); ++row) {
					ValidateFinite(vector[row], "GramSchmidt input");
					scale = std::max(scale, std::abs(vector[row]));
				}
				if (scale == 0)
					return false;
				for (int row = 0; row < vector.size(); ++row)
					vector[row] /= scale;
				const Real inputNorm = vector.NormL2();
				// Two-pass MGS; see Bjorck, Numerical Methods for Least Squares Problems.
				for (int pass = 0; pass < 2; ++pass)
					for (int previous = 0; previous < basisCount; ++previous) {
						Real projection = 0;
						for (int row = 0; row < vector.size(); ++row)
							projection += vector[row] * basisComponent(previous, row);
						for (int row = 0; row < vector.size(); ++row)
							vector[row] -= projection * basisComponent(previous, row);
					}
				const Real norm = vector.NormL2();
				if (IsNegligibleRelativeTo(norm, inputNorm, tolerance))
					return false;
				for (int row = 0; row < vector.size(); ++row)
					vector[row] /= norm;
				return true;
			}
		}
	} // namespace Utils

	inline Matrix<Real> GramSchmidt(const Matrix<Real>& matrix, Real tolerance = Precision::LinearDependenceTolerance)
	{
		ValidateTolerance(tolerance, "GramSchmidt");
		const int rows = matrix.rows();
		const int cols = matrix.cols();
		Matrix<Real> result(rows, cols);
		for (int col = 0; col < cols; ++col) {
			Vector<Real> vector(rows);
			for (int row = 0; row < rows; ++row)
				vector[row] = matrix(row, col);
			if (Utils::Detail::OrthonormalizeColumn(vector, col,
				[&](int previous, int row) { return result(row, previous); }, tolerance))
				for (int row = 0; row < rows; ++row)
					result(row, col) = vector[row];
		}
		return result;
	}

	template<int N>
	std::vector<VectorN<Real, N>> GramSchmidtVectors(
			const std::vector<VectorN<Real, N>>& vectors, Real tolerance = Precision::LinearDependenceTolerance)
	{
		ValidateTolerance(tolerance, "GramSchmidtVectors");
		if (vectors.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
			throw ArgumentError("GramSchmidtVectors - vector count exceeds integer dimensions");
		std::vector<VectorN<Real, N>> result;
		result.reserve(vectors.size());
		for (const auto& input : vectors) {
			VectorN<Real, N> vector = input;
			if (Utils::Detail::OrthonormalizeColumn(vector, static_cast<int>(result.size()),
				[&](int previous, int row) { return result[previous][row]; }, tolerance))
				result.push_back(vector);
		}
		return result;
	}

	template<int N>
	bool IsOrthonormal(
			const std::vector<VectorN<Real, N>>& vectors,
			Real tolerance = Precision::OrthogonalityTolerance)
	{
		if (!std::isfinite(tolerance) || tolerance < 0)
			return false;
		for (size_t left = 0; left < vectors.size(); ++left) {
			for (int component = 0; component < N; ++component)
				if (!std::isfinite(vectors[left][component]))
					return false;
			if (!(std::abs(vectors[left].NormL2() - static_cast<MML::Real>(1.0)) <= tolerance))
				return false;
			for (size_t right = left + 1; right < vectors.size(); ++right)
				if (!(std::abs(Utils::ScalarProduct<N>(vectors[left], vectors[right])) <= tolerance))
					return false;
		}
		return true;
	}
} // namespace MML

#endif // MML_MATRIX_OPS_H
