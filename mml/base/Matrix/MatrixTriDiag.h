///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixTriDiag.h                                                     ///
///  Description: Tridiagonal matrix with O(n) storage and specialized solvers        ///
///               Thomas algorithm implementation for efficient solving               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_TRIDIAG_H
#define MML_MATRIX_TRIDIAG_H

#include <mml/base/PrintFormatting.h>
#include <mml/MMLBase.h>
#include <algorithm>
#include <mml/MMLExceptions.h>
#include <mml/base/Matrix/MatrixPrintFormat.h>
#include <mml/base/Matrix/Matrix.h>

#include <mml/base/Vector/Vector.h>

#include <limits>

// Standard headers - include what we use
#include <iomanip>
#include <initializer_list>
#include <iostream>

namespace MML
{
	template<class Type>
	class TridiagonalMatrix
	{
	private:
		int _dim;
		Vector<Type> _belowDiag;
		Vector<Type> _diag;
		Vector<Type> _aboveDiag;

	public:
		TridiagonalMatrix(int dim, const Vector<Type>& a, const Vector<Type>& diag, const Vector<Type>& c) : _dim(dim), _belowDiag(a), _diag(diag), _aboveDiag(c)
		{
			if( dim < 1)
				throw MatrixDimensionError("TridiagonalMatrix constructor: dimension must be >= 1");

			if (a.size() != dim || diag.size() != dim || c.size() != dim)
				throw MatrixDimensionError("TridiagonalMatrix constructor: vector dimensions must match matrix dimension");
		}

		TridiagonalMatrix(int dim, std::initializer_list<Type> values) : _dim(dim), _belowDiag(std::max(dim, 0)), _diag(std::max(dim, 0)), _aboveDiag(std::max(dim, 0))
		{
			if (dim < 1)
				throw MatrixDimensionError("TridiagonalMatrix constructor: dimension must be >= 1");
			if (static_cast<int>(values.size()) != dim * 3 - 2)
				throw MatrixDimensionError("TridiagonalMatrix constructor: initializer list size must be 3*dim-2");

			auto val = values.begin();
			if (dim == 1) {
				_belowDiag[0] = Type{};
				_diag[0] = *val;
				_aboveDiag[0] = Type{};
				return;
			}
			_belowDiag[0] = Type{};
			_diag[0] = *val++;
			_aboveDiag[0] = *val++;
			for (int i = 1; i < dim - 1; ++i)
			{
				_belowDiag[i] = *val++;
				_diag[i] = *val++;
				_aboveDiag[i] = *val++;
			}
			_belowDiag[dim - 1] = *val++;
			_diag[dim - 1] = *val++;
			_aboveDiag[dim - 1] = Type{};
		}

		int rows() const noexcept { return _dim; }
		int cols() const noexcept { return _dim; }

		Type  operator()(int i, int j) const {
			if (i < 0 || i >= _dim || j < 0 || j >= _dim)
				throw MatrixAccessBoundsError("TridiagonalMatrix::operator() const", i, j, _dim, _dim);
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _aboveDiag[i];
			else if (i == j + 1)
				return _belowDiag[i];
			else
				return 0.0;
		}
		Type& operator()(int i, int j) {
			if (i < 0 || i >= _dim || j < 0 || j >= _dim)
				throw MatrixAccessBoundsError("TridiagonalMatrix::operator()", i, j, _dim, _dim);
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _aboveDiag[i];
			else if (i == j + 1)
				return _belowDiag[i];
			else
				throw MatrixAccessBoundsError("TridiagonalMatrix::operator() - element outside the three diagonals is not writable", i, j, _dim, _dim);
		}

		/// @brief Inclusive absolute-tolerance comparison of all stored band entries.
		/// @return False for mismatched dimensions, nonfinite entries or invalid tolerance.
		bool IsEqualTo(const TridiagonalMatrix& b, Real eps = Precision::MatrixIsEqualTolerance) const
		{
			if (_dim != b._dim || !std::isfinite(eps) || eps < 0)
				return false;

			const auto equal = [eps](const Type& first, const Type& second) {
				if constexpr (std::is_integral_v<Type>) {
					const auto difference = first >= second
						? static_cast<unsigned long long>(first) - static_cast<unsigned long long>(second)
						: static_cast<unsigned long long>(second) - static_cast<unsigned long long>(first);
					// Compare in integer space: converting a large difference to Real can round it down to eps.
					return eps >= static_cast<Real>(std::numeric_limits<unsigned long long>::max()) ||
						difference <= static_cast<unsigned long long>(eps);
				} else {
					return IsFiniteValue(first) && IsFiniteValue(second) && Abs(first - second) <= eps;
				}
			};
			for (int i = 0; i < _dim; i++) {
				if (!equal(_diag[i], b._diag[i]) || !equal(_belowDiag[i], b._belowDiag[i]) ||
					!equal(_aboveDiag[i], b._aboveDiag[i]))
					return false;
			}
			return true;
		}

		TridiagonalMatrix operator+(const TridiagonalMatrix& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("TridiagonalMatrix::operator+() - must be same dim", _dim, _dim, b._dim, b._dim);

			TridiagonalMatrix temp(*this);
			for (int i = 0; i < _belowDiag.size(); i++)
				temp._belowDiag[i] += b._belowDiag[i];
			for (int i = 0; i < _diag.size(); i++)
				temp._diag[i] += b._diag[i];
			for (int i = 0; i < _aboveDiag.size(); i++)
				temp._aboveDiag[i] += b._aboveDiag[i];

			return temp;
		}
		TridiagonalMatrix operator-(const TridiagonalMatrix& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("TridiagonalMatrix::operator+() - must be same dim", _dim, _dim, b._dim, b._dim);

			TridiagonalMatrix temp(*this);
			for (int i = 0; i < _belowDiag.size(); i++)
				temp._belowDiag[i] -= b._belowDiag[i];
			for (int i = 0; i < _diag.size(); i++)
				temp._diag[i] -= b._diag[i];
			for (int i = 0; i < _aboveDiag.size(); i++)
				temp._aboveDiag[i] -= b._aboveDiag[i];

			return temp;
		}

		friend TridiagonalMatrix operator*(const TridiagonalMatrix& a, Type b)
		{
			TridiagonalMatrix	ret(a);

			for (int i = 0; i < ret._belowDiag.size(); i++)
				ret._belowDiag[i] *= b;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] *= b;
			for (int i = 0; i < ret._aboveDiag.size(); i++)
				ret._aboveDiag[i] *= b;

			return ret;
		}
		friend TridiagonalMatrix operator*(Type a, const TridiagonalMatrix& b)
		{
			TridiagonalMatrix	ret(b);

			for (int i = 0; i < ret._belowDiag.size(); i++)
				ret._belowDiag[i] *= a;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] *= a;
			for (int i = 0; i < ret._aboveDiag.size(); i++)
				ret._aboveDiag[i] *= a;

			return ret;
		}
		friend TridiagonalMatrix operator/(const TridiagonalMatrix& a, Type b)
		{
			TridiagonalMatrix	ret(a);

			for (int i = 0; i < ret._belowDiag.size(); i++)
				ret._belowDiag[i] /= b;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] /= b;
			for (int i = 0; i < ret._aboveDiag.size(); i++)
				ret._aboveDiag[i] /= b;

			return ret;
		}

		/// @brief Returns the transpose (copy). Swaps and shifts the above/below diagonals.
		[[nodiscard]] TridiagonalMatrix Transposed() const
		{
			// For a tridiagonal matrix, transpose swaps and shifts the above and below diagonals
			// belowDiag_new[i] = aboveDiag_old[i-1] (with belowDiag_new[0] = 0)
			// aboveDiag_new[i] = belowDiag_old[i+1] (with aboveDiag_new[dim-1] = 0)
			Vector<Type> newBelowDiag(_dim);
			Vector<Type> newAboveDiag(_dim);

			newBelowDiag[0] = 0.0;
			for (int i = 1; i < _dim; i++)
				newBelowDiag[i] = _aboveDiag[i - 1];

			for (int i = 0; i < _dim - 1; i++)
				newAboveDiag[i] = _belowDiag[i + 1];
			newAboveDiag[_dim - 1] = 0.0;

			return TridiagonalMatrix(_dim, newBelowDiag, _diag, newAboveDiag);
		}

		/// @deprecated Use Transposed() (CONTRACT_POLICY §11: noun/participle returns a copy).
		[[deprecated("use Transposed()")]]
		TridiagonalMatrix GetTranspose() const { return Transposed(); }

		/// @brief Converts to full matrix representation (dense n×n copy).
		Matrix<Type> ToFullMatrix() const
		{
			Matrix<Type> full(_dim, _dim);
			for (int i = 0; i < _dim; i++)
				for (int j = 0; j < _dim; j++)
					full(i, j) = (*this)(i, j);
			return full;
		}


		void Solve(const Vector<Type>& rhs, Vector<Type>& sol) const
		{
			int j, n = _belowDiag.size();
			if (rhs.size() != n)
				throw VectorDimensionError("TridiagonalMatrix::Solve - rhs size must equal matrix dimension", rhs.size(), n);
			if (sol.size() != n)
				throw VectorDimensionError("TridiagonalMatrix::Solve - solution vector size must equal matrix dimension", sol.size(), n);
			if (n == 0)
				return;

			Type bet;  // Must be Type to work with complex matrices
			Vector<Type> gam(n);

			Real matrixScale = Real{0};
			const auto coefficientMagnitude = [](const Type& value) {
				if constexpr (MMLComplex<Type>)
					return std::max(Abs(value.real()), Abs(value.imag()));
				else
					return Abs(value);
			};
			for (int i = 0; i < n; i++) {
				if (!IsFiniteValue(_diag[i]) || !IsFiniteValue(rhs[i]) ||
					(i > 0 && !IsFiniteValue(_belowDiag[i])) ||
					(i + 1 < n && !IsFiniteValue(_aboveDiag[i])))
					throw NumericInputError("TridiagonalMatrix::Solve - matrix and rhs must be finite");
				matrixScale = std::max(matrixScale, coefficientMagnitude(_diag[i]));
				if (i > 0) matrixScale = std::max(matrixScale, coefficientMagnitude(_belowDiag[i]));
				if (i + 1 < n) matrixScale = std::max(matrixScale, coefficientMagnitude(_aboveDiag[i]));
			}
			if (matrixScale == 0)
				throw SingularMatrixError("TridiagonalMatrix::Solve - zero pivot at element 0", 0, 0);
			// Normalize the Thomas recurrence before elimination. Both epsilon*scale
			// and unscaled elimination products can underflow in subnormal systems.
			// Integral instantiations retain their existing integer-division arithmetic.
			// A real divisor keeps complex normalization componentwise, avoiding complex-division overflow.
			using CoefficientScale = std::conditional_t<MMLComplex<Type>, Real, Type>;
			const CoefficientScale coefficientScale = std::is_integral_v<Type> ? CoefficientScale{1} : CoefficientScale(matrixScale);
			const Real pivotScale = std::is_integral_v<Type> ? matrixScale : Real{1};
			const auto unsafePivot = [pivotScale](const Type& pivot) {
				const Real magnitude = Abs(pivot);
				return !std::isfinite(magnitude) || magnitude == 0 ||
					magnitude / pivotScale < Precision::Eps;
			};

			bet = _diag[0] / coefficientScale;
			if (unsafePivot(bet))
				throw SingularMatrixError("TridiagonalMatrix::Solve - zero pivot at element 0", 0, 0);

			sol[0] = (rhs[0] / coefficientScale) / bet;

			for (j = 1; j < n; j++) {
				gam[j] = (_aboveDiag[j - 1] / coefficientScale) / bet;
				bet = (_diag[j] / coefficientScale) - (_belowDiag[j] / coefficientScale) * gam[j];

				if (unsafePivot(bet))
					throw SingularMatrixError("TridiagonalMatrix::Solve - zero pivot during elimination", Real(j), j);

				sol[j] = ((rhs[j] / coefficientScale) - (_belowDiag[j] / coefficientScale) * sol[j - 1]) / bet;
			}
			for (j = (n - 2); j >= 0; j--)
				sol[j] -= gam[j + 1] * sol[j + 1];
		}
		Vector<Type> Solve(const Vector<Type>& rhs) const
		{
			Vector<Type> sol(rhs.size());
			Solve(rhs, sol);
			return sol;
		}

		// New formatted print with MatrixPrintFormat
		void Print(std::ostream& stream, const MatrixPrintFormat& fmt = MatrixPrintFormat::Default()) const
		{
			PrintFormatting::Matrix(stream, rows(), cols(), [&](int i, int j) { return (*this)(i, j); }, fmt, " (Tridiagonal)");
		}

		/// @brief Print with explicit width and precision.
		void   Print(std::ostream& stream, int width, int precision) const
		{
			MatrixPrintFormat fmt = MatrixPrintFormat::Default();
			fmt.width = width;
			fmt.precision = precision;
			Print(stream, fmt);
		}
	};
}

#endif