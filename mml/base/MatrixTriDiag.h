///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixTriDiag.h                                                     ///
///  Description: Tridiagonal matrix with O(n) storage and specialized solvers        ///
///               Thomas algorithm implementation for efficient solving               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_TRIDIAG_H
#define MML_MATRIX_TRIDIAG_H

#include "MMLBase.h"

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
		TridiagonalMatrix(int dim, const Vector<Type>& a, const Vector<Type>& diag, const Vector<Type>& c) : _belowDiag(a), _diag(diag), _aboveDiag(c), _dim(dim)
		{
			if( dim < 3)
				throw("Error in TridiagonalMatrix constructor: dim < 3");

			if (a.size() != dim || diag.size() != dim || c.size() != dim)
				throw("Error in TridiagonalMatrix constructor: wrong dimensions");
		}

		TridiagonalMatrix(int dim, std::initializer_list<Type> values) : _dim(dim), _belowDiag(dim), _diag(dim), _aboveDiag(dim)
		{
			if (values.size() != dim * 3 - 2)
				throw("Error in TridiagonalMatrix constructor: wrong dimensions");

			auto val = values.begin();
			_belowDiag[0] = 0.0;
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
			_aboveDiag[dim - 1] = 0.0;
		}

		int RowNum() const { return _dim; }
		int ColNum() const { return _dim; }

		Type  operator()(int i, int j) const {
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _aboveDiag[i];
			else if (i == j + 1 && j < _dim - 1)
				return _belowDiag[i];
			else
				return 0.0;
		}
		Type& operator()(int i, int j) {
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _aboveDiag[i];
			else if (i == j + 1 && j < _dim - 1)
				return _belowDiag[i];
			else
				throw MatrixAccessBoundsError("TridiagonalMatrix::operator()", i, j, _dim, _dim);
		}

		bool IsEqual(const TridiagonalMatrix& b, Type eps = Defaults::MatrixIsEqualTolerance) const
		{
			if (_dim != b._dim)
				return false;

			for (int i = 0; i < _dim; i++) {
				if (std::abs(_diag[i] - b._diag[i]) > eps)
					return false;
				if (std::abs(_belowDiag[i] - b._belowDiag[i]) > eps)
					return false;
				if (std::abs(_aboveDiag[i] - b._aboveDiag[i]) > eps)
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
			int	i, j;
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
			int	i, j;
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
			int	i, j;
			TridiagonalMatrix	ret(a);

			for (int i = 0; i < ret._belowDiag.size(); i++)
				ret._belowDiag[i] /= b;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] /= b;
			for (int i = 0; i < ret._aboveDiag.size(); i++)
				ret._aboveDiag[i] /= b;

			return ret;
		}

		TridiagonalMatrix GetTranspose() const
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

		TridiagonalMatrix GetInverse() const
		{
			// For tridiagonal matrices, computing the exact inverse is complex
			// and often results in a full matrix (not tridiagonal)
			// This would require returning a full Matrix<Type> instead
			throw std::runtime_error("TridiagonalMatrix::GetInverse() - "
				"inverse of tridiagonal matrix is generally not tridiagonal. "
				"Use Matrix<Type> conversion or solve linear systems directly with Solve()");
		}

		void Solve(const Vector<Type>& rhs, Vector<Type>& sol)
		{
			int j, n = _belowDiag.size();
			Real bet;
			Vector<Type> gam(n);

			if (_diag[0] == 0.0)
				throw("Error 1 in tridag");

			sol[0] = rhs[0] / (bet = _diag[0]);

			for (j = 1; j < n; j++) {
				gam[j] = _aboveDiag[j - 1] / bet;
				bet = _diag[j] - _belowDiag[j] * gam[j];

				if (bet == 0.0)
					throw("Error 2 in tridag");

				sol[j] = (rhs[j] - _belowDiag[j] * sol[j - 1]) / bet;
			}
			for (j = (n - 2); j >= 0; j--)
				sol[j] -= gam[j + 1] * sol[j + 1];
		}
		Vector<Type> Solve(const Vector<Type>& rhs)
		{
			Vector<Type> sol(rhs.size());
			Solve(rhs, sol);
			return sol;
		}

		// New formatted print with MatrixPrintFormat
		void Print(std::ostream& stream, const MatrixPrintFormat& fmt = MatrixPrintFormat::Default()) const
		{
			// Show header if requested
			if (fmt.showHeader) {
				stream << "Rows: " << _dim << " Cols: " << _dim << " (Tridiagonal)" << std::endl;
			}

			// Set formatting flags
			std::ios_base::fmtflags oldFlags = stream.flags();
			if (fmt.scientific)
				stream << std::scientific;
			else if (fmt.fixed)
				stream << std::fixed;

			// Compact mode - print on one line for small matrices
			if (fmt.compactMode && _dim <= 3) {
				if (fmt.showBrackets)
					stream << "[";
				for (int i = 0; i < _dim; i++) {
					if (i > 0) stream << "; ";
					for (int j = 0; j < _dim; j++) {
						if (j > 0) stream << fmt.delimiter;
						stream << std::setw(fmt.width) << std::setprecision(fmt.precision) << (*this)(i, j);
					}
				}
				if (fmt.showBrackets)
					stream << "]";
				stream << std::endl;
			}
			else {
				// Normal multi-line mode
				for (int i = 0; i < _dim; i++)
				{
					if (fmt.showBrackets)
						stream << "[ ";
					for (int j = 0; j < _dim; j++)
					{
						stream << std::setw(fmt.width) << std::setprecision(fmt.precision) << (*this)(i, j);
						if (j < _dim - 1)
							stream << fmt.delimiter;
					}
					if (fmt.showBrackets)
						stream << " ]";
					if (i < _dim - 1)
						stream << std::endl;
				}
			}

			// Restore original flags
			stream.flags(oldFlags);
		}

		// Legacy print method (backward compatibility)
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