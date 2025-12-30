///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Matrix3D.h                                                          ///
///  Description: Three-dimensional matrix (tensor) class for volumetric data         ///
///               Indexed access to 3D arrays with arithmetic operations              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_MATRIX_3D_H
#define MML_MATRIX_3D_H

#include "MMLBase.h"

namespace MML
{
	template <class Type>
	class Matrix3D {
	private:
		int _n;
		int _m;
		int _k;
		Type*** _v;

	public:
		Matrix3D() : _n(0), _m(0), _k(0), _v(nullptr) {}

		Matrix3D(int n, int m, int k) : _n(n), _m(m), _k(k), _v(new Type** [n])
		{
			int i, j;

			_v[0] = new Type * [n * m];
			_v[0][0] = new Type[n * m * k];

			for (j = 1; j < m; j++)
				_v[0][j] = _v[0][j - 1] + k;

			for (i = 1; i < n; i++) {
				_v[i] = _v[i - 1] + m;
				_v[i][0] = _v[i - 1][0] + m * k;

				for (j = 1; j < m; j++)
					_v[i][j] = _v[i][j - 1] + k;
			}
		}

		~Matrix3D()
		{
			if (_v != NULL) {
				delete[](_v[0][0]);
				delete[](_v[0]);
				delete[](_v);
			}
		}

		//subscripting: pointer to row i
		inline Type** operator[](const int i) { return _v[i]; }
		inline const Type* const* operator[](const int i) const { return _v[i]; }

		Type  operator()(int i, int j, int k) const { return _v[i][j][k]; }
		Type& operator()(int i, int j, int k) { return _v[i][j][k]; }

		inline int dim1() const { return _n; }
		inline int dim2() const { return _m; }
		inline int dim3() const { return _k; }
	};
}

#endif
