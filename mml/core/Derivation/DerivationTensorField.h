///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationTensorField.h                                             ///
///  Description: Derivatives of tensor fields                                        ///
///               Covariant derivatives, Lie derivatives, tensor calculus             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_TENSOR_FIELD_H
#define MML_DERIVATION_TENSOR_FIELD_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"
#include "interfaces/ITensorField.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

#include "DerivationBase.h"

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer1Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
			}
			return diff / h;
		}

		template <int N>
		static Real NDer1Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, k, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, k, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, k, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, k, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
			}
			return diff / h;
		}

		template <int N>
		static Real NDer1Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, k, l, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, k, l, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, k, l, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, k, l, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
			}
			return diff / h;
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer2Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, x);

				*error = Constants::Eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, k, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, k, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, k, x);

				*error = Constants::Eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, k, l, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, l, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, k, l, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, k, l, x);

				*error = Constants::Eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer4Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, k, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, k, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, k, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, k, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}


		template <int N>
		static Real NDer4Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, k, l, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, k, l, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, k, l, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, k, l, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
	}
}

#endif // MML_DERIVATION_TENSOR_FIELD_H