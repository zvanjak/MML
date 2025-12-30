///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationVectorFunction.h                                          ///
///  Description: Derivatives of vector-valued functions f:R^n->R^m                   ///
///               Jacobian matrices, divergence, curl                                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_VECTOR_FUNCTION_H
#define MML_DERIVATION_VECTOR_FUNCTION_H

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer1Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, 
														 const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f(x)[func_index];

			x[deriv_index] = x_orig + h;
			Real yh = f(x)[func_index];

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f(x)[func_index];
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
			}
			return diff / h;
		}

		template <int N>
		static Real NDer1Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, 
														 const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, func_index, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N>& f, int func_index, 
																						  const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer1Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N>& f, int func_index, 
																							const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer1PartialByAll(f, func_index, point, NDer1_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer1PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, 
																										 Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer1Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer1Partial(f, i, j, point, h);
				}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer1PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, 
																										 MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer1PartialAllByAll(f, point, NDer1_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer2Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, 
														 const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f(x)[func_index];

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f(x)[func_index];

				*error = Constants::Eps * (std::abs(yh) +  std::abs(ymh)) / (2 * h) + 
								 std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		template <int N>
		static Real NDer2Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, func_index, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer2Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer2PartialByAll(f, func_index, point, NDer2_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer2PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer2Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer2Partial(f, i, j, point, h);
				}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer2PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer2PartialAllByAll(f, point, NDer2_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer4Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, 
														 Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f(x)[func_index];

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f(x)[func_index];

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		template <int N>
		static Real NDer4Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, 
														 Real* error = nullptr)
		{
			return NDer4Partial(f, func_index, deriv_index, point, NDer4_h, error);
		}
		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, 
																							VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer4Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, func_index, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, 
																							VectorN<Real, N>* error = nullptr)
		{
			return NDer4PartialByAll(f, func_index, point, NDer4_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer4PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, 
																										 MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer4Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer4Partial(f, i, j, point, h);
				}

			return ret;
		}
		template <int N>
		static MatrixNM<Real, N, N> NDer4PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, 
																										 MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer4PartialAllByAll(f, point, NDer4_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer6Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x)[func_index];

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x)[func_index];

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;

			if (error)
			{
				x[deriv_index] = orig_x + 4 * h;
				Real y4h = f(x)[func_index];

				x[deriv_index] = orig_x - 4 * h;
				Real ym4h = f(x)[func_index];

				Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		template <int N>
		static Real NDer6Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, func_index, deriv_index, point, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer6Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, func_index, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer6PartialByAll(f, func_index, point, NDer6_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer6PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer6Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer6Partial(f, i, j, point, h);
				}

			return ret;
		}
		template <int N>
		static MatrixNM<Real, N, N> NDer6PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer6PartialAllByAll(f, point, NDer6_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer8Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x)[func_index];

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x)[func_index];

			x[deriv_index] = orig_x + 4 * h;
			Real y4h = f(x)[func_index];

			x[deriv_index] = orig_x - 4 * h;
			Real ym4h = f(x)[func_index];

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;
			Real y4 = ym4h - y4h;

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				x[deriv_index] = orig_x + 5 * h;
				Real y5h = f(x)[func_index];

				x[deriv_index] = orig_x - 5 * h;
				Real ym5h = f(x)[func_index];

				Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}
		template <int N>
		static Real NDer8Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, func_index, deriv_index, point, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer8Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer8PartialByAll(f, func_index, point, NDer8_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer8PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer8Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer8Partial(f, i, j, point, h);
				}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer8PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer8PartialAllByAll(f, point, NDer8_h, error);
		}

		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		template<int N>
		static inline Real(*DeriveVecPartial)(const IVectorFunction<N>& f, int func_index, int deriv_index, 
																					const VectorN<Real, N>& point, Real* error) = Derivation::NDer4Partial;
		template<int N>
		static inline VectorN<Real, N>(*DeriveVecPartialAll)(const IVectorFunction<N>& f, int func_index, 
																												 const VectorN<Real, N>& point, 
																												 VectorN<Real, N>* error) = Derivation::NDer4PartialByAll;
		template<int N>
		static inline MatrixNM<Real, N, N>(*DeriveVecPartialAllByAll)(const IVectorFunction<N>& f, const VectorN<Real, N>& point, 
																																	MatrixNM<Real, N, N>* error) = Derivation::NDer4PartialAllByAll;

	}
}
#endif // MML_DERIVATION_VECTOR_FUNCTION_H