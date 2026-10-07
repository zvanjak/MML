///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationVectorFunction.h                                          ///
///  Description: Derivatives of vector-valued functions f:R^n->R^m                   ///
///               Jacobian matrices, divergence, curl                                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_VECTOR_FUNCTION_H
#define MML_DERIVATION_VECTOR_FUNCTION_H

#include <mml/MMLBase.h>

#include "DerivationBase.h"
#include "FirstDerivativeStencil.h"

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>

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
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::One>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x)[func_index]; }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}

		template <int N>
		static Real NDer1Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, 
														 const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, func_index, deriv_index, point, ScaleStep(NDer1_h, point[deriv_index]), error);
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
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer1Partial(f, func_index, i, point, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, func_index, i, point);
			}
			return ret;
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
			MatrixNM<Real, N, N> ret;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer1Partial(f, i, j, point, &((*error)(i, j)));
					else
						ret(i, j) = NDer1Partial(f, i, j, point);
				}
			return ret;
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer2Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, 
														 const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::Two>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x)[func_index]; }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer2Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, func_index, deriv_index, point, ScaleStep(NDer2_h, point[deriv_index]), error);
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
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer2Partial(f, func_index, i, point, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, func_index, i, point);
			}
			return ret;
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
			MatrixNM<Real, N, N> ret;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer2Partial(f, i, j, point, &((*error)(i, j)));
					else
						ret(i, j) = NDer2Partial(f, i, j, point);
				}
			return ret;
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer4Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, 
														 Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::Four>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x)[func_index]; }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer4Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, 
														 Real* error = nullptr)
		{
			return NDer4Partial(f, func_index, deriv_index, point, ScaleStep(NDer4_h, point[deriv_index]), error);
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
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer4Partial(f, func_index, i, point, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, func_index, i, point);
			}
			return ret;
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
			MatrixNM<Real, N, N> ret;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer4Partial(f, i, j, point, &((*error)(i, j)));
					else
						ret(i, j) = NDer4Partial(f, i, j, point);
				}
			return ret;
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer6Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::Six>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x)[func_index]; }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer6Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, func_index, deriv_index, point, ScaleStep(NDer6_h, point[deriv_index]), error);
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
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer6Partial(f, func_index, i, point, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, func_index, i, point);
			}
			return ret;
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
			MatrixNM<Real, N, N> ret;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer6Partial(f, i, j, point, &((*error)(i, j)));
					else
						ret(i, j) = NDer6Partial(f, i, j, point);
				}
			return ret;
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer8Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::Eight>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x)[func_index]; }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer8Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, func_index, deriv_index, point, ScaleStep(NDer8_h, point[deriv_index]), error);
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
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer8Partial(f, func_index, i, point, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, func_index, i, point);
			}
			return ret;
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
			MatrixNM<Real, N, N> ret;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer8Partial(f, i, j, point, &((*error)(i, j)));
					else
						ret(i, j) = NDer8Partial(f, i, j, point);
				}
			return ret;
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