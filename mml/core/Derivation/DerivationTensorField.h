///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationTensorField.h                                             ///
///  Description: Derivatives of tensor fields                                        ///
///               Covariant derivatives, Lie derivatives, tensor calculus             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_TENSOR_FIELD_H
#define MML_DERIVATION_TENSOR_FIELD_H

#include <mml/MMLBase.h>

#include <mml/interfaces/IFunction.h>
#include <mml/interfaces/ITensorField.h>

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>

#include "DerivationBase.h"
#include "FirstDerivativeStencil.h"

namespace MML
{
	namespace Derivation
	{
		namespace Detail
		{
			template<FirstDerivativeOrder Order, int N, typename Component>
			Real EvaluateTensorComponentPartial(Component&& component, int deriv_index,
			                                    const VectorN<Real, N>& point, Real h, Real* error)
			{
				auto result = EvaluateScalarPartialFirstDerivativeStencil<Order>(
					[&](int offset) { auto x = point; x[deriv_index] += offset * h; return component(x); }, h, error != nullptr);
				if (error) *error = result.error;
				return result.value;
			}
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer1Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::One>(
				[&](const auto& x) { return f.Component(i, j, x); }, deriv_index, point, h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, deriv_index, point, ScaleStep(NDer1_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, k, deriv_index, point, ScaleStep(NDer1_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::One>(
				[&](const auto& x) { return f.Component(i, j, k, x); }, deriv_index, point, h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, k, l, deriv_index, point, ScaleStep(NDer1_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::One>(
				[&](const auto& x) { return f.Component(i, j, k, l, x); }, deriv_index, point, h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer2Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::Two>(
				[&](const auto& x) { return f.Component(i, j, x); }, deriv_index, point, h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, deriv_index, point, ScaleStep(NDer2_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, k, deriv_index, point, ScaleStep(NDer2_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::Two>(
				[&](const auto& x) { return f.Component(i, j, k, x); }, deriv_index, point, h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, k, l, deriv_index, point, ScaleStep(NDer2_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::Two>(
				[&](const auto& x) { return f.Component(i, j, k, l, x); }, deriv_index, point, h, error);
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer4Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::Four>(
				[&](const auto& x) { return f.Component(i, j, x); }, deriv_index, point, h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, deriv_index, point, ScaleStep(NDer4_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, k, deriv_index, point, ScaleStep(NDer4_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::Four>(
				[&](const auto& x) { return f.Component(i, j, k, x); }, deriv_index, point, h, error);
		}


		template <int N>
		static Real NDer4Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, k, l, deriv_index, point, ScaleStep(NDer4_h, point[deriv_index]), error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			return Detail::EvaluateTensorComponentPartial<Detail::FirstDerivativeOrder::Four>(
				[&](const auto& x) { return f.Component(i, j, k, l, x); }, deriv_index, point, h, error);
		}
	}
}

#endif // MML_DERIVATION_TENSOR_FIELD_H