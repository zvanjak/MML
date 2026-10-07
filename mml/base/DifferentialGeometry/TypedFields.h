///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/TypedFields.h                                  ///
///  Description: Typed scalar, vector, and differential form field interfaces        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DIFFERENTIAL_GEOMETRY_TYPED_FIELDS_H
#define MML_DIFFERENTIAL_GEOMETRY_TYPED_FIELDS_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/base/DifferentialGeometry/Point.h>
#include <mml/interfaces/IFunction.h>

namespace MML
{
	template<int N, class FrameTag>
	class IScalarField
	{
	public:
		virtual Real operator()(const Point<Real, N, FrameTag>& point) const = 0;
		virtual ~IScalarField() { }
	};

	template<int N, class FrameTag>
	class IVectorField
	{
	public:
		virtual TangentVector<N, FrameTag> operator()(const Point<Real, N, FrameTag>& point) const = 0;
		virtual ~IVectorField() { }
	};

	template<int N, int K, class FrameTag>
	class IFormField
	{
	public:
		virtual DifferentialForm<Real, N, K, FrameTag> operator()(const Point<Real, N, FrameTag>& point) const = 0;
		virtual ~IFormField() { }
	};

	template<int N, class FrameTag>
	class ScalarFunctionFieldAdapter : public IScalarField<N, FrameTag>
	{
		const IScalarFunction<N>& _func;

	public:
		explicit ScalarFunctionFieldAdapter(const IScalarFunction<N>& func) : _func(func) { }

		Real operator()(const Point<Real, N, FrameTag>& point) const override
		{
			return _func(point.coordinates());
		}
	};

	template<int N, class FrameTag>
	class VectorFunctionFieldAdapter : public IVectorField<N, FrameTag>
	{
		const IVectorFunction<N>& _func;

	public:
		explicit VectorFunctionFieldAdapter(const IVectorFunction<N>& func) : _func(func) { }

		TangentVector<N, FrameTag> operator()(const Point<Real, N, FrameTag>& point) const override
		{
			return TangentVector<N, FrameTag>(_func(point.coordinates()));
		}
	};
}

#endif // MML_DIFFERENTIAL_GEOMETRY_TYPED_FIELDS_H