///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/FieldOperations.h                              ///
///  Description: Differential operators for typed fields                             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_FIELD_OPERATIONS_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_FIELD_OPERATIONS_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/Metric.h>
#include <mml/base/DifferentialGeometry/TypedFields.h>
#include <mml/core/Derivation.h>
#include <mml/core/Fields/FieldOperations.h>

namespace MML
{
	namespace Detail
	{
		template<int N, class FrameTag>
		class ScalarFieldFunctionAdapter : public IScalarFunction<N>
		{
			const IScalarField<N, FrameTag>& _field;

		public:
			explicit ScalarFieldFunctionAdapter(const IScalarField<N, FrameTag>& field) : _field(field) { }

			Real operator()(const VectorN<Real, N>& x) const override
			{
				return _field(Point<Real, N, FrameTag>(x));
			}
		};

		template<int N, class FrameTag>
		class VectorFieldFunctionAdapter : public IVectorFunction<N>
		{
			const IVectorField<N, FrameTag>& _field;

		public:
			explicit VectorFieldFunctionAdapter(const IVectorField<N, FrameTag>& field) : _field(field) { }

			VectorN<Real, N> operator()(const VectorN<Real, N>& x) const override
			{
				return _field(Point<Real, N, FrameTag>(x)).components();
			}
		};

		template<int K>
		std::array<int, K == 0 ? 1 : K> OmitFormIndex(const std::array<int, K + 1>& indices, int omitted)
		{
			std::array<int, K == 0 ? 1 : K> result{};
			for (int src = 0, dst = 0; src < K + 1; src++) {
				if (src == omitted)
					continue;
				result[dst++] = indices[src];
			}
			return result;
		}
	}

	template<int N, class FrameTag>
	class ExteriorDerivativeScalarField : public IFormField<N, 1, FrameTag>
	{
		const IScalarField<N, FrameTag>& _field;
		Real _step;

	public:
		explicit ExteriorDerivativeScalarField(const IScalarField<N, FrameTag>& field, Real step = REAL(0.0))
			: _field(field), _step(step)
		{
		}

		Form1<N, FrameTag> operator()(const Point<Real, N, FrameTag>& point) const override
		{
			Detail::ScalarFieldFunctionAdapter<N, FrameTag> adapter(_field);
			VectorN<Real, N> derivatives = _step > REAL(0.0)
				? Derivation::NDer4PartialByAll(adapter, point.coordinates(), _step)
				: Derivation::NDer4PartialByAll(adapter, point.coordinates());

			Covector<N, FrameTag> covector(derivatives);
			return ToForm(covector);
		}
	};

	template<int N, class FrameTag>
	ExteriorDerivativeScalarField<N, FrameTag> ExteriorDerivative(const IScalarField<N, FrameTag>& field, Real step = REAL(0.0))
	{
		return ExteriorDerivativeScalarField<N, FrameTag>(field, step);
	}

	template<int N, class FrameTag>
	ExteriorDerivativeScalarField<N, FrameTag> exterior_derivative(const IScalarField<N, FrameTag>& field, Real step = REAL(0.0))
	{
		return ExteriorDerivative(field, step);
	}

	template<int N, int K, class FrameTag>
	class ExteriorDerivativeFormField : public IFormField<N, K + 1, FrameTag>
	{
		static_assert(K >= 0, "Form degree cannot be negative");
		static_assert(K < N, "Exterior derivative of a top-degree form is not defined in this dimension");

		const IFormField<N, K, FrameTag>& _field;
		Real _step;

		class ComponentFunction : public IScalarFunction<N>
		{
			const IFormField<N, K, FrameTag>& _field;
			std::array<int, K == 0 ? 1 : K> _componentIndices;

		public:
			ComponentFunction(const IFormField<N, K, FrameTag>& field,
				const std::array<int, K == 0 ? 1 : K>& componentIndices)
				: _field(field), _componentIndices(componentIndices)
			{
			}

			Real operator()(const VectorN<Real, N>& x) const override
			{
				return _field(Point<Real, N, FrameTag>(x)).ComponentAt(_componentIndices);
			}
		};

		Real componentDerivative(const Point<Real, N, FrameTag>& point,
			const std::array<int, K == 0 ? 1 : K>& sourceIndices,
			int derivativeIndex) const
		{
			ComponentFunction component(_field, sourceIndices);
			if (_step > REAL(0.0))
				return Derivation::NDer4Partial(component, derivativeIndex, point.coordinates(), _step);
			return Derivation::NDer4Partial(component, derivativeIndex, point.coordinates());
		}

	public:
		explicit ExteriorDerivativeFormField(const IFormField<N, K, FrameTag>& field, Real step = REAL(0.0))
			: _field(field), _step(step)
		{
		}

		DifferentialForm<Real, N, K + 1, FrameTag> operator()(const Point<Real, N, FrameTag>& point) const override
		{
			DifferentialForm<Real, N, K + 1, FrameTag> result;
			constexpr int OutIndexStorageSize = K + 1;

			for (int flat = 0; flat < result.ComponentCount; flat++) {
				std::array<int, OutIndexStorageSize> indices{};
				int remaining = flat;
				for (int i = K; i >= 0; i--) {
					indices[i] = remaining % N;
					remaining /= N;
				}

				Real value = REAL(0.0);
				for (int omitted = 0; omitted < K + 1; omitted++) {
					auto sourceIndices = Detail::OmitFormIndex<K>(indices, omitted);
					Real term = componentDerivative(point, sourceIndices, indices[omitted]);
					value += (omitted % 2 == 0) ? term : -term;
				}
				result.ComponentAt(indices) = value;
			}

			return result;
		}
	};

	template<int N, int K, class FrameTag>
	ExteriorDerivativeFormField<N, K, FrameTag> ExteriorDerivative(const IFormField<N, K, FrameTag>& field,
		Real step = REAL(0.0))
		requires (K < N)
	{
		return ExteriorDerivativeFormField<N, K, FrameTag>(field, step);
	}

	template<int N, int K, class FrameTag>
	ExteriorDerivativeFormField<N, K, FrameTag> exterior_derivative(const IFormField<N, K, FrameTag>& field,
		Real step = REAL(0.0))
		requires (K < N)
	{
		return ExteriorDerivative(field, step);
	}

	template<int N, class FrameTag>
	Real DirectionalDerivative(const IScalarField<N, FrameTag>& field,
		const Point<Real, N, FrameTag>& point,
		const TangentVector<N, FrameTag>& direction,
		Real step = REAL(0.0))
	{
		return ExteriorDerivative(field, step)(point)(direction);
	}

	template<int N, class FrameTag>
	Real directional_derivative(const IScalarField<N, FrameTag>& field,
		const Point<Real, N, FrameTag>& point,
		const TangentVector<N, FrameTag>& direction,
		Real step = REAL(0.0))
	{
		return DirectionalDerivative(field, point, direction, step);
	}

	template<int N, class FrameTag>
	class MetricGradientField : public IVectorField<N, FrameTag>
	{
		const IScalarField<N, FrameTag>& _field;
		Metric<N, FrameTag> _metric;
		Real _step;

	public:
		MetricGradientField(const IScalarField<N, FrameTag>& field, const Metric<N, FrameTag>& metric, Real step = REAL(0.0))
			: _field(field), _metric(metric), _step(step)
		{
		}

		TangentVector<N, FrameTag> operator()(const Point<Real, N, FrameTag>& point) const override
		{
			Form1<N, FrameTag> df = ExteriorDerivative(_field, _step)(point);
			return Sharp(_metric, ToCovector(df));
		}
	};

	template<int N, class FrameTag>
	MetricGradientField<N, FrameTag> Gradient(const IScalarField<N, FrameTag>& field, const Metric<N, FrameTag>& metric, Real step = REAL(0.0))
	{
		return MetricGradientField<N, FrameTag>(field, metric, step);
	}

	template<int N, class FrameTag>
	MetricGradientField<N, FrameTag> gradient(const IScalarField<N, FrameTag>& field, const Metric<N, FrameTag>& metric, Real step = REAL(0.0))
	{
		return Gradient(field, metric, step);
	}

	template<int N, class FrameTag>
	TangentVector<N, FrameTag> GradientCart(const IScalarField<N, FrameTag>& field,
		const Point<Real, N, FrameTag>& point)
	{
		Detail::ScalarFieldFunctionAdapter<N, FrameTag> adapter(field);
		return TangentVector<N, FrameTag>(ScalarFieldOperations::GradientCart<N>(adapter, point.coordinates()));
	}

	template<int N, class FrameTag>
	TangentVector<N, FrameTag> GradientCart(const IScalarField<N, FrameTag>& field,
		const Point<Real, N, FrameTag>& point,
		int derivativeOrder)
	{
		Detail::ScalarFieldFunctionAdapter<N, FrameTag> adapter(field);
		return TangentVector<N, FrameTag>(ScalarFieldOperations::GradientCart<N>(adapter, point.coordinates(), derivativeOrder));
	}

	template<int N, class FrameTag>
	Real LaplacianCart(const IScalarField<N, FrameTag>& field, const Point<Real, N, FrameTag>& point)
	{
		Detail::ScalarFieldFunctionAdapter<N, FrameTag> adapter(field);
		return ScalarFieldOperations::LaplacianCart<N>(adapter, point.coordinates());
	}

	template<int N, class FrameTag>
	Real DivergenceCart(const IVectorField<N, FrameTag>& field, const Point<Real, N, FrameTag>& point)
	{
		Detail::VectorFieldFunctionAdapter<N, FrameTag> adapter(field);
		return VectorFieldOperations::DivCart<N>(adapter, point.coordinates());
	}

	template<class FrameTag>
	TangentVector<3, FrameTag> CurlCart(const IVectorField<3, FrameTag>& field,
		const Point<Real, 3, FrameTag>& point)
	{
		Detail::VectorFieldFunctionAdapter<3, FrameTag> adapter(field);
		return TangentVector<3, FrameTag>(VectorFieldOperations::CurlCart(adapter, point.coordinates()));
	}
}

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_FIELD_OPERATIONS_H