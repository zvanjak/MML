///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/Chart.h                                        ///
///  Description: Single-chart embeddings and pullback/push-forward helpers          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_CHART_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_CHART_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/base/DifferentialGeometry/Metric.h>
#include <mml/base/DifferentialGeometry/Point.h>
#include <mml/core/Derivation.h>

namespace MML::DifferentialGeometry
{
	namespace Detail
	{
		template<int Dimension, int K>
		std::array<int, K == 0 ? 1 : K> ChartUnflattenIndex(int flat)
		{
			std::array<int, K == 0 ? 1 : K> indices{};
			for (int i = K - 1; i >= 0; i--) {
				indices[i] = flat % Dimension;
				flat /= Dimension;
			}
			return indices;
		}
	}

	template<int DomainN, int AmbientN, class DomainFrame, class AmbientFrame>
	class Chart
	{
	public:
		virtual Point<Real, AmbientN, AmbientFrame> map_point(const Point<Real, DomainN, DomainFrame>& point) const = 0;
		virtual MatrixNM<Real, AmbientN, DomainN> jacobian(const Point<Real, DomainN, DomainFrame>& point) const = 0;

		TangentVector<AmbientN, AmbientFrame> push_forward(
			const TangentVector<DomainN, DomainFrame>& vector,
			const Point<Real, DomainN, DomainFrame>& atPoint) const
		{
			MatrixNM<Real, AmbientN, DomainN> jac = jacobian(atPoint);
			TangentVector<AmbientN, AmbientFrame> result;
			for (int i = 0; i < AmbientN; i++) {
				result[i] = REAL(0.0);
				for (int a = 0; a < DomainN; a++)
					result[i] += jac(i, a) * vector[a];
			}
			return result;
		}

		Covector<DomainN, DomainFrame> pull_back(
			const Covector<AmbientN, AmbientFrame>& covector,
			const Point<Real, DomainN, DomainFrame>& atPoint) const
		{
			MatrixNM<Real, AmbientN, DomainN> jac = jacobian(atPoint);
			Covector<DomainN, DomainFrame> result;
			for (int a = 0; a < DomainN; a++) {
				result[a] = REAL(0.0);
				for (int i = 0; i < AmbientN; i++)
					result[a] += covector[i] * jac(i, a);
			}
			return result;
		}

		Form1<DomainN, DomainFrame> pull_back(
			const Form1<AmbientN, AmbientFrame>& form,
			const Point<Real, DomainN, DomainFrame>& atPoint) const
		{
			return ToForm(pull_back(ToCovector(form), atPoint));
		}

		template<int K>
		DifferentialForm<Real, DomainN, K, DomainFrame> pull_back(
			const DifferentialForm<Real, AmbientN, K, AmbientFrame>& form,
			const Point<Real, DomainN, DomainFrame>& atPoint) const
			requires (K <= DomainN && K <= AmbientN)
		{
			DifferentialForm<Real, DomainN, K, DomainFrame> result;
			MatrixNM<Real, AmbientN, DomainN> jac = jacobian(atPoint);

			for (int sourceFlat = 0; sourceFlat < result.ComponentCount; sourceFlat++) {
				auto sourceIndices = Detail::ChartUnflattenIndex<DomainN, K>(sourceFlat);
				Real component = REAL(0.0);

				for (int targetFlat = 0; targetFlat < form.ComponentCount; targetFlat++) {
					auto targetIndices = Detail::ChartUnflattenIndex<AmbientN, K>(targetFlat);
					Real term = form.ComponentAt(targetIndices);
					for (int index = 0; index < K; index++)
						term *= jac(targetIndices[index], sourceIndices[index]);
					component += term;
				}

				result.ComponentAt(sourceIndices) = component;
			}

			return result;
		}

		Metric<DomainN, DomainFrame> induced_metric(const Point<Real, DomainN, DomainFrame>& atPoint) const
		{
			MatrixNM<Real, AmbientN, DomainN> jac = jacobian(atPoint);
			MatrixNM<Real, DomainN, DomainN> components;

			for (int a = 0; a < DomainN; a++) {
				for (int b = 0; b < DomainN; b++) {
					Real value = REAL(0.0);
					for (int i = 0; i < AmbientN; i++)
						value += jac(i, a) * jac(i, b);
					components(a, b) = value;
				}
			}

			return Metric<DomainN, DomainFrame>(components);
		}

		virtual ~Chart() { }
	};

	template<int DomainN, int AmbientN, class DomainFrame, class AmbientFrame>
	class FunctionChart : public Chart<DomainN, AmbientN, DomainFrame, AmbientFrame>
	{
		const IVectorFunctionNM<DomainN, AmbientN>& _embedding;
		Real _jacobianStep;

	public:
		FunctionChart(const IVectorFunctionNM<DomainN, AmbientN>& embedding, Real jacobianStep = REAL(0.0))
			: _embedding(embedding), _jacobianStep(jacobianStep)
		{
		}

		Point<Real, AmbientN, AmbientFrame> map_point(const Point<Real, DomainN, DomainFrame>& point) const override
		{
			return Point<Real, AmbientN, AmbientFrame>(_embedding(point.coordinates()));
		}

		MatrixNM<Real, AmbientN, DomainN> jacobian(const Point<Real, DomainN, DomainFrame>& point) const override
		{
			return Derivation::calcJacobian<DomainN, AmbientN>(_embedding, point.coordinates(), _jacobianStep);
		}
	};

	template<int DomainN, int AmbientN, class DomainFrame, class AmbientFrame>
	Point<Real, AmbientN, AmbientFrame> map_point(
		const Chart<DomainN, AmbientN, DomainFrame, AmbientFrame>& chart,
		const Point<Real, DomainN, DomainFrame>& point)
	{
		return chart.map_point(point);
	}

	template<int DomainN, int AmbientN, class DomainFrame, class AmbientFrame>
	TangentVector<AmbientN, AmbientFrame> push_forward(
		const Chart<DomainN, AmbientN, DomainFrame, AmbientFrame>& chart,
		const TangentVector<DomainN, DomainFrame>& vector,
		const Point<Real, DomainN, DomainFrame>& atPoint)
	{
		return chart.push_forward(vector, atPoint);
	}

	template<int DomainN, int AmbientN, class DomainFrame, class AmbientFrame>
	Covector<DomainN, DomainFrame> pull_back(
		const Chart<DomainN, AmbientN, DomainFrame, AmbientFrame>& chart,
		const Covector<AmbientN, AmbientFrame>& covector,
		const Point<Real, DomainN, DomainFrame>& atPoint)
	{
		return chart.pull_back(covector, atPoint);
	}

	template<int DomainN, int AmbientN, int K, class DomainFrame, class AmbientFrame>
	DifferentialForm<Real, DomainN, K, DomainFrame> pull_back(
		const Chart<DomainN, AmbientN, DomainFrame, AmbientFrame>& chart,
		const DifferentialForm<Real, AmbientN, K, AmbientFrame>& form,
		const Point<Real, DomainN, DomainFrame>& atPoint)
		requires (K <= DomainN && K <= AmbientN)
	{
		return chart.template pull_back<K>(form, atPoint);
	}
} // namespace MML::DifferentialGeometry

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_CHART_H