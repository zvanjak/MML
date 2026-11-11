///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/InducedMetric.h                                ///
///  Description: Induced metric bridge from Cartesian surfaces to metric fields      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_INDUCED_METRIC_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_INDUCED_METRIC_H

#include <mml/MMLBase.h>

#include <mml/core/MetricTensor.h>
#include <mml/core/Surfaces.h>

namespace MML::DifferentialGeometry
{
	class InducedMetric2D : public MetricTensorField<2>
	{
		const Surfaces::ISurfaceCartesian& _surface;

	public:
		explicit InducedMetric2D(const Surfaces::ISurfaceCartesian& surface)
			: MetricTensorField<2>(0, 2), _surface(surface)
		{
		}

		const Surfaces::ISurfaceCartesian& surface() const noexcept { return _surface; }

		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, 2>& pos) const override {
			ValidateChristoffelArguments(i, j, k, pos);
			if (typeid(*this) == typeid(InducedMetric2D))
				if (const auto* provider = dynamic_cast<const IMetricJetProvider<2>*>(&_surface))
					return ConnectionFromJet(provider->EvaluateMetricJet(pos), i, j, k);
			return MetricTensorField<2>::GetChristoffelSymbolSecondKind(i, j, k, pos);
		}

		Real component(int i, int j, const VectorN<Real, 2>& pos) const override
		{
			Real E, F, G;
			_surface.GetFirstNormalFormCoefficients(pos[0], pos[1], E, F, G);

			if (i == 0 && j == 0)
				return E;
			if ((i == 0 && j == 1) || (i == 1 && j == 0))
				return F;
			if (i == 1 && j == 1)
				return G;

			throw IndexError("InducedMetric2D::Component - index out of range");
		}

	protected:
		Real DeriveChristoffelSymbolSecondKind(int i, int j, int k, int axis, const VectorN<Real, 2>& pos) const override {
			ValidateChristoffelArguments(i, j, k, pos);
			if (axis < 0 || axis >= 2)
				throw IndexError("Connection derivative index outside metric dimension");
			if (typeid(*this) == typeid(InducedMetric2D))
				if (const auto* provider = dynamic_cast<const IMetricJetProvider<2>*>(&_surface))
					return ConnectionFromJet(provider->EvaluateMetricJet(pos), i, j, k, axis);
			return MetricTensorField<2>::DeriveChristoffelSymbolSecondKind(i, j, k, axis, pos);
		}
	};

	inline MatrixNM<Real, 2, 2> FirstFundamentalFormMatrix(const Surfaces::ISurfaceCartesian& surface, Real u, Real w)
	{
		Real E, F, G;
		surface.GetFirstNormalFormCoefficients(u, w, E, F, G);

		MatrixNM<Real, 2, 2> metric;
		metric(0, 0) = E;
		metric(0, 1) = F;
		metric(1, 0) = F;
		metric(1, 1) = G;
		return metric;
	}

	inline Real GaussianCurvatureIntrinsic(const MetricTensorField<2>& metric, const VectorN<Real, 2>& pos)
	{
		return static_cast<MML::Real>(0.5) * metric.GetRicciScalar(pos);
	}

	inline Real GaussianCurvatureIntrinsic(const Surfaces::ISurfaceCartesian& surface, Real u, Real w)
	{
		InducedMetric2D metric(surface);
		return GaussianCurvatureIntrinsic(metric, VectorN<Real, 2>{u, w});
	}

} // namespace MML::DifferentialGeometry

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_INDUCED_METRIC_H
