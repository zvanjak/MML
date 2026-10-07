///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/MetricAdapters.h                               ///
///  Description: Adapters from runtime metric tensor fields to typed metrics         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_METRIC_ADAPTERS_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_METRIC_ADAPTERS_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/Metric.h>
#include <mml/core/MetricTensor.h>

namespace MML
{
	template<int N, class FrameTag>
	Metric<N, FrameTag> MetricAt(const MetricTensorField<N>& metricField, const VectorN<Real, N>& pos)
	{
		return Metric<N, FrameTag>(metricField.GetCovariantMetric(pos));
	}
}

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_METRIC_ADAPTERS_H