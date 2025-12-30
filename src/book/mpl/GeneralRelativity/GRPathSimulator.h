#if !defined MPL_GRPATH_SIMULATOR_H
#define MPL_GRPATH_SIMULATOR_H

#include "MMLBase.h"

#include "base/Matrix.h"
#include "base/VectorN.h"
#include "base/Function.h"
#include "base/InterpolatedFunction.h"

#include "core/Derivation.h"
#include "core/MetricTensor.h"

#include "core/Integration/PathIntegration.h"

#include "tools/Visualizer.h"


using namespace MML;

namespace MPL
{
	// given a path in 4D spacetime, and ITensorField<4> representing spacetime metric
	// calculate proper time along given path

	class HelperCurveProperTimeGeneral : public IRealFunction
	{
		const MetricTensorField<4>& _metric;
		const IParametricCurve<4>& _curve;

	public:
		HelperCurveProperTimeGeneral(const IParametricCurve<4>& curve, const MetricTensorField<4>& metric) : _curve(curve), _metric(metric) {}

		Real operator()(Real t) const
		{
			auto pos			= _curve(t);
			auto vec_tang = Derivation::DeriveCurve<4>(_curve, t, nullptr);

			Tensor2<4> metricTensor = _metric(pos);			// get metric tensor at curve point

			Real intVal = -metricTensor(vec_tang, vec_tang);		// apply it to the tangent vector

			if (intVal < 0.0)
				throw std::runtime_error("HelperCurveProperTime: negative value under square root, spacelike point/vector");

			return sqrt(intVal);
		}
	};

	static Real CalcProperTime(const IParametricCurve<4>& curve, Real t1, Real t2, const MetricTensorField<4>& metric)
	{
		return IntegrateTrap(HelperCurveProperTimeGeneral(curve, metric), t1, t2);
	}
}

#endif // MPL_GRPATH_SIMULATOR_H