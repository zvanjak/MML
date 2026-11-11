///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  Common configuration and result types for interpolation                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_INTERPOLATION_TYPES_H
#define MML_INTERPOLATION_TYPES_H

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Vector/Vector.h>
#include <mml/MMLExceptions.h>

#include <cmath>
#include <string>

namespace MML {
	namespace detail {
		inline void ValidateInterpolationNodes(const Vector<Real>& nodes, int minimum, const char* context) {
			if (nodes.size() < minimum)
				throw RealFuncInterpInitError(std::string(context) + ": insufficient nodes");
			for (int index = 0; index < nodes.size(); ++index)
				if (!std::isfinite(nodes[index]))
					throw RealFuncInterpInitError(std::string(context) + ": nodes must be finite");
			if (nodes.size() < 2)
				return;
			const bool ascending = nodes[1] > nodes[0];
			for (int index = 1; index < nodes.size(); ++index)
				if (ascending ? nodes[index] <= nodes[index - 1] : nodes[index] >= nodes[index - 1])
					throw RealFuncInterpInitError(std::string(context) + ": nodes must be strictly monotonic");
		}
	}

	enum class ExtrapolationPolicy {
		Allow,
		Throw,
		Clamp
	};

	enum class DuplicateHandlingPolicy {
		Reject,
		KeepFirst,
		KeepLast
	};

	enum class InterpolationStatus {
		Success,
		Extrapolated,
		Clamped,
		OutOfRange,
		SingularStencil,
		InvalidInput
	};

	struct InterpolationConfig : public EvaluationConfigBase {
		ExtrapolationPolicy extrapolation_policy = ExtrapolationPolicy::Allow;
		DuplicateHandlingPolicy duplicate_policy = DuplicateHandlingPolicy::Reject;
		Real tolerance = PrecisionValues<Real>::DefaultTolerance;
		bool validate_monotonic_x = true;
	};

	struct InterpolationResult : public EvaluationResultBase {
		Real value = 0.0;
		Real error_estimate = NoEstimate;
		Real query_point = 0.0;
		Real evaluated_point = 0.0;
		int interval_index = -1;
		InterpolationStatus interpolation_status = InterpolationStatus::Success;

		[[nodiscard]] bool WasExtrapolated() const noexcept {
			return interpolation_status == InterpolationStatus::Extrapolated;
		}

		[[nodiscard]] bool WasClamped() const noexcept {
			return interpolation_status == InterpolationStatus::Clamped;
		}
	};

} // namespace MML

#endif // MML_INTERPOLATION_TYPES_H
