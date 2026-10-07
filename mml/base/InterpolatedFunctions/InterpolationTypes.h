///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  Common configuration and result types for interpolation                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_INTERPOLATION_TYPES_H
#define MML_INTERPOLATION_TYPES_H

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>

namespace MML {

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
		Real error_estimate = 0.0;
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
