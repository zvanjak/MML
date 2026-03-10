#ifndef MML_INTEGRATION_BASE_H
#define MML_INTEGRATION_BASE_H

#include "MMLBase.h"

namespace MML
{
	/// Result structure for numerical integration algorithms
	/// Provides convergence status and diagnostics
	/// @note For production code, check error_estimate and converged fields!
	struct IntegrationResult {
		Real value;					 ///< Computed integral value
		Real error_estimate; ///< Estimated absolute error
		int iterations;			 ///< Number of iterations/refinements performed
		bool converged;			 ///< True if convergence criteria met

		/// Implicit conversion to Real for convenience
		/// @note Discards error_estimate, iterations, and converged fields
		operator Real() const { return value; }

		/// Constructor for easy initialization
		IntegrationResult(Real val = 0.0, Real err = 0.0, int iter = 0, bool conv = true)
				: value(val)
				, error_estimate(err)
				, iterations(iter)
				, converged(conv) {}
	};
}

#endif  // MML_INTEGRATION_BASE_H
