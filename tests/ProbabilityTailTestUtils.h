#ifndef MML_PROBABILITY_TAIL_TEST_UTILS_H
#define MML_PROBABILITY_TAIL_TEST_UTILS_H

#include "TestPrecision.h"
#include <algorithm>
#include <cmath>

namespace MML::Testing {
	inline Real GammaLogRoundoff(Real a, Real x) {
		// Roundoff in a*log(x)-x-lgamma(a), including the final log-tail addition.
		return 4 * Precision::Eps * (1 + std::abs(a * std::log(x)) + x + std::abs(std::lgamma(a)));
	}

	inline Real GammaTailTolerance(Real a, Real x) {
		return Tol(Precision::EqualityTolerance,
			std::max(Precision::EqualityTolerance, GammaLogRoundoff(a, x)));
	}
}

#endif
