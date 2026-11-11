#if !defined MML_IMPLICIT_NEWTON_MATRIX_H
#define MML_IMPLICIT_NEWTON_MATRIX_H

#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <cmath>
#include <optional>

namespace MML::Detail {

	// Step-local simplified Newton matrix (Hairer-Wanner II, IV.8).
	class ImplicitNewtonMatrix {
		std::optional<LUSolver<Real>> _lu;
		Real _previousNorm = 0;
	public:
		bool NeedsRefresh(Real norm, Real tolerance, int remainingIterations, bool reuse) {
			bool refresh = !_lu || !reuse;
			if (!refresh) {
				Real rate = norm / _previousNorm;
				// Refresh on poor contraction or predicted exhaustion of the Newton budget.
				refresh = !std::isfinite(rate) || rate >= Precision::NewtonRefreshRatio ||
					norm * std::pow(rate, remainingIterations) >= tolerance;
			}
			_previousNorm = norm;
			return refresh;
		}

		void Factor(const Matrix<Real>& matrix) {
			for (int i = 0; i < matrix.rows(); ++i)
				for (int j = 0; j < matrix.cols(); ++j)
					if (!std::isfinite(matrix(i, j)))
						throw MatrixNumericalError("Non-finite implicit Newton iteration matrix");
			_lu.emplace(matrix);
		}

		Vector<Real> Solve(const Vector<Real>& rhs) {
			return _lu.value().Solve(rhs);
		}
	};
}

#endif // MML_IMPLICIT_NEWTON_MATRIX_H
