///////////////////////////////////////////////////////////////////////////////////////////
///  MinimalMathLibrary - internal scaled Lyapunov-frame orthonormalization
///  Copyright (c) 2024-2026 Zvonimir Vanjak. MIT License (see LICENSE.md).
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LYAPUNOV_QR_H
#define MML_LYAPUNOV_QR_H

#include <mml/systems/DynamicalSystem/DynamicalSystemValidation.h>
#include <algorithm>
#include <cmath>
#include <string>

namespace MML::Systems::LyapunovDetail
{
	// Benettin QR accumulation, with column scaling and twice-applied modified
	// Gram-Schmidt (DGKS reorthogonalization) to retain small independent stretches.
	inline Vector<Real> Orthonormalize(Matrix<Real>& frame)
	{
		const int dimension = frame.rows();
		Vector<Real> logarithms(dimension);
		for (int column = 0; column < dimension; ++column) {
			Real scale = 0;
			for (int row = 0; row < dimension; ++row) {
				if (!std::isfinite(frame(row, column)))
					throw NumericalMethodError("Lyapunov: non-finite tangent frame");
				scale = std::max(scale, std::abs(frame(row, column)));
			}
			if (scale == 0)
				throw SingularMatrixError("Lyapunov: zero tangent stretch; full spectrum is unavailable");

			Vector<Real> vector(dimension);
			StableNormL2Accumulator originalNorm;
			for (int row = 0; row < dimension; ++row) {
				vector[row] = frame(row, column) / scale;
				originalNorm.add(std::abs(vector[row]));
			}
			for (int pass = 0; pass < 2; ++pass)
				for (int previous = 0; previous < column; ++previous) {
					Real projection = 0;
					for (int row = 0; row < dimension; ++row)
						projection += vector[row] * frame(row, previous);
					for (int row = 0; row < dimension; ++row)
						vector[row] -= projection * frame(row, previous);
				}

			StableNormL2Accumulator residualNorm;
			for (Real value : vector)
				residualNorm.add(std::abs(value));
			const Real norm = residualNorm.value();
			if (!std::isfinite(norm))
				throw NumericalMethodError("Lyapunov: non-finite QR residual");
			if (IsNegligibleRelativeTo(norm, originalNorm.value(), Precision::Eps * dimension))
				throw SingularMatrixError("Lyapunov: numerically dependent tangent columns; full spectrum is unavailable");

			logarithms[column] = std::log(scale) + std::log(norm);
			for (int row = 0; row < dimension; ++row)
				frame(row, column) = vector[row] / norm;
		}
		return logarithms;
	}
}
#endif // MML_LYAPUNOV_QR_H
