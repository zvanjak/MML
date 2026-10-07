///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/BVPDiagnostics1D.h                                   ///
///  Description: Diagnostics and sampling helpers for one-dimensional BVP results    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BVP_DIAGNOSTICS_1D_H
#define MML_BVP_DIAGNOSTICS_1D_H

#include <mml/core/FunctionSpaces/DenseBVPSolver1D.h>

namespace MML::FunctionSpaces
{
	inline Vector<Real> SampleExpansion(const FunctionExpansion1D& expansion, const Vector<Real>& points)
	{
		Vector<Real> values(points.size());
		for (int index = 0; index < points.size(); ++index)
			values[index] = expansion.evaluate(points[index]);
		return values;
	}

	inline Vector<Real> UniformSamplePoints(const FunctionSpace1D& space, int pointCount)
	{
		if (pointCount <= 0)
			throw FunctionSpaceInputError("UniformSamplePoints: point count must be positive");

		Vector<Real> points(pointCount);
		if (pointCount == 1) {
			points[0] = (space.domainMin() + space.domainMax()) / REAL(2.0);
			return points;
		}

		for (int index = 0; index < pointCount; ++index)
			points[index] = space.domainMin() + space.domainLength() * static_cast<Real>(index) / static_cast<Real>(pointCount - 1);
		return points;
	}

	inline Real CollocationResidualMaxNorm(const FunctionExpansion1D& expansion,
		const LinearDifferentialOperator1D& differentialOperator,
		const IRealFunction& rhsFunction,
		const TrialSpace1D& trialSpace)
	{
		Matrix<Real> operatorMatrix = BuildCollocationMatrix(differentialOperator, trialSpace);
		Real maxResidual = REAL(0.0);
		for (int row = 0; row < trialSpace.dimension(); ++row) {
			Real applied = REAL(0.0);
			for (int column = 0; column < trialSpace.dimension(); ++column)
				applied += operatorMatrix(row, column) * expansion.coefficients()[column];
			Real residual = applied - rhsFunction(trialSpace.node(row));
			maxResidual = std::max(maxResidual, std::abs(residual));
		}
		return maxResidual;
	}
}

#endif // MML_BVP_DIAGNOSTICS_1D_H
