///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLPrecision.h                                                      ///
///  Description: Precision constants and tolerances for float/double types           ///
///               Configurable thresholds for numerical comparisons                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PRECISION_H
#define MML_PRECISION_H

namespace MML
{
	// Template struct for precision values
	template<typename T>
	struct PrecisionValues;

	// Specialization for float
	template<>
	struct PrecisionValues<float>
	{
		static constexpr float ComplexAreEqualTolerance = 1e-6f;
		static constexpr float ComplexAreEqualAbsTolerance = 1e-6f;

		static constexpr float MatrixIsEqualTolerance = 1e-6f;
		static constexpr float VectorIsEqualTolerance = 1e-6f;

		static constexpr float Pnt2CartIsEqualTolerance = 1e-6f;
		static constexpr float Pnt2PolarIsEqualTolerance = 1e-6f;
		static constexpr float Pnt3CartIsEqualTolerance = 1e-6f;
		static constexpr float Pnt3SphIsEqualTolerance = 1e-6f;
		static constexpr float Pnt3CylIsEqualTolerance = 1e-6f;

		static constexpr float Vec2CartIsEqualTolerance = 1e-6f;
		static constexpr float Vec3CartIsEqualTolerance = 1e-6f;
		static constexpr float Vec3CartIsParallelTolerance = 1e-6f;

		static constexpr float Vec3SphIsEqualTolerance = 1e-6f;

		static constexpr float Line3DAreEqualTolerance = 1e-6f;
		static constexpr float Line3DIsPointOnLineTolerance = 1e-6f;
		static constexpr float Line3DIsPerpendicularTolerance = 1e-6f;
		static constexpr float Line3DIsParallelTolerance = 1e-6f;
		static constexpr float Line3DIntersectionTolerance = 1e-6f;

		static constexpr float Plane3DIsPointOnPlaneTolerance = 1e-6f;

		static constexpr float Triangle3DIsPointInsideTolerance = 1e-6f;	static constexpr float Triangle3DIsRightTolerance = 1e-6f;
	static constexpr float Triangle3DIsIsoscelesTolerance = 1e-6f;
	static constexpr float Triangle3DIsEquilateralTolerance = 1e-6f;
		static constexpr float IsMatrixSymmetricTolerance = 1e-6f;
		static constexpr float IsMatrixDiagonalTolerance = 1e-6f;
		static constexpr float IsMatrixUnitTolerance = 1e-6f;
		static constexpr float IsMatrixOrthogonalTolerance = 1e-6f;

		static constexpr float RankAlgEPS = 1e-10f;

		// Numerical computation thresholds
		static constexpr float NumericalZeroThreshold = 1e-12f;      // For checking near-zero values
		static constexpr float QuaternionZeroThreshold = 1e-12f;     // Quaternion singularity detection
		static constexpr float SurfaceNormalThreshold = 1e-12f;      // Surface normal magnitude threshold
		static constexpr float DerivativeStepSize = 1e-6f;           // Default step for numerical derivatives
		static constexpr float DefaultTolerance = 1e-6f;             // General purpose tolerance
		static constexpr float DefaultToleranceRelaxed = 1e-8f;      // Relaxed tolerance for geometric tests

		// Algorithm-specific thresholds
		static constexpr float EigenSolverZeroThreshold = 1e-12f;    // Eigenvalue solver zero detection
		static constexpr float PolynomialCoeffZeroThreshold = 1e-10f; // Polynomial coefficient zero check
		static constexpr float MatrixElementZeroThreshold = 1e-12f;  // Matrix element near-zero check
		static constexpr float DeterminantZeroThreshold = 1e-12f;    // Determinant zero detection
		static constexpr float DivisionSafetyThreshold = 1e-25f;     // Safe division threshold
		static constexpr float OrthogonalityTolerance = 1e-8f;       // Orthogonal vector check tolerance
		static constexpr float LinearDependenceTolerance = 1e-10f;   // Linear dependence check tolerance
	};

	// Specialization for double
	template<>
	struct PrecisionValues<double>
	{
		static constexpr double ComplexAreEqualTolerance = 1e-10;
		static constexpr double ComplexAreEqualAbsTolerance = 1e-10;

		static constexpr double MatrixIsEqualTolerance = 1e-10;
		static constexpr double VectorIsEqualTolerance = 1e-10;

		static constexpr double Pnt2CartIsEqualTolerance = 1e-10;
		static constexpr double Pnt2PolarIsEqualTolerance = 1e-10;
		static constexpr double Pnt3CartIsEqualTolerance = 1e-10;
		static constexpr double Pnt3SphIsEqualTolerance = 1e-10;
		static constexpr double Pnt3CylIsEqualTolerance = 1e-10;

		static constexpr double Vec2CartIsEqualTolerance = 1e-10;
		static constexpr double Vec3CartIsEqualTolerance = 1e-10;
		static constexpr double Vec3CartIsParallelTolerance = 1e-10;

		static constexpr double Vec3SphIsEqualTolerance = 1e-10;

		static constexpr double Line3DAreEqualTolerance = 1e-10;
		static constexpr double Line3DIsPointOnLineTolerance = 1e-10;
		static constexpr double Line3DIsPerpendicularTolerance = 1e-10;
		static constexpr double Line3DIsParallelTolerance = 1e-10;
		static constexpr double Line3DIntersectionTolerance = 1e-10;

		static constexpr double Plane3DIsPointOnPlaneTolerance = 1e-10;

		static constexpr double Triangle3DIsPointInsideTolerance = 1e-10f;	static constexpr double Triangle3DIsRightTolerance = 1e-10;
	static constexpr double Triangle3DIsIsoscelesTolerance = 1e-10;
	static constexpr double Triangle3DIsEquilateralTolerance = 1e-10;
		static constexpr double IsMatrixSymmetricTolerance = 1e-10;
		static constexpr double IsMatrixDiagonalTolerance = 1e-10;
		static constexpr double IsMatrixUnitTolerance = 1e-10;
		static constexpr double IsMatrixOrthogonalTolerance = 1e-10;

		static constexpr double RankAlgEPS = 1e-12f;

		// Numerical computation thresholds
		static constexpr double NumericalZeroThreshold = 1e-12;      // For checking near-zero values
		static constexpr double QuaternionZeroThreshold = 1e-12;     // Quaternion singularity detection
		static constexpr double SurfaceNormalThreshold = 1e-12;      // Surface normal magnitude threshold
		static constexpr double DerivativeStepSize = 1e-6;           // Default step for numerical derivatives
		static constexpr double DefaultTolerance = 1e-6;             // General purpose tolerance
		static constexpr double DefaultToleranceRelaxed = 1e-8;      // Relaxed tolerance for geometric tests

		// Algorithm-specific thresholds
		static constexpr double EigenSolverZeroThreshold = 1e-15;    // Eigenvalue solver zero detection
		static constexpr double PolynomialCoeffZeroThreshold = 1e-12; // Polynomial coefficient zero check
		static constexpr double MatrixElementZeroThreshold = 1e-15;  // Matrix element near-zero check
		static constexpr double DeterminantZeroThreshold = 1e-15;    // Determinant zero detection
		static constexpr double DivisionSafetyThreshold = 1e-30;     // Safe division threshold
		static constexpr double OrthogonalityTolerance = 1e-10;      // Orthogonal vector check tolerance
		static constexpr double LinearDependenceTolerance = 1e-12;   // Linear dependence check tolerance
	};

	// Specialization for long double
	template<>
	struct PrecisionValues<long double>
	{
		static constexpr long double ComplexAreEqualTolerance = 1e-15L;
		static constexpr long double ComplexAreEqualAbsTolerance = 1e-15L;

		static constexpr long double MatrixIsEqualTolerance = 1e-15L;
		static constexpr long double VectorIsEqualTolerance = 1e-15L;

		static constexpr long double Pnt2CartIsEqualTolerance = 1e-15;
		static constexpr long double Pnt2PolarIsEqualTolerance = 1e-15;
		static constexpr long double Pnt3CartIsEqualTolerance = 1e-15;
		static constexpr long double Pnt3SphIsEqualTolerance = 1e-15;
		static constexpr long double Pnt3CylIsEqualTolerance = 1e-15;

		static constexpr long double Vec2CartIsEqualTolerance = 1e-15;
		static constexpr long double Vec3CartIsEqualTolerance = 1e-15;
		static constexpr long double Vec3CartIsParallelTolerance = 1e-15;

		static constexpr long double Vec3SphIsEqualTolerance = 1e-15;

		static constexpr long double Line3DAreEqualTolerance = 1e-15;
		static constexpr long double Line3DIsPointOnLineTolerance = 1e-15;
		static constexpr long double Line3DIsPerpendicularTolerance = 1e-15;
		static constexpr long double Line3DIsParallelTolerance = 1e-15;
		static constexpr long double Line3DIntersectionTolerance = 1e-15;

		static constexpr long double Plane3DIsPointOnPlaneTolerance = 1e-15;

		static constexpr long double Triangle3DIsPointInsideTolerance = 1e-15f;	static constexpr long double Triangle3DIsRightTolerance = 1e-12L;
	static constexpr long double Triangle3DIsIsoscelesTolerance = 1e-12L;
	static constexpr long double Triangle3DIsEquilateralTolerance = 1e-12L;
		static constexpr long double IsMatrixSymmetricTolerance = 1e-15;
		static constexpr long double IsMatrixDiagonalTolerance = 1e-15;
		static constexpr long double IsMatrixUnitTolerance = 1e-15;
		static constexpr long double IsMatrixOrthogonalTolerance = 1e-15;

		static constexpr long double RankAlgEPS = 1e-13f;

		// Numerical computation thresholds
		static constexpr long double NumericalZeroThreshold = 1e-12L;      // For checking near-zero values
		static constexpr long double QuaternionZeroThreshold = 1e-12L;     // Quaternion singularity detection
		static constexpr long double SurfaceNormalThreshold = 1e-12L;      // Surface normal magnitude threshold
		static constexpr long double DerivativeStepSize = 1e-6L;           // Default step for numerical derivatives
		static constexpr long double DefaultTolerance = 1e-6L;             // General purpose tolerance
		static constexpr long double DefaultToleranceRelaxed = 1e-8L;      // Relaxed tolerance for geometric tests

		// Algorithm-specific thresholds
		static constexpr long double EigenSolverZeroThreshold = 1e-18L;    // Eigenvalue solver zero detection
		static constexpr long double PolynomialCoeffZeroThreshold = 1e-15L; // Polynomial coefficient zero check
		static constexpr long double MatrixElementZeroThreshold = 1e-18L;  // Matrix element near-zero check
		static constexpr long double DeterminantZeroThreshold = 1e-18L;    // Determinant zero detection
		static constexpr long double DivisionSafetyThreshold = 1e-35L;     // Safe division threshold
		static constexpr long double OrthogonalityTolerance = 1e-12L;      // Orthogonal vector check tolerance
		static constexpr long double LinearDependenceTolerance = 1e-15L;   // Linear dependence check tolerance
	};
}


#endif

