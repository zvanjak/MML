///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLPrecision.h                                                      ///
///  Description: Precision constants and tolerances for float/double types           ///
///               Configurable thresholds for numerical comparisons                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PRECISION_H
#define MML_PRECISION_H

#include "MMLTypeDefs.h"
#include <limits>

namespace MML
{
	// True when long double genuinely has more precision than double.
	// GCC/Linux: long double is 80-bit (64 mantissa bits) vs double's 53 -> true
	// MSVC/Windows: long double == double (both 53 mantissa bits) -> false
	constexpr bool LongDoubleHasExtraPrecision =
		std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits;

	// Template struct for precision values
	template<typename T>
	struct PrecisionValues;

	// Specialization for float
	template<>
	struct PrecisionValues<float>
	{
		static constexpr float GeometryEpsilon = 1e-6f;               // Floating-point comparison epsilon for geometry

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

		// Angle comparison tolerance (for wrap-aware angle equality)
		static constexpr float AngleIsEqualTolerance = 1e-6f;

		// Shape property tolerance (for geometric shape classification: IsRight, IsSquare, etc.)
		static constexpr float ShapePropertyTolerance = 1e-6f;

		static constexpr float Line3DAreEqualTolerance = 1e-6f;
		static constexpr float Line3DIsPointOnLineTolerance = 1e-6f;
		static constexpr float Line3DIsPerpendicularTolerance = 1e-6f;
		static constexpr float Line3DIsParallelTolerance = 1e-6f;
		static constexpr float Line3DIntersectionTolerance = 1e-6f;

		static constexpr float Plane3DIsPointOnPlaneTolerance = 1e-4f;

		static constexpr float Triangle3DIsPointInsideTolerance = 1e-6f;
		static constexpr float Triangle3DIsRightTolerance = 1e-6f;
		static constexpr float Triangle3DIsIsoscelesTolerance = 1e-6f;
		static constexpr float Triangle3DIsEquilateralTolerance = 1e-6f;
		static constexpr float IsMatrixSymmetricTolerance = 1e-6f;
		static constexpr float IsMatrixDiagonalTolerance = 1e-6f;
		static constexpr float IsMatrixUnitTolerance = 1e-6f;
		static constexpr float IsMatrixZeroTolerance = 1e-6f;
		static constexpr float IsMatrixOrthogonalTolerance = 1e-6f;

		static constexpr float RankAlgEPS = 1e-5f;

		// Numerical computation thresholds
		static constexpr float NumericalZeroThreshold = 1e-6f;       // For checking near-zero values
		static constexpr float QuaternionZeroThreshold = 1e-6f;      // Quaternion singularity detection
		static constexpr float SurfaceNormalThreshold = 1e-6f;       // Surface normal magnitude threshold
		static constexpr float DerivativeStepSize = 5e-3f;           // Default step for numerical derivatives
		static constexpr float DefaultTolerance = 1e-6f;             // General purpose tolerance
		static constexpr float DefaultToleranceStrict = 1e-6f;       // Strict tolerance for precise geometric tests

		// Algorithm-specific thresholds
		static constexpr float EigenSolverConvergenceTolerance = 5e-6f; // Eigensolver convergence tolerance
		static constexpr float IterativeSolverTolerance = 5e-6f;      // Iterative linear solver convergence tolerance
		static constexpr float EigenSolverZeroThreshold = 1e-6f;     // Eigenvalue solver zero detection
		static constexpr float PolynomialCoeffZeroThreshold = 1e-5f; // Polynomial coefficient zero check
		static constexpr float MatrixElementZeroThreshold = 1e-6f;   // Matrix element near-zero check
		static constexpr float DeterminantZeroThreshold = 1e-6f;     // Determinant zero detection
		static constexpr float DivisionSafetyThreshold = 1e-20f;     // Safe division threshold
		static constexpr float OrthogonalityTolerance = 1e-5f;       // Orthogonal vector check tolerance
		static constexpr float LinearDependenceTolerance = 1e-5f;    // Linear dependence check tolerance

		// Optimization algorithm thresholds
		static constexpr float OptimizationTolerance = 3e-4f;       // Brent/GoldenSection convergence tolerance (≈ sqrt(eps))
		static constexpr float OptimizationGradientTolerance = 1e-4f; // Gradient-based optimizer convergence tolerance

		// Integration thresholds
		static constexpr float IntegrationTolerance = 1e-5f;        // Adaptive integration convergence tolerance

		// ODE solver thresholds
		static constexpr float ODEDefaultTolerance = 1e-5f;          // ODE default error tolerance
		static constexpr float ODEMinStepSize = 1e-7f;              // ODE minimum step size
		static constexpr float NewtonTolerance = 1e-4f;             // Newton iteration convergence (stiff/BVP)
		static constexpr float BVPTolerance = 1e-5f;               // BVP shooting method tolerance
		static constexpr float JacobianDelta = 1e-4f;              // Numerical Jacobian perturbation
		static constexpr float EventTolerance = 1e-6f;             // Event detection tolerance
		static constexpr float StiffMinStepSize = 1e-6f;           // Stiff solver minimum step size
		static constexpr float ErrorNormFloor = 1e-5f;             // Error norm floor for step control
		static constexpr float InitialStepScale = 1e-5f;           // Initial step estimation scale

		// Complex analysis thresholds
		static constexpr float ComplexAnalysisTolerance = 1e-5f;     // Default tolerance for complex contour integration
		static constexpr float ComplexRootFindingTolerance = 1e-5f;  // Default tolerance for complex root finding

		// Complex-step derivative step sizes
		static constexpr float ComplexStepH = 1e-20f;               // Step for 1st-order complex-step derivative
		static constexpr float ComplexStepH2 = 1e-10f;              // Step for 2nd-order complex-step derivative
	};

	// Specialization for double
	template<>
	struct PrecisionValues<double>
	{
		static constexpr double GeometryEpsilon = 1e-10;                // Floating-point comparison epsilon for geometry

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

		// Angle comparison tolerance (for wrap-aware angle equality)
		static constexpr double AngleIsEqualTolerance = 1e-10;

		// Shape property tolerance (for geometric shape classification: IsRight, IsSquare, etc.)
		static constexpr double ShapePropertyTolerance = 1e-10;

		static constexpr double Line3DAreEqualTolerance = 1e-10;
		static constexpr double Line3DIsPointOnLineTolerance = 1e-10;
		static constexpr double Line3DIsPerpendicularTolerance = 1e-10;
		static constexpr double Line3DIsParallelTolerance = 1e-10;
		static constexpr double Line3DIntersectionTolerance = 1e-10;

		static constexpr double Plane3DIsPointOnPlaneTolerance = 1e-10;

		static constexpr double Triangle3DIsPointInsideTolerance = 1e-10;
		static constexpr double Triangle3DIsRightTolerance = 1e-10;
		static constexpr double Triangle3DIsIsoscelesTolerance = 1e-10;
		static constexpr double Triangle3DIsEquilateralTolerance = 1e-10;
		static constexpr double IsMatrixSymmetricTolerance = 1e-10;
		static constexpr double IsMatrixDiagonalTolerance = 1e-10;
		static constexpr double IsMatrixUnitTolerance = 1e-10;
		static constexpr double IsMatrixZeroTolerance = 1e-10;
		static constexpr double IsMatrixOrthogonalTolerance = 1e-10;

		static constexpr double RankAlgEPS = 1e-12;

		// Numerical computation thresholds
		static constexpr double NumericalZeroThreshold = 1e-12;      // For checking near-zero values
		static constexpr double QuaternionZeroThreshold = 1e-12;     // Quaternion singularity detection
		static constexpr double SurfaceNormalThreshold = 1e-12;      // Surface normal magnitude threshold
		static constexpr double DerivativeStepSize = 1e-6;           // Default step for numerical derivatives
		static constexpr double DefaultTolerance = 1e-6;             // General purpose tolerance
		static constexpr double DefaultToleranceStrict = 1e-8;        // Strict tolerance for precise geometric tests

		// Algorithm-specific thresholds
		static constexpr double EigenSolverConvergenceTolerance = 1e-10; // Eigensolver convergence tolerance
		static constexpr double IterativeSolverTolerance = 1e-10;     // Iterative linear solver convergence tolerance
		static constexpr double EigenSolverZeroThreshold = 1e-15;    // Eigenvalue solver zero detection
		static constexpr double PolynomialCoeffZeroThreshold = 1e-12; // Polynomial coefficient zero check
		static constexpr double MatrixElementZeroThreshold = 1e-15;  // Matrix element near-zero check
		static constexpr double DeterminantZeroThreshold = 1e-15;    // Determinant zero detection
		static constexpr double DivisionSafetyThreshold = 1e-30;     // Safe division threshold
		static constexpr double OrthogonalityTolerance = 1e-10;      // Orthogonal vector check tolerance
		static constexpr double LinearDependenceTolerance = 1e-12;   // Linear dependence check tolerance

		// Optimization algorithm thresholds
		static constexpr double OptimizationTolerance = 3e-8;       // Brent/GoldenSection convergence tolerance (≈ sqrt(eps))
		static constexpr double OptimizationGradientTolerance = 1e-8; // Gradient-based optimizer convergence tolerance

		// Integration thresholds
		static constexpr double IntegrationTolerance = 1e-10;       // Adaptive integration convergence tolerance

		// ODE solver thresholds
		static constexpr double ODEDefaultTolerance = 1e-10;         // ODE default error tolerance
		static constexpr double ODEMinStepSize = 1e-15;             // ODE minimum step size
		static constexpr double NewtonTolerance = 1e-8;             // Newton iteration convergence (stiff/BVP)
		static constexpr double BVPTolerance = 1e-10;              // BVP shooting method tolerance
		static constexpr double JacobianDelta = 1e-8;              // Numerical Jacobian perturbation
		static constexpr double EventTolerance = 1e-12;            // Event detection tolerance
		static constexpr double StiffMinStepSize = 1e-12;          // Stiff solver minimum step size
		static constexpr double ErrorNormFloor = 1e-10;            // Error norm floor for step control
		static constexpr double InitialStepScale = 1e-10;          // Initial step estimation scale

		// Complex analysis thresholds
		static constexpr double ComplexAnalysisTolerance = 1e-10;    // Default tolerance for complex contour integration
		static constexpr double ComplexRootFindingTolerance = 1e-10; // Default tolerance for complex root finding

		// Complex-step derivative step sizes
		static constexpr double ComplexStepH = 1e-100;              // Step for 1st-order complex-step derivative
		static constexpr double ComplexStepH2 = 1e-50;              // Step for 2nd-order complex-step derivative
	};

	// Specialization for long double
	// When long double has genuine extra precision (GCC/Linux 80-bit), use tighter thresholds.
	// When long double == double (MSVC), fall back to double-equivalent values.
	template<>
	struct PrecisionValues<long double>
	{
		static constexpr long double GeometryEpsilon = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;

		static constexpr long double ComplexAreEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double ComplexAreEqualAbsTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		static constexpr long double MatrixIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double VectorIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		static constexpr long double Pnt2CartIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Pnt2PolarIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Pnt3CartIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Pnt3SphIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Pnt3CylIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		static constexpr long double Vec2CartIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Vec3CartIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Vec3CartIsParallelTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		static constexpr long double Vec3SphIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		// Angle comparison tolerance (for wrap-aware angle equality)
		static constexpr long double AngleIsEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		// Shape property tolerance (for geometric shape classification: IsRight, IsSquare, etc.)
		static constexpr long double ShapePropertyTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		static constexpr long double Line3DAreEqualTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Line3DIsPointOnLineTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Line3DIsPerpendicularTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Line3DIsParallelTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Line3DIntersectionTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		static constexpr long double Plane3DIsPointOnPlaneTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		static constexpr long double Triangle3DIsPointInsideTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double Triangle3DIsRightTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double Triangle3DIsIsoscelesTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double Triangle3DIsEquilateralTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double IsMatrixSymmetricTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double IsMatrixDiagonalTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double IsMatrixUnitTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double IsMatrixZeroTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;
		static constexpr long double IsMatrixOrthogonalTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;

		static constexpr long double RankAlgEPS = LongDoubleHasExtraPrecision ? 1e-13L : 1e-12L;

		// Numerical computation thresholds
		static constexpr long double NumericalZeroThreshold = LongDoubleHasExtraPrecision ? 1e-15L : 1e-12L;
		static constexpr long double QuaternionZeroThreshold = LongDoubleHasExtraPrecision ? 1e-15L : 1e-12L;
		static constexpr long double SurfaceNormalThreshold = LongDoubleHasExtraPrecision ? 1e-15L : 1e-12L;
		static constexpr long double DerivativeStepSize = LongDoubleHasExtraPrecision ? 1e-8L : 1e-6L;
		static constexpr long double DefaultTolerance = 1e-6L;             // General purpose tolerance
		static constexpr long double DefaultToleranceStrict = 1e-8L;

		// Algorithm-specific thresholds
		static constexpr long double EigenSolverConvergenceTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double IterativeSolverTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double EigenSolverZeroThreshold = LongDoubleHasExtraPrecision ? 1e-18L : 1e-15L;
		static constexpr long double PolynomialCoeffZeroThreshold = LongDoubleHasExtraPrecision ? 1e-15L : 1e-12L;
		static constexpr long double MatrixElementZeroThreshold = LongDoubleHasExtraPrecision ? 1e-18L : 1e-15L;
		static constexpr long double DeterminantZeroThreshold = LongDoubleHasExtraPrecision ? 1e-18L : 1e-15L;
		static constexpr long double DivisionSafetyThreshold = LongDoubleHasExtraPrecision ? 1e-35L : 1e-30L;
		static constexpr long double OrthogonalityTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double LinearDependenceTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-12L;

		// Optimization algorithm thresholds
		static constexpr long double OptimizationTolerance = LongDoubleHasExtraPrecision ? 3e-10L : 3e-8L;
		static constexpr long double OptimizationGradientTolerance = LongDoubleHasExtraPrecision ? 1e-10L : 1e-8L;

		// Integration thresholds
		static constexpr long double IntegrationTolerance = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;

		// ODE solver thresholds
		static constexpr long double ODEDefaultTolerance = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;
		static constexpr long double ODEMinStepSize = LongDoubleHasExtraPrecision ? 1e-18L : 1e-15L;
		static constexpr long double NewtonTolerance = LongDoubleHasExtraPrecision ? 1e-10L : 1e-8L;
		static constexpr long double BVPTolerance = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;
		static constexpr long double JacobianDelta = LongDoubleHasExtraPrecision ? 1e-10L : 1e-8L;
		static constexpr long double EventTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-12L;
		static constexpr long double StiffMinStepSize = LongDoubleHasExtraPrecision ? 1e-15L : 1e-12L;
		static constexpr long double ErrorNormFloor = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;
		static constexpr long double InitialStepScale = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;

		// Complex analysis thresholds
		static constexpr long double ComplexAnalysisTolerance = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;
		static constexpr long double ComplexRootFindingTolerance = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;

		// Complex-step derivative step sizes
		static constexpr long double ComplexStepH = LongDoubleHasExtraPrecision ? 1e-200L : 1e-100L;
		static constexpr long double ComplexStepH2 = LongDoubleHasExtraPrecision ? 1e-100L : 1e-50L;
	};

	// Convenience alias: PrecisionValues for the current Real type
	// Usage: MML::Precision::MatrixIsEqualTolerance instead of MML::PrecisionValues<Real>::...
	using Precision = PrecisionValues<Real>;

	// NOTE: PrecisionValues<__float128> does not exist yet.
	// __float128 support is under active development — see MMLTypeDefs.h and
	// analysis/2026-01-21/float128 support.md for the implementation roadmap.
}

#endif

