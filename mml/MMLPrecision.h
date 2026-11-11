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

#include <mml/MMLTypeDefs.h>
#include <limits>
#include <cstddef>

namespace MML
{
	namespace ParallelDefaults {
		inline constexpr std::size_t MinTasks = 2;
		inline constexpr std::size_t MonteCarloBlockSize = 4096;
		inline constexpr std::size_t MonteCarloMinSamples = 8192;
		inline constexpr std::size_t PartitionsPerWorker = 4;
	}

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
		// ---- Tiers (docs/TOLERANCE_POLICY.md section 2). Every comparison-family constant below is one of these. ----
		static constexpr float Eps = std::numeric_limits<float>::epsilon();
		static constexpr float SqrtEps = 0.000345267f;  // eps^(1/2): optimal step for forward differences, Brent/golden-section floor
		static constexpr float CbrtEps = 0.00492157f;  // eps^(1/3): optimal step for 2nd-order central differences
		static constexpr float FourthRootEps = 0.0185814f;  // eps^(1/4): optimal step for 2nd-order second derivatives
		static constexpr float FifthRootEps = 0.0412346f;  // eps^(1/5): optimal step for 4th-order central differences (CalcJacobian)
		static constexpr float CoarseTolerance = 1e-4f;  // fast/preview solver presets and sampling heuristics over a region
		static constexpr float LooseTolerance = 1e-6f;  // exploratory / sampling checks, user-facing default tolerance
		static constexpr float ModerateTolerance = 1e-6f;  // convergence of iterative solvers on well-conditioned problems
		static constexpr float EqualityTolerance = 1e-6f;  // two computed values of the same quantity agree
		static constexpr float ZeroTolerance = 1e-6f;  // a computed quantity is negligible (IsNearlyZero default)
		static constexpr float StrictTolerance = 1e-6f;  // structural zero detection: pivots, determinants, eigenvalues
		// ---- Scaled roundoff and pragmatic limits ----
		static constexpr float RoundoffZeroThreshold = 100 * Eps;  // zero up to roundoff of O(100) flops at unit scale: norm before normalizing, integer-valuedness, trailing polynomial coefficients
		static constexpr float UnboundedIntervalSamplingLimit = 1e6f;  // finite stand-in for +-infinity when sampling an unbounded interval
		static constexpr float IntegrationRoundoffFloorScale = 50 * Eps;   // QUADPACK: convergence floor 50*eps*integral(|f|) below which a step error is roundoff
		static constexpr float KrylovBreakdownThreshold = 1000 * Eps;      // Arnoldi/Lanczos: h_(j+1,j) below this fraction of the start norm = (lucky) breakdown
		static constexpr float SVDCompletionTolerancePerDimension = 32 * Eps;   // complex SVD: a column is already spanned when its residual norm <= this * dimension
		static constexpr float MonteCarloTargetRelativeError = 1e-3f;   // default statistical target for adaptive Monte Carlo
		// ---- Algorithm constants that used to be literals in mml/algorithms and mml/systems (TOLERANCE_POLICY section 1) ----
		static constexpr float OverflowDetectionThreshold = 1e30f;   // |f(x)| beyond this is treated as overflow by the function analyzers
		static constexpr float InfiniteLimitThreshold = 1e15f;       // one-sided limit magnitude beyond this is reported as infinite
		static constexpr float DivergenceMagnitudeThreshold = 1e6f;  // |value| beyond this counts as divergent / unbounded (Lyapunov time, jump size)
		static constexpr float PenaltyParameterMax = 1e6f;           // augmented Lagrangian: cap on the penalty parameter
		static constexpr float LagrangeMultiplierBound = 1e6f;    // augmented Lagrangian: |multiplier| cap
		static constexpr float BigMPenalty = 1e6f;                      // LP big-M value when not using two-phase simplex
		static constexpr float ArmijoSufficientDecrease = 1e-4f;       // c1 of the Armijo/Wolfe sufficient-decrease condition (Nocedal-Wright 3.4)
		static constexpr float MinimumStepRatio = 1e-6f;            // adaptive step below this fraction of the nominal step counts as collapsed
		static constexpr float DAEMinStepSize = 1e-5f;                // DAE drivers: smallest admissible step
		static constexpr float RODASMinStepFactor = 0.2f;             // RODAS4 controller shrink bound
		static constexpr float RODASMaxStepFactor = 6.0f;             // RODAS4 controller growth bound
		static constexpr float RODASRejectMaxFactor = 0.5f;           // A rejected RODAS4 attempt must shrink
		static constexpr float StepControllerErrorFloor = 1e-4f;     // PI step controller: error norm floor so a zero error cannot give an infinite growth factor
		static constexpr float WellConditionedLimit = 1e2f;             // condition number below this: well conditioned
		static constexpr float IllConditionedLimit = 1e4f;               // condition number above this: ill conditioned
		static constexpr float BrentZeroEpsilon = 1e-3f * SqrtEps;             // Brent/golden-section minimization ZEPS: absolute floor protecting x = 0 (NR 10.3)

		static constexpr float GeometryEpsilon = LooseTolerance;               // Floating-point comparison epsilon for geometry

		static constexpr float ComplexAreEqualTolerance = LooseTolerance;
		static constexpr float ComplexAreEqualAbsTolerance = LooseTolerance;

		static constexpr float MatrixIsEqualTolerance = LooseTolerance;
		static constexpr float VectorIsEqualTolerance = LooseTolerance;

		static constexpr float Pnt2CartIsEqualTolerance = LooseTolerance;
		static constexpr float Pnt2PolarIsEqualTolerance = LooseTolerance;
		static constexpr float Pnt3CartIsEqualTolerance = LooseTolerance;
		static constexpr float Pnt3SphIsEqualTolerance = LooseTolerance;
		static constexpr float Pnt3CylIsEqualTolerance = LooseTolerance;

		static constexpr float Vec2CartIsEqualTolerance = LooseTolerance;
		static constexpr float Vec3CartIsEqualTolerance = LooseTolerance;
		static constexpr float Vec3CartIsParallelTolerance = LooseTolerance;

		static constexpr float Vec3SphIsEqualTolerance = LooseTolerance;

		// Angle comparison tolerance (for wrap-aware angle equality)
		static constexpr float AngleIsEqualTolerance = LooseTolerance;

		// Shape property tolerance (for geometric shape classification: IsRight, IsSquare, etc.)
		static constexpr float ShapePropertyTolerance = LooseTolerance;

		static constexpr float Line3DAreEqualTolerance = LooseTolerance;
		static constexpr float Line3DIsPointOnLineTolerance = LooseTolerance;
		static constexpr float Line3DIsPerpendicularTolerance = LooseTolerance;
		static constexpr float Line3DIsParallelTolerance = LooseTolerance;
		static constexpr float Line3DIntersectionTolerance = LooseTolerance;

		static constexpr float Plane3DIsPointOnPlaneTolerance = 1e-4f;

		static constexpr float Triangle3DIsPointInsideTolerance = LooseTolerance;
		static constexpr float Triangle3DIsRightTolerance = LooseTolerance;
		static constexpr float Triangle3DIsIsoscelesTolerance = LooseTolerance;
		static constexpr float Triangle3DIsEquilateralTolerance = LooseTolerance;
		static constexpr float IsMatrixSymmetricTolerance = LooseTolerance;
		static constexpr float IsMatrixDiagonalTolerance = LooseTolerance;
		static constexpr float IsMatrixUnitTolerance = LooseTolerance;
		static constexpr float IsMatrixZeroTolerance = LooseTolerance;
		static constexpr float IsMatrixOrthogonalTolerance = LooseTolerance;

		static constexpr float RankAlgEPS = 1e-5f;

		// Numerical computation thresholds
		static constexpr float NumericalZeroThreshold = LooseTolerance;       // For checking near-zero values
		static constexpr float QuaternionZeroThreshold = LooseTolerance;      // Quaternion singularity detection
		static constexpr float SurfaceNormalThreshold = LooseTolerance;       // Surface normal magnitude threshold
		static constexpr float DerivativeStepSize = 5e-3f;           // Default step for numerical derivatives
		static constexpr float DefaultTolerance = LooseTolerance;             // General purpose tolerance
		static constexpr float DefaultToleranceStrict = LooseTolerance;       // Strict tolerance for precise geometric tests

		// Algorithm-specific thresholds
		static constexpr float EigenSolverConvergenceTolerance = 5e-6f; // Eigensolver convergence tolerance
		static constexpr float IterativeSolverTolerance = 5e-6f;      // Iterative linear solver convergence tolerance
		static constexpr float EigenSolverZeroThreshold = LooseTolerance;     // Eigenvalue solver zero detection
		static constexpr float PolynomialCoeffZeroThreshold = 1e-5f; // Polynomial coefficient zero check
		static constexpr float MatrixElementZeroThreshold = LooseTolerance;   // Matrix element near-zero check
		static constexpr float DeterminantZeroThreshold = LooseTolerance;     // Determinant zero detection
		static constexpr float DivisionSafetyThreshold = 1e-20f;     // Safe division threshold
		static constexpr float OrthogonalityTolerance = 1e-5f;       // Orthogonal vector check tolerance
		static constexpr float LinearDependenceTolerance = 1e-5f;    // Linear dependence check tolerance

		// Optimization algorithm thresholds
		static constexpr float OptimizationTolerance = 3e-4f;       // Brent/GoldenSection convergence tolerance (≈ sqrt(eps))
		static constexpr float OptimizationGradientTolerance = 1e-4f; // Gradient-based optimizer convergence tolerance

		// Integration thresholds
		static constexpr float IntegrationTolerance = 1e-5f;        // Adaptive integration convergence tolerance
		static constexpr float TrapezoidIntegrationTolerance = 1e-4f;   // Legacy Defaults::TrapezoidIntegrationEPS
		static constexpr float SimpsonIntegrationTolerance = 1e-5f;     // Legacy Defaults::SimpsonIntegrationEPS
		static constexpr float RombergIntegrationTolerance = LooseTolerance; // Legacy Defaults::RombergIntegrationEPS
		static constexpr float LineIntegralTolerance = 1e-5f;           // Legacy Defaults::LineIntegralPrecision / WorkIntegralPrecision

		// Special-function and dynamical-systems thresholds
		static constexpr float PolynomialRootTolerance = 1e-6f;      // Polynomial root refinement tolerance
		static constexpr float AGMConvergenceTolerance = 1e-7f;      // AGM iteration convergence (elliptic integrals)
		static constexpr float LyapunovChaosThreshold = 1e-4f;       // Max Lyapunov exponent above this = chaotic
		static constexpr float FixedPointTolerance = 1e-6f;          // Fixed-point Newton convergence
		static constexpr float FixedPointUniquenessTolerance = 1e-4f; // Distance below which fixed points are duplicates

		// ODE solver thresholds
		static constexpr float ODEDefaultTolerance = 1e-5f;          // ODE default error tolerance
		static constexpr float ODEMinStepSize = 1e-7f;              // ODE minimum step size
		static constexpr float NewtonTolerance = 1e-4f;             // Newton iteration convergence (stiff/BVP)
		static constexpr float NewtonRefreshRatio = 0.5f;           // Refresh simplified Newton matrix on slow contraction
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
		// ---- Tiers (docs/TOLERANCE_POLICY.md section 2). Every comparison-family constant below is one of these. ----
		static constexpr double Eps = std::numeric_limits<double>::epsilon();
		static constexpr double SqrtEps = 1.49012e-08;  // eps^(1/2): optimal step for forward differences, Brent/golden-section floor
		static constexpr double CbrtEps = 6.05545e-06;  // eps^(1/3): optimal step for 2nd-order central differences
		static constexpr double FourthRootEps = 0.00012207;  // eps^(1/4): optimal step for 2nd-order second derivatives
		static constexpr double FifthRootEps = 0.000740096;  // eps^(1/5): optimal step for 4th-order central differences (CalcJacobian)
		static constexpr double CoarseTolerance = 1e-4;  // fast/preview solver presets and sampling heuristics over a region
		static constexpr double LooseTolerance = 1e-6;  // exploratory / sampling checks, user-facing default tolerance
		static constexpr double ModerateTolerance = 1e-8;  // convergence of iterative solvers on well-conditioned problems
		static constexpr double EqualityTolerance = 1e-10;  // two computed values of the same quantity agree
		static constexpr double ZeroTolerance = 1e-12;  // a computed quantity is negligible (IsNearlyZero default)
		static constexpr double StrictTolerance = 1e-15;  // structural zero detection: pivots, determinants, eigenvalues
		// ---- Scaled roundoff and pragmatic limits ----
		static constexpr double RoundoffZeroThreshold = 100 * Eps;  // zero up to roundoff of O(100) flops at unit scale: norm before normalizing, integer-valuedness, trailing polynomial coefficients
		static constexpr double UnboundedIntervalSamplingLimit = 1e10;  // finite stand-in for +-infinity when sampling an unbounded interval
		static constexpr double IntegrationRoundoffFloorScale = 50 * Eps;   // QUADPACK: convergence floor 50*eps*integral(|f|) below which a step error is roundoff
		static constexpr double KrylovBreakdownThreshold = 1000 * Eps;      // Arnoldi/Lanczos: h_(j+1,j) below this fraction of the start norm = (lucky) breakdown
		static constexpr double SVDCompletionTolerancePerDimension = 32 * Eps;   // complex SVD: a column is already spanned when its residual norm <= this * dimension
		static constexpr double MonteCarloTargetRelativeError = 1e-3;   // default statistical target for adaptive Monte Carlo
		// ---- Algorithm constants that used to be literals in mml/algorithms and mml/systems (TOLERANCE_POLICY section 1) ----
		static constexpr double OverflowDetectionThreshold = 1e300;   // |f(x)| beyond this is treated as overflow by the function analyzers
		static constexpr double InfiniteLimitThreshold = 1e150;       // one-sided limit magnitude beyond this is reported as infinite
		static constexpr double DivergenceMagnitudeThreshold = 1e10;  // |value| beyond this counts as divergent / unbounded (Lyapunov time, jump size)
		static constexpr double PenaltyParameterMax = 1e12;           // augmented Lagrangian: cap on the penalty parameter
		static constexpr double LagrangeMultiplierBound = 1e10;    // augmented Lagrangian: |multiplier| cap
		static constexpr double BigMPenalty = 1e9;                      // LP big-M value when not using two-phase simplex
		static constexpr double ArmijoSufficientDecrease = 1e-4;       // c1 of the Armijo/Wolfe sufficient-decrease condition (Nocedal-Wright 3.4)
		static constexpr double MinimumStepRatio = 1e-10;            // adaptive step below this fraction of the nominal step counts as collapsed
		static constexpr double DAEMinStepSize = 1e-8;                // DAE drivers: smallest admissible step
		static constexpr double RODASMinStepFactor = 0.2;             // RODAS4 controller shrink bound
		static constexpr double RODASMaxStepFactor = 6.0;             // RODAS4 controller growth bound
		static constexpr double RODASRejectMaxFactor = 0.5;           // A rejected RODAS4 attempt must shrink
		static constexpr double StepControllerErrorFloor = 1e-4;     // PI step controller: error norm floor so a zero error cannot give an infinite growth factor
		static constexpr double WellConditionedLimit = 1e4;             // condition number below this: well conditioned
		static constexpr double IllConditionedLimit = 1e8;               // condition number above this: ill conditioned
		static constexpr double BrentZeroEpsilon = 1e-3 * SqrtEps;             // Brent/golden-section minimization ZEPS: absolute floor protecting x = 0 (NR 10.3)

		static constexpr double GeometryEpsilon = EqualityTolerance;                // Floating-point comparison epsilon for geometry

		static constexpr double ComplexAreEqualTolerance = EqualityTolerance;
		static constexpr double ComplexAreEqualAbsTolerance = EqualityTolerance;

		static constexpr double MatrixIsEqualTolerance = EqualityTolerance;
		static constexpr double VectorIsEqualTolerance = EqualityTolerance;

		static constexpr double Pnt2CartIsEqualTolerance = EqualityTolerance;
		static constexpr double Pnt2PolarIsEqualTolerance = EqualityTolerance;
		static constexpr double Pnt3CartIsEqualTolerance = EqualityTolerance;
		static constexpr double Pnt3SphIsEqualTolerance = EqualityTolerance;
		static constexpr double Pnt3CylIsEqualTolerance = EqualityTolerance;

		static constexpr double Vec2CartIsEqualTolerance = EqualityTolerance;
		static constexpr double Vec3CartIsEqualTolerance = EqualityTolerance;
		static constexpr double Vec3CartIsParallelTolerance = EqualityTolerance;

		static constexpr double Vec3SphIsEqualTolerance = EqualityTolerance;

		// Angle comparison tolerance (for wrap-aware angle equality)
		static constexpr double AngleIsEqualTolerance = EqualityTolerance;

		// Shape property tolerance (for geometric shape classification: IsRight, IsSquare, etc.)
		static constexpr double ShapePropertyTolerance = EqualityTolerance;

		static constexpr double Line3DAreEqualTolerance = EqualityTolerance;
		static constexpr double Line3DIsPointOnLineTolerance = EqualityTolerance;
		static constexpr double Line3DIsPerpendicularTolerance = EqualityTolerance;
		static constexpr double Line3DIsParallelTolerance = EqualityTolerance;
		static constexpr double Line3DIntersectionTolerance = EqualityTolerance;

		static constexpr double Plane3DIsPointOnPlaneTolerance = EqualityTolerance;

		static constexpr double Triangle3DIsPointInsideTolerance = EqualityTolerance;
		static constexpr double Triangle3DIsRightTolerance = EqualityTolerance;
		static constexpr double Triangle3DIsIsoscelesTolerance = EqualityTolerance;
		static constexpr double Triangle3DIsEquilateralTolerance = EqualityTolerance;
		static constexpr double IsMatrixSymmetricTolerance = EqualityTolerance;
		static constexpr double IsMatrixDiagonalTolerance = EqualityTolerance;
		static constexpr double IsMatrixUnitTolerance = EqualityTolerance;
		static constexpr double IsMatrixZeroTolerance = EqualityTolerance;
		static constexpr double IsMatrixOrthogonalTolerance = EqualityTolerance;

		static constexpr double RankAlgEPS = 1e-12;

		// Numerical computation thresholds
		static constexpr double NumericalZeroThreshold = ZeroTolerance;      // For checking near-zero values
		static constexpr double QuaternionZeroThreshold = ZeroTolerance;     // Quaternion singularity detection
		static constexpr double SurfaceNormalThreshold = ZeroTolerance;      // Surface normal magnitude threshold
		static constexpr double DerivativeStepSize = 1e-6;           // Default step for numerical derivatives
		static constexpr double DefaultTolerance = LooseTolerance;             // General purpose tolerance
		static constexpr double DefaultToleranceStrict = ModerateTolerance;        // Strict tolerance for precise geometric tests

		// Algorithm-specific thresholds
		static constexpr double EigenSolverConvergenceTolerance = 1e-10; // Eigensolver convergence tolerance
		static constexpr double IterativeSolverTolerance = 1e-10;     // Iterative linear solver convergence tolerance
		static constexpr double EigenSolverZeroThreshold = StrictTolerance;    // Eigenvalue solver zero detection
		static constexpr double PolynomialCoeffZeroThreshold = ZeroTolerance; // Polynomial coefficient zero check
		static constexpr double MatrixElementZeroThreshold = StrictTolerance;  // Matrix element near-zero check
		static constexpr double DeterminantZeroThreshold = StrictTolerance;    // Determinant zero detection
		static constexpr double DivisionSafetyThreshold = 1e-30;     // Safe division threshold
		static constexpr double OrthogonalityTolerance = EqualityTolerance;      // Orthogonal vector check tolerance
		static constexpr double LinearDependenceTolerance = ZeroTolerance;   // Linear dependence check tolerance

		// Optimization algorithm thresholds
		static constexpr double OptimizationTolerance = 3e-8;       // Brent/GoldenSection convergence tolerance (≈ sqrt(eps))
		static constexpr double OptimizationGradientTolerance = 1e-8; // Gradient-based optimizer convergence tolerance

		// Integration thresholds
		static constexpr double IntegrationTolerance = 1e-10;       // Adaptive integration convergence tolerance
		static constexpr double TrapezoidIntegrationTolerance = 1e-4;   // Legacy Defaults::TrapezoidIntegrationEPS
		static constexpr double SimpsonIntegrationTolerance = 1e-5;     // Legacy Defaults::SimpsonIntegrationEPS
		static constexpr double RombergIntegrationTolerance = LooseTolerance; // Legacy Defaults::RombergIntegrationEPS
		static constexpr double LineIntegralTolerance = 1e-5;           // Legacy Defaults::LineIntegralPrecision / WorkIntegralPrecision

		// Special-function and dynamical-systems thresholds
		static constexpr double PolynomialRootTolerance = 1e-10;     // Polynomial root refinement tolerance
		static constexpr double AGMConvergenceTolerance = 1e-15;     // AGM iteration convergence (elliptic integrals)
		static constexpr double LyapunovChaosThreshold = 1e-6;       // Max Lyapunov exponent above this = chaotic
		static constexpr double FixedPointTolerance = 1e-10;         // Fixed-point Newton convergence
		static constexpr double FixedPointUniquenessTolerance = 1e-6; // Distance below which fixed points are duplicates

		// ODE solver thresholds
		static constexpr double ODEDefaultTolerance = EqualityTolerance;         // ODE default error tolerance
		static constexpr double ODEMinStepSize = 1e-15;             // ODE minimum step size
		static constexpr double NewtonTolerance = 1e-8;             // Newton iteration convergence (stiff/BVP)
		static constexpr double NewtonRefreshRatio = 0.5;
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
		// ---- Tiers (docs/TOLERANCE_POLICY.md section 2). Every comparison-family constant below is one of these. ----
		static constexpr long double Eps = std::numeric_limits<long double>::epsilon();
		static constexpr long double SqrtEps = LongDoubleHasExtraPrecision ? 3.29272e-10L : 1.49012e-08L;  // eps^(1/2): optimal step for forward differences, Brent/golden-section floor
		static constexpr long double CbrtEps = LongDoubleHasExtraPrecision ? 4.76837e-07L : 6.05545e-06L;  // eps^(1/3): optimal step for 2nd-order central differences
		static constexpr long double FourthRootEps = LongDoubleHasExtraPrecision ? 1.81459e-05L : 0.00012207L;  // eps^(1/4): optimal step for 2nd-order second derivatives
		static constexpr long double FifthRootEps = LongDoubleHasExtraPrecision ? 0.000161073L : 0.000740096L;  // eps^(1/5): optimal step for 4th-order central differences (CalcJacobian)
		static constexpr long double CoarseTolerance = 1e-4L;  // fast/preview solver presets and sampling heuristics over a region
		static constexpr long double LooseTolerance = 1e-6L;  // exploratory / sampling checks, user-facing default tolerance
		static constexpr long double ModerateTolerance = 1e-8L;  // convergence of iterative solvers on well-conditioned problems
		static constexpr long double EqualityTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-10L;  // two computed values of the same quantity agree
		static constexpr long double ZeroTolerance = LongDoubleHasExtraPrecision ? 1e-15L : 1e-12L;  // a computed quantity is negligible (IsNearlyZero default)
		static constexpr long double StrictTolerance = LongDoubleHasExtraPrecision ? 1e-18L : 1e-15L;  // structural zero detection: pivots, determinants, eigenvalues
		// ---- Scaled roundoff and pragmatic limits ----
		static constexpr long double RoundoffZeroThreshold = 100 * Eps;  // zero up to roundoff of O(100) flops at unit scale: norm before normalizing, integer-valuedness, trailing polynomial coefficients
		static constexpr long double UnboundedIntervalSamplingLimit = 1e10L;  // finite stand-in for +-infinity when sampling an unbounded interval
		static constexpr long double IntegrationRoundoffFloorScale = 50 * Eps;   // QUADPACK: convergence floor 50*eps*integral(|f|) below which a step error is roundoff
		static constexpr long double KrylovBreakdownThreshold = 1000 * Eps;      // Arnoldi/Lanczos: h_(j+1,j) below this fraction of the start norm = (lucky) breakdown
		static constexpr long double SVDCompletionTolerancePerDimension = 32 * Eps;   // complex SVD: a column is already spanned when its residual norm <= this * dimension
		static constexpr long double MonteCarloTargetRelativeError = 1e-3L;   // default statistical target for adaptive Monte Carlo
		// ---- Algorithm constants that used to be literals in mml/algorithms and mml/systems (TOLERANCE_POLICY section 1) ----
		static constexpr long double OverflowDetectionThreshold = 1e300L;   // |f(x)| beyond this is treated as overflow by the function analyzers
		static constexpr long double InfiniteLimitThreshold = 1e150L;       // one-sided limit magnitude beyond this is reported as infinite
		static constexpr long double DivergenceMagnitudeThreshold = 1e10L;  // |value| beyond this counts as divergent / unbounded (Lyapunov time, jump size)
		static constexpr long double PenaltyParameterMax = 1e12L;           // augmented Lagrangian: cap on the penalty parameter
		static constexpr long double LagrangeMultiplierBound = 1e10L;    // augmented Lagrangian: |multiplier| cap
		static constexpr long double BigMPenalty = 1e9L;                      // LP big-M value when not using two-phase simplex
		static constexpr long double ArmijoSufficientDecrease = 1e-4L;       // c1 of the Armijo/Wolfe sufficient-decrease condition (Nocedal-Wright 3.4)
		static constexpr long double MinimumStepRatio = 1e-10L;            // adaptive step below this fraction of the nominal step counts as collapsed
		static constexpr long double DAEMinStepSize = 1e-8L;                // DAE drivers: smallest admissible step
		static constexpr long double RODASMinStepFactor = 0.2L;             // RODAS4 controller shrink bound
		static constexpr long double RODASMaxStepFactor = 6.0L;             // RODAS4 controller growth bound
		static constexpr long double RODASRejectMaxFactor = 0.5L;           // A rejected RODAS4 attempt must shrink
		static constexpr long double StepControllerErrorFloor = 1e-4L;     // PI step controller: error norm floor so a zero error cannot give an infinite growth factor
		static constexpr long double WellConditionedLimit = 1e4L;             // condition number below this: well conditioned
		static constexpr long double IllConditionedLimit = 1e8L;               // condition number above this: ill conditioned
		static constexpr long double BrentZeroEpsilon = 1e-3L * SqrtEps;             // Brent/golden-section minimization ZEPS: absolute floor protecting x = 0 (NR 10.3)

		static constexpr long double GeometryEpsilon = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;

		static constexpr long double ComplexAreEqualTolerance = EqualityTolerance;
		static constexpr long double ComplexAreEqualAbsTolerance = EqualityTolerance;

		static constexpr long double MatrixIsEqualTolerance = EqualityTolerance;
		static constexpr long double VectorIsEqualTolerance = EqualityTolerance;

		static constexpr long double Pnt2CartIsEqualTolerance = EqualityTolerance;
		static constexpr long double Pnt2PolarIsEqualTolerance = EqualityTolerance;
		static constexpr long double Pnt3CartIsEqualTolerance = EqualityTolerance;
		static constexpr long double Pnt3SphIsEqualTolerance = EqualityTolerance;
		static constexpr long double Pnt3CylIsEqualTolerance = EqualityTolerance;

		static constexpr long double Vec2CartIsEqualTolerance = EqualityTolerance;
		static constexpr long double Vec3CartIsEqualTolerance = EqualityTolerance;
		static constexpr long double Vec3CartIsParallelTolerance = EqualityTolerance;

		static constexpr long double Vec3SphIsEqualTolerance = EqualityTolerance;

		// Angle comparison tolerance (for wrap-aware angle equality)
		static constexpr long double AngleIsEqualTolerance = EqualityTolerance;

		// Shape property tolerance (for geometric shape classification: IsRight, IsSquare, etc.)
		static constexpr long double ShapePropertyTolerance = EqualityTolerance;

		static constexpr long double Line3DAreEqualTolerance = EqualityTolerance;
		static constexpr long double Line3DIsPointOnLineTolerance = EqualityTolerance;
		static constexpr long double Line3DIsPerpendicularTolerance = EqualityTolerance;
		static constexpr long double Line3DIsParallelTolerance = EqualityTolerance;
		static constexpr long double Line3DIntersectionTolerance = EqualityTolerance;

		static constexpr long double Plane3DIsPointOnPlaneTolerance = EqualityTolerance;

		static constexpr long double Triangle3DIsPointInsideTolerance = EqualityTolerance;
		static constexpr long double Triangle3DIsRightTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double Triangle3DIsIsoscelesTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double Triangle3DIsEquilateralTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double IsMatrixSymmetricTolerance = EqualityTolerance;
		static constexpr long double IsMatrixDiagonalTolerance = EqualityTolerance;
		static constexpr long double IsMatrixUnitTolerance = EqualityTolerance;
		static constexpr long double IsMatrixZeroTolerance = EqualityTolerance;
		static constexpr long double IsMatrixOrthogonalTolerance = EqualityTolerance;

		static constexpr long double RankAlgEPS = LongDoubleHasExtraPrecision ? 1e-13L : 1e-12L;

		// Numerical computation thresholds
		static constexpr long double NumericalZeroThreshold = ZeroTolerance;
		static constexpr long double QuaternionZeroThreshold = ZeroTolerance;
		static constexpr long double SurfaceNormalThreshold = ZeroTolerance;
		static constexpr long double DerivativeStepSize = LongDoubleHasExtraPrecision ? 1e-8L : 1e-6L;
		static constexpr long double DefaultTolerance = 1e-6L;             // General purpose tolerance
		static constexpr long double DefaultToleranceStrict = 1e-8L;

		// Algorithm-specific thresholds
		static constexpr long double EigenSolverConvergenceTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double IterativeSolverTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double EigenSolverZeroThreshold = StrictTolerance;
		static constexpr long double PolynomialCoeffZeroThreshold = ZeroTolerance;
		static constexpr long double MatrixElementZeroThreshold = StrictTolerance;
		static constexpr long double DeterminantZeroThreshold = StrictTolerance;
		static constexpr long double DivisionSafetyThreshold = LongDoubleHasExtraPrecision ? 1e-35L : 1e-30L;
		static constexpr long double OrthogonalityTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double LinearDependenceTolerance = ZeroTolerance;

		// Optimization algorithm thresholds
		static constexpr long double OptimizationTolerance = LongDoubleHasExtraPrecision ? 3e-10L : 3e-8L;
		static constexpr long double OptimizationGradientTolerance = LongDoubleHasExtraPrecision ? 1e-10L : 1e-8L;

		// Integration thresholds
		static constexpr long double IntegrationTolerance = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;
		static constexpr long double TrapezoidIntegrationTolerance = 1e-4L;   // Legacy Defaults::TrapezoidIntegrationEPS
		static constexpr long double SimpsonIntegrationTolerance = 1e-5L;     // Legacy Defaults::SimpsonIntegrationEPS
		static constexpr long double RombergIntegrationTolerance = LooseTolerance; // Legacy Defaults::RombergIntegrationEPS
		static constexpr long double LineIntegralTolerance = 1e-5L;           // Legacy Defaults::LineIntegralPrecision / WorkIntegralPrecision

		// Special-function and dynamical-systems thresholds
		static constexpr long double PolynomialRootTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double AGMConvergenceTolerance = LongDoubleHasExtraPrecision ? 1e-18L : 1e-15L;
		static constexpr long double LyapunovChaosThreshold = 1e-6L;
		static constexpr long double FixedPointTolerance = LongDoubleHasExtraPrecision ? 1e-12L : 1e-10L;
		static constexpr long double FixedPointUniquenessTolerance = 1e-6L;

		// ODE solver thresholds
		static constexpr long double ODEDefaultTolerance = LongDoubleHasExtraPrecision ? 1e-14L : 1e-10L;
		static constexpr long double ODEMinStepSize = LongDoubleHasExtraPrecision ? 1e-18L : 1e-15L;
		static constexpr long double NewtonTolerance = LongDoubleHasExtraPrecision ? 1e-10L : 1e-8L;
		static constexpr long double NewtonRefreshRatio = 0.5L;
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
