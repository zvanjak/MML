///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ScalarFieldOperations.h                                             ///
///  Description: Gradient and Laplacian operations for scalar fields                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_SCALAR_FIELD_OPERATIONS_H
#define MML_SCALAR_FIELD_OPERATIONS_H

#include <mml/core/Fields/FieldOperationsCommon.h>
#include <mml/core/DifferentialGeometry/CoordinateSingularities.h>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////////////////
	/// ScalarFieldOperations - Differential operations on scalar fields
	///
	/// Provides gradient and Laplacian operations for scalar-valued functions of N variables.
	/// Supports:
	///   - Cartesian coordinates (any dimension N)
	///   - Spherical coordinates (3D): (r, θ, φ) where θ=polar, φ=azimuthal
	///   - Cylindrical coordinates (3D): (r, φ, z)
	///   - General curvilinear coordinates via metric tensor
	///
	/// Derivative accuracy can be controlled via der_order parameter (1, 2, 4, 6, or 8 points).
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace ScalarFieldOperations
	{
		///////////////////////////////////////////////////////////////////////////////////////////
		//                           GENERAL COORDINATE OPERATIONS
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes gradient of scalar field in general curvilinear coordinates.
		///
		/// Uses the metric tensor to properly handle non-Cartesian coordinate systems.
		/// The gradient is computed as: ∇f = gⁱʲ ∂ⱼf where gⁱʲ is the contravariant metric.
		///
		/// @tparam N          Number of dimensions
		/// @param scalarField Scalar function f: ℝᴺ → ℝ
		/// @param pos         Position vector in the coordinate system
		/// @param metricTensorField  Metric tensor field defining the coordinate geometry
		/// @return            Gradient vector ∇f with contravariant components
		///
		/// @note For Cartesian coordinates, use GradientCart() which is more efficient.
		/// @see MetricTensorField, GradientCart
		template<int N>
		static VectorN<Real, N> Gradient(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos,
		                                                                 const MetricTensorField<N>& metricTensorField)
		{
			// Gradient in general coordinates: ∇ᶠ = gⁱʲ ∂ᵢf
			// The partial derivatives ∂ᵢf give covariant components
			// We need the contravariant metric gⁱʲ to raise indices
			VectorN<Real, N> covar_derivs = Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);

			MatrixNM<Real, N, N> g_contravar = metricTensorField.GetContravariantMetric(pos);

			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0.0;
				for (int j = 0; j < N; j++)
					ret[i] += g_contravar[i][j] * covar_derivs[j];
			}

			return ret;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                              GRADIENT - CARTESIAN
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes gradient of scalar field in Cartesian coordinates.
		///
		/// For f(x₁, x₂, ..., xₙ), returns ∇f = (∂f/∂x₁, ∂f/∂x₂, ..., ∂f/∂xₙ)
		/// Uses default numerical differentiation (4-point formula).
		///
		/// @tparam N          Number of dimensions
		/// @param scalarField Scalar function f: ℝᴺ → ℝ
		/// @param pos         Position vector (x₁, x₂, ..., xₙ)
		/// @return            Gradient vector ∇f
		///
		/// @example
		///   // f(x,y,z) = x² + 2y + 3z
		///   ScalarFunction<3> f([](auto& p){ return p[0]*p[0] + 2*p[1] + 3*p[2]; });
		///   auto grad = GradientCart(f, {1,1,1});  // Returns (2, 2, 3)
		template<int N>
		static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos)
		{
			return Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);
		}

		/// Computes gradient in Cartesian coordinates with specified derivative accuracy.
		///
		/// @tparam N          Number of dimensions
		/// @param scalarField Scalar function f: ℝᴺ → ℝ
		/// @param pos         Position vector (x₁, x₂, ..., xₙ)
		/// @param der_order   Derivative accuracy: 1, 2, 4, 6, or 8 (higher = more accurate but slower)
		/// @return            Gradient vector ∇f
		/// @throws std::invalid_argument if der_order is not in {1, 2, 4, 6, 8}
		template<int N>
		static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos,
		                                                                         int der_order)
		{
			switch (der_order)
			{
			case 1: return Derivation::NDer1PartialByAll<N>(scalarField, pos, nullptr);
			case 2: return Derivation::NDer2PartialByAll<N>(scalarField, pos, nullptr);
			case 4: return Derivation::NDer4PartialByAll<N>(scalarField, pos, nullptr);
			case 6: return Derivation::NDer6PartialByAll<N>(scalarField, pos, nullptr);
			case 8: return Derivation::NDer8PartialByAll<N>(scalarField, pos, nullptr);
			default:
				throw ArgumentError("GradientCart: der_order must be in 1, 2, 4, 6 or 8");
			}
		}

		/// Computes gradient in Cartesian coordinates with structured result.
		///
		/// Returns an EvaluationResult containing the gradient vector, optional
		/// per-component error estimates, timing, and AlgorithmStatus.
		///
		/// @tparam N          Number of dimensions
		/// @param scalarField Scalar function f: ℝᴺ → ℝ
		/// @param pos         Position vector (x₁, x₂, ..., xₙ)
		/// @param config      Field operation configuration
		/// @return            EvaluationResult with gradient value and diagnostics
		template<int N>
		static EvaluationResult<VectorN<Real, N>, VectorN<Real, N>>
		GradientCartDetailed(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos,
		                     const FieldOperationConfig& config = {})
		{
			using ResultType = EvaluationResult<VectorN<Real, N>, VectorN<Real, N>>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"GradientCart", config,
				[&](ResultType& result, int& func_evals) {
					VectorN<Real, N> error_vec{};
					result.value = FieldOperationDetail::DispatchGradient<N>(
						scalarField, pos, config.derivative_order,
						config.estimate_error, &error_vec, func_evals);
					if (config.estimate_error)
						result.error = error_vec;
				});
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                              GRADIENT - SPHERICAL
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes gradient in spherical coordinates (r, θ, φ).
		///
		/// For spherical coordinates where:
		///   - r ∈ [0, ∞)     : radial distance from origin
		///   - θ ∈ [0, π]     : polar angle from z-axis
		///   - φ ∈ [0, 2π)    : azimuthal angle in xy-plane
		///
		/// Returns: ∇f = (∂f/∂r, (1/r)∂f/∂θ, (1/r·sinθ)∂f/∂φ)
		///
		/// @param scalarField Scalar function f(r, θ, φ)
		/// @param pos         Position vector (r, θ, φ) in spherical coordinates
		/// @return            Gradient vector in spherical basis (êᵣ, êθ, êφ)
		///
		/// @warning Position must have r > 0 and θ ≠ 0, π to avoid singularities.
		/// @param policy Singularity handling policy (default: Throw)
		static VectorN<Real, 3> GradientSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos,
		                                      SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			Vector3Spherical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

			const Real r = pos[0];
			const Real theta = pos[1];

			// θ-component: (1/r) ∂f/∂θ - singular at r=0
			ret[1] = ret[1] * Singularity::SafeInverseR(r, policy, "GradientSpher θ-component");

			// φ-component: (1/(r·sinθ)) ∂f/∂φ - singular at r=0 or poles
			ret[2] = ret[2] * Singularity::SafeInverseRSinTheta(r, theta, policy, "GradientSpher φ-component");

			return ret;
		}

		/// Computes gradient in spherical coordinates with specified derivative accuracy.
		///
		/// @param scalarField Scalar function f(r, θ, φ)
		/// @param pos         Position vector (r, θ, φ)
		/// @param der_order   Derivative accuracy: 1, 2, 4, 6, or 8
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            Gradient vector in spherical basis
		/// @throws std::invalid_argument if der_order is not in {1, 2, 4, 6, 8}
		/// @throws DomainError if at singularity and policy is Throw
		static Vec3Sph GradientSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos,
		                             int der_order,
		                             SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			Vector3Spherical ret;

			switch (der_order)
			{
			case 1: ret = Derivation::NDer1PartialByAll<3>(scalarField, pos, nullptr); break;
			case 2: ret = Derivation::NDer2PartialByAll<3>(scalarField, pos, nullptr); break;
			case 4: ret = Derivation::NDer4PartialByAll<3>(scalarField, pos, nullptr); break;
			case 6: ret = Derivation::NDer6PartialByAll<3>(scalarField, pos, nullptr); break;
			case 8: ret = Derivation::NDer8PartialByAll<3>(scalarField, pos, nullptr); break;
			default:
				throw ArgumentError("GradientSpher: der_order must be in 1, 2, 4, 6 or 8");
			}

			const Real r = pos[0];
			const Real theta = pos[1];

			// θ-component: (1/r) ∂f/∂θ - singular at r=0
			ret[1] = ret[1] * Singularity::SafeInverseR(r, policy, "GradientSpher θ-component");

			// φ-component: (1/(r·sinθ)) ∂f/∂φ - singular at r=0 or poles
			ret[2] = ret[2] * Singularity::SafeInverseRSinTheta(r, theta, policy, "GradientSpher φ-component");

			return ret;
		}

		/// Computes gradient in spherical coordinates with structured result.
		///
		/// @param scalarField Scalar function f(r, θ, φ)
		/// @param pos         Position (r, θ, φ) with r > 0, 0 < θ < π
		/// @param config      Field operation configuration
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            EvaluationResult with gradient vector and diagnostics
		static EvaluationResult<Vec3Sph, Vec3Sph>
		GradientSpherDetailed(const IScalarFunction<3>& scalarField, const Vec3Sph& pos,
		                      const FieldOperationConfig& config = {},
		                      SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			using ResultType = EvaluationResult<Vec3Sph, Vec3Sph>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"GradientSpher", config,
				[&](ResultType& result, int& func_evals) {
					VectorN<Real, 3> error_vec{};
					result.value = FieldOperationDetail::DispatchGradient<3>(
						scalarField, pos, config.derivative_order,
						config.estimate_error, &error_vec, func_evals);
					if (config.estimate_error)
						result.error = error_vec;

					const Real r = pos[0];
					const Real theta = pos[1];
					result.value[1] *= Singularity::SafeInverseR(r, policy, "GradientSpher θ-component");
					result.value[2] *= Singularity::SafeInverseRSinTheta(r, theta, policy, "GradientSpher φ-component");
					// Scale error estimates by the same factors
					if (config.estimate_error) {
						result.error[1] *= std::abs(Singularity::SafeInverseR(r, policy, "GradientSpher θ-error"));
						result.error[2] *= std::abs(Singularity::SafeInverseRSinTheta(r, theta, policy, "GradientSpher φ-error"));
					}
				});
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                              GRADIENT - CYLINDRICAL
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes gradient in cylindrical coordinates (r, φ, z).
		///
		/// For cylindrical coordinates where:
		///   - r ∈ [0, ∞)     : radial distance from z-axis
		///   - φ ∈ [0, 2π)    : azimuthal angle in xy-plane
		///   - z ∈ (-∞, ∞)   : height along z-axis
		///
		/// Returns: ∇f = (∂f/∂r, (1/r)∂f/∂φ, ∂f/∂z)
		///
		/// @param scalarField Scalar function f(r, φ, z)
		/// @param pos         Position vector (r, φ, z) in cylindrical coordinates
		/// @return            Gradient vector in cylindrical basis (êᵣ, êφ, êz)
		///
		/// @warning Position must have r > 0 to avoid singularity at z-axis.
		/// @param policy Singularity handling policy (default: Throw)
		static Vec3Cyl GradientCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos,
		                           SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			Vector3Cylindrical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

			const Real r = pos[0];

			// φ-component: (1/r) ∂f/∂φ - singular at r=0
			ret[1] = ret[1] * Singularity::SafeInverseR(r, policy, "GradientCyl φ-component");

			return ret;
		}

		/// Computes gradient in cylindrical coordinates with specified derivative accuracy.
		///
		/// @param scalarField Scalar function f(r, φ, z)
		/// @param pos         Position vector (r, φ, z)
		/// @param der_order   Derivative accuracy: 1, 2, 4, 6, or 8
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            Gradient vector in cylindrical basis
		/// @throws std::invalid_argument if der_order is not in {1, 2, 4, 6, 8}
		/// @throws DomainError if at singularity and policy is Throw
		static Vec3Cyl GradientCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos,
		                           int der_order,
		                           SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			Vector3Cylindrical ret;

			switch (der_order)
			{
			case 1: ret = Derivation::NDer1PartialByAll<3>(scalarField, pos, nullptr); break;
			case 2: ret = Derivation::NDer2PartialByAll<3>(scalarField, pos, nullptr); break;
			case 4: ret = Derivation::NDer4PartialByAll<3>(scalarField, pos, nullptr); break;
			case 6: ret = Derivation::NDer6PartialByAll<3>(scalarField, pos, nullptr); break;
			case 8: ret = Derivation::NDer8PartialByAll<3>(scalarField, pos, nullptr); break;
			default:
				throw ArgumentError("GradientCyl: der_order must be in 1, 2, 4, 6 or 8");
			}

			const Real r = pos[0];

			// φ-component: (1/r) ∂f/∂φ - singular at r=0
			ret[1] = ret[1] * Singularity::SafeInverseR(r, policy, "GradientCyl φ-component");

			return ret;
		}

		/// Computes gradient in cylindrical coordinates with structured result.
		///
		/// @param scalarField Scalar function f(r, φ, z)
		/// @param pos         Position (r, φ, z) with r > 0
		/// @param config      Field operation configuration
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            EvaluationResult with gradient vector and diagnostics
		static EvaluationResult<Vec3Cyl, Vec3Cyl>
		GradientCylDetailed(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos,
		                    const FieldOperationConfig& config = {},
		                    SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			using ResultType = EvaluationResult<Vec3Cyl, Vec3Cyl>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"GradientCyl", config,
				[&](ResultType& result, int& func_evals) {
					VectorN<Real, 3> error_vec{};
					result.value = FieldOperationDetail::DispatchGradient<3>(
						scalarField, pos, config.derivative_order,
						config.estimate_error, &error_vec, func_evals);
					if (config.estimate_error)
						result.error = error_vec;

					const Real r = pos[0];
					result.value[1] *= Singularity::SafeInverseR(r, policy, "GradientCyl φ-component");
					if (config.estimate_error)
						result.error[1] *= std::abs(Singularity::SafeInverseR(r, policy, "GradientCyl φ-error"));
				});
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                                    LAPLACIAN
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes Laplacian of scalar field in Cartesian coordinates.
		///
		/// The Laplacian ∇²f = Σᵢ ∂²f/∂xᵢ² measures the "curvature" of the scalar field,
		/// or equivalently, how much the value at a point differs from the average
		/// of surrounding points. It appears in:
		///   - Heat equation: ∂T/∂t = α∇²T
		///   - Wave equation: ∂²u/∂t² = c²∇²u
		///   - Poisson equation: ∇²φ = -ρ/ε₀
		///
		/// @tparam N          Number of dimensions
		/// @param scalarField Scalar function f: ℝᴺ → ℝ
		/// @param pos         Position vector
		/// @return            Scalar Laplacian value ∇²f
		///
		/// @example
		///   // f(x,y,z) = x² + y² + z² (paraboloid)
		///   ScalarFunction<3> f([](auto& p){ return p[0]*p[0] + p[1]*p[1] + p[2]*p[2]; });
		///   Real lapl = LaplacianCart(f, {1,1,1});  // Returns 6.0 (constant curvature)
		template<int N>
		static Real LaplacianCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos)
		{
			Real lapl = 0.0;
			for (int i = 0; i < N; i++)
				lapl += Derivation::DeriveSecPartial<N>(scalarField, i, i, pos, nullptr);

			return lapl;
		}

		/// Computes Laplacian in Cartesian coordinates with structured result.
		///
		/// @tparam N          Number of dimensions
		/// @param scalarField Scalar function f: ℝᴺ → ℝ
		/// @param pos         Position vector
		/// @param config      Field operation configuration
		/// @return            EvaluationResult with Laplacian value and diagnostics
		template<int N>
		static EvaluationResult<Real>
		LaplacianCartDetailed(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos,
		                      const FieldOperationConfig& config = {})
		{
			using ResultType = EvaluationResult<Real>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"LaplacianCart", config,
				[&](ResultType& result, int& func_evals) {
					Real lapl = 0.0;
					Real total_error = 0.0;
					for (int i = 0; i < N; i++) {
						Real err = 0.0;
						lapl += Derivation::DeriveSecPartial<N>(scalarField, i, i, pos,
						                                        config.estimate_error ? &err : nullptr);
						if (config.estimate_error)
							total_error += std::abs(err);
					}
					result.value = lapl;
					if (config.estimate_error)
						result.error = total_error;
					func_evals = N * 3; // second derivative uses ~3 evaluations per dimension
				});
		}

		/// Computes Laplacian of scalar field in spherical coordinates.
		///
		/// In spherical coordinates (r, θ, φ):
		///   ∇²f = (1/r²)∂/∂r(r²∂f/∂r) + (1/r²sinθ)∂/∂θ(sinθ·∂f/∂θ) + (1/r²sin²θ)∂²f/∂φ²
		///
		/// @param scalarField Scalar function f(r, θ, φ)
		/// @param pos         Position (r, θ, φ) with r > 0, 0 < θ < π
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            Scalar Laplacian value ∇²f
		///
		/// @warning Singular at r = 0 and θ = 0, π (coordinate singularities)
		/// @throws DomainError if at singularity and policy is Throw
		static Real LaplacianSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos,
		                           SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			const Real r = pos.R();
			const Real theta = pos.Theta();
			// Note: phi = pos.Phi() is not needed directly in the Laplacian formula

			// ∇²f = ∂²f/∂r² + (2/r)∂f/∂r + (cotθ/r²)∂f/∂θ + (1/r²)∂²f/∂θ² + (1/r²sin²θ)∂²f/∂φ²

			// Radial term: ∂²f/∂r²
			Real d2f_dr2 = Derivation::DeriveSecPartial<3>(scalarField, 0, 0, pos, nullptr);

			// Radial correction: (2/r)∂f/∂r - singular at r=0
			Real df_dr = Derivation::DerivePartial<3>(scalarField, 0, pos, nullptr);
			Real inv_r = Singularity::SafeInverseR(r, policy, "LaplacianSpher radial term");
			Real radial_correction = Real(2) * inv_r * df_dr;

			// Polar (θ) terms - singular at r=0 and poles
			Real df_dtheta = Derivation::DerivePartial<3>(scalarField, 1, pos, nullptr);
			Real d2f_dtheta2 = Derivation::DeriveSecPartial<3>(scalarField, 1, 1, pos, nullptr);
			Real inv_r2 = Singularity::SafeInverseR2(r, policy, "LaplacianSpher 1/r²");
			Real cot_over_r2 = Singularity::SafeCotThetaOverR2(r, theta, policy, "LaplacianSpher cotθ/r²");
			Real theta_term = cot_over_r2 * df_dtheta + inv_r2 * d2f_dtheta2;

			// Azimuthal (φ) term: (1/r²sin²θ)∂²f/∂φ² - singular at r=0 and poles
			Real d2f_dphi2 = Derivation::DeriveSecPartial<3>(scalarField, 2, 2, pos, nullptr);
			Real inv_r2_sin2 = Singularity::SafeInverseR2Sin2Theta(r, theta, policy, "LaplacianSpher φ-term");
			Real phi_term = inv_r2_sin2 * d2f_dphi2;

			return d2f_dr2 + radial_correction + theta_term + phi_term;
		}

		/// Computes Laplacian of scalar field in cylindrical coordinates.
		///
		/// In cylindrical coordinates (r, φ, z):
		///   ∇²f = (1/r)∂/∂r(r·∂f/∂r) + (1/r²)∂²f/∂φ² + ∂²f/∂z²
		///
		/// @param scalarField Scalar function f(r, φ, z)
		/// @param pos         Position (r, φ, z) with r > 0
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            Scalar Laplacian value ∇²f
		///
		/// @warning Singular at r = 0 (z-axis singularity)
		/// @throws DomainError if at singularity and policy is Throw
		static Real LaplacianCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos,
		                         SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			const Real r = pos[0];

			// ∇²f = (1/r)∂/∂r(r·∂f/∂r) + (1/r²)∂²f/∂φ² + ∂²f/∂z²
			//     = (1/r)∂f/∂r + ∂²f/∂r² + (1/r²)∂²f/∂φ² + ∂²f/∂z²

			// Radial term: (1/r)∂f/∂r + ∂²f/∂r² - singular at r=0
			Real df_dr = Derivation::DerivePartial<3>(scalarField, 0, pos, nullptr);
			Real d2f_dr2 = Derivation::DeriveSecPartial<3>(scalarField, 0, 0, pos, nullptr);
			Real inv_r = Singularity::SafeInverseR(r, policy, "LaplacianCyl radial term");
			Real radial_term = inv_r * df_dr + d2f_dr2;

			// Azimuthal term: (1/r²)∂²f/∂φ² - singular at r=0
			Real d2f_dphi2 = Derivation::DeriveSecPartial<3>(scalarField, 1, 1, pos, nullptr);
			Real inv_r2 = Singularity::SafeInverseR2(r, policy, "LaplacianCyl φ-term");
			Real phi_term = inv_r2 * d2f_dphi2;

			// Axial term: ∂²f/∂z² (no singularity)
			Real d2f_dz2 = Derivation::DeriveSecPartial<3>(scalarField, 2, 2, pos, nullptr);

			return radial_term + phi_term + d2f_dz2;
		}

		/// Computes Laplacian in spherical coordinates with structured result.
		///
		/// @param scalarField Scalar function f(r, θ, φ)
		/// @param pos         Position (r, θ, φ) with r > 0, 0 < θ < π
		/// @param config      Field operation configuration
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            EvaluationResult with Laplacian value and diagnostics
		static EvaluationResult<Real>
		LaplacianSpherDetailed(const IScalarFunction<3>& scalarField, const Vec3Sph& pos,
		                       const FieldOperationConfig& config = {},
		                       SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			using ResultType = EvaluationResult<Real>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"LaplacianSpher", config,
				[&](ResultType& result, int& func_evals) {
					result.value = LaplacianSpher(scalarField, pos, policy);
					func_evals = 7; // 2 first partials + 3 second partials + field evals
				});
		}

		/// Computes Laplacian in cylindrical coordinates with structured result.
		///
		/// @param scalarField Scalar function f(r, φ, z)
		/// @param pos         Position (r, φ, z) with r > 0
		/// @param config      Field operation configuration
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            EvaluationResult with Laplacian value and diagnostics
		static EvaluationResult<Real>
		LaplacianCylDetailed(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos,
		                     const FieldOperationConfig& config = {},
		                     SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			using ResultType = EvaluationResult<Real>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"LaplacianCyl", config,
				[&](ResultType& result, int& func_evals) {
					result.value = LaplacianCyl(scalarField, pos, policy);
					func_evals = 5; // 1 first partial + 3 second partials
				});
		}
	};
}

#endif