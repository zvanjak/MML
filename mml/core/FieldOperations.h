///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FieldOperations.h                                                   ///
///  Description: Vector field differential operations for mathematical physics       ///
///               - Gradient: scalar field → vector field (∇f)                        ///
///               - Divergence: vector field → scalar field (∇·F)                     ///
///               - Curl: vector field → vector field (∇×F)                           ///
///               - Laplacian: scalar field → scalar field (∇²f)                      ///
///               Supports Cartesian, spherical, and cylindrical coordinates          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
///
/// MATHEMATICAL BACKGROUND
/// =======================
///
/// Gradient (∇f)
/// -------------
/// The gradient of a scalar field f points in the direction of steepest increase
/// with magnitude equal to the rate of change in that direction.
///
///   Cartesian (x,y,z):    ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z)
///
///   Spherical (r,θ,φ):    ∇f = (∂f/∂r, (1/r)∂f/∂θ, (1/r·sinθ)∂f/∂φ)
///
///   Cylindrical (r,φ,z):  ∇f = (∂f/∂r, (1/r)∂f/∂φ, ∂f/∂z)
///
/// Divergence (∇·F)
/// ----------------
/// The divergence measures the "outflow" of a vector field at a point.
/// Positive divergence indicates a source, negative indicates a sink.
///
///   Cartesian:    ∇·F = ∂Fx/∂x + ∂Fy/∂y + ∂Fz/∂z
///
///   Spherical:    ∇·F = (1/r²)∂(r²Fr)/∂r + (1/r·sinθ)∂(sinθ·Fθ)/∂θ + (1/r·sinθ)∂Fφ/∂φ
///
///   Cylindrical:  ∇·F = (1/r)∂(r·Fr)/∂r + (1/r)∂Fφ/∂φ + ∂Fz/∂z
///
/// Curl (∇×F)
/// ----------
/// The curl measures the rotation or circulation of a vector field.
/// In 3D, it produces a vector perpendicular to the plane of rotation.
///
///   Cartesian:    ∇×F = (∂Fz/∂y - ∂Fy/∂z, ∂Fx/∂z - ∂Fz/∂x, ∂Fy/∂x - ∂Fx/∂y)
///
/// Laplacian (∇²f)
/// ---------------
/// The Laplacian measures the difference between the value at a point and
/// the average of surrounding values. Important in heat/wave equations.
///
///   Cartesian:    ∇²f = ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z²
///
///   Spherical:    ∇²f = (1/r²)∂/∂r(r²∂f/∂r) + (1/r²sinθ)∂/∂θ(sinθ·∂f/∂θ) + (1/r²sin²θ)∂²f/∂φ²
///
///   Cylindrical:  ∇²f = (1/r)∂/∂r(r·∂f/∂r) + (1/r²)∂²f/∂φ² + ∂²f/∂z²
///
/// General Coordinates
/// -------------------
/// For arbitrary curvilinear coordinates, operations use the metric tensor gᵢⱼ
/// and Christoffel symbols Γⁱⱼₖ to handle coordinate curvature properly.
///
/// USAGE EXAMPLE
/// =============
///   // Scalar field: f(x,y,z) = x² + y² + z²
///   ScalarFunction<3> f([](const VectorN<Real,3>& p) {
///       return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
///   });
///   
///   VectorN<Real,3> pos{1.0, 2.0, 3.0};
///   auto grad = ScalarFieldOperations::GradientCart<3>(f, pos);  // (2, 4, 6)
///   Real lapl = ScalarFieldOperations::LaplacianCart<3>(f, pos); // 6.0
///
///   // Vector field: F(x,y,z) = (y, -x, 0)  (rotation about z-axis)
///   VectorFunction<3> F([](const VectorN<Real,3>& p) {
///       return VectorN<Real,3>{p[1], -p[0], 0.0};
///   });
///   
///   Real div = VectorFieldOperations::DivCart<3>(F, pos);       // 0 (incompressible)
///   auto curl = VectorFieldOperations::CurlCart(F, pos);        // (0, 0, -2)
///
/// SEE ALSO
/// ========
///   - Fields.h: Predefined physical fields (gravity, electric, magnetic)
///   - MetricTensor.h: Metric tensors for curvilinear coordinates
///   - Derivation.h: Numerical differentiation routines
///   - CoordTransf.h: Coordinate transformations
///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_FIELD_OPERATIONS_H
#define MML_FIELD_OPERATIONS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "mml/base/Geometry/Geometry.h"
#include "base/Vector/VectorN.h"
#include "base/Vector/VectorTypes.h"
#include "base/Matrix/MatrixNM.h"
#include "base/Matrix/Matrix.h"
#include "base/Function.h"
#include "base/Tensor.h"

#include "core/Derivation.h"
#include "core/MetricTensor.h"
#include "core/MatrixUtils.h"
#include "core/SingularityHandling.h"
#include "core/AlgorithmTypes.h"

namespace MML
{
	/******************************************************************************/
	/*****               Field Operation Configuration                       *****/
	/******************************************************************************/

	/// Configuration for field differential operations (gradient, divergence, curl, Laplacian).
	///
	/// Extends EvaluationConfigBase with derivative order selection.
	/// Field operations are non-iterative evaluations that internally use numerical
	/// differentiation, so derivative accuracy is the primary tunable parameter.
	struct FieldOperationConfig : public EvaluationConfigBase {
		/// Derivative order for internal numerical differentiation.
		/// 0 = use library default (NDer4), or 1, 2, 4, 6, 8.
		int derivative_order = 0;
	};

	/******************************************************************************/
	/*****               Field Operation Detail Helpers                      *****/
	/******************************************************************************/

	namespace FieldOperationDetail
	{
		/// Execute a field operation with structured result reporting.
		///
		/// Wraps the raw computation in timing, finiteness checking, and
		/// exception handling according to the configured policy.
		///
		/// @tparam ResultType  EvaluationResult specialization for the output
		/// @tparam ComputeFn   Lambda (int& func_evals) -> void that populates result.value/error
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteFieldDetailed(const char* algorithm_name,
		                                const FieldOperationConfig& config,
		                                ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;

				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);

				int func_evals = 0;
				compute(result, func_evals);
				result.function_evaluations = func_evals;

				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const DomainError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const NumericInputError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const NumericalMethodError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::NumericalInstability, ex.what(), algorithm_name);
			}
			catch (const std::invalid_argument& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const std::exception& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
			}
		}

		/// Dispatch derivative-order selection for DerivePartialAll.
		/// Returns the gradient and optionally populates per-component error estimates.
		template<int N>
		VectorN<Real, N> DispatchGradient(const IScalarFunction<N>& f, const VectorN<Real, N>& pos,
		                                  int order, bool want_error, VectorN<Real, N>* error,
		                                  int& func_evals)
		{
			// Function pointer type for NDerKPartialByAll(f, pos, error*)
			using DeriveFn = VectorN<Real, N>(*)(const IScalarFunction<N>&,
			                                     const VectorN<Real, N>&,
			                                     VectorN<Real, N>*);

			DeriveFn deriveAll = nullptr;
			int stencil = 0;

			switch (order) {
			case 1:  deriveAll = &Derivation::template NDer1PartialByAll<N>; stencil = 2; break;
			case 2:  deriveAll = &Derivation::template NDer2PartialByAll<N>; stencil = 3; break;
			case 0:  // fall through to default (NDer4)
			case 4:  deriveAll = &Derivation::template NDer4PartialByAll<N>; stencil = 5; break;
			case 6:  deriveAll = &Derivation::template NDer6PartialByAll<N>; stencil = 7; break;
			case 8:  deriveAll = &Derivation::template NDer8PartialByAll<N>; stencil = 9; break;
			default:
				throw std::invalid_argument("FieldOperation: derivative_order must be 0, 1, 2, 4, 6, or 8");
			}

			func_evals = N * (want_error ? stencil + 1 : stencil);
			return deriveAll(f, pos, want_error ? error : nullptr);
		}
	} // namespace FieldOperationDetail

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
		static VectorN<Real, N> Gradient(IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos, 
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
				throw std::invalid_argument("GradientCart: der_order must be in 1, 2, 4, 6 or 8");
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
				throw std::invalid_argument("GradientSpher: der_order must be in 1, 2, 4, 6 or 8");
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
				throw std::invalid_argument("GradientCyl: der_order must be in 1, 2, 4, 6 or 8");
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

	///////////////////////////////////////////////////////////////////////////////////////////
	/// VectorFieldOperations - Differential operations on vector fields
	///
	/// Provides divergence and curl operations for vector-valued functions.
	/// Supports:
	///   - Cartesian coordinates (any dimension N for divergence, 3D for curl)
	///   - Spherical coordinates (3D): (r, θ, φ)
	///   - Cylindrical coordinates (3D): (r, φ, z)
	///   - General curvilinear coordinates via metric tensor and Christoffel symbols
	///
	/// Physical interpretation:
	///   - Divergence: measures local "expansion" or "compression" of a flow field
	///   - Curl: measures local rotation or circulation of a flow field
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace VectorFieldOperations
	{
		///////////////////////////////////////////////////////////////////////////////////////////
		//                           GENERAL COORDINATE OPERATIONS
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes divergence of vector field in general curvilinear coordinates.
		///
		/// Uses the metric tensor and Christoffel symbols to properly handle
		/// non-Cartesian coordinate systems. The divergence in general coordinates is:
		///   ∇·F = ∂ᵢFⁱ + Γⁱᵢₖ Fᵏ
		///
		/// where Γⁱⱼₖ are the Christoffel symbols of the second kind.
		///
		/// @tparam N                Number of dimensions
		/// @param vectorField       Vector function F: ℝᴺ → ℝᴺ
		/// @param pos               Position vector in the coordinate system
		/// @param metricTensorField Metric tensor field with Christoffel symbols
		/// @return                  Scalar divergence value ∇·F
		///
		/// @note For Cartesian coordinates, use DivCart() which is more efficient.
		/// @see MetricTensorField, DivCart
		template<int N>
		static Real Divergence(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos,
													const MetricTensorField<N>& metricTensorField)
		{
			Real div = 0.0;
			VectorN<Real, N> vec_val = vectorField(pos);

			for (int i = 0; i < N; i++)
			{
				// Standard partial derivative term: ∂Fⁱ/∂xⁱ
				div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, nullptr);

				// Christoffel symbol correction for curved coordinates: Γⁱᵢₖ Fᵏ
				for (int k = 0; k < N; k++)
				{
					div += vec_val[k] * metricTensorField.GetChristoffelSymbolSecondKind(i, i, k, pos);
				}
			}
			return div;
		}

		/// Computes divergence via the metric determinant formula (alternative method).
		///
		/// Uses the identity:
		///   ∇·F = (1/√g) ∂ᵢ(√g Fⁱ)
		///
		/// where g = det(gᵢⱼ) is the determinant of the covariant metric tensor.
		/// This is mathematically equivalent to the Christoffel-based formula but:
		///   - Avoids computing N² Christoffel symbols
		///   - Uses only the metric determinant (single scalar per point)
		///   - Can be more numerically stable for some coordinate systems
		///
		/// @tparam N                Number of dimensions
		/// @param vectorField       Vector function F: ℝᴺ → ℝᴺ (contravariant components)
		/// @param pos               Position vector in the coordinate system
		/// @param metricTensorField Metric tensor field defining the coordinate geometry
		/// @return                  Scalar divergence value ∇·F
		///
		/// @note For Cartesian coordinates, √g = 1 everywhere and this reduces to DivCart.
		/// @note For spherical coordinates, √g = r²sinθ.
		/// @note For cylindrical coordinates, √g = r.
		/// @see Divergence, DivCart
		template<int N>
		static Real DivergenceViaDet(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos,
		                             const MetricTensorField<N>& metricTensorField)
		{
			// Compute √g at the evaluation point
			Real g = Utils::Det(metricTensorField.GetCovariantMetric(pos));
			Real sqrt_g = std::sqrt(std::abs(g));

			if (sqrt_g < Defaults::DefaultTolerance)
				return Real(0);  // Degenerate metric — divergence undefined

			Real div = 0.0;
			for (int i = 0; i < N; i++)
			{
				// Differentiate the scalar function √g(x) · Fⁱ(x) with respect to xⁱ
				ScalarFunctionFromStdFunc<N> sqrt_g_Fi(
					[&vectorField, &metricTensorField, i](const VectorN<Real, N>& x) -> Real {
						Real gx = Utils::Det(metricTensorField.GetCovariantMetric(x));
						return std::sqrt(std::abs(gx)) * vectorField(x)[i];
					}
				);
				div += Derivation::DerivePartial<N>(sqrt_g_Fi, i, pos, nullptr);
			}

			return div / sqrt_g;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                              DIVERGENCE - CARTESIAN
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes divergence of vector field in Cartesian coordinates.
		///
		/// For F = (F₁, F₂, ..., Fₙ), returns ∇·F = ∂F₁/∂x₁ + ∂F₂/∂x₂ + ... + ∂Fₙ/∂xₙ
		///
		/// Physical interpretation:
		///   - div > 0: source (field "emanates" from the point)
		///   - div < 0: sink (field "converges" to the point)
		///   - div = 0: incompressible flow (volume-preserving)
		///
		/// @tparam N          Number of dimensions
		/// @param vectorField Vector function F: ℝᴺ → ℝᴺ
		/// @param pos         Position vector
		/// @return            Scalar divergence value ∇·F
		///
		/// @example
		///   // F(x,y,z) = (x, y, z) - radial outward field
		///   VectorFunction<3> F([](auto& p){ return p; });
		///   Real div = DivCart(F, {1,1,1});  // Returns 3.0 (uniform expansion)
		///
		///   // F(x,y,z) = (y, -x, 0) - rotation about z-axis
		///   VectorFunction<3> G([](auto& p){ return VectorN<Real,3>{p[1], -p[0], 0}; });
		///   Real div = DivCart(G, {1,1,1});  // Returns 0.0 (incompressible)
		template<int N>
		static Real DivCart(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos)
		{
			Real div = 0.0;
			for (int i = 0; i < N; i++)
				div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, nullptr);

			return div;
		}

		/// Computes divergence in Cartesian coordinates with structured result.
		///
		/// @tparam N          Number of dimensions
		/// @param vectorField Vector function F: ℝᴺ → ℝᴺ
		/// @param pos         Position vector
		/// @param config      Field operation configuration
		/// @return            EvaluationResult with divergence value and diagnostics
		template<int N>
		static EvaluationResult<Real>
		DivCartDetailed(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos,
		                const FieldOperationConfig& config = {})
		{
			using ResultType = EvaluationResult<Real>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"DivCart", config,
				[&](ResultType& result, int& func_evals) {
					Real div = 0.0;
					Real total_error = 0.0;
					for (int i = 0; i < N; i++) {
						Real err = 0.0;
						div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos,
						                                       config.estimate_error ? &err : nullptr);
						if (config.estimate_error)
							total_error += std::abs(err);
					}
					result.value = div;
					if (config.estimate_error)
						result.error = total_error;
					func_evals = N * 5; // NDer4 default: ~5 evaluations per component
				});
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                              DIVERGENCE - SPHERICAL
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes divergence in spherical coordinates (r, θ, φ).
		///
		/// In spherical coordinates:
		///   ∇·F = (1/r²)∂(r²Fᵣ)/∂r + (1/r·sinθ)∂(sinθ·Fθ)/∂θ + (1/r·sinθ)∂Fφ/∂φ
		///
		/// @param vectorField Vector function F(r, θ, φ) = (Fᵣ, Fθ, Fφ)
		/// @param x           Position (r, θ, φ) with r > 0, 0 < θ < π
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            Scalar divergence value ∇·F
		///
		/// @warning Singular at r = 0 and θ = 0, π (coordinate singularities)
		/// @throws DomainError if at singularity and policy is Throw
		static Real DivSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x,
		                     SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			const Real r = x[0];
			const Real theta = x[1];
			
			VectorN<Real, 3> vals = vectorField(x);

			VectorN<Real, 3> derivs;
			for (int i = 0; i < 3; i++)
				derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);

			// Calculate inverse factors with singularity handling
			Real inv_r2 = Singularity::SafeInverseR2(r, policy, "DivSpher 1/r²");
			Real inv_r_sin = Singularity::SafeInverseRSinTheta(r, theta, policy, "DivSpher 1/(r·sinθ)");

			Real div = 0.0;
			// r-component: (1/r²)∂(r²Fᵣ)/∂r = (2/r)Fᵣ + ∂Fᵣ/∂r
			div += inv_r2 * (2 * r * vals[0] + r * r * derivs[0]);
			// θ-component: (1/r·sinθ)∂(sinθ·Fθ)/∂θ = (cotθ/r)Fθ + (1/r)∂Fθ/∂θ
			div += inv_r_sin * (cos(theta) * vals[1] + sin(theta) * derivs[1]);
			// φ-component: (1/r·sinθ)∂Fφ/∂φ
			div += inv_r_sin * derivs[2];

			return div;
		}

		/// Computes divergence in spherical coordinates with structured result.
		///
		/// @param vectorField Vector function F(r, θ, φ) = (Fᵣ, Fθ, Fφ)
		/// @param pos         Position (r, θ, φ) with r > 0, 0 < θ < π
		/// @param config      Field operation configuration
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            EvaluationResult with divergence value and diagnostics
		static EvaluationResult<Real>
		DivSpherDetailed(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos,
		                 const FieldOperationConfig& config = {},
		                 SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			using ResultType = EvaluationResult<Real>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"DivSpher", config,
				[&](ResultType& result, int& func_evals) {
					result.value = DivSpher(vectorField, pos, policy);
					func_evals = 3 * 5 + 1; // 3 partial derivs + 1 function eval
				});
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                              DIVERGENCE - CYLINDRICAL
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes divergence in cylindrical coordinates (r, φ, z).
		///
		/// In cylindrical coordinates:
		///   ∇·F = (1/r)∂(r·Fᵣ)/∂r + (1/r)∂Fφ/∂φ + ∂Fz/∂z
		///
		/// @param vectorField Vector function F(r, φ, z) = (Fᵣ, Fφ, Fz)
		/// @param x           Position (r, φ, z) with r > 0
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            Scalar divergence value ∇·F
		///
		/// @warning Singular at r = 0 (z-axis singularity)
		/// @throws DomainError if at singularity and policy is Throw
		static Real DivCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x,
		                   SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			const Real r = x[0];
			
			VectorN<Real, 3> vals = vectorField(x);

			VectorN<Real, 3> derivs;
			for (int i = 0; i < 3; i++)
				derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);

			// Calculate inverse factor with singularity handling
			Real inv_r = Singularity::SafeInverseR(r, policy, "DivCyl 1/r");

			Real div = 0.0;
			// r-component: (1/r)∂(r·Fᵣ)/∂r = (1/r)Fᵣ + ∂Fᵣ/∂r
			div += inv_r * (vals[0] + r * derivs[0]);
			// φ-component: (1/r)∂Fφ/∂φ
			div += inv_r * derivs[1];
			// z-component: ∂Fz/∂z
			div += derivs[2];

			return div;
		}

		/// Computes divergence in cylindrical coordinates with structured result.
		///
		/// @param vectorField Vector function F(r, φ, z) = (Fᵣ, Fφ, Fz)
		/// @param pos         Position (r, φ, z) with r > 0
		/// @param config      Field operation configuration
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            EvaluationResult with divergence value and diagnostics
		static EvaluationResult<Real>
		DivCylDetailed(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos,
		               const FieldOperationConfig& config = {},
		               SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			using ResultType = EvaluationResult<Real>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"DivCyl", config,
				[&](ResultType& result, int& func_evals) {
					result.value = DivCyl(vectorField, pos, policy);
					func_evals = 3 * 5 + 1; // 3 partial derivs + 1 function eval
				});
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                                CURL - CARTESIAN
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes curl of vector field in 3D Cartesian coordinates.
		///
		/// The curl (∇×F) measures the local rotation or circulation of a vector field.
		///
		/// In Cartesian coordinates:
		///   ∇×F = (∂Fz/∂y - ∂Fy/∂z, ∂Fx/∂z - ∂Fz/∂x, ∂Fy/∂x - ∂Fx/∂y)
		///
		/// Physical interpretation:
		///   - |curl| > 0: field has rotational component at that point
		///   - curl = 0: field is irrotational (conservative)
		///   - Direction of curl indicates axis of rotation (right-hand rule)
		///
		/// @param vectorField Vector function F(x, y, z) = (Fx, Fy, Fz)
		/// @param pos         Position vector (x, y, z)
		/// @return            Curl vector ∇×F in Cartesian basis
		///
		/// @example
		///   // F(x,y,z) = (y, -x, 0) - rotation about z-axis
		///   VectorFunction<3> F([](auto& p){ return VectorN<Real,3>{p[1], -p[0], 0}; });
		///   auto curl = CurlCart(F, {1,1,1});  // Returns (0, 0, -2)
		///
		///   // Uniform field F = (1, 0, 0)
		///   VectorFunction<3> G([](auto& p){ return VectorN<Real,3>{1, 0, 0}; });
		///   auto curl = CurlCart(G, {1,1,1});  // Returns (0, 0, 0) - irrotational
		static Vec3Cart CurlCart(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			// ∂Fz/∂y and ∂Fy/∂z for x-component of curl
			Real dzdy = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
			Real dydz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

			// ∂Fx/∂z and ∂Fz/∂x for y-component of curl
			Real dxdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
			Real dzdx = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

			// ∂Fy/∂x and ∂Fx/∂y for z-component of curl
			Real dydx = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
			Real dxdy = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

			Vector3Cartesian curl{ dzdy - dydz, dxdz - dzdx, dydx - dxdy };

			return curl;
		}

		/// Computes curl of vector field in 3D Cartesian coordinates with structured result.
		///
		/// @param vectorField Vector function F(x, y, z) = (Fx, Fy, Fz)
		/// @param pos         Position vector (x, y, z)
		/// @param config      Field operation configuration
		/// @return            EvaluationResult with curl vector and per-component error estimates
		static EvaluationResult<Vec3Cart, Vec3Cart>
		CurlCartDetailed(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos,
		                 const FieldOperationConfig& config = {})
		{
			using ResultType = EvaluationResult<Vec3Cart, Vec3Cart>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"CurlCart", config,
				[&](ResultType& result, int& func_evals) {
					Real e_dzdy = 0, e_dydz = 0, e_dxdz = 0, e_dzdx = 0, e_dydx = 0, e_dxdy = 0;
					Real* ep = config.estimate_error ? &e_dzdy : nullptr;

					Real dzdy = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, ep);
					ep = config.estimate_error ? &e_dydz : nullptr;
					Real dydz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, ep);

					ep = config.estimate_error ? &e_dxdz : nullptr;
					Real dxdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, ep);
					ep = config.estimate_error ? &e_dzdx : nullptr;
					Real dzdx = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, ep);

					ep = config.estimate_error ? &e_dydx : nullptr;
					Real dydx = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, ep);
					ep = config.estimate_error ? &e_dxdy : nullptr;
					Real dxdy = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, ep);

					result.value = Vec3Cart{dzdy - dydz, dxdz - dzdx, dydx - dxdy};
					if (config.estimate_error) {
						result.error = Vec3Cart{
							std::abs(e_dzdy) + std::abs(e_dydz),
							std::abs(e_dxdz) + std::abs(e_dzdx),
							std::abs(e_dydx) + std::abs(e_dxdy)
						};
					}
					func_evals = 6 * 5; // 6 partial derivatives, ~5 evals each (NDer4)
				});
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                                CURL - SPHERICAL
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes curl of vector field in 3D spherical coordinates.
		///
		/// In spherical coordinates (r, θ, φ), the curl components are:
		///   (∇×F)ᵣ = (1/r·sinθ)[∂(sinθ·Fφ)/∂θ - ∂Fθ/∂φ]
		///   (∇×F)θ = (1/r)[(1/sinθ)∂Fᵣ/∂φ - ∂(r·Fφ)/∂r]
		///   (∇×F)φ = (1/r)[∂(r·Fθ)/∂r - ∂Fᵣ/∂θ]
		///
		/// @param vectorField Vector function F(r, θ, φ) = (Fᵣ, Fθ, Fφ)
		/// @param pos         Position (r, θ, φ) with r > 0, 0 < θ < π
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            Curl vector ∇×F in spherical basis (êᵣ, êθ, êφ)
		///
		/// @warning Singular at r = 0 and θ = 0, π
		/// @throws DomainError if at singularity and policy is Throw
		static Vec3Sph CurlSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos,
		                         SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			VectorN<Real, 3> vals = vectorField(pos);

			// Partial derivatives of each component
			Real dphidtheta = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);  // ∂Fφ/∂θ
			Real dthetadphi = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);  // ∂Fθ/∂φ

			Real drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);      // ∂Fᵣ/∂φ
			Real dphidr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);      // ∂Fφ/∂r

			Real dthetadr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);    // ∂Fθ/∂r
			Real drdtheta = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);    // ∂Fᵣ/∂θ

			const Real& r = pos[0];
			const Real& theta = pos[1];

			// Calculate inverse factors with singularity handling
			Real inv_r = Singularity::SafeInverseR(r, policy, "CurlSpher 1/r");
			Real inv_r_sin = Singularity::SafeInverseRSinTheta(r, theta, policy, "CurlSpher 1/(r·sinθ)");
			Real inv_sin = Singularity::SafeDivide(1.0, sin(theta), policy, "CurlSpher 1/sinθ");

			Vector3Spherical ret;
			// r-component: (1/r·sinθ)[cosθ·Fφ + sinθ·∂Fφ/∂θ - ∂Fθ/∂φ]
			ret[0] = inv_r_sin * (cos(theta) * vals[2] + sin(theta) * dphidtheta - dthetadphi);
			// θ-component: (1/r)[(1/sinθ)∂Fᵣ/∂φ - Fφ - r·∂Fφ/∂r]
			ret[1] = inv_r * (inv_sin * drdphi - vals[2] - r * dphidr);
			// φ-component: (1/r)[Fθ + r·∂Fθ/∂r - ∂Fᵣ/∂θ]
			ret[2] = inv_r * (vals[1] + r * dthetadr - drdtheta);

			return ret;
		}

		/// Computes curl of vector field in 3D spherical coordinates with structured result.
		///
		/// @param vectorField Vector function F(r, θ, φ) = (Fᵣ, Fθ, Fφ)
		/// @param pos         Position (r, θ, φ) with r > 0, 0 < θ < π
		/// @param config      Field operation configuration
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            EvaluationResult with curl vector and diagnostics
		static EvaluationResult<Vec3Sph, Vec3Sph>
		CurlSpherDetailed(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos,
		                  const FieldOperationConfig& config = {},
		                  SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			using ResultType = EvaluationResult<Vec3Sph, Vec3Sph>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"CurlSpher", config,
				[&](ResultType& result, int& func_evals) {
					result.value = CurlSpher(vectorField, pos, policy);
					func_evals = 6 * 5 + 1; // 6 partial derivs + 1 function eval
				});
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		//                                CURL - CYLINDRICAL
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Computes curl of vector field in 3D cylindrical coordinates.
		///
		/// In cylindrical coordinates (r, φ, z), the curl components are:
		///   (∇×F)ᵣ = (1/r)∂Fz/∂φ - ∂Fφ/∂z
		///   (∇×F)φ = ∂Fᵣ/∂z - ∂Fz/∂r
		///   (∇×F)z = (1/r)[Fφ + r·∂Fφ/∂r - ∂Fᵣ/∂φ]
		///
		/// @param vectorField Vector function F(r, φ, z) = (Fᵣ, Fφ, Fz)
		/// @param pos         Position (r, φ, z) with r > 0
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            Curl vector ∇×F in cylindrical basis (êᵣ, êφ, êz)
		///
		/// @warning Singular at r = 0 (z-axis singularity)
		/// @throws DomainError if at singularity and policy is Throw
		static Vec3Cyl CurlCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos,
		                       SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			const Real r = pos[0];
			
			VectorN<Real, 3> vals = vectorField(pos);

			// Partial derivatives
			Real dzdphi = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);  // ∂Fz/∂φ
			Real dphidz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);  // ∂Fφ/∂z

			Real drdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);    // ∂Fᵣ/∂z
			Real dzdr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);    // ∂Fz/∂r

			Real dphidr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);  // ∂Fφ/∂r
			Real drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);  // ∂Fᵣ/∂φ

			// Calculate inverse factor with singularity handling
			Real inv_r = Singularity::SafeInverseR(r, policy, "CurlCyl 1/r");

			// r-component: (1/r)∂Fz/∂φ - ∂Fφ/∂z
			// φ-component: ∂Fᵣ/∂z - ∂Fz/∂r
			// z-component: (1/r)[Fφ + r·∂Fφ/∂r - ∂Fᵣ/∂φ]
			Vector3Cylindrical ret{
				(inv_r * dzdphi - dphidz), 
				drdz - dzdr, 
				inv_r * (vals[1] + r * dphidr - drdphi)
			};

			return ret;
		}

		/// Computes curl of vector field in 3D cylindrical coordinates with structured result.
		///
		/// @param vectorField Vector function F(r, φ, z) = (Fᵣ, Fφ, Fz)
		/// @param pos         Position (r, φ, z) with r > 0
		/// @param config      Field operation configuration
		/// @param policy      Singularity handling policy (default: Throw)
		/// @return            EvaluationResult with curl vector and diagnostics
		static EvaluationResult<Vec3Cyl, Vec3Cyl>
		CurlCylDetailed(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos,
		                const FieldOperationConfig& config = {},
		                SingularityPolicy policy = Singularity::DEFAULT_POLICY)
		{
			using ResultType = EvaluationResult<Vec3Cyl, Vec3Cyl>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>(
				"CurlCyl", config,
				[&](ResultType& result, int& func_evals) {
					result.value = CurlCyl(vectorField, pos, policy);
					func_evals = 6 * 5 + 1; // 6 partial derivs + 1 function eval
				});
		}
	};
}
#endif
