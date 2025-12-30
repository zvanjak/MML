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
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
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

#include "base/Geometry.h"
#include "base/VectorN.h"
#include "base/VectorTypes.h"
#include "base/MatrixNM.h"
#include "base/Tensor.h"

#include "core/Derivation.h"
#include "core/MetricTensor.h"

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
		static VectorN<Real, 3> GradientSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos)
		{
			Vector3Spherical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

			ret[1] = ret[1] / pos[0];                    // θ-component: (1/r) ∂f/∂θ
			ret[2] = ret[2] / (pos[0] * sin(pos[1]));    // φ-component: (1/r·sinθ) ∂f/∂φ

			return ret;
		}

		/// Computes gradient in spherical coordinates with specified derivative accuracy.
		///
		/// @param scalarField Scalar function f(r, θ, φ)
		/// @param pos         Position vector (r, θ, φ)
		/// @param der_order   Derivative accuracy: 1, 2, 4, 6, or 8
		/// @return            Gradient vector in spherical basis
		/// @throws std::invalid_argument if der_order is not in {1, 2, 4, 6, 8}
		static Vec3Sph GradientSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos, 
																 int der_order)
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

			ret[1] = ret[1] / pos[0];                    // θ-component
			ret[2] = ret[2] / (pos[0] * sin(pos[1]));    // φ-component

			return ret;
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
		static Vec3Cyl GradientCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos)
		{
			Vector3Cylindrical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

			ret[1] = ret[1] / pos[0];    // φ-component: (1/r) ∂f/∂φ

			return ret;
		}

		/// Computes gradient in cylindrical coordinates with specified derivative accuracy.
		///
		/// @param scalarField Scalar function f(r, φ, z)
		/// @param pos         Position vector (r, φ, z)
		/// @param der_order   Derivative accuracy: 1, 2, 4, 6, or 8
		/// @return            Gradient vector in cylindrical basis
		/// @throws std::invalid_argument if der_order is not in {1, 2, 4, 6, 8}
		static Vec3Cyl GradientCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos, 
															 int der_order)
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
			ret[1] = ret[1] / pos[0];    // φ-component

			return ret;
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

		/// Computes Laplacian of scalar field in spherical coordinates.
		///
		/// In spherical coordinates (r, θ, φ):
		///   ∇²f = (1/r²)∂/∂r(r²∂f/∂r) + (1/r²sinθ)∂/∂θ(sinθ·∂f/∂θ) + (1/r²sin²θ)∂²f/∂φ²
		///
		/// @param scalarField Scalar function f(r, θ, φ)
		/// @param pos         Position (r, θ, φ) with r > 0, 0 < θ < π
		/// @return            Scalar Laplacian value ∇²f
		///
		/// @warning Singular at r = 0 and θ = 0, π (coordinate singularities)
		static Real LaplacianSpher(const IScalarFunction<3>& scalarField, const Vec3Sph& pos)
		{
			const Real r = pos.R();
			const Real theta = pos.Theta();
			// Note: phi = pos.Phi() is not needed directly in the Laplacian formula

			// ∇²f = ∂²f/∂r² + (2/r)∂f/∂r + (1/r²sinθ)[cosθ·∂f/∂θ + sinθ·∂²f/∂θ²] + (1/r²sin²θ)∂²f/∂φ²
			//     = ∂²f/∂r² + (2/r)∂f/∂r + (cotθ/r²)∂f/∂θ + (1/r²)∂²f/∂θ² + (1/r²sin²θ)∂²f/∂φ²
			
			// Radial term: ∂²f/∂r²
			Real d2f_dr2 = Derivation::DeriveSecPartial<3>(scalarField, 0, 0, pos, nullptr);
			
			// Radial correction: (2/r)∂f/∂r
			Real df_dr = Derivation::DerivePartial<3>(scalarField, 0, pos, nullptr);
			Real radial_correction = (2.0 / r) * df_dr;
			
			// Polar (θ) term: (cotθ/r²)∂f/∂θ + (1/r²)∂²f/∂θ²
			Real df_dtheta = Derivation::DerivePartial<3>(scalarField, 1, pos, nullptr);
			Real d2f_dtheta2 = Derivation::DeriveSecPartial<3>(scalarField, 1, 1, pos, nullptr);
			Real theta_term = (cos(theta) / (r * r * sin(theta))) * df_dtheta + (1.0 / (r * r)) * d2f_dtheta2;
			
			// Azimuthal (φ) term: (1/r²sin²θ)∂²f/∂φ²
			Real d2f_dphi2 = Derivation::DeriveSecPartial<3>(scalarField, 2, 2, pos, nullptr);
			Real phi_term = (1.0 / (r * r * sin(theta) * sin(theta))) * d2f_dphi2;

			return d2f_dr2 + radial_correction + theta_term + phi_term;
		}

		/// Computes Laplacian of scalar field in cylindrical coordinates.
		///
		/// In cylindrical coordinates (r, φ, z):
		///   ∇²f = (1/r)∂/∂r(r·∂f/∂r) + (1/r²)∂²f/∂φ² + ∂²f/∂z²
		///
		/// @param scalarField Scalar function f(r, φ, z)
		/// @param pos         Position (r, φ, z) with r > 0
		/// @return            Scalar Laplacian value ∇²f
		///
		/// @warning Singular at r = 0 (z-axis singularity)
		static Real LaplacianCyl(const IScalarFunction<3>& scalarField, const Vec3Cyl& pos)
		{
			const Real r = pos[0];

			// ∇²f = (1/r)∂/∂r(r·∂f/∂r) + (1/r²)∂²f/∂φ² + ∂²f/∂z²
			//     = (1/r)∂f/∂r + ∂²f/∂r² + (1/r²)∂²f/∂φ² + ∂²f/∂z²
			
			// Radial term: (1/r)∂f/∂r + ∂²f/∂r²
			Real df_dr = Derivation::DerivePartial<3>(scalarField, 0, pos, nullptr);
			Real d2f_dr2 = Derivation::DeriveSecPartial<3>(scalarField, 0, 0, pos, nullptr);
			Real radial_term = (1.0 / r) * df_dr + d2f_dr2;
			
			// Azimuthal term: (1/r²)∂²f/∂φ²
			Real d2f_dphi2 = Derivation::DeriveSecPartial<3>(scalarField, 1, 1, pos, nullptr);
			Real phi_term = (1.0 / (r * r)) * d2f_dphi2;
			
			// Axial term: ∂²f/∂z²
			Real d2f_dz2 = Derivation::DeriveSecPartial<3>(scalarField, 2, 2, pos, nullptr);

			return radial_term + phi_term + d2f_dz2;
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
		/// @return            Scalar divergence value ∇·F
		///
		/// @warning Singular at r = 0 and θ = 0, π (coordinate singularities)
		static Real DivSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x)
		{
			VectorN<Real, 3> vals = vectorField(x);

			VectorN<Real, 3> derivs;
			for (int i = 0; i < 3; i++)
				derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);

			Real div = 0.0;
			// r-component: (1/r²)∂(r²Fᵣ)/∂r = (2/r)Fᵣ + ∂Fᵣ/∂r
			div += 1 / (x[0] * x[0]) * (2 * x[0] * vals[0] + x[0] * x[0] * derivs[0]);
			// θ-component: (1/r·sinθ)∂(sinθ·Fθ)/∂θ = (cotθ/r)Fθ + (1/r)∂Fθ/∂θ
			div += 1 / (x[0] * sin(x[1])) * (cos(x[1]) * vals[1] + sin(x[1]) * derivs[1]);
			// φ-component: (1/r·sinθ)∂Fφ/∂φ
			div += 1 / (x[0] * sin(x[1])) * derivs[2];

			return div;
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
		/// @return            Scalar divergence value ∇·F
		///
		/// @warning Singular at r = 0 (z-axis singularity)
		static Real DivCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x)
		{
			VectorN<Real, 3> vals = vectorField(x);

			VectorN<Real, 3> derivs;
			for (int i = 0; i < 3; i++)
				derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);

			Real div = 0.0;
			// r-component: (1/r)∂(r·Fᵣ)/∂r = (1/r)Fᵣ + ∂Fᵣ/∂r
			div += 1 / x[0] * (vals[0] + x[0] * derivs[0]);
			// φ-component: (1/r)∂Fφ/∂φ
			div += 1 / x[0] * derivs[1];
			// z-component: ∂Fz/∂z
			div += derivs[2];

			return div;
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
		/// @return            Curl vector ∇×F in spherical basis (êᵣ, êθ, êφ)
		///
		/// @warning Singular at r = 0 and θ = 0, π
		static Vec3Sph CurlSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			VectorN<Real, 3> vals = vectorField(pos);

			// Partial derivatives of each component
			Real dphidtheta = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);  // ∂Fφ/∂θ
			Real dthetadphi = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);  // ∂Fθ/∂φ

			Real drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);      // ∂Fᵣ/∂φ
			Real dphidr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);      // ∂Fφ/∂r

			Real dthetadr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);    // ∂Fθ/∂r
			Real drdtheta = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);    // ∂Fᵣ/∂θ

			Vector3Spherical ret;
			const Real& r = pos[0];
			const Real& theta = pos[1];
			const Real& phi = pos[2];

			// r-component: (1/r·sinθ)[cosθ·Fφ + sinθ·∂Fφ/∂θ - ∂Fθ/∂φ]
			ret[0] = 1 / (r * sin(theta)) * (cos(theta) * vals[2] + sin(theta) * dphidtheta - dthetadphi);
			// θ-component: (1/r)[(1/sinθ)∂Fᵣ/∂φ - Fφ - r·∂Fφ/∂r]
			ret[1] = 1 / r * (1 / sin(theta) * drdphi - vals[2] - r * dphidr);
			// φ-component: (1/r)[Fθ + r·∂Fθ/∂r - ∂Fᵣ/∂θ]
			ret[2] = 1 / r * (vals[1] + r * dthetadr - drdtheta);

			return ret;
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
		/// @return            Curl vector ∇×F in cylindrical basis (êᵣ, êφ, êz)
		///
		/// @warning Singular at r = 0 (z-axis singularity)
		static Vec3Cyl CurlCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			VectorN<Real, 3> vals = vectorField(pos);

			// Partial derivatives
			Real dzdphi = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);  // ∂Fz/∂φ
			Real dphidz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);  // ∂Fφ/∂z

			Real drdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);    // ∂Fᵣ/∂z
			Real dzdr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);    // ∂Fz/∂r

			Real dphidr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);  // ∂Fφ/∂r
			Real drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);  // ∂Fᵣ/∂φ

			// r-component: (1/r)∂Fz/∂φ - ∂Fφ/∂z
			// φ-component: ∂Fᵣ/∂z - ∂Fz/∂r
			// z-component: (1/r)[Fφ + r·∂Fφ/∂r - ∂Fᵣ/∂φ]
			Vector3Cylindrical ret{
				(1 / pos[0] * dzdphi - dphidz), 
				drdz - dzdr, 
				1 / pos[0] * (vals[1] + pos[0] * dphidr - drdphi)
			};

			return ret;
		}
	};
}
#endif
