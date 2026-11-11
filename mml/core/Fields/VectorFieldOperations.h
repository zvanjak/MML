///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorFieldOperations.h                                             ///
///  Description: Divergence and curl operations for vector fields                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_VECTOR_FIELD_OPERATIONS_H
#define MML_VECTOR_FIELD_OPERATIONS_H

#include <mml/algorithms/MatrixSingularity.h>
#include <mml/core/Fields/FieldOperationsCommon.h>
#include <mml/core/DifferentialGeometry/CoordinateSingularities.h>

namespace MML
{
	namespace FieldOperationDetail {
		template<int N>
		VectorN<Real, N> ConvertOrthogonalComponents(const VectorN<Real, N>& components, const VectorN<Real, N>& pos,
																								 const MetricTensorField<N>& metric, bool toPhysical) {
			const char* operation = toPhysical ? "ToPhysical" : "ToContravariant";
			if (!metric.IsRiemannian() || metric.getNumCovar() != 2 || metric.getNumContravar() != 0)
				throw DomainError(std::string(operation) + ": a covariant Riemannian metric is required");
			for (int i = 0; i < N; ++i)
				if (!std::isfinite(pos[i]) || !std::isfinite(components[i]))
					throw DomainError(std::string(operation) + ": position and vector components must be finite");

			const auto g = metric.GetCovariantMetric(pos);
			VectorN<Real, N> scales;
			for (int i = 0; i < N; ++i) {
				if (!std::isfinite(g(i, i)) || g(i, i) <= 0)
					throw DomainError(std::string(operation) + ": metric diagonal must be finite and positive");
				scales[i] = std::sqrt(g(i, i));
				if (IsUnsafeDenominator(scales[i]))
					throw DomainError(std::string(operation) + ": metric scale factor is too small for safe conversion");
			}
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j) {
					if (!std::isfinite(g(i, j)))
						throw DomainError(std::string(operation) + ": metric components must be finite");
					// Normalize by both basis lengths without forming their potentially overflowing product.
					if (i != j && !IsNearlyZero(g(i, j) / scales[i] / scales[j], Precision::OrthogonalityTolerance))
						throw DomainError(std::string(operation) + ": nonorthogonal coordinates require an explicit orthonormal frame");
				}

			VectorN<Real, N> result;
			for (int i = 0; i < N; ++i) {
				result[i] = toPhysical ? components[i] * scales[i] : components[i] / scales[i];
				if (!std::isfinite(result[i]))
					throw DomainError(std::string(operation) + ": converted component is not finite");
			}
			return result;
		}
	} // namespace FieldOperationDetail

	///////////////////////////////////////////////////////////////////////////////////////////
	/// VectorFieldOperations - Differential operations and vector component conversions
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
		/// Convert contravariant coordinate-basis components v^i to physical components h_i v^i.
		/// Requires a finite, positive diagonal covariant Riemannian metric, with h_i = sqrt(g_ii).
		/// Off-diagonal correlations must be below Precision::OrthogonalityTolerance.
		/// Throws DomainError for invalid inputs, unsafe scales or nonfinite output. This does not lower indices.
		template<int N>
		VectorN<Real, N> ToPhysical(const VectorN<Real, N>& components, const VectorN<Real, N>& pos, const MetricTensorField<N>& metric) {
			return FieldOperationDetail::ConvertOrthogonalComponents(components, pos, metric, true);
		}

		/// Convert physical orthonormal-basis components to contravariant components v^i = v^hat_i / h_i.
		/// Has the same metric and finiteness requirements as ToPhysical; this does not raise indices.
		template<int N>
		VectorN<Real, N> ToContravariant(const VectorN<Real, N>& components, const VectorN<Real, N>& pos, const MetricTensorField<N>& metric) {
			return FieldOperationDetail::ConvertOrthogonalComponents(components, pos, metric, false);
		}

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
		/// Input components are contravariant coordinate-basis components, not physical components.
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
				div += Differential::VectorPartial<N>(vectorField, i, i, pos);

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
		/// @return                  Divergence value with status and derivative diagnostics
		///
		/// @note For Cartesian coordinates, √g = 1 everywhere and this reduces to DivCart.
		/// @note For spherical coordinates, √g = r²sinθ.
		/// @note For cylindrical coordinates, √g = r.
		/// @see Divergence, DivCart
		/// Determinant divergence with status reporting and derivative diagnostics.
		/// Degenerate center/stencil metrics are InvalidInput with ConvertToStatus, otherwise DomainError.
		/// Normalize volume ratios in log space to avoid determinant underflow/overflow.
		template<int N>
		static EvaluationResult<Real> DivergenceViaDetDetailed(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos,
																													 const MetricTensorField<N>& metricTensorField,
																													 const FieldOperationConfig& config = {}) {
			using ResultType = EvaluationResult<Real>;
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>("DivergenceViaDet", config, [&](ResultType& result, int& func_evals) {
				if (metricTensorField.getNumContravar() != 0 || metricTensorField.getNumCovar() != 2)
					throw DomainError("DivergenceViaDet requires covariant metric storage");
				const auto* vectorDerivatives = dynamic_cast<const IVectorDerivativeProvider<N>*>(&vectorField);
				const auto* metricDerivatives = dynamic_cast<const IMetricJetProvider<N>*>(&metricTensorField);
				if (config.derivative_order == 0 && vectorDerivatives && metricDerivatives) {
					Singularity::RequireFiniteCoordinates<N>(pos, "DivergenceViaDet");
					const auto jet = metricDerivatives->EvaluateMetricJet(pos);
					MatrixNM<Real, N, N> inverse;
					try {
						inverse = Singularity::InverseMatrix(jet.value, SingularityPolicy::Throw);
					} catch (const SingularMatrixError&) {
						throw DomainError("DivergenceViaDet: degenerate metric");
					}
					const auto field = vectorField(pos);
					++func_evals;
					Real div = 0;
					for (int i = 0; i < N; ++i) {
						Singularity::RequireFiniteMatrix(jet.first[i], "DivergenceViaDet metric derivative");
						Real logVolumePartial = 0;
						for (int a = 0; a < N; ++a)
							for (int b = 0; b < N; ++b)
								logVolumePartial += Real(0.5) * inverse(a, b) * jet.first[i](b, a);
						div += vectorDerivatives->EvaluateComponentPartial(pos, i, i) + field[i] * logVolumePartial;
						++func_evals;
					}
					if (!std::isfinite(div))
						throw NumericalMethodError("DivergenceViaDet: nonfinite provider result");
					result.value = div;
					FieldOperationDetail::MarkProviderResult(result, func_evals, func_evals);
					return;
				}
				auto logVolume = [&](const VectorN<Real, N>& x) {
					Singularity::RequireFiniteCoordinates<N>(x, "DivergenceViaDet");
					try {
						return Real(0.5) * Singularity::LogAbsDeterminant(metricTensorField.GetCovariantMetric(x));
					} catch (const SingularMatrixError& ex) {
						throw DomainError(std::string("DivergenceViaDet: degenerate metric: ") + ex.what());
					}
				};
				const Real reference = logVolume(pos);
				Real div = 0, error = 0;
				for (int i = 0; i < N; ++i) {
					ScalarFunctionFromStdFunc<N> weightedComponent([&](const VectorN<Real, N>& x) {
						const Real weight = std::exp(logVolume(x) - reference);
						if (!std::isfinite(weight) || weight == 0)
							throw NumericalMethodError("DivergenceViaDet: volume ratio underflow or overflow");
						++func_evals;
						const Real value = weight * vectorField(x)[i];
						if (config.check_finite && !std::isfinite(value))
							throw NumericalMethodError("DivergenceViaDet: nonfinite volume-weighted field value");
						return value;
					});
					Real componentError = 0;
					Real* estimate = config.estimate_error ? &componentError : nullptr;
					switch (config.derivative_order) {
					case 0:
					case 4:
						div += Derivation::NDer4Partial(weightedComponent, i, pos, estimate);
						break;
					case 1:
						div += Derivation::NDer1Partial(weightedComponent, i, pos, estimate);
						break;
					case 2:
						div += Derivation::NDer2Partial(weightedComponent, i, pos, estimate);
						break;
					case 6:
						div += Derivation::NDer6Partial(weightedComponent, i, pos, estimate);
						break;
					case 8:
						div += Derivation::NDer8Partial(weightedComponent, i, pos, estimate);
						break;
					default:
						throw ArgumentError("DivergenceViaDet: derivative_order must be 0, 1, 2, 4, 6, or 8");
					}
					error += std::abs(componentError);
				}
				if (!std::isfinite(div) || !std::isfinite(error))
					throw NumericalMethodError("DivergenceViaDet: nonfinite derivative result");
				result.value = div;
				if (config.estimate_error)
					result.error = error;
			});
		}

		/// Scalar determinant divergence; throws on degenerate metrics instead of returning zero.
		template<int N>
		static Real DivergenceViaDet(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos,
																 const MetricTensorField<N>& metricTensorField) {
			FieldOperationConfig config;
			config.estimate_error = false;
			return DivergenceViaDetDetailed(vectorField, pos, metricTensorField, config).value;
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
				div += Differential::VectorPartial<N>(vectorField, i, i, pos);

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
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>("DivCart", config, [&](ResultType& result, int& func_evals) {
				if (FieldOperationDetail::UsesProvidedDerivatives(vectorField, config)) {
					result.value = DivCart(vectorField, pos);
					FieldOperationDetail::MarkProviderResult(result, N, func_evals);
					return;
				}
				Real div = 0.0;
				Real total_error = 0.0;
				for (int i = 0; i < N; i++) {
					Real err = 0.0;
					div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, config.estimate_error ? &err : nullptr);
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
		/// @note F uses physical components in the local orthonormal basis, as returned by GradientSpher.
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
				derivs[i] = Differential::VectorPartial<3>(vectorField, i, i, x);

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
					func_evals = dynamic_cast<const IVectorDerivativeProvider<3>*>(&vectorField) ? 4 : 3 * 5 + 1;
					if (dynamic_cast<const IVectorDerivativeProvider<3>*>(&vectorField))
						FieldOperationDetail::MarkProviderResult(result, func_evals, func_evals);
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
		/// @note F uses physical components in the local orthonormal basis, as returned by GradientCyl.
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
				derivs[i] = Differential::VectorPartial<3>(vectorField, i, i, x);

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
					func_evals = dynamic_cast<const IVectorDerivativeProvider<3>*>(&vectorField) ? 4 : 3 * 5 + 1;
					if (dynamic_cast<const IVectorDerivativeProvider<3>*>(&vectorField))
						FieldOperationDetail::MarkProviderResult(result, func_evals, func_evals);
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
			Real dzdy = Differential::VectorPartial<3>(vectorField, 2, 1, pos);
			Real dydz = Differential::VectorPartial<3>(vectorField, 1, 2, pos);

			// ∂Fx/∂z and ∂Fz/∂x for y-component of curl
			Real dxdz = Differential::VectorPartial<3>(vectorField, 0, 2, pos);
			Real dzdx = Differential::VectorPartial<3>(vectorField, 2, 0, pos);

			// ∂Fy/∂x and ∂Fx/∂y for z-component of curl
			Real dydx = Differential::VectorPartial<3>(vectorField, 1, 0, pos);
			Real dxdy = Differential::VectorPartial<3>(vectorField, 0, 1, pos);

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
			return FieldOperationDetail::ExecuteFieldDetailed<ResultType>("CurlCart", config, [&](ResultType& result, int& func_evals) {
				if (FieldOperationDetail::UsesProvidedDerivatives(vectorField, config)) {
					result.value = CurlCart(vectorField, pos);
					FieldOperationDetail::MarkProviderResult(result, 6, func_evals);
					return;
				}
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
					result.error =
						Vec3Cart{std::abs(e_dzdy) + std::abs(e_dydz), std::abs(e_dxdz) + std::abs(e_dzdx), std::abs(e_dydx) + std::abs(e_dxdy)};
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
			Real dphidtheta = Differential::VectorPartial<3>(vectorField, 2, 1, pos); // ∂Fφ/∂θ
			Real dthetadphi = Differential::VectorPartial<3>(vectorField, 1, 2, pos); // ∂Fθ/∂φ

			Real drdphi = Differential::VectorPartial<3>(vectorField, 0, 2, pos); // ∂Fᵣ/∂φ
			Real dphidr = Differential::VectorPartial<3>(vectorField, 2, 0, pos); // ∂Fφ/∂r

			Real dthetadr = Differential::VectorPartial<3>(vectorField, 1, 0, pos); // ∂Fθ/∂r
			Real drdtheta = Differential::VectorPartial<3>(vectorField, 0, 1, pos); // ∂Fᵣ/∂θ

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
					func_evals = dynamic_cast<const IVectorDerivativeProvider<3>*>(&vectorField) ? 7 : 6 * 5 + 1;
					if (dynamic_cast<const IVectorDerivativeProvider<3>*>(&vectorField))
						FieldOperationDetail::MarkProviderResult(result, func_evals, func_evals);
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
			Real dzdphi = Differential::VectorPartial<3>(vectorField, 2, 1, pos); // ∂Fz/∂φ
			Real dphidz = Differential::VectorPartial<3>(vectorField, 1, 2, pos); // ∂Fφ/∂z

			Real drdz = Differential::VectorPartial<3>(vectorField, 0, 2, pos); // ∂Fᵣ/∂z
			Real dzdr = Differential::VectorPartial<3>(vectorField, 2, 0, pos); // ∂Fz/∂r

			Real dphidr = Differential::VectorPartial<3>(vectorField, 1, 0, pos); // ∂Fφ/∂r
			Real drdphi = Differential::VectorPartial<3>(vectorField, 0, 1, pos); // ∂Fᵣ/∂φ

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
					func_evals = dynamic_cast<const IVectorDerivativeProvider<3>*>(&vectorField) ? 7 : 6 * 5 + 1;
					if (dynamic_cast<const IVectorDerivativeProvider<3>*>(&vectorField))
						FieldOperationDetail::MarkProviderResult(result, func_evals, func_evals);
				});
		}
	};
}

#endif