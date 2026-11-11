///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FieldAnalyzers.h                                                    ///
///  Description: Vector field analysis tools (streamlines, potential, flux)          ///
///               Divergence, curl, gradient visualizations                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FIELDS_ANALYZERS_H
#define MML_FIELDS_ANALYZERS_H

#include <mml/MMLBase.h>
#include <mml/MMLNumericValidation.h>
#include <algorithm>
#include <limits>
#include <numeric>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/Integration.h>
#include <mml/core/Integration/PathIntegration.h>
#include <mml/core/FunctionHelpers.h>

namespace MML {
	namespace FieldAnalyzerDetail {
		inline void ValidatePoint(const Vec3Cart& point) {
			for (int axis = 0; axis < 3; ++axis)
				ValidateFinite(point[axis], "Field analyzer coordinate");
		}

		inline FieldOperationConfig CheckedConfig() {
			FieldOperationConfig config;
			config.estimate_error = false;
			return config;
		}

		inline Real Magnitude(const Vec3Cart& value) {
			const Real magnitude = std::hypot(value.X(), value.Y(), value.Z());
			FieldOperationDetail::CheckFinite(magnitude, true);
			return magnitude;
		}

		// n subintervals produce (n + 1)^3 endpoint-inclusive samples; false stops traversal.
		template<class Visitor>
		bool SampleGrid(const Vec3Cart& lower, const Vec3Cart& upper, int subintervals, Visitor visit) {
			ValidatePoint(lower);
			ValidatePoint(upper);
			if (subintervals <= 0)
				throw ArgumentError("Field analyzer samplesPerDim must be positive");
			for (int axis = 0; axis < 3; ++axis)
				if (lower[axis] > upper[axis])
					throw ArgumentError("Field analyzer bounds must be ordered");
			const std::size_t nodes = static_cast<std::size_t>(subintervals) + 1;
			const auto limit = std::numeric_limits<std::size_t>::max();
			if (nodes > limit / nodes || nodes * nodes > limit / nodes)
				throw ResourceLimitError("Field analyzer grid sample count is unrepresentable");
			for (std::size_t i = 0; i < nodes; ++i)
				for (std::size_t j = 0; j < nodes; ++j)
					for (std::size_t k = 0; k < nodes; ++k) {
						const Vec3Cart pos{
							std::lerp(lower.X(), upper.X(), Real(i) / Real(subintervals)),
							std::lerp(lower.Y(), upper.Y(), Real(j) / Real(subintervals)),
							std::lerp(lower.Z(), upper.Z(), Real(k) / Real(subintervals))};
						if (!visit(pos)) return false;
					}
			return true;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                           SCALAR FIELD ANALYZER                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	// Analyzes properties of scalar fields in 3D space
	//
	// Borrows its field; the field must outlive the analyzer.
	// Region work is O((samplesPerDim + 1)^3 * derivative_cost).
	class ScalarFieldAnalyzer {
		const IScalarFunction<3>& _field;

	public:
		ScalarFieldAnalyzer(const IScalarFunction<3>& f)
			: _field(f) {}
		ScalarFieldAnalyzer(const IScalarFunction<3>&&) = delete;

		// Get the gradient at a point
		Vec3Cart GradientAt(const Vec3Cart& pos) const {
			FieldAnalyzerDetail::ValidatePoint(pos);
			return ScalarFieldOperations::GradientCartDetailed(_field, pos, FieldAnalyzerDetail::CheckedConfig()).value;
		}

		// Check if gradient is zero at a point (critical point)
		bool IsCriticalPoint(const Vec3Cart& pos, Real tol = Precision::LooseTolerance) const {
			ValidateTolerance(tol, "IsCriticalPoint");
			return FieldAnalyzerDetail::Magnitude(GradientAt(pos)) < tol;
		}

		// Find approximate critical points by sampling a region
		std::vector<Vec3Cart> FindCriticalPoints(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, int samplesPerDim = 10,
												 Real tol = Precision::CoarseTolerance) const {
			std::vector<Vec3Cart> criticalPoints;
			ValidateTolerance(tol, "FindCriticalPoints");
			FieldAnalyzerDetail::SampleGrid(minCorner, maxCorner, samplesPerDim, [&](const Vec3Cart& pos) {
				if (IsCriticalPoint(pos, tol)) criticalPoints.push_back(pos);
				return true;
			});
			return criticalPoints;
		}

		// Compute Laplacian at a point (∇²f = ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z²)
		Real LaplacianAt(const Vec3Cart& pos) const {
			FieldAnalyzerDetail::ValidatePoint(pos);
			return ScalarFieldOperations::LaplacianCartDetailed(_field, pos, FieldAnalyzerDetail::CheckedConfig()).value;
		}

		// Check if field satisfies Laplace equation (harmonic) in a region
		bool IsHarmonic(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, int samplesPerDim = 5, Real tol = Precision::CoarseTolerance) const {
			ValidateTolerance(tol, "IsHarmonic");
			return FieldAnalyzerDetail::SampleGrid(minCorner, maxCorner, samplesPerDim, [&](const Vec3Cart& pos) {
				return std::abs(LaplacianAt(pos)) <= tol;
			});
		}
	};

	/////////////////////////////////////////////////////////////////////////////////////
	///                           VECTOR FIELD ANALYZER                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	// Analyzes properties of vector fields in 3D space
	//
	// Borrows its field; the field must outlive the analyzer.
	// Region work is O((samplesPerDim + 1)^3 * derivative_cost).
	class VectorFieldAnalyzer {
		const IVectorFunction<3>& _field;

	public:
		VectorFieldAnalyzer(const IVectorFunction<3>& f)
			: _field(f) {}
		VectorFieldAnalyzer(const IVectorFunction<3>&&) = delete;

		// Get divergence at a point
		Real DivergenceAt(const Vec3Cart& pos) const {
			FieldAnalyzerDetail::ValidatePoint(pos);
			return VectorFieldOperations::DivCartDetailed(_field, pos, FieldAnalyzerDetail::CheckedConfig()).value;
		}

		// Get curl at a point
		Vec3Cart CurlAt(const Vec3Cart& pos) const {
			FieldAnalyzerDetail::ValidatePoint(pos);
			return VectorFieldOperations::CurlCartDetailed(_field, pos, FieldAnalyzerDetail::CheckedConfig()).value;
		}

		// Check if field is solenoidal (divergence-free) at a point
		bool IsSolenoidalAt(const Vec3Cart& pos, Real tol = Precision::LooseTolerance) const {
			ValidateTolerance(tol, "IsSolenoidalAt");
			return std::abs(DivergenceAt(pos)) < tol;
		}

		// Check if field is irrotational (curl-free) at a point
		bool IsIrrotationalAt(const Vec3Cart& pos, Real tol = Precision::LooseTolerance) const {
			ValidateTolerance(tol, "IsIrrotationalAt");
			return FieldAnalyzerDetail::Magnitude(CurlAt(pos)) < tol;
		}

		// Check if field is solenoidal throughout a region (sampling-based)
		// Solenoidal fields have ∇·F = 0 (no sources or sinks)
		bool IsSolenoidal(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, int samplesPerDim = 5, Real tol = Precision::CoarseTolerance) const {
			ValidateTolerance(tol, "IsSolenoidal");
			return FieldAnalyzerDetail::SampleGrid(minCorner, maxCorner, samplesPerDim, [&](const Vec3Cart& pos) {
				return IsSolenoidalAt(pos, tol);
			});
		}

		// Check if field is irrotational throughout a region (sampling-based)
		// Irrotational fields have ∇×F = 0
		bool IsIrrotational(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, int samplesPerDim = 5, Real tol = Precision::CoarseTolerance) const {
			ValidateTolerance(tol, "IsIrrotational");
			return FieldAnalyzerDetail::SampleGrid(minCorner, maxCorner, samplesPerDim, [&](const Vec3Cart& pos) {
				return IsIrrotationalAt(pos, tol);
			});
		}

		// Compatibility spelling for sampled irrotationality, not a global proof.
		// Global conservativeness additionally requires smoothness and suitable domain topology.
		bool IsConservative(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, int samplesPerDim = 5, Real tol = Precision::CoarseTolerance) const {
			return IsIrrotational(minCorner, maxCorner, samplesPerDim, tol);
		}

		// Get maximum divergence magnitude in a region
		Real MaxDivergence(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, int samplesPerDim = 5) const {
			Real maxDiv = 0.0;
			FieldAnalyzerDetail::SampleGrid(minCorner, maxCorner, samplesPerDim, [&](const Vec3Cart& pos) {
				maxDiv = std::max(maxDiv, std::abs(DivergenceAt(pos)));
				return true;
			});
			return maxDiv;
		}

		// Get maximum curl magnitude in a region
		Real MaxCurlMagnitude(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, int samplesPerDim = 5) const {
			Real maxCurl = 0.0;
			FieldAnalyzerDetail::SampleGrid(minCorner, maxCorner, samplesPerDim, [&](const Vec3Cart& pos) {
				maxCurl = std::max(maxCurl, FieldAnalyzerDetail::Magnitude(CurlAt(pos)));
				return true;
			});
			return maxCurl;
		}

		// Compute circulation around a closed circular path
		Real CirculationAroundCircle(const Vec3Cart& center, const Vec3Cart& normal, Real radius, int numPoints = 100) const {
			FieldAnalyzerDetail::ValidatePoint(center);
			FieldAnalyzerDetail::ValidatePoint(normal);
			ValidateFinite(radius, "CirculationAroundCircle radius");
			if (radius < 0 || numPoints <= 0)
				throw ArgumentError("CirculationAroundCircle requires nonnegative radius and positive numPoints");
			const Real normalScale = std::max({std::abs(normal.X()), std::abs(normal.Y()), std::abs(normal.Z())});
			if (normalScale == 0)
				throw ArgumentError("CirculationAroundCircle requires a nonzero normal");
			// Create orthonormal basis for the plane
			Vec3Cart n = (normal / normalScale).Normalized();
			Vec3Cart u, v;

			// Find a vector not parallel to n
			if (std::abs(n.X()) < 0.9) {
				u = VectorProduct(n, Vec3Cart(1, 0, 0)).Normalized();
			} else {
				u = VectorProduct(n, Vec3Cart(0, 1, 0)).Normalized();
			}
			v = VectorProduct(n, u);

			// Integrate F·dr around the circle
			Real circulation = 0.0;
			Real dtheta = MML::Real(2.0) * Constants::PI / numPoints;

			for (int i = 0; i < numPoints; ++i) {
				Real theta = i * dtheta;
				Real theta_next = (i + 1) * dtheta;

				// Current and next points on circle
				Vec3Cart p1 = center + u * (radius * std::cos(theta)) + v * (radius * std::sin(theta));
				Vec3Cart p2 = center + u * (radius * std::cos(theta_next)) + v * (radius * std::sin(theta_next));
				FieldOperationDetail::CheckFinite(p1, true);
				FieldOperationDetail::CheckFinite(p2, true);

				// Tangent direction (dr)
				Vec3Cart dr = p2 - p1;
				FieldOperationDetail::CheckFinite(dr, true);

				// Midpoint for field evaluation
				Vec3Cart midpoint{std::midpoint(p1.X(), p2.X()), std::midpoint(p1.Y(), p2.Y()), std::midpoint(p1.Z(), p2.Z())};
				VectorN<Real, 3> F = _field(midpoint);
				FieldOperationDetail::CheckFinite(F, true);

				// Add contribution F·dr
				circulation += F[0] * dr.X() + F[1] * dr.Y() + F[2] * dr.Z();
				FieldOperationDetail::CheckFinite(circulation, true);
			}

			return circulation;
		}
	};

} // end namespace MML

#endif // MML_FIELDS_ANALYZERS_H
