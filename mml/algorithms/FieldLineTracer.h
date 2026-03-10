///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FieldLineTracer.h                                                   ///
///  Description: Field line tracing for vector fields                                ///
///               Traces paths following field direction (not magnitude)              ///
///               Arc-length parameterized for uniform visual spacing                 ///
///                                                                                   ///
///  Applications:                                                                    ///
///    - Electric/magnetic field line visualization                                   ///
///    - Gravitational field lines                                                    ///
///    - Fluid flow streamlines (direction-only variant)                              ///
///    - Heat flow visualization                                                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FIELD_LINE_TRACER_H
#define MML_FIELD_LINE_TRACER_H

#include "MMLBase.h"
#include "interfaces/IFunction.h"
#include "base/Vector/VectorN.h"

#include <vector>
#include <cmath>
#include <functional>

namespace MML {

	/////////////////////////////////////////////////////////////////////////////////////
	///                              FIELD LINE TYPES                                 ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Result of tracing a single field line
	/// @tparam N Dimension (2 or 3)
	template<int N>
	struct FieldLine {
		std::vector<VectorN<Real, N>> points;  ///< Points along the field line
		Real totalLength = 0.0;                 ///< Total arc length of the line
		bool hitBoundary = false;               ///< True if terminated at boundary
		bool hitSingularity = false;            ///< True if terminated due to weak/zero field
		
		/// Number of points in this field line
		size_t size() const { return points.size(); }
		
		/// Check if field line is empty
		bool empty() const { return points.empty(); }
	};

	/// Convenience aliases
	using FieldLine2D = FieldLine<2>;
	using FieldLine3D = FieldLine<3>;

	/// @brief Bounding box for termination conditions
	template<int N>
	struct BoundingBox {
		VectorN<Real, N> minCorner;
		VectorN<Real, N> maxCorner;

		/// Check if a point is inside the bounding box
		bool Contains(const VectorN<Real, N>& point) const {
			for (int i = 0; i < N; ++i) {
				if (point[i] < minCorner[i] || point[i] > maxCorner[i]) {
					return false;
				}
			}
			return true;
		}

		/// Create a 2D bounding box
		static BoundingBox<2> Create2D(Real xMin, Real xMax, Real yMin, Real yMax) {
			return BoundingBox<2>{{xMin, yMin}, {xMax, yMax}};
		}

		/// Create a 3D bounding box
		static BoundingBox<3> Create3D(Real xMin, Real xMax, Real yMin, Real yMax, Real zMin, Real zMax) {
			return BoundingBox<3>{{xMin, yMin, zMin}, {xMax, yMax, zMax}};
		}
	};

	using BoundingBox2D = BoundingBox<2>;
	using BoundingBox3D = BoundingBox<3>;

	/////////////////////////////////////////////////////////////////////////////////////
	///                           FIELD LINE TRACER                                   ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Traces field lines through a vector field
	/// 
	/// Field lines follow the DIRECTION of the vector field, not the magnitude.
	/// Uses arc-length parameterization: dx/ds = F(x)/|F(x)|
	/// 
	/// This is different from streamlines which follow dx/dt = F(x) (velocity field).
	/// 
	/// Complexity:
	///   Trace:             O(maxPoints × 4) field evaluations (RK4, 4 f-evals/step)
	///   TraceFromSource:   O(numLines × maxPoints × 4) field evaluations
	///   TraceFromSeedPoints: O(numSeeds × maxPoints × 4) field evaluations
	/// 
	/// @tparam N Dimension (2 or 3)
	template<int N>
	class FieldLineTracer {
	public:
		/// @brief Configuration for field line tracing
		struct Config {
			Real stepSize = 0.05;           ///< Arc-length step size (ds)
			Real maxLength = 100.0;         ///< Maximum arc length before termination
			int maxPoints = 10000;          ///< Maximum number of points (safety limit)
			Real minFieldMagnitude = 1e-10; ///< Terminate if field gets too weak (singularity)
			bool traceForward = true;       ///< Trace in field direction
			bool traceBackward = false;     ///< Trace against field direction
		};

	private:
		Config _config;

	public:
		/// Default constructor with default configuration
		FieldLineTracer() = default;

		/// Constructor with custom configuration
		explicit FieldLineTracer(const Config& config) : _config(config) {}

		/// Set configuration
		void SetConfig(const Config& config) { _config = config; }

		/// Get current configuration
		const Config& GetConfig() const { return _config; }

		//===================================================================================
		// Single field line tracing
		//===================================================================================

		/// @brief Trace a single field line from a starting point
		/// @param field The vector field to trace through
		/// @param startPoint Starting position
		/// @param bounds Bounding box for termination
		/// @return FieldLine containing the traced points
		FieldLine<N> Trace(const IVectorFunction<N>& field,
		                   const VectorN<Real, N>& startPoint,
		                   const BoundingBox<N>& bounds) const
		{
			FieldLine<N> result;
			
			// Check if start point is in bounds
			if (!bounds.Contains(startPoint)) {
				return result;
			}

			// Trace forward
			if (_config.traceForward) {
				auto forwardLine = TraceDirection(field, startPoint, bounds, 1.0);
				result.points = std::move(forwardLine.points);
				result.totalLength = forwardLine.totalLength;
				result.hitBoundary = forwardLine.hitBoundary;
				result.hitSingularity = forwardLine.hitSingularity;
			}

			// Trace backward (prepend points in reverse order)
			if (_config.traceBackward) {
				auto backwardLine = TraceDirection(field, startPoint, bounds, -1.0);
				
				if (!backwardLine.empty()) {
					// Reverse backward points and prepend (skip first point to avoid duplicate)
					std::vector<VectorN<Real, N>> combined;
					combined.reserve(backwardLine.size() + result.size());
					
					for (auto it = backwardLine.points.rbegin(); it != backwardLine.points.rend(); ++it) {
						combined.push_back(*it);
					}
					
					// Skip first point of forward trace if we have backward points (it's the start point)
					if (!result.empty() && !combined.empty()) {
						for (size_t i = 1; i < result.size(); ++i) {
							combined.push_back(result.points[i]);
						}
					} else {
						for (const auto& p : result.points) {
							combined.push_back(p);
						}
					}
					
					result.points = std::move(combined);
					result.totalLength += backwardLine.totalLength;
					result.hitBoundary = result.hitBoundary || backwardLine.hitBoundary;
					result.hitSingularity = result.hitSingularity || backwardLine.hitSingularity;
				}
			}

			return result;
		}

		//===================================================================================
		// Multiple field lines from source point
		//===================================================================================

		/// @brief Trace multiple field lines radiating from a source point
		/// 
		/// Creates seed points uniformly distributed on a circle (2D) or sphere (3D)
		/// around the source point, then traces field lines from each seed.
		/// 
		/// @param field The vector field to trace through
		/// @param sourcePoint Center point (e.g., charge or mass location)
		/// @param numLines Number of field lines to generate
		/// @param seedRadius Distance from source to start tracing (avoid singularity)
		/// @param bounds Bounding box for termination
		/// @return Vector of field lines
		std::vector<FieldLine<N>> TraceFromSource(const IVectorFunction<N>& field,
		                                           const VectorN<Real, N>& sourcePoint,
		                                           int numLines,
		                                           Real seedRadius,
		                                           const BoundingBox<N>& bounds) const
		{
			std::vector<VectorN<Real, N>> seeds = GenerateSeedPoints(sourcePoint, numLines, seedRadius);
			return TraceFromSeedPoints(field, seeds, bounds);
		}

		/// @brief Trace field lines from a list of seed points
		/// @param field The vector field to trace through
		/// @param seedPoints Starting positions for each field line
		/// @param bounds Bounding box for termination
		/// @return Vector of field lines (one per seed point)
		std::vector<FieldLine<N>> TraceFromSeedPoints(const IVectorFunction<N>& field,
		                                               const std::vector<VectorN<Real, N>>& seedPoints,
		                                               const BoundingBox<N>& bounds) const
		{
			std::vector<FieldLine<N>> lines;
			lines.reserve(seedPoints.size());

			for (const auto& seed : seedPoints) {
				lines.push_back(Trace(field, seed, bounds));
			}

			return lines;
		}

	private:
		//===================================================================================
		// Internal tracing implementation
		//===================================================================================

		/// @brief Trace in one direction using RK4
		/// @param direction +1.0 for forward, -1.0 for backward
		FieldLine<N> TraceDirection(const IVectorFunction<N>& field,
		                            const VectorN<Real, N>& startPoint,
		                            const BoundingBox<N>& bounds,
		                            Real direction) const
		{
			FieldLine<N> result;
			result.points.push_back(startPoint);

			VectorN<Real, N> current = startPoint;
			Real arcLength = 0.0;
			const Real ds = _config.stepSize * direction;

			while (arcLength < _config.maxLength && 
			       static_cast<int>(result.points.size()) < _config.maxPoints) {
				
				// RK4 step with normalized field direction
				VectorN<Real, N> next;
				bool stepOk = RK4Step(field, current, ds, next);

				if (!stepOk) {
					result.hitSingularity = true;
					break;
				}

				// Check bounds
				if (!bounds.Contains(next)) {
					result.hitBoundary = true;
					break;
				}

				result.points.push_back(next);
				arcLength += std::abs(ds);
				current = next;
			}

			result.totalLength = arcLength;
			return result;
		}

		/// @brief Single RK4 step with normalized field direction
		/// @return false if field is too weak (singularity)
		bool RK4Step(const IVectorFunction<N>& field,
		             const VectorN<Real, N>& x,
		             Real ds,
		             VectorN<Real, N>& xNext) const
		{
			// k1 = normalized field at x
			VectorN<Real, N> k1;
			if (!GetNormalizedField(field, x, k1)) return false;

			// k2 = normalized field at x + ds/2 * k1
			VectorN<Real, N> x2 = x + k1 * (ds * 0.5);
			VectorN<Real, N> k2;
			if (!GetNormalizedField(field, x2, k2)) return false;

			// k3 = normalized field at x + ds/2 * k2
			VectorN<Real, N> x3 = x + k2 * (ds * 0.5);
			VectorN<Real, N> k3;
			if (!GetNormalizedField(field, x3, k3)) return false;

			// k4 = normalized field at x + ds * k3
			VectorN<Real, N> x4 = x + k3 * ds;
			VectorN<Real, N> k4;
			if (!GetNormalizedField(field, x4, k4)) return false;

			// Combine: x_next = x + ds * (k1 + 2*k2 + 2*k3 + k4) / 6
			xNext = x + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (ds / 6.0);
			return true;
		}

		/// @brief Get normalized field direction at a point
		/// @return false if field magnitude is below threshold
		bool GetNormalizedField(const IVectorFunction<N>& field,
		                        const VectorN<Real, N>& x,
		                        VectorN<Real, N>& normalized) const
		{
			VectorN<Real, N> F = field(x);
			Real mag = F.NormL2();

			if (mag < _config.minFieldMagnitude) {
				return false;  // Singularity or zero field
			}

			normalized = F / mag;
			return true;
		}

		//===================================================================================
		// Seed point generation
		//===================================================================================

		/// @brief Generate seed points uniformly distributed around a center
		std::vector<VectorN<Real, N>> GenerateSeedPoints(const VectorN<Real, N>& center,
		                                                  int numPoints,
		                                                  Real radius) const
		{
			if constexpr (N == 2) {
				return GenerateSeedPoints2D(center, numPoints, radius);
			} else {
				return GenerateSeedPoints3D(center, numPoints, radius);
			}
		}

		/// @brief Generate points on a circle (2D)
		std::vector<VectorN<Real, 2>> GenerateSeedPoints2D(const VectorN<Real, 2>& center,
		                                                    int numPoints,
		                                                    Real radius) const
		{
			std::vector<VectorN<Real, 2>> seeds;
			seeds.reserve(numPoints);

			const Real dTheta = 2.0 * Constants::PI / numPoints;
			for (int i = 0; i < numPoints; ++i) {
				Real theta = i * dTheta;
				VectorN<Real, 2> seed{
					center[0] + radius * std::cos(theta),
					center[1] + radius * std::sin(theta)
				};
				seeds.push_back(seed);
			}

			return seeds;
		}

		/// @brief Generate points approximately uniformly on a sphere (3D)
		/// Uses Fibonacci sphere algorithm for good distribution
		std::vector<VectorN<Real, 3>> GenerateSeedPoints3D(const VectorN<Real, 3>& center,
		                                                    int numPoints,
		                                                    Real radius) const
		{
			std::vector<VectorN<Real, 3>> seeds;
			seeds.reserve(numPoints);

			// Fibonacci sphere for uniform distribution
			const Real goldenRatio = (1.0 + std::sqrt(5.0)) / 2.0;
			const Real angleIncrement = 2.0 * Constants::PI / goldenRatio;

			for (int i = 0; i < numPoints; ++i) {
				Real t = static_cast<Real>(i) / static_cast<Real>(numPoints - 1);
				Real phi = std::acos(1.0 - 2.0 * t);  // Polar angle [0, π]
				Real theta = angleIncrement * i;      // Azimuthal angle

				VectorN<Real, 3> seed{
					center[0] + radius * std::sin(phi) * std::cos(theta),
					center[1] + radius * std::sin(phi) * std::sin(theta),
					center[2] + radius * std::cos(phi)
				};
				seeds.push_back(seed);
			}

			return seeds;
		}
	};

	/////////////////////////////////////////////////////////////////////////////////////
	///                          CONVENIENCE ALIASES                                  ///
	/////////////////////////////////////////////////////////////////////////////////////

	using FieldLineTracer2D = FieldLineTracer<2>;
	using FieldLineTracer3D = FieldLineTracer<3>;

} // namespace MML

#endif // MML_FIELD_LINE_TRACER_H
