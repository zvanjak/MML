///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DBodyBase.h                                                ///
///  Description: Base interfaces and classes for 3D solid bodies                     ///
///               IBody, ISolidBodyWithBoundary, mesh body base classes               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry3DBodyBase.h
/// @brief Base interfaces and abstract classes for 3D solid bodies.
/// - IBody: Abstract interface for all solid bodies
/// - ISolidBodyWithBoundary: Bodies defined by boundary functions
/// - BodyWithTriangleSurfaces: Base for triangulated mesh bodies
/// - BodyWithRectSurfaces: Base for quad mesh bodies
/// - ComposedSolidSurfaces3D: Composite solid from multiple bodies
/// @see Geometry3DBodies.h for the aggregate header

#if !defined MML_GEOMETRY_3D_BODY_BASE_H
#define MML_GEOMETRY_3D_BODY_BASE_H

#include <cmath>
#include <algorithm>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DBounding.h>

namespace MML {
	/// Resolution controls for estimates, not certified quadrature error bounds.
	struct BoundaryGeometrySettings {
		int integrationIntervals = Precision::BoundaryIntegrationIntervals;
		Real derivativeStep = Precision::CbrtEps; ///< Step in normalized [0,1] parameters.
	};

	namespace Geometry3DBodyDetail {
		inline Real FiniteValue(Real value) {
			if (!std::isfinite(value))
				throw GeometryError("Function-bounded geometry requires finite values and representable arithmetic");
			return value;
		}

		inline void ValidateBounds(Real lower, Real upper) {
			FiniteValue(lower);
			FiniteValue(upper);
			if (lower > upper)
				throw GeometryError("Function-bounded geometry requires ordered boundary values");
		}

		inline Real IntervalWidth(Real lower, Real upper) {
			ValidateBounds(lower, upper);
			return FiniteValue(upper - lower);
		}

		inline void ValidateSettings(const BoundaryGeometrySettings& settings) {
			if (settings.integrationIntervals < 2 || settings.integrationIntervals % 2 != 0)
				throw GeometryError("Boundary integration requires a positive even interval count");
			if (!std::isfinite(settings.derivativeStep) || settings.derivativeStep <= 0
				|| settings.derivativeStep > Precision::BoundaryDerivativeMaximumStep)
				throw GeometryError("Boundary derivative step must be finite and in (0, 1/2]");
		}

		template<typename Function>
		Real IntegrateSimpson(Function&& function, Real lower, Real upper, int intervals) {
			const Real width = IntervalWidth(lower, upper);
			if (lower == upper)
				return 0;
			Real left = FiniteValue(function(lower));
			Real average = 0;
			for (int panel = 0; panel < intervals / 2; ++panel) {
				const Real middle = FiniteValue(function(std::lerp(lower, upper, Real(2 * panel + 1) / intervals)));
				const Real right = FiniteValue(function(std::lerp(lower, upper, Real(2 * panel + 2) / intervals)));
				// Average Simpson panels before scaling, avoiding a large weighted sum.
				const Real panelAverage = std::lerp(std::lerp(left, right, Real(1) / 2), middle, Real(2) / 3);
				average = panel == 0 ? panelAverage : std::lerp(average, panelAverage, Real(1) / (panel + 1));
				left = right;
			}
			return FiniteValue(average * width);
		}

		template<typename Function>
		Real UnitDerivative(Function&& function, Real parameter, Real step) {
			Real first, second;
			if (parameter >= step && parameter <= 1 - step) {
				first = parameter - step;
				second = parameter + step;
			} else if (parameter <= Real(1) / 2) {
				step = std::min(step, (1 - parameter) / 2);
				first = parameter + step;
				second = std::min(Real(1), parameter + 2 * step);
			} else {
				step = std::min(step, parameter / 2);
				first = parameter - step;
				second = std::max(Real(0), parameter - 2 * step);
			}
			first = std::clamp(first, Real(0), Real(1));
			second = std::clamp(second, Real(0), Real(1));
			const Real delta1 = first - parameter, delta2 = second - parameter;
			if (delta1 == 0 || delta2 == 0 || delta1 == delta2)
				throw GeometryError("Boundary derivative step does not resolve distinct parameter nodes");
			const Real center = FiniteValue(function(parameter));
			const Real slope1 = FiniteValue((FiniteValue(function(first)) - center) / delta1);
			const Real slope2 = FiniteValue((FiniteValue(function(second)) - center) / delta2);
			// Three-point Lagrange derivative using actual, possibly asymmetric node spacing.
			return FiniteValue(std::lerp(slope1, slope2, -delta1 / (delta2 - delta1)));
		}

		inline std::vector<Triangle3D> Triangulate(const std::vector<RectSurface3D>& surfaces) {
			std::vector<Triangle3D> triangles;
			triangles.reserve(2 * surfaces.size());
			for (const auto& surface : surfaces) {
				triangles.emplace_back(surface._pnt1, surface._pnt2, surface._pnt3);
				triangles.emplace_back(surface._pnt1, surface._pnt3, surface._pnt4);
			}
			return triangles;
		}

		inline std::vector<Pnt3Cart> Vertices(const std::vector<Triangle3D>& triangles) {
			std::vector<Pnt3Cart> vertices;
			vertices.reserve(3 * triangles.size());
			for (const auto& triangle : triangles) {
				vertices.push_back(triangle.Pnt1());
				vertices.push_back(triangle.Pnt2());
				vertices.push_back(triangle.Pnt3());
			}
			return vertices;
		}

		inline Real SignedTetrahedronVolume(const Pnt3Cart& reference, const Triangle3D& triangle) {
			const Vec3Cart a(reference, triangle.Pnt1());
			const Vec3Cart b(reference, triangle.Pnt2());
			const Vec3Cart c(reference, triangle.Pnt3());
			return ScalarProduct(a, VectorProduct(b, c)) / MML::Real(6.0);
		}

		inline Real Volume(const std::vector<Triangle3D>& triangles) {
			if (triangles.empty())
				return 0.0;

			const Pnt3Cart reference = triangles.front().Pnt1();
			Real signedVolume = 0.0;
			for (const auto& triangle : triangles)
				signedVolume += SignedTetrahedronVolume(reference, triangle);
			return std::abs(signedVolume);
		}

		inline Real SurfaceArea(const std::vector<Triangle3D>& triangles) {
			Real area = 0.0;
			for (const auto& triangle : triangles)
				area += triangle.Area();
			return area;
		}

		inline Pnt3Cart Center(const std::vector<Triangle3D>& triangles) {
			if (triangles.empty())
				throw GeometryError("Cannot calculate the center of an empty surface mesh");

			const Pnt3Cart reference = triangles.front().Pnt1();
			Real signedVolume = 0.0;
			Real momentX = 0.0;
			Real momentY = 0.0;
			Real momentZ = 0.0;
			for (const auto& triangle : triangles) {
				const Real tetrahedronVolume = SignedTetrahedronVolume(reference, triangle);
				const Vec3Cart centroid = (Vec3Cart(reference, triangle.Pnt1()) + Vec3Cart(reference, triangle.Pnt2())
					+ Vec3Cart(reference, triangle.Pnt3())) / 4;
				signedVolume += tetrahedronVolume;
				momentX += tetrahedronVolume * centroid.X();
				momentY += tetrahedronVolume * centroid.Y();
				momentZ += tetrahedronVolume * centroid.Z();
			}

			if (std::abs(signedVolume) <= Precision::GeometryEpsilon)
				throw GeometryError("Cannot calculate the center of a zero-volume surface mesh");
			return Pnt3Cart(reference.X() + momentX / signedVolume, reference.Y() + momentY / signedVolume,
				reference.Z() + momentZ / signedVolume);
		}

		inline Box3D BoundingBox(const std::vector<Triangle3D>& triangles) {
			return Box3D::FromPoints(Vertices(triangles));
		}

		inline BoundingSphere3D BoundingSphere(const std::vector<Triangle3D>& triangles) {
			const Pnt3Cart center = Center(triangles);
			Real radius = 0.0;
			for (const auto& vertex : Vertices(triangles))
				radius = std::max(radius, center.Dist(vertex));
			return BoundingSphere3D(center, radius);
		}

		inline bool IsPointOnTriangle(const Pnt3Cart& point, const Triangle3D& triangle) {
			const Vec3Cart edge1(triangle.Pnt1(), triangle.Pnt2());
			const Vec3Cart edge2(triangle.Pnt1(), triangle.Pnt3());
			const Vec3Cart normal = VectorProduct(edge1, edge2);
			const Real normalLength = normal.NormL2();
			if (normalLength <= Precision::GeometryEpsilon)
				return false;
			const Real distance = std::abs(ScalarProduct(Vec3Cart(triangle.Pnt1(), point), normal)) / normalLength;
			return distance <= Precision::GeometryEpsilon && triangle.IsPointInside(point);
		}

		inline bool IsInside(const Pnt3Cart& point, const std::vector<Triangle3D>& triangles) {
			Real solidAngle = 0.0;
			for (const auto& triangle : triangles) {
				if (IsPointOnTriangle(point, triangle))
					return true;

				const Vec3Cart a(point, triangle.Pnt1());
				const Vec3Cart b(point, triangle.Pnt2());
				const Vec3Cart c(point, triangle.Pnt3());
				const Real denominator = a.NormL2() * b.NormL2() * c.NormL2()
					+ ScalarProduct(a, b) * c.NormL2()
					+ ScalarProduct(b, c) * a.NormL2()
					+ ScalarProduct(c, a) * b.NormL2();
				const Real determinant = ScalarProduct(a, VectorProduct(b, c));
				solidAngle += MML::Real(2.0) * std::atan2(determinant, denominator);
			}
			return std::abs(solidAngle) > Constants::PI;
		}
	}

	/// @brief Abstract interface for 3D solid bodies.
	/// Defines the contract for all solid body representations in MML.
	/// Implementations provide geometric properties (volume, surface area),
	/// spatial queries (containment), and bounding volumes.

	class IBody {
	public:
		/// @name Geometric Properties
		/// @{
		virtual Real Volume() const = 0;		///< Total enclosed volume
		virtual Real SurfaceArea() const = 0;	///< Total surface area
		virtual Pnt3Cart getCenter() const = 0; ///< Geometric centroid
		/// @}

		/// @name Bounding Volumes
		/// @{
		virtual Box3D getBoundingBox() const = 0;				///< Axis-aligned bounding box
		virtual BoundingSphere3D getBoundingSphere() const = 0; ///< Bounding sphere
		/// @}

		/// @name Spatial Queries
		/// @{

		/// @brief Test if a point lies inside the body.
		/// @param pnt Point to test
		/// @return true if point is strictly inside

		virtual bool IsInside(const Pnt3Cart& pnt) const = 0;
		/// @}

		virtual std::string to_string() const = 0;
		/// @deprecated Use to_string() (CONTRACT_POLICY section 11: STL spelling). Non-virtual forwarder.
		[[deprecated("use to_string()")]]
		std::string ToString() const { return to_string(); }

		/// @deprecated Use getCenter() (CONTRACT_POLICY §11: accessors are getX()). Non-virtual forwarder;
		/// implementers override the lowerCamel name.
		[[deprecated("use getCenter()")]]
		Pnt3Cart GetCenter() const { return getCenter(); }
		/// @deprecated Use getBoundingBox().
		[[deprecated("use getBoundingBox()")]]
		Box3D GetBoundingBox() const { return getBoundingBox(); }
		/// @deprecated Use getBoundingSphere().
		[[deprecated("use getBoundingSphere()")]]
		BoundingSphere3D GetBoundingSphere() const { return getBoundingSphere(); }

		virtual ~IBody() = default;
	};

	/// @brief Solid body defined by boundary functions for numerical integration.
	/// This class represents a 3D solid defined by functional boundaries:
	/// - x ∈ [x1, x2]
	/// - y ∈ [y1(x), y2(x)]
	/// - z ∈ [z1(x,y), z2(x,y)]
	/// This representation is ideal for:
	/// - Volume integration (mass, center of mass)
	/// - Moment of inertia calculations
	/// - Bodies with variable density
	/// @note Geometric measurements are numerical estimates. Guaranteed bounding
	/// volumes require a caller-supplied conservative box through SetBoundingBox.
	/// @see SolidBodyWithBoundary for variable density implementation
	/// @see SolidBodyWithBoundaryConstDensity for constant density

	class ISolidBodyWithBoundary : public IBody {
		using YBoundary = Real (*)(Real);
		using ZBoundary = Real (*)(Real, Real);
		struct GuaranteedBounds {
			Box3D box;
			Real x1, x2;
			YBoundary y1, y2;
			ZBoundary z1, z2;
		};
		std::optional<GuaranteedBounds> _guaranteedBounds;

		void ValidateDefinition() const {
			Geometry3DBodyDetail::ValidateBounds(_x1, _x2);
			if (!_y1 || !_y2 || !_z1 || !_z2)
				throw GeometryError("Function-bounded geometry requires nonnull boundary callbacks");
		}

		void ValidateEstimate(const BoundaryGeometrySettings& settings) const {
			ValidateDefinition();
			Geometry3DBodyDetail::ValidateSettings(settings);
		}

		std::pair<Real, Real> YBounds(Real x) const {
			const Real lower = _y1(x), upper = _y2(x);
			Geometry3DBodyDetail::ValidateBounds(lower, upper);
			return {lower, upper};
		}

		std::pair<Real, Real> ZBounds(Real x, Real y) const {
			const Real lower = _z1(x, y), upper = _z2(x, y);
			Geometry3DBodyDetail::ValidateBounds(lower, upper);
			return {lower, upper};
		}

		Real XAt(Real s) const { return std::lerp(_x1, _x2, s); }

		Real YAt(Real s, Real t) const {
			const auto bounds = YBounds(XAt(s));
			return std::lerp(bounds.first, bounds.second, t);
		}

		Real ZAt(bool upper, Real s, Real t) const {
			const Real x = XAt(s);
			const auto y = YBounds(x);
			const auto z = ZBounds(x, std::lerp(y.first, y.second, t));
			return upper ? z.second : z.first;
		}

		template<typename Function>
		Real IntegrateDomain(Function&& function, const BoundaryGeometrySettings& settings) const {
			const Real widthX = Geometry3DBodyDetail::IntervalWidth(_x1, _x2);
			const Real average = Geometry3DBodyDetail::IntegrateSimpson([&](Real s) {
				const Real x = XAt(s);
				const auto bounds = YBounds(x);
				const Real widthY = Geometry3DBodyDetail::IntervalWidth(bounds.first, bounds.second);
				return Geometry3DBodyDetail::IntegrateSimpson([&](Real t) {
					const Real y = std::lerp(bounds.first, bounds.second, t);
					const auto z = ZBounds(x, y);
					return function(x, y, z.first, z.second) * widthY;
				}, 0, 1, settings.integrationIntervals);
			}, 0, 1, settings.integrationIntervals);
			return Geometry3DBodyDetail::FiniteValue(widthX * average);
		}

		static BoundingSphere3D SphereFromBox(const Box3D& box) {
			const auto lower = box.Min(), upper = box.Max();
			Geometry3DBodyDetail::ValidateBounds(lower.X(), upper.X());
			Geometry3DBodyDetail::ValidateBounds(lower.Y(), upper.Y());
			Geometry3DBodyDetail::ValidateBounds(lower.Z(), upper.Z());
			const Pnt3Cart center(std::lerp(lower.X(), upper.X(), Real(1) / 2),
				std::lerp(lower.Y(), upper.Y(), Real(1) / 2), std::lerp(lower.Z(), upper.Z(), Real(1) / 2));
			const auto roundUp = [](Real value) {
				return Geometry3DBodyDetail::FiniteValue(std::nextafter(value, std::numeric_limits<Real>::infinity()));
			};
			const auto extent = [&](Real lo, Real mid, Real hi) {
				const Real difference = std::max(mid - lo, hi - mid);
				return difference == 0 ? Real(0) : roundUp(difference);
			};
			const Real x = extent(lower.X(), center.X(), upper.X());
			const Real y = extent(lower.Y(), center.Y(), upper.Y());
			const Real z = extent(lower.Z(), center.Z(), upper.Z());
			const Real scale = std::max({x, y, z});
			if (scale == 0) return BoundingSphere3D(center, 0);
			Real squaredNorm = 0;
			// Outward-round each scaled operation so the sphere encloses the supplied box.
			for (Real component : {x, y, z}) {
				if (component == 0) continue;
				const Real ratio = roundUp(component / scale);
				squaredNorm = roundUp(squaredNorm + roundUp(ratio * ratio));
			}
			return BoundingSphere3D(center, roundUp(scale * roundUp(std::sqrt(squaredNorm))));
		}

	public:
		Real _x1, _x2; ///< X-axis bounds

		Real (*_y1)(Real); ///< Lower Y bound as function of x
		Real (*_y2)(Real); ///< Upper Y bound as function of x

		Real (*_z1)(Real, Real); ///< Lower Z bound as function of (x, y)
		Real (*_z2)(Real, Real); ///< Upper Z bound as function of (x, y)

	public:
		/// @brief Construct body from boundary functions.
		/// @param x1,x2 X-axis bounds
		/// @param y1,y2 Y bounds as functions of x
		/// @param z1,z2 Z bounds as functions of (x, y)

		ISolidBodyWithBoundary(Real x1, Real x2, Real (*y1)(Real), Real (*y2)(Real), Real (*z1)(Real, Real), Real (*z2)(Real, Real))
			: _x1(x1)
			, _x2(x2)
			, _y1(y1)
			, _y2(y2)
			, _z1(z1)
			, _z2(z2) { ValidateDefinition(); }

		/// @brief Get density at a point inside the body.
		/// @param x Position vector (3D)
		/// @return Density at that position

		virtual Real getDensity(const VectorN<Real, 3>& x) const = 0;

		virtual bool IsInside(const Pnt3Cart& pnt) const override {
			ValidateDefinition();
			detail::ValidateGeometryCoordinates({pnt.X(), pnt.Y(), pnt.Z()}, "ISolidBodyWithBoundary::IsInside");
			const Real x = pnt.X();
			const Real y = pnt.Y();
			const Real z = pnt.Z();

			if (_x1 < x && x < _x2) {
				const auto boundsY = YBounds(x);
				if (boundsY.first < y && y < boundsY.second) {
					const auto boundsZ = ZBounds(x, y);
					if (boundsZ.first < z && z < boundsZ.second) {
						return true; // point is inside the solid body
					}
				}
			}
			return false;
		}

		virtual Real Volume() const override {
			return EstimateVolume();
		}

		/// Fixed-grid Simpson estimate; resolution is not a certified error bound.
		Real EstimateVolume(const BoundaryGeometrySettings& settings = {}) const {
			ValidateEstimate(settings);
			return IntegrateDomain([](Real, Real, Real lower, Real upper) { return upper - lower; }, settings);
		}

		virtual Real SurfaceArea() const override { return EstimateSurfaceArea(); }

		Real EstimateSurfaceArea(const BoundaryGeometrySettings& settings = {}) const {
			ValidateEstimate(settings);
			const Real widthX = Geometry3DBodyDetail::IntervalWidth(_x1, _x2);
			const auto graphArea = [&](bool upper) {
				return Geometry3DBodyDetail::IntegrateSimpson([&](Real s) {
					const auto bounds = YBounds(XAt(s));
					const Real widthY = Geometry3DBodyDetail::IntervalWidth(bounds.first, bounds.second);
					return Geometry3DBodyDetail::IntegrateSimpson([&](Real t) {
						const Real ys = Geometry3DBodyDetail::UnitDerivative([&](Real u) { return YAt(u, t); }, s, settings.derivativeStep);
						const Real zs = Geometry3DBodyDetail::UnitDerivative([&](Real u) { return ZAt(upper, u, t); }, s, settings.derivativeStep);
						const Real zt = Geometry3DBodyDetail::UnitDerivative([&](Real u) { return ZAt(upper, s, u); }, t, settings.derivativeStep);
						// Norm of r_s cross r_t for r(s,t)=(x(s),y(s,t),z(s,t)).
						return std::hypot(ys * zt - zs * widthY, widthX * zt, widthX * widthY);
					}, 0, 1, settings.integrationIntervals);
				}, 0, 1, settings.integrationIntervals);
			};

			const auto yWallArea = [&](YBoundary boundary) {
				return Geometry3DBodyDetail::IntegrateSimpson([&](Real s) {
					const Real x = XAt(s), y = Geometry3DBodyDetail::FiniteValue(boundary(x));
					const auto z = ZBounds(x, y);
					const Real derivative = Geometry3DBodyDetail::UnitDerivative(
						[&](Real u) { return boundary(XAt(u)); }, s, settings.derivativeStep);
					return (z.second - z.first) * std::hypot(widthX, derivative);
				}, 0, 1, settings.integrationIntervals);
			};

			const auto xWallArea = [&](Real x) {
				const auto bounds = YBounds(x);
				return Geometry3DBodyDetail::IntegrateSimpson([&](Real y) {
					const auto z = ZBounds(x, y);
					return z.second - z.first;
				}, bounds.first, bounds.second, settings.integrationIntervals);
			};

			return Geometry3DBodyDetail::FiniteValue(graphArea(false) + graphArea(true)
				+ yWallArea(_y1) + yWallArea(_y2) + xWallArea(_x1) + xWallArea(_x2));
		}

		virtual Pnt3Cart getCenter() const override { return EstimateCenter(); }

		Pnt3Cart EstimateCenter(const BoundaryGeometrySettings& settings = {}) const {
			ValidateEstimate(settings);
			const Real volume = EstimateVolume(settings);
			if (volume == 0)
				throw GeometryError("Cannot calculate the center of a zero-volume bounded solid");

			const Real momentX = IntegrateDomain([](Real x, Real, Real lower, Real upper) { return x * (upper - lower); }, settings);
			const Real momentY = IntegrateDomain([](Real, Real y, Real lower, Real upper) { return y * (upper - lower); }, settings);
			const Real momentZ = IntegrateDomain([](Real, Real, Real lower, Real upper) {
				return std::lerp(lower, upper, Real(1) / 2) * (upper - lower);
			}, settings);
			return Pnt3Cart(Geometry3DBodyDetail::FiniteValue(momentX / volume),
				Geometry3DBodyDetail::FiniteValue(momentY / volume), Geometry3DBodyDetail::FiniteValue(momentZ / volume));
		}

		/// Caller asserts this box encloses the complete current boundary functions.
		void SetBoundingBox(const Box3D& box) {
			ValidateDefinition();
			Geometry3DBodyDetail::ValidateBounds(box.Min().X(), box.Max().X());
			Geometry3DBodyDetail::ValidateBounds(box.Min().Y(), box.Max().Y());
			Geometry3DBodyDetail::ValidateBounds(box.Min().Z(), box.Max().Z());
			if (box.Min().X() > _x1 || box.Max().X() < _x2)
				throw GeometryError("Supplied box does not enclose the body's x interval");
			_guaranteedBounds = GuaranteedBounds{box, _x1, _x2, _y1, _y2, _z1, _z2};
		}

		void ClearBoundingBox() { _guaranteedBounds.reset(); }

		virtual Box3D getBoundingBox() const override {
			ValidateDefinition();
			if (!_guaranteedBounds)
				throw GeometryError("A caller-guaranteed box is required; use EstimateBoundingBox for sampled extrema");
			const auto& bound = *_guaranteedBounds;
			if (bound.x1 != _x1 || bound.x2 != _x2 || bound.y1 != _y1 || bound.y2 != _y2
				|| bound.z1 != _z1 || bound.z2 != _z2)
				throw GeometryError("Boundary definition changed; supply a new conservative bounding box");
			return bound.box;
		}

		/// Sampled extrema only; this box need not enclose unsampled portions of the body.
		Box3D EstimateBoundingBox(const BoundaryGeometrySettings& settings = {}) const {
			ValidateEstimate(settings);
			Real minY = std::numeric_limits<Real>::infinity(), maxY = -minY;
			Real minZ = minY, maxZ = maxY;
			for (int xIndex = 0; xIndex <= settings.integrationIntervals; ++xIndex) {
				const Real x = XAt(Real(xIndex) / settings.integrationIntervals);
				const auto bounds = YBounds(x);
				const Real lowerY = bounds.first, upperY = bounds.second;
				minY = std::min(minY, lowerY);
				maxY = std::max(maxY, upperY);
				for (int yIndex = 0; yIndex <= settings.integrationIntervals; ++yIndex) {
					const Real y = std::lerp(lowerY, upperY, Real(yIndex) / settings.integrationIntervals);
					const auto z = ZBounds(x, y);
					minZ = std::min(minZ, z.first);
					maxZ = std::max(maxZ, z.second);
				}
			}
			return Box3D(Pnt3Cart(_x1, minY, minZ), Pnt3Cart(_x2, maxY, maxZ));
		}

		virtual BoundingSphere3D getBoundingSphere() const override {
			return SphereFromBox(getBoundingBox());
		}

		BoundingSphere3D EstimateBoundingSphere(const BoundaryGeometrySettings& settings = {}) const {
			return SphereFromBox(EstimateBoundingBox(settings));
		}

		virtual std::string to_string() const override {
			std::ostringstream oss;
			oss << "ISolidBodyWithBoundary[x∈[" << _x1 << "," << _x2 << "]]";
			return oss.str();
		}
	};

	/// @brief Solid body with variable density distribution.
	/// Extension of ISolidBodyWithBoundary where density varies spatially
	/// according to a user-provided function.

	class SolidBodyWithBoundary : public ISolidBodyWithBoundary {
		Real (*_density)(const VectorN<Real, 3>& x); ///< Density function
	public:
		/// @brief Construct with boundary functions and density function.

		SolidBodyWithBoundary(Real x1, Real x2, Real (*y1)(Real), Real (*y2)(Real), Real (*z1)(Real, Real), Real (*z2)(Real, Real),
							  Real (*density)(const VectorN<Real, 3>& x))
			: ISolidBodyWithBoundary(x1, x2, y1, y2, z1, z2)
			, _density(density) {}

		virtual Real getDensity(const VectorN<Real, 3>& x) const { return _density(x); }
	};

	/// @brief Solid body with uniform constant density.
	/// Simplified ISolidBodyWithBoundary for homogeneous materials.

	class SolidBodyWithBoundaryConstDensity : public ISolidBodyWithBoundary {
		Real _density; ///< Constant density value
	public:
		SolidBodyWithBoundaryConstDensity(Real x1, Real x2, Real (*y1)(Real), Real (*y2)(Real), Real (*z1)(Real, Real),
										  Real (*z2)(Real, Real), Real density)
			: ISolidBodyWithBoundary(x1, x2, y1, y2, z1, z2)
			, _density(density) {}

		virtual Real getDensity(const VectorN<Real, 3>& x) const { return _density; }
	};

	/// @brief Solid body defined by triangular surface mesh.
	/// Base class for bodies represented as a collection of Triangle3D faces.
	/// Used for mesh-based geometric representations that support:
	/// - Surface rendering
	/// - Surface integration
	/// - Ray casting (with appropriate algorithms)
	/// @note Faces must form a closed, consistently oriented mesh for volume,
	/// centroid, and containment calculations.

	class BodyWithTriangleSurfaces : public IBody {
	protected:
		std::vector<Triangle3D> _surfaces; ///< Triangle faces
	public:
		/// /** @brief Get the number of triangular faces. */

		int getSurfaceCount() const {
			if (_surfaces.size() > static_cast<size_t>((std::numeric_limits<int>::max)()))
				throw std::length_error("BodyWithTriangleSurfaces: face count exceeds the int index range");
			return static_cast<int>(_surfaces.size());
		}
		/// @deprecated Use getSurfaceCount() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getSurfaceCount()")]]
		int GetSurfaceCount() const { return getSurfaceCount(); }

		/// @brief Access a triangle face by index.
		/// @param index Face index (0 to getSurfaceCount()-1)
		/// @throws IndexError if index invalid

		const Triangle3D& getSurface(int index) const {
			if (index < 0 || index >= getSurfaceCount())
				throw IndexError("BodyWithTriangleSurfaces::getSurface - index out of range");
			return _surfaces[index];
		}
		/// @deprecated Use getSurface() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getSurface()")]]
		const Triangle3D& GetSurface(int index) const { return getSurface(index); }

		Real Volume() const override { return Geometry3DBodyDetail::Volume(_surfaces); }
		Real SurfaceArea() const override { return Geometry3DBodyDetail::SurfaceArea(_surfaces); }
		Pnt3Cart getCenter() const override { return Geometry3DBodyDetail::Center(_surfaces); }
		Box3D getBoundingBox() const override { return Geometry3DBodyDetail::BoundingBox(_surfaces); }
		BoundingSphere3D getBoundingSphere() const override { return Geometry3DBodyDetail::BoundingSphere(_surfaces); }
		bool IsInside(const Pnt3Cart& pnt) const override { return Geometry3DBodyDetail::IsInside(pnt, _surfaces); }
		std::string to_string() const override {
			std::ostringstream oss;
			oss << "BodyWithTriangleSurfaces{Triangles=" << _surfaces.size() << ", Volume=" << Volume()
				<< ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}
	};

	/// @brief Solid body defined by rectangular (quad) surface mesh.
	/// Base class for bodies represented as RectSurface3D faces.
	/// Commonly used for:
	/// - Box-shaped objects
	/// - Parametric surfaces (torus, etc.)
	/// - CAD-style representations
	/// @note Faces must form a closed, consistently oriented mesh for volume,
	/// centroid, and containment calculations.

	class BodyWithRectSurfaces : public IBody {
	protected:
		std::vector<RectSurface3D> _surfaces; ///< Rectangular faces

	public:
		int getSurfaceCount() const { return static_cast<int>(_surfaces.size()); }
		/// @deprecated Use getSurfaceCount() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getSurfaceCount()")]]
		int GetSurfaceCount() const { return getSurfaceCount(); }
		const RectSurface3D& getSurface(int index) const {
			if (index < 0 || index >= getSurfaceCount())
				throw IndexError("BodyWithRectSurfaces::getSurface - index out of range");
			return _surfaces[index];
		}
		/// @deprecated Use getSurface() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getSurface()")]]
		const RectSurface3D& GetSurface(int index) const { return getSurface(index); }

		Real Volume() const override { return Geometry3DBodyDetail::Volume(Geometry3DBodyDetail::Triangulate(_surfaces)); }
		Real SurfaceArea() const override { return Geometry3DBodyDetail::SurfaceArea(Geometry3DBodyDetail::Triangulate(_surfaces)); }
		Pnt3Cart getCenter() const override { return Geometry3DBodyDetail::Center(Geometry3DBodyDetail::Triangulate(_surfaces)); }
		Box3D getBoundingBox() const override { return Geometry3DBodyDetail::BoundingBox(Geometry3DBodyDetail::Triangulate(_surfaces)); }
		BoundingSphere3D getBoundingSphere() const override {
			return Geometry3DBodyDetail::BoundingSphere(Geometry3DBodyDetail::Triangulate(_surfaces));
		}
		bool IsInside(const Pnt3Cart& pnt) const override {
			return Geometry3DBodyDetail::IsInside(pnt, Geometry3DBodyDetail::Triangulate(_surfaces));
		}
		std::string to_string() const override {
			std::ostringstream oss;
			oss << "BodyWithRectSurfaces{Quads=" << _surfaces.size() << ", Volume=" << Volume()
				<< ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}
	};

	/// @brief Composite solid made of multiple IBody components.
	/// Represents a solid that is the union of other solids.
	/// Point containment returns true if point is inside ANY component.
	/// @note For tight/watertight validation, compute flux through all surfaces.

	class ComposedSolidSurfaces3D {
	private:
		std::vector<IBody*> _solids; ///< Component bodies (not owned)

	public:
		ComposedSolidSurfaces3D() = default;
		ComposedSolidSurfaces3D(const std::vector<IBody*>& solids)
			: _solids(solids) {}

		void AddSolid(IBody* solid) { _solids.push_back(solid); }
		size_t getSolidCount() const { return _solids.size(); }
		/// @deprecated Use getSolidCount() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getSolidCount()")]]
		size_t GetSolidCount() const { return getSolidCount(); }
		IBody* getSolid(size_t index) { return _solids[index]; }
		/// @deprecated Use getSolid() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getSolid()")]]
		IBody* GetSolid(size_t index)  { return getSolid(index); }
		const IBody* getSolid(size_t index) const { return _solids[index]; }
		/// @deprecated Use getSolid() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getSolid()")]]
		const IBody* GetSolid(size_t index) const { return getSolid(index); }

		/// /** @brief Test if point is inside any component solid. */

		bool IsInside(const Pnt3Cart& pnt) const {
			for (const auto* solid : _solids) {
				if (solid && solid->IsInside(pnt))
					return true;
			}
			return false;
		}
	};

} // namespace MML

#endif
