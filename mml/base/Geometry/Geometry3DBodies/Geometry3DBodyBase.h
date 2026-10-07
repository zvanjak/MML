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
#include <sstream>
#include <string>
#include <vector>

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DBounding.h>

namespace MML {
	namespace Geometry3DBodyDetail {
		constexpr int BoundaryIntegrationIntervals = 64;

		template<typename Function>
		Real IntegrateSimpson(Function&& function, Real lower, Real upper) {
			if (lower == upper)
				return 0.0;
			const Real step = (upper - lower) / BoundaryIntegrationIntervals;
			Real sum = function(lower) + function(upper);
			for (int index = 1; index < BoundaryIntegrationIntervals; ++index)
				sum += (index % 2 == 0 ? 2.0 : 4.0) * function(lower + index * step);
			return sum * step / 3.0;
		}

		template<typename Function, typename LowerBound, typename UpperBound>
		Real IntegrateBoundaryDomain(Function&& function, Real x1, Real x2, LowerBound&& y1, UpperBound&& y2) {
			return IntegrateSimpson(
				[&](Real x) {
					return IntegrateSimpson([&](Real y) { return function(x, y); }, y1(x), y2(x));
				},
				x1, x2);
		}

		inline Real DerivativeX(Real (*function)(Real, Real), Real x, Real y) {
			const Real step = 1.0e-5 * (1.0 + std::abs(x));
			return (function(x + step, y) - function(x - step, y)) / (2.0 * step);
		}

		inline Real DerivativeY(Real (*function)(Real, Real), Real x, Real y) {
			const Real step = 1.0e-5 * (1.0 + std::abs(y));
			return (function(x, y + step) - function(x, y - step)) / (2.0 * step);
		}

		inline Real Derivative(Real (*function)(Real), Real x) {
			const Real step = 1.0e-5 * (1.0 + std::abs(x));
			return (function(x + step) - function(x - step)) / (2.0 * step);
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
			return ScalarProduct(a, VectorProduct(b, c)) / 6.0;
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
				const Pnt3Cart centroid = (reference + triangle.Pnt1() + triangle.Pnt2() + triangle.Pnt3()) / 4.0;
				signedVolume += tetrahedronVolume;
				momentX += tetrahedronVolume * centroid.X();
				momentY += tetrahedronVolume * centroid.Y();
				momentZ += tetrahedronVolume * centroid.Z();
			}

			if (std::abs(signedVolume) <= Constants::GEOMETRY_EPSILON)
				throw GeometryError("Cannot calculate the center of a zero-volume surface mesh");
			return Pnt3Cart(momentX / signedVolume, momentY / signedVolume, momentZ / signedVolume);
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
			if (normalLength <= Constants::GEOMETRY_EPSILON)
				return false;
			const Real distance = std::abs(ScalarProduct(Vec3Cart(triangle.Pnt1(), point), normal)) / normalLength;
			return distance <= Constants::GEOMETRY_EPSILON && triangle.IsPointInside(point);
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
				solidAngle += 2.0 * std::atan2(determinant, denominator);
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
		virtual Pnt3Cart GetCenter() const = 0; ///< Geometric centroid
		/// @}

		/// @name Bounding Volumes
		/// @{
		virtual Box3D GetBoundingBox() const = 0;				///< Axis-aligned bounding box
		virtual BoundingSphere3D GetBoundingSphere() const = 0; ///< Bounding sphere
		/// @}

		/// @name Spatial Queries
		/// @{

		/// @brief Test if a point lies inside the body.
		/// @param pnt Point to test
		/// @return true if point is strictly inside

		virtual bool IsInside(const Pnt3Cart& pnt) const = 0;
		/// @}

		virtual std::string ToString() const = 0;

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
	/// @note Volume/SurfaceArea require numerical integration via
	/// ContinuousMassMomentOfInertiaTensorCalculator.
	/// @see SolidBodyWithBoundary for variable density implementation
	/// @see SolidBodyWithBoundaryConstDensity for constant density

	class ISolidBodyWithBoundary : public IBody {
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
			, _z2(z2) {}

		/// @brief Get density at a point inside the body.
		/// @param x Position vector (3D)
		/// @return Density at that position

		virtual Real getDensity(const VectorN<Real, 3>& x) const = 0;

		virtual bool IsInside(const Pnt3Cart& pnt) const override {
			const Real x = pnt.X();
			const Real y = pnt.Y();
			const Real z = pnt.Z();

			if (_x1 < x && x < _x2) {
				// check y bounds
				if (_y1(x) < y && y < _y2(x)) {
					// check z bounds
					if (_z1(x, y) < z && z < _z2(x, y)) {
						return true; // point is inside the solid body
					}
				}
			}
			return false;
		}

		virtual Real Volume() const override {
			return Geometry3DBodyDetail::IntegrateBoundaryDomain(
				[&](Real x, Real y) { return _z2(x, y) - _z1(x, y); }, _x1, _x2, _y1, _y2);
		}

		virtual Real SurfaceArea() const override {
			const auto graphArea = [&](Real (*surface)(Real, Real)) {
				return Geometry3DBodyDetail::IntegrateBoundaryDomain(
					[&](Real x, Real y) {
						const Real dx = Geometry3DBodyDetail::DerivativeX(surface, x, y);
						const Real dy = Geometry3DBodyDetail::DerivativeY(surface, x, y);
						return std::sqrt(1.0 + dx * dx + dy * dy);
					},
					_x1, _x2, _y1, _y2);
			};

			const auto yWallArea = [&](Real (*boundary)(Real)) {
				return Geometry3DBodyDetail::IntegrateSimpson(
					[&](Real x) {
						const Real y = boundary(x);
						const Real derivative = Geometry3DBodyDetail::Derivative(boundary, x);
						return (_z2(x, y) - _z1(x, y)) * std::sqrt(1.0 + derivative * derivative);
					},
					_x1, _x2);
			};

			const auto xWallArea = [&](Real x) {
				return Geometry3DBodyDetail::IntegrateSimpson(
					[&](Real y) { return _z2(x, y) - _z1(x, y); }, _y1(x), _y2(x));
			};

			return graphArea(_z1) + graphArea(_z2) + yWallArea(_y1) + yWallArea(_y2)
				+ xWallArea(_x1) + xWallArea(_x2);
		}

		virtual Pnt3Cart GetCenter() const override {
			const Real volume = Volume();
			if (std::abs(volume) <= Constants::GEOMETRY_EPSILON)
				throw GeometryError("Cannot calculate the center of a zero-volume bounded solid");

			const auto integrate = [&](const auto& integrand) {
				return Geometry3DBodyDetail::IntegrateBoundaryDomain(integrand, _x1, _x2, _y1, _y2);
			};
			const Real momentX = integrate([&](Real x, Real y) { return x * (_z2(x, y) - _z1(x, y)); });
			const Real momentY = integrate([&](Real x, Real y) { return y * (_z2(x, y) - _z1(x, y)); });
			const Real momentZ = integrate([&](Real x, Real y) {
				const Real lower = _z1(x, y);
				const Real upper = _z2(x, y);
				return 0.5 * (upper * upper - lower * lower);
			});
			return Pnt3Cart(momentX / volume, momentY / volume, momentZ / volume);
		}

		virtual Box3D GetBoundingBox() const override {
			Real minY = _y1(_x1);
			Real maxY = _y2(_x1);
			Real minZ = _z1(_x1, minY);
			Real maxZ = _z2(_x1, maxY);
			for (int xIndex = 0; xIndex <= Geometry3DBodyDetail::BoundaryIntegrationIntervals; ++xIndex) {
				const Real x = _x1 + (_x2 - _x1) * xIndex / Geometry3DBodyDetail::BoundaryIntegrationIntervals;
				const Real lowerY = _y1(x);
				const Real upperY = _y2(x);
				minY = std::min(minY, lowerY);
				maxY = std::max(maxY, upperY);
				for (int yIndex = 0; yIndex <= Geometry3DBodyDetail::BoundaryIntegrationIntervals; ++yIndex) {
					const Real y = lowerY + (upperY - lowerY) * yIndex / Geometry3DBodyDetail::BoundaryIntegrationIntervals;
					minZ = std::min(minZ, _z1(x, y));
					maxZ = std::max(maxZ, _z2(x, y));
				}
			}
			return Box3D(Pnt3Cart(_x1, minY, minZ), Pnt3Cart(_x2, maxY, maxZ));
		}

		virtual BoundingSphere3D GetBoundingSphere() const override {
			Pnt3Cart center = GetCenter();
			Box3D box = GetBoundingBox();

			Real radius = 0.0;
			for (Real x : {box.Min().X(), box.Max().X()})
				for (Real y : {box.Min().Y(), box.Max().Y()})
					for (Real z : {box.Min().Z(), box.Max().Z()})
						radius = std::max(radius, center.Dist(Pnt3Cart(x, y, z)));
			return BoundingSphere3D(center, radius);
		}

		virtual std::string ToString() const override {
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

		int GetSurfaceCount() const { return _surfaces.size(); }

		/// @brief Access a triangle face by index.
		/// @param index Face index (0 to GetSurfaceCount()-1)
		/// @throws IndexError if index invalid

		const Triangle3D& GetSurface(int index) const {
			if (index < 0 || index >= GetSurfaceCount())
				throw IndexError("BodyWithTriangleSurfaces::GetSurface - index out of range");
			return _surfaces[index];
		}

		Real Volume() const override { return Geometry3DBodyDetail::Volume(_surfaces); }
		Real SurfaceArea() const override { return Geometry3DBodyDetail::SurfaceArea(_surfaces); }
		Pnt3Cart GetCenter() const override { return Geometry3DBodyDetail::Center(_surfaces); }
		Box3D GetBoundingBox() const override { return Geometry3DBodyDetail::BoundingBox(_surfaces); }
		BoundingSphere3D GetBoundingSphere() const override { return Geometry3DBodyDetail::BoundingSphere(_surfaces); }
		bool IsInside(const Pnt3Cart& pnt) const override { return Geometry3DBodyDetail::IsInside(pnt, _surfaces); }
		std::string ToString() const override {
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
		int GetSurfaceCount() const { return static_cast<int>(_surfaces.size()); }
		const RectSurface3D& GetSurface(int index) const {
			if (index < 0 || index >= GetSurfaceCount())
				throw IndexError("BodyWithRectSurfaces::GetSurface - index out of range");
			return _surfaces[index];
		}

		Real Volume() const override { return Geometry3DBodyDetail::Volume(Geometry3DBodyDetail::Triangulate(_surfaces)); }
		Real SurfaceArea() const override { return Geometry3DBodyDetail::SurfaceArea(Geometry3DBodyDetail::Triangulate(_surfaces)); }
		Pnt3Cart GetCenter() const override { return Geometry3DBodyDetail::Center(Geometry3DBodyDetail::Triangulate(_surfaces)); }
		Box3D GetBoundingBox() const override { return Geometry3DBodyDetail::BoundingBox(Geometry3DBodyDetail::Triangulate(_surfaces)); }
		BoundingSphere3D GetBoundingSphere() const override {
			return Geometry3DBodyDetail::BoundingSphere(Geometry3DBodyDetail::Triangulate(_surfaces));
		}
		bool IsInside(const Pnt3Cart& pnt) const override {
			return Geometry3DBodyDetail::IsInside(pnt, Geometry3DBodyDetail::Triangulate(_surfaces));
		}
		std::string ToString() const override {
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
		size_t GetSolidCount() const { return _solids.size(); }
		IBody* GetSolid(size_t index) { return _solids[index]; }
		const IBody* GetSolid(size_t index) const { return _solids[index]; }

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
