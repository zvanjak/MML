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

#include <sstream>
#include <string>
#include <vector>

#include "mml/MMLBase.h"
#include "mml/MMLExceptions.h"
#include "mml/base/Vector/VectorN.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry3D.h"
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBounding.h"

namespace MML {

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

		// IBody interface - provide default implementations
		// These are primarily used for physics calculations, not geometric queries
		virtual Real Volume() const override {
			// Volume would require integration - not implemented for general case
			throw NotImplementedError("Volume() requires integration - use ContinuousMassMomentOfInertiaTensorCalculator for mass properties");
		}

		virtual Real SurfaceArea() const override {
			// Surface area would require integration - not implemented for general case
			throw NotImplementedError("SurfaceArea() requires integration - not implemented for general ISolidBodyWithBoundary");
		}

		virtual Pnt3Cart GetCenter() const override {
			// Return geometric center of the bounding box as approximation
			Real cx = (_x1 + _x2) / 2.0;
			Real cy = (_y1(cx) + _y2(cx)) / 2.0;
			Real cz = (_z1(cx, cy) + _z2(cx, cy)) / 2.0;
			return Pnt3Cart(cx, cy, cz);
		}

		virtual Box3D GetBoundingBox() const override {
			// For general case, we need to find min/max across the entire domain
			// This is an approximation - assumes boundaries are monotonic
			Real minY = std::min(_y1(_x1), _y1(_x2));
			Real maxY = std::max(_y2(_x1), _y2(_x2));

			Real midX = (_x1 + _x2) / 2.0;
			Real midY = (_y1(midX) + _y2(midX)) / 2.0;
			Real minZ = std::min(_z1(_x1, midY), std::min(_z1(_x2, midY), _z1(midX, minY)));
			Real maxZ = std::max(_z2(_x1, midY), std::max(_z2(_x2, midY), _z2(midX, maxY)));

			return Box3D(Pnt3Cart(_x1, minY, minZ), Pnt3Cart(_x2, maxY, maxZ));
		}

		virtual BoundingSphere3D GetBoundingSphere() const override {
			Pnt3Cart center = GetCenter();
			Box3D box = GetBoundingBox();

			// Use half the diagonal of bounding box as radius
			Real dx = box.Max().X() - box.Min().X();
			Real dy = box.Max().Y() - box.Min().Y();
			Real dz = box.Max().Z() - box.Min().Z();
			Real radius = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);

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
	/// @note IBody methods throw by default; override in derived classes.

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

		// IBody interface - to be implemented by derived classes
		Real Volume() const override { throw NotImplementedError("Volume() not implemented"); }
		Real SurfaceArea() const override { throw NotImplementedError("SurfaceArea() not implemented"); }
		Pnt3Cart GetCenter() const override { throw NotImplementedError("GetCenter() not implemented"); }
		Box3D GetBoundingBox() const override { throw NotImplementedError("GetBoundingBox() not implemented"); }
		BoundingSphere3D GetBoundingSphere() const override { throw NotImplementedError("GetBoundingSphere() not implemented"); }
		bool IsInside(const Pnt3Cart& pnt) const override { throw NotImplementedError("IsInside() not implemented"); }
		std::string ToString() const override { throw NotImplementedError("ToString() not implemented"); }
	};

	/// @brief Solid body defined by rectangular (quad) surface mesh.
	/// Base class for bodies represented as RectSurface3D faces.
	/// Commonly used for:
	/// - Box-shaped objects
	/// - Parametric surfaces (torus, etc.)
	/// - CAD-style representations
	/// @note IBody methods throw by default; override in derived classes.

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

		// IBody interface - to be implemented by derived classes
		Real Volume() const override { throw NotImplementedError("Volume() not implemented"); }
		Real SurfaceArea() const override { throw NotImplementedError("SurfaceArea() not implemented"); }
		Pnt3Cart GetCenter() const override { throw NotImplementedError("GetCenter() not implemented"); }
		Box3D GetBoundingBox() const override { throw NotImplementedError("GetBoundingBox() not implemented"); }
		BoundingSphere3D GetBoundingSphere() const override { throw NotImplementedError("GetBoundingSphere() not implemented"); }
		bool IsInside(const Pnt3Cart& pnt) const override { throw NotImplementedError("IsInside() not implemented"); }
		std::string ToString() const override { throw NotImplementedError("ToString() not implemented"); }

		// TODO: isClosed() - cast rays from center of mass in all directions
		// and check if they hit at least one surface
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
