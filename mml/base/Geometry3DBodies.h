///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DBodies.h                                                  ///
///  Description: 3D solid bodies (Sphere, Box, Cylinder, Cone, Torus)                ///
///               Volume, surface area, intersection calculations                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_3D_BODIES_H
#define MML_GEOMETRY_3D_BODIES_H

#include <sstream>

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/BaseUtils.h"

#include "base/VectorN.h"
#include "base/VectorTypes.h"
#include "base/Geometry3D.h"

using namespace MML::Utils;

namespace MML
{
	// ============================================================================
	// Bounding Volumes for Spatial Queries and Collision Detection
	// ============================================================================

	/// <summary>
	/// Axis-Aligned Bounding Box (AABB) for fast spatial queries and collision broad-phase.
	/// Defined by minimum and maximum corner points.
	/// </summary>
	class BoundingBox3D
	{
	private:
		Pnt3Cart _min;
		Pnt3Cart _max;

	public:
		// Constructors
		BoundingBox3D() : _min(0, 0, 0), _max(0, 0, 0) {}
		
		BoundingBox3D(const Pnt3Cart& min, const Pnt3Cart& max) 
			: _min(min), _max(max)
		{
			// Validate that min is actually less than max in all dimensions
			if (min.X() > max.X() || min.Y() > max.Y() || min.Z() > max.Z())
				throw std::invalid_argument("BoundingBox3D: min point must be <= max point in all dimensions");
		}

		// Accessors
		Pnt3Cart Min() const { return _min; }
		Pnt3Cart Max() const { return _max; }

		// Dimensions
		Real Width() const { return _max.X() - _min.X(); }
		Real Height() const { return _max.Y() - _min.Y(); }
		Real Depth() const { return _max.Z() - _min.Z(); }
		Real Volume() const { return Width() * Height() * Depth(); }

		// Center point
		Pnt3Cart Center() const 
		{ 
			return Pnt3Cart(
				(_min.X() + _max.X()) / 2.0,
				(_min.Y() + _max.Y()) / 2.0,
				(_min.Z() + _max.Z()) / 2.0
			);
		}

		// Point containment test
		bool Contains(const Pnt3Cart& pnt) const 
		{
			return pnt.X() >= _min.X() && pnt.X() <= _max.X() &&
			       pnt.Y() >= _min.Y() && pnt.Y() <= _max.Y() &&
			       pnt.Z() >= _min.Z() && pnt.Z() <= _max.Z();
		}

		// Box-box intersection test
		bool Intersects(const BoundingBox3D& other) const 
		{
			// Boxes don't intersect if separated along any axis
			return !(_max.X() < other._min.X() || _min.X() > other._max.X() ||
			         _max.Y() < other._min.Y() || _min.Y() > other._max.Y() ||
			         _max.Z() < other._min.Z() || _min.Z() > other._max.Z());
		}

		// Expand to include a point
		void ExpandToInclude(const Pnt3Cart& pnt)
		{
			if (pnt.X() < _min.X()) _min.X() = pnt.X();
			if (pnt.Y() < _min.Y()) _min.Y() = pnt.Y();
			if (pnt.Z() < _min.Z()) _min.Z() = pnt.Z();
			if (pnt.X() > _max.X()) _max.X() = pnt.X();
			if (pnt.Y() > _max.Y()) _max.Y() = pnt.Y();
			if (pnt.Z() > _max.Z()) _max.Z() = pnt.Z();
		}

		// Expand to include another box
		void ExpandToInclude(const BoundingBox3D& box)
		{
			ExpandToInclude(box._min);
			ExpandToInclude(box._max);
		}

		// String representation
		std::string ToString() const
		{
			std::ostringstream oss;
			oss << "BoundingBox3D[Min=(" << _min.X() << "," << _min.Y() << "," << _min.Z() << ")"
			    << ", Max=(" << _max.X() << "," << _max.Y() << "," << _max.Z() << ")"
			    << ", Size=(" << Width() << "×" << Height() << "×" << Depth() << ")]";
			return oss.str();
		}
	};

	/// <summary>
	/// Bounding Sphere for fast spatial queries and collision detection.
	/// Often faster than AABB for intersection tests but less tight fitting.
	/// </summary>
	class BoundingSphere3D
	{
	private:
		Pnt3Cart _center;
		Real _radius;

	public:
		// Constructors
		BoundingSphere3D() : _center(0, 0, 0), _radius(0) {}
		
		BoundingSphere3D(const Pnt3Cart& center, Real radius) 
			: _center(center), _radius(radius)
		{
			if (radius < 0.0)
				throw std::invalid_argument("BoundingSphere3D: radius must be non-negative");
		}

		// Accessors
		Pnt3Cart Center() const { return _center; }
		Real Radius() const { return _radius; }

		// Volume (4/3 * π * r³)
		Real Volume() const 
		{ 
			return (4.0 / 3.0) * Constants::PI * _radius * _radius * _radius; 
		}

		// Point containment test
		bool Contains(const Pnt3Cart& pnt) const 
		{
			return _center.Dist(pnt) <= _radius;
		}

		// Sphere-sphere intersection test
		bool Intersects(const BoundingSphere3D& other) const 
		{
			Real centerDist = _center.Dist(other._center);
			return centerDist <= (_radius + other._radius);
		}

		// Expand to include a point
		void ExpandToInclude(const Pnt3Cart& pnt)
		{
			Real dist = _center.Dist(pnt);
			if (dist > _radius)
				_radius = dist;
		}

		// String representation
		std::string ToString() const
		{
			std::ostringstream oss;
			oss << "BoundingSphere3D[Center=(" << _center.X() << "," << _center.Y() << "," << _center.Z() << ")"
			    << ", Radius=" << _radius << "]";
			return oss.str();
		}
	};

	class IBody
	{
	public:
		// Essential geometric properties
		virtual Real Volume() const = 0;
		virtual Real SurfaceArea() const = 0;
		virtual Pnt3Cart GetCenter() const = 0;
		
		// Bounding volumes
		virtual BoundingBox3D GetBoundingBox() const = 0;
		virtual BoundingSphere3D GetBoundingSphere() const = 0;
		
		// Spatial queries
		virtual bool IsInside(const Pnt3Cart& pnt) const = 0;
		
		// String representation
		virtual std::string ToString() const = 0;
		
		// Virtual destructor for proper cleanup
		virtual ~IBody() = default;
	};

	class ISolidBodyWithBoundary : public IBody
	{
	public:
		Real _x1, _x2;

		Real(*_y1)(Real);
		Real(*_y2)(Real);

		Real(*_z1)(Real, Real);
		Real(*_z2)(Real, Real);

	public:
		ISolidBodyWithBoundary(Real x1, Real x2,
			Real(*y1)(Real), Real(*y2)(Real),
			Real(*z1)(Real, Real), Real(*z2)(Real, Real))
			: _x1(x1), _x2(x2), _y1(y1), _y2(y2), _z1(z1), _z2(z2)
		{	}

		virtual Real getDensity(const VectorN<Real, 3>& x) const = 0;

		virtual bool IsInside(const Pnt3Cart& pnt) const override
		{
			const Real x = pnt.X();
			const Real y = pnt.Y();
			const Real z = pnt.Z();

			if (_x1 < x && x < _x2)
			{
				// check y bounds
				if (_y1(x) < y && y < _y2(x))
				{
					// check z bounds
					if (_z1(x, y) < z && z < _z2(x, y))
					{
						return true;	// point is inside the solid body
					}
				}
			}
			return false;
		}

		// IBody interface - provide default implementations
		// These are primarily used for physics calculations, not geometric queries
		virtual Real Volume() const override 
		{ 
			// Volume would require integration - not implemented for general case
			throw std::logic_error("Volume() requires integration - use ContinuousMassMomentOfInertiaTensorCalculator for mass properties"); 
		}
		
		virtual Real SurfaceArea() const override 
		{ 
			// Surface area would require integration - not implemented for general case
			throw std::logic_error("SurfaceArea() requires integration - not implemented for general ISolidBodyWithBoundary"); 
		}
		
		virtual Pnt3Cart GetCenter() const override 
		{ 
			// Return geometric center of the bounding box as approximation
			Real cx = (_x1 + _x2) / 2.0;
			Real cy = (_y1(cx) + _y2(cx)) / 2.0;
			Real cz = (_z1(cx, cy) + _z2(cx, cy)) / 2.0;
			return Pnt3Cart(cx, cy, cz);
		}
		
		virtual BoundingBox3D GetBoundingBox() const override 
		{ 
			// For general case, we need to find min/max across the entire domain
			// This is an approximation - assumes boundaries are monotonic
			Real minY = std::min(_y1(_x1), _y1(_x2));
			Real maxY = std::max(_y2(_x1), _y2(_x2));
			
			Real midX = (_x1 + _x2) / 2.0;
			Real midY = (_y1(midX) + _y2(midX)) / 2.0;
			Real minZ = std::min(_z1(_x1, midY), std::min(_z1(_x2, midY), _z1(midX, minY)));
			Real maxZ = std::max(_z2(_x1, midY), std::max(_z2(_x2, midY), _z2(midX, maxY)));
			
			return BoundingBox3D(Pnt3Cart(_x1, minY, minZ), Pnt3Cart(_x2, maxY, maxZ));
		}
		
		virtual BoundingSphere3D GetBoundingSphere() const override 
		{ 
			Pnt3Cart center = GetCenter();
			BoundingBox3D box = GetBoundingBox();
			
			// Use half the diagonal of bounding box as radius
			Real dx = box.Max().X() - box.Min().X();
			Real dy = box.Max().Y() - box.Min().Y();
			Real dz = box.Max().Z() - box.Min().Z();
			Real radius = 0.5 * std::sqrt(dx*dx + dy*dy + dz*dz);
			
			return BoundingSphere3D(center, radius);
		}
		
		virtual std::string ToString() const override 
		{ 
			std::ostringstream oss;
			oss << "ISolidBodyWithBoundary[x∈[" << _x1 << "," << _x2 << "]]";
			return oss.str();
		}

	};

	class SolidBodyWithBoundary : public ISolidBodyWithBoundary
	{
		Real(*_density)(const VectorN<Real, 3>& x);
	public:
		SolidBodyWithBoundary(Real x1, Real x2, Real(*y1)(Real), Real(*y2)(Real),
			Real(*z1)(Real, Real), Real(*z2)(Real, Real),
			Real(*density)(const VectorN<Real, 3>& x))
			: ISolidBodyWithBoundary(x1, x2, y1, y2, z1, z2), _density(density)
		{	}

		virtual Real getDensity(const VectorN<Real, 3>& x) const
		{
			return _density(x);
		}
	};

	class SolidBodyWithBoundaryConstDensity : public ISolidBodyWithBoundary
	{
		Real _density;
	public:
		SolidBodyWithBoundaryConstDensity(Real x1, Real x2, Real(*y1)(Real), Real(*y2)(Real),
			Real(*z1)(Real, Real), Real(*z2)(Real, Real),
			Real density)
			: ISolidBodyWithBoundary(x1, x2, y1, y2, z1, z2), _density(density)
		{	}

		virtual Real getDensity(const VectorN<Real, 3>& x) const
		{
			return _density;
		}
	};

	// solid body in 3D defined by triangular surfaces
	class BodyWithTriangleSurfaces : public IBody
	{
	protected:
		std::vector<Triangle3D> _surfaces;
	public:
		int GetSurfaceCount() const
		{
			return _surfaces.size();
		}
		const Triangle3D& GetSurface(int index) const
		{
			if (index < 0 || index >= GetSurfaceCount())
				throw std::out_of_range("BodyWithTriangleSurfaces::GetSurface - index out of range");
			return _surfaces[index];
		}
		
		// IBody interface - to be implemented by derived classes
		Real Volume() const override { throw std::logic_error("Volume() not implemented"); }
		Real SurfaceArea() const override { throw std::logic_error("SurfaceArea() not implemented"); }
		Pnt3Cart GetCenter() const override { throw std::logic_error("GetCenter() not implemented"); }
		BoundingBox3D GetBoundingBox() const override { throw std::logic_error("GetBoundingBox() not implemented"); }
		BoundingSphere3D GetBoundingSphere() const override { throw std::logic_error("GetBoundingSphere() not implemented"); }
		bool IsInside(const Pnt3Cart& pnt) const override { throw std::logic_error("IsInside() not implemented"); }
		std::string ToString() const override { throw std::logic_error("ToString() not implemented"); }
	};

	// solid body in 3D defined by rectangular surfaces
	class BodyWithRectSurfaces : public IBody
	{
	protected:
		std::vector<RectSurface3D> _surfaces;

	public:
		int GetSurfaceCount() const
		{
			return static_cast<int>(_surfaces.size());
		}
		const RectSurface3D& GetSurface(int index) const
		{
			if (index < 0 || index >= GetSurfaceCount())
				throw std::out_of_range("BodyWithRectSurfaces::GetSurface - index out of range");
			return _surfaces[index];
		}

		// IBody interface - to be implemented by derived classes
		Real Volume() const override { throw std::logic_error("Volume() not implemented"); }
		Real SurfaceArea() const override { throw std::logic_error("SurfaceArea() not implemented"); }
		Pnt3Cart GetCenter() const override { throw std::logic_error("GetCenter() not implemented"); }
		BoundingBox3D GetBoundingBox() const override { throw std::logic_error("GetBoundingBox() not implemented"); }
		BoundingSphere3D GetBoundingSphere() const override { throw std::logic_error("GetBoundingSphere() not implemented"); }
		bool IsInside(const Pnt3Cart& pnt) const override { throw std::logic_error("IsInside() not implemented"); }
		std::string ToString() const override { throw std::logic_error("ToString() not implemented"); }
		
		// isClosed() - iz centra mase (?) odasilje zrake u svim smje
		// rovima i gleda da li je pogodio iti jednu povrsinu
		// vraca listu povrsina, koje bi TREBALE omedjivati tijelo!

		// i po tome se moze obaviti surface integracija!!!
	};

	// represents solid that is composed of other solids
	// dobar nacin za provjeriti je li composed solid "tight" je izracunati fluks kroz sve povrsine
	class ComposedSolidSurfaces3D
	{
	private:
		std::vector<IBody*> _solids;

	public:
		ComposedSolidSurfaces3D() = default;
		ComposedSolidSurfaces3D(const std::vector<IBody*>& solids) : _solids(solids) {}

		void AddSolid(IBody* solid) { _solids.push_back(solid); }
		size_t GetSolidCount() const { return _solids.size(); }
		IBody* GetSolid(size_t index) { return _solids[index]; }
		const IBody* GetSolid(size_t index) const { return _solids[index]; }

		bool IsInside(const Pnt3Cart& pnt) const
		{
			// For a composed solid, point is inside if it's inside ANY of the component solids
			for (const auto* solid : _solids)
			{
				if (solid && solid->IsInside(pnt))
					return true;
			}
			return false;
		}
	};

	class Cube3D : public BodyWithRectSurfaces
	{
		// kocka
		Real _a;
		Pnt3Cart _center;
	public:
		Cube3D(Real a) : _a(a), _center(0, 0, 0)
		{
			Pnt3Cart pnt1( a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt2( a / 2,  a / 2, -a / 2);
			Pnt3Cart pnt3(-a / 2,  a / 2, -a / 2);
			Pnt3Cart pnt4(-a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt5( a / 2, -a / 2,  a / 2);
			Pnt3Cart pnt6( a / 2,  a / 2,  a / 2);
			Pnt3Cart pnt7(-a / 2,  a / 2,  a / 2);
			Pnt3Cart pnt8(-a / 2, -a / 2,  a / 2);

			// dodati svih 6 stranica u popis povrsina
			_surfaces.push_back(RectSurface3D(pnt1, pnt4, pnt3, pnt2));     // lower side in xy plane
			_surfaces.push_back(RectSurface3D(pnt5, pnt6, pnt7, pnt8));     // upper side in xy plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt2, pnt6, pnt5));     // front side in yz plane
			_surfaces.push_back(RectSurface3D(pnt4, pnt8, pnt7, pnt3));     // back side in yz plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt5, pnt8, pnt4));     // left side in xz plane
			_surfaces.push_back(RectSurface3D(pnt2, pnt3, pnt7, pnt6));     // right side in xz plane
		}
		Cube3D(Real a, const Pnt3Cart& center) : _a(a), _center(center)
		{
			// create 8 points for corners of the cube, centered at the given center point
			Pnt3Cart pnt1(center.X() + a / 2, center.Y() - a / 2, center.Z() - a / 2);
			Pnt3Cart pnt2(center.X() + a / 2, center.Y() + a / 2, center.Z() - a / 2);
			Pnt3Cart pnt3(center.X() - a / 2, center.Y() + a / 2, center.Z() - a / 2);
			Pnt3Cart pnt4(center.X() - a / 2, center.Y() - a / 2, center.Z() - a / 2);
			Pnt3Cart pnt5(center.X() + a / 2, center.Y() - a / 2, center.Z() + a / 2);
			Pnt3Cart pnt6(center.X() + a / 2, center.Y() + a / 2, center.Z() + a / 2);
			Pnt3Cart pnt7(center.X() - a / 2, center.Y() + a / 2, center.Z() + a / 2);
			Pnt3Cart pnt8(center.X() - a / 2, center.Y() - a / 2, center.Z() + a / 2);

			// add all 6 faces of the cube
			_surfaces.push_back(RectSurface3D(pnt1, pnt4, pnt3, pnt2));     // lower side in xy plane
			_surfaces.push_back(RectSurface3D(pnt5, pnt6, pnt7, pnt8));     // upper side in xy plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt2, pnt6, pnt5));     // front side in yz plane
			_surfaces.push_back(RectSurface3D(pnt4, pnt8, pnt7, pnt3));     // back side in yz plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt5, pnt8, pnt4));     // left side in xz plane
			_surfaces.push_back(RectSurface3D(pnt2, pnt3, pnt7, pnt6));     // right side in xz plane
		}

		bool IsInside(const Pnt3Cart& pnt) const override
		{
			// Check if point is inside the cube
			Real half_a = _a / 2.0;
			return (pnt.X() >= _center.X() - half_a && pnt.X() <= _center.X() + half_a &&
							pnt.Y() >= _center.Y() - half_a && pnt.Y() <= _center.Y() + half_a &&
							pnt.Z() >= _center.Z() - half_a && pnt.Z() <= _center.Z() + half_a);
		}

		// IBody interface implementation
		Real Volume() const override
		{
			// Volume = a³
			return _a * _a * _a;
		}

		Real SurfaceArea() const override
		{
			// Surface area = 6a²
			return 6.0 * _a * _a;
		}

		Pnt3Cart GetCenter() const override
		{
			return _center;
		}

		BoundingBox3D GetBoundingBox() const override
		{
			// AABB with min/max at center ± half-side in all dimensions
			Real half_a = _a / 2.0;
			Pnt3Cart min(_center.X() - half_a, _center.Y() - half_a, _center.Z() - half_a);
			Pnt3Cart max(_center.X() + half_a, _center.Y() + half_a, _center.Z() + half_a);
			return BoundingBox3D(min, max);
		}

		BoundingSphere3D GetBoundingSphere() const override
		{
			// Bounding sphere has radius = half the space diagonal = (a√3)/2
			Real radius = _a * std::sqrt(3.0) / 2.0;
			return BoundingSphere3D(_center, radius);
		}

		std::string ToString() const override
		{
			std::ostringstream oss;
			oss << "Cube3D{Center=(" << _center.X() << ", " << _center.Y() << ", " << _center.Z()
				<< "), Side=" << _a << ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}

		// Getters
		Real GetSide() const { return _a; }
	};

	class CubeWithTriangles3D : public BodyWithTriangleSurfaces
	{
		Real _a;
		Pnt3Cart _center;
	public:
		// Cube centered at origin
		CubeWithTriangles3D(Real a) : _a(a), _center(0, 0, 0)
		{
			Real h = a / 2.0;
			// 8 vertices of the cube
			Pnt3Cart p1(h, -h, -h);
			Pnt3Cart p2(h, h, -h);
			Pnt3Cart p3(-h, h, -h);
			Pnt3Cart p4(-h, -h, -h);
			Pnt3Cart p5(h, -h, h);
			Pnt3Cart p6(h, h, h);
			Pnt3Cart p7(-h, h, h);
			Pnt3Cart p8(-h, -h, h);

			// Each face: two triangles (CCW order for outward normals)
			// Bottom face (z = -h)
			_surfaces.emplace_back(p1, p4, p3);
			_surfaces.emplace_back(p1, p3, p2);

			// Top face (z = +h)
			_surfaces.emplace_back(p5, p6, p7);
			_surfaces.emplace_back(p5, p7, p8);

			// Front face (y = +h)
			_surfaces.emplace_back(p2, p3, p7);
			_surfaces.emplace_back(p2, p7, p6);

			// Back face (y = -h)
			_surfaces.emplace_back(p1, p5, p8);
			_surfaces.emplace_back(p1, p8, p4);

			// Left face (x = -h)
			_surfaces.emplace_back(p4, p8, p7);
			_surfaces.emplace_back(p4, p7, p3);

			// Right face (x = +h)
			_surfaces.emplace_back(p1, p2, p6);
			_surfaces.emplace_back(p1, p6, p5);
		}

		// Cube centered at arbitrary point
		CubeWithTriangles3D(Real a, const Pnt3Cart& center) : _a(a), _center(center)
		{
			Real h = a / 2.0;
			// 8 vertices of the cube, centered at _center
			Pnt3Cart p1(center.X() + h, center.Y() - h, center.Z() - h);
			Pnt3Cart p2(center.X() + h, center.Y() + h, center.Z() - h);
			Pnt3Cart p3(center.X() - h, center.Y() + h, center.Z() - h);
			Pnt3Cart p4(center.X() - h, center.Y() - h, center.Z() - h);
			Pnt3Cart p5(center.X() + h, center.Y() - h, center.Z() + h);
			Pnt3Cart p6(center.X() + h, center.Y() + h, center.Z() + h);
			Pnt3Cart p7(center.X() - h, center.Y() + h, center.Z() + h);
			Pnt3Cart p8(center.X() - h, center.Y() - h, center.Z() + h);

			// Each face: two triangles (CCW order for outward normals)
			// Bottom face (z = -h)
			_surfaces.emplace_back(p1, p4, p3);
			_surfaces.emplace_back(p1, p3, p2);

			// Top face (z = +h)
			_surfaces.emplace_back(p5, p6, p7);
			_surfaces.emplace_back(p5, p7, p8);

			// Front face (y = +h)
			_surfaces.emplace_back(p2, p3, p7);
			_surfaces.emplace_back(p2, p7, p6);

			// Back face (y = -h)
			_surfaces.emplace_back(p1, p5, p8);
			_surfaces.emplace_back(p1, p8, p4);

			// Left face (x = -h)
			_surfaces.emplace_back(p4, p8, p7);
			_surfaces.emplace_back(p4, p7, p3);

			// Right face (x = +h)
			_surfaces.emplace_back(p1, p2, p6);
			_surfaces.emplace_back(p1, p6, p5);
		}

		bool IsInside(const Pnt3Cart& pnt) const override
		{
			Real h = _a / 2.0;
			return (pnt.X() >= _center.X() - h && pnt.X() <= _center.X() + h &&
				pnt.Y() >= _center.Y() - h && pnt.Y() <= _center.Y() + h &&
				pnt.Z() >= _center.Z() - h && pnt.Z() <= _center.Z() + h);
		}

		// IBody interface implementation
		Real Volume() const override { return _a * _a * _a; }
		
		Real SurfaceArea() const override { return 6.0 * _a * _a; }
		
		Pnt3Cart GetCenter() const override { return _center; }
		
		BoundingBox3D GetBoundingBox() const override
		{
			Real h = _a / 2.0;
			Pnt3Cart min(_center.X() - h, _center.Y() - h, _center.Z() - h);
			Pnt3Cart max(_center.X() + h, _center.Y() + h, _center.Z() + h);
			return BoundingBox3D(min, max);
		}
		
		BoundingSphere3D GetBoundingSphere() const override
		{
			// Sphere centered at cube center with radius to corner
			// Distance to corner = (a/2)*sqrt(3)
			Real radius = (_a / 2.0) * std::sqrt(3.0);
			return BoundingSphere3D(_center, radius);
		}
		
		std::string ToString() const override
		{
			std::ostringstream oss;
			oss << "CubeWithTriangles3D: Center=(" << _center.X() << ", " << _center.Y() << ", " << _center.Z() << ")"
				<< ", Side=" << _a
				<< ", Volume=" << Volume()
				<< ", SurfaceArea=" << SurfaceArea()
				<< ", Triangles=" << _surfaces.size();
			return oss.str();
		}

		Real GetSide() const { return _a; }
	};

	class Torus3D : public BodyWithRectSurfaces
	{
		Real _R;  // Major radius (distance from center to tube center)
		Real _r;  // Minor radius (tube radius)
		Pnt3Cart _center;
		int _numU;  // Number of divisions around major circle
		int _numV;  // Number of divisions around tube
	public:
		// Torus centered at origin with specified major and minor radii
		Torus3D(Real R, Real r, int numU = 20, int numV = 12)
			: _R(R), _r(r), _center(0, 0, 0), _numU(numU), _numV(numV)
		{
			constructSurfaces();
		}

		// Torus centered at specified point
		Torus3D(Real R, Real r, const Pnt3Cart& center, int numU = 20, int numV = 12)
			: _R(R), _r(r), _center(center), _numU(numU), _numV(numV)
		{
			constructSurfaces();
		}

		void constructSurfaces()
		{
			_surfaces.clear();
			Real du = 2.0 * Constants::PI / _numU;
			Real dv = 2.0 * Constants::PI / _numV;

			// Generate rectangular patches for the torus surface
			for (int i = 0; i < _numU; ++i)
			{
				for (int j = 0; j < _numV; ++j)
				{
					Real u1 = i * du;
					Real u2 = (i + 1) * du;
					Real v1 = j * dv;
					Real v2 = (j + 1) * dv;

					// Calculate the four corners of this patch
					// Torus parametrization: x = (R + r*cos(v))*cos(u), y = (R + r*cos(v))*sin(u), z = r*sin(v)
					Pnt3Cart p1 = torusPoint(u1, v1);
					Pnt3Cart p2 = torusPoint(u2, v1);
					Pnt3Cart p3 = torusPoint(u2, v2);
					Pnt3Cart p4 = torusPoint(u1, v2);

					_surfaces.push_back(RectSurface3D(p1, p2, p3, p4));
				}
			}
		}

		Pnt3Cart torusPoint(Real u, Real v) const
		{
			Real x = (_R + _r * std::cos(v)) * std::cos(u);
			Real y = (_R + _r * std::cos(v)) * std::sin(u);
			Real z = _r * std::sin(v);
			return Pnt3Cart(_center.X() + x, _center.Y() + y, _center.Z() + z);
		}

		bool IsInside(const Pnt3Cart& pnt) const override
		{
			// Translate point to torus-centered coordinates
			Real x = pnt.X() - _center.X();
			Real y = pnt.Y() - _center.Y();
			Real z = pnt.Z() - _center.Z();

			// Distance from point to Z-axis
			Real d = std::sqrt(x * x + y * y);

			// Distance from point to the tube center circle
			Real dist = std::sqrt((d - _R) * (d - _R) + z * z);

			return dist <= _r;
		}

		// IBody interface implementation
		Real Volume() const override
		{
			// Volume = 2π²Rr²
			return 2.0 * Constants::PI * Constants::PI * _R * _r * _r;
		}

		Real SurfaceArea() const override
		{
			// Surface area = 4π²Rr
			return 4.0 * Constants::PI * Constants::PI * _R * _r;
		}

		Pnt3Cart GetCenter() const override
		{
			return _center;
		}

		BoundingBox3D GetBoundingBox() const override
		{
			// AABB with extent based on major and minor radii
			Real outerRadius = _R + _r;  // Outer radius of torus
			Pnt3Cart min(_center.X() - outerRadius, _center.Y() - outerRadius, _center.Z() - _r);
			Pnt3Cart max(_center.X() + outerRadius, _center.Y() + outerRadius, _center.Z() + _r);
			return BoundingBox3D(min, max);
		}

		BoundingSphere3D GetBoundingSphere() const override
		{
			// Bounding sphere with radius = R + r (outer radius)
			Real radius = _R + _r;
			return BoundingSphere3D(_center, radius);
		}

		std::string ToString() const override
		{
			std::ostringstream oss;
			oss << "Torus3D{Center=(" << _center.X() << ", " << _center.Y() << ", " << _center.Z()
				<< "), MajorRadius=" << _R << ", MinorRadius=" << _r 
				<< ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}

		// Getters
		Real GetMajorRadius() const { return _R; }
		Real GetMinorRadius() const { return _r; }
		int GetNumU() const { return _numU; }
		int GetNumV() const { return _numV; }
	};

	class Cylinder3D : public BodyWithTriangleSurfaces
	{
		Real _R;      // Radius
		Real _H;      // Height
		Pnt3Cart _center;
		int _numSegments;  // Number of segments around the circumference
	public:
		// Cylinder centered at origin with radius R and height H
		Cylinder3D(Real R, Real H, int numSegments = 20)
			: _R(R), _H(H), _center(0, 0, -H/2), _numSegments(numSegments)
		{
			constructSurfaces();
		}

		// Cylinder centered at specified point (center is at bottom)
		Cylinder3D(Real R, Real H, const Pnt3Cart& center, int numSegments = 20)
			: _R(R), _H(H), _center(center), _numSegments(numSegments)
		{
			constructSurfaces();
		}

		void constructSurfaces()
		{
			_surfaces.clear();
			Real angleStep = 2.0 * Constants::PI / _numSegments;

			// Bottom and top centers
			Pnt3Cart bottomCenter(_center.X(), _center.Y(), _center.Z());
			Pnt3Cart topCenter(_center.X(), _center.Y(), _center.Z() + _H);

			// Generate triangles for the cylinder
			for (int i = 0; i < _numSegments; ++i)
			{
				Real angle1 = i * angleStep;
				Real angle2 = (i + 1) * angleStep;

				// Points on bottom circle
				Pnt3Cart b1(_center.X() + _R * std::cos(angle1),
										_center.Y() + _R * std::sin(angle1),
										_center.Z());
				Pnt3Cart b2(_center.X() + _R * std::cos(angle2),
										_center.Y() + _R * std::sin(angle2),
										_center.Z());

				// Points on top circle
				Pnt3Cart t1(_center.X() + _R * std::cos(angle1),
										_center.Y() + _R * std::sin(angle1),
										_center.Z() + _H);
				Pnt3Cart t2(_center.X() + _R * std::cos(angle2),
										_center.Y() + _R * std::sin(angle2),
										_center.Z() + _H);

				// Bottom cap triangle (looking down from above, CCW)
				_surfaces.emplace_back(bottomCenter, b2, b1);

				// Top cap triangle (looking up from below, CCW)
				_surfaces.emplace_back(topCenter, t1, t2);

				// Side surface: two triangles forming a rectangle
				_surfaces.emplace_back(b1, b2, t2);
				_surfaces.emplace_back(b1, t2, t1);
			}
		}

		bool IsInside(const Pnt3Cart& pnt) const override
		{
			// Translate to cylinder-centered coordinates
			Real x = pnt.X() - _center.X();
			Real y = pnt.Y() - _center.Y();
			Real z = pnt.Z() - _center.Z();

			// Check if point is within height range
			if (z < 0.0 || z > _H)
				return false;

			// Check if point is within circular cross-section
			Real distSq = x * x + y * y;
			return distSq <= _R * _R;
		}

		// IBody interface implementation
		Real Volume() const override
		{
			// Volume = π * R² * H
			return Constants::PI * _R * _R * _H;
		}

		Real SurfaceArea() const override
		{
			// Surface area = 2πR² (caps) + 2πRH (lateral surface)
			return 2.0 * Constants::PI * _R * _R + 2.0 * Constants::PI * _R * _H;
		}

		Pnt3Cart GetCenter() const override
		{
			// Return geometric center (middle of cylinder)
			return Pnt3Cart(_center.X(), _center.Y(), _center.Z() + _H / 2.0);
		}

		BoundingBox3D GetBoundingBox() const override
		{
			// AABB with min/max based on radius and height
			Pnt3Cart min(_center.X() - _R, _center.Y() - _R, _center.Z());
			Pnt3Cart max(_center.X() + _R, _center.Y() + _R, _center.Z() + _H);
			return BoundingBox3D(min, max);
		}

		BoundingSphere3D GetBoundingSphere() const override
		{
			// Bounding sphere centered at geometric center
			// Radius is distance from center to corner of bounding cylinder
			Pnt3Cart geomCenter = GetCenter();
			Real halfH = _H / 2.0;
			Real radius = std::sqrt(_R * _R + halfH * halfH);
			return BoundingSphere3D(geomCenter, radius);
		}

		std::string ToString() const override
		{
			std::ostringstream oss;
			Pnt3Cart geomCenter = GetCenter();
			oss << "Cylinder3D{Center=(" << geomCenter.X() << ", " << geomCenter.Y() << ", " << geomCenter.Z()
				<< "), Radius=" << _R << ", Height=" << _H 
				<< ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}

		// Getters
		Real GetRadius() const { return _R; }
		Real GetHeight() const { return _H; }
		Pnt3Cart GetBaseCenter() const { return _center; }
		int GetNumSegments() const { return _numSegments; }
	};

	class Sphere3D : public BodyWithTriangleSurfaces
	{
		Real _R;      // Radius
		Pnt3Cart _center;
		int _numLatitude;   // Number of latitude divisions
		int _numLongitude;  // Number of longitude divisions
	public:
		// Sphere centered at origin with radius R
		Sphere3D(Real R, int numLatitude = 16, int numLongitude = 20)
			: _R(R), _center(0, 0, 0), _numLatitude(numLatitude), _numLongitude(numLongitude)
		{
			constructSurfaces();
		}

		// Sphere centered at specified point
		Sphere3D(Real R, const Pnt3Cart& center, int numLatitude = 16, int numLongitude = 20)
			: _R(R), _center(center), _numLatitude(numLatitude), _numLongitude(numLongitude)
		{
			constructSurfaces();
		}

		void constructSurfaces()
		{
			_surfaces.clear();
			Real dTheta = Constants::PI / _numLatitude;        // Latitude step (0 to π)
			Real dPhi = 2.0 * Constants::PI / _numLongitude;   // Longitude step (0 to 2π)

			// Generate triangles for the sphere surface
			for (int i = 0; i < _numLatitude; ++i)
			{
				Real theta1 = i * dTheta;
				Real theta2 = (i + 1) * dTheta;

				for (int j = 0; j < _numLongitude; ++j)
				{
					Real phi1 = j * dPhi;
					Real phi2 = (j + 1) * dPhi;

					// Four corners of the current patch
					Pnt3Cart p1 = spherePoint(theta1, phi1);
					Pnt3Cart p2 = spherePoint(theta1, phi2);
					Pnt3Cart p3 = spherePoint(theta2, phi2);
					Pnt3Cart p4 = spherePoint(theta2, phi1);

					// Special handling for poles (avoid degenerate triangles)
					if (i == 0)
					{
						// Top pole: only one triangle
						_surfaces.emplace_back(p1, p4, p3);
					}
					else if (i == _numLatitude - 1)
					{
						// Bottom pole: only one triangle
						_surfaces.emplace_back(p1, p2, p3);
					}
					else
					{
						// Regular patch: two triangles
						_surfaces.emplace_back(p1, p2, p3);
						_surfaces.emplace_back(p1, p3, p4);
					}
				}
			}
		}

		Pnt3Cart spherePoint(Real theta, Real phi) const
		{
			// Spherical coordinates: theta (latitude, 0 to π), phi (longitude, 0 to 2π)
			Real x = _R * std::sin(theta) * std::cos(phi);
			Real y = _R * std::sin(theta) * std::sin(phi);
			Real z = _R * std::cos(theta);
			return Pnt3Cart(_center.X() + x, _center.Y() + y, _center.Z() + z);
		}

		bool IsInside(const Pnt3Cart& pnt) const override
		{
			// Check if point is within spherical radius
			Real dx = pnt.X() - _center.X();
			Real dy = pnt.Y() - _center.Y();
			Real dz = pnt.Z() - _center.Z();
			Real distSq = dx * dx + dy * dy + dz * dz;
			return distSq <= _R * _R;
		}

	// IBody interface implementation
	Real Volume() const override
	{
		// Volume = (4/3) * π * R³
		return (4.0 / 3.0) * Constants::PI * _R * _R * _R;
	}

	Real SurfaceArea() const override
	{
		// Surface area = 4 * π * R²
		return 4.0 * Constants::PI * _R * _R;
	}

	Pnt3Cart GetCenter() const override
	{
		return _center;
	}

	BoundingBox3D GetBoundingBox() const override
	{
		// AABB with min/max at center ± radius in all dimensions
		Pnt3Cart min(_center.X() - _R, _center.Y() - _R, _center.Z() - _R);
		Pnt3Cart max(_center.X() + _R, _center.Y() + _R, _center.Z() + _R);
		return BoundingBox3D(min, max);
	}

	BoundingSphere3D GetBoundingSphere() const override
	{
		// Trivial - the sphere itself is the bounding sphere
		return BoundingSphere3D(_center, _R);
	}

	std::string ToString() const override
	{
		std::ostringstream oss;
		oss << "Sphere3D{Center=(" << _center.X() << ", " << _center.Y() << ", " << _center.Z() 
			<< "), Radius=" << _R << ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
		return oss.str();
	}

	// Getters
	Real GetRadius() const { return _R; }
	int GetNumLatitude() const { return _numLatitude; }
	int GetNumLongitude() const { return _numLongitude; }
	};

	class Pyramid3D : public BodyWithTriangleSurfaces
	{
		Real _a;			// base side length
		Real _h;			// pyramid height
		Vec3Cart _center;
	public:
		Pyramid3D(Real a, Real h) : _a(a), _h(h), _center(0, 0, 0)
		{
			Pnt3Cart pnt1( a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt2( a / 2,  a / 2, -a / 2);
			Pnt3Cart pnt3(-a / 2,  a / 2, -a / 2);
			Pnt3Cart pnt4(-a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt5(     0,      0,      h);	// apex of the pyramid
			
			_surfaces.push_back(TriangleSurface3D(pnt1, pnt4, pnt3));     // lower side in xy plane
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt1, pnt2));     // front side in yz plane
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt2, pnt3));     // right side in yz plane
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt3, pnt4));     // back side in yz plane
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt4, pnt1));     // left side in xz plane
		}
		Pyramid3D(Real a, Real h, const Vec3Cart& center) : Pyramid3D(a, h)
		{
			_center = center;
		}

		bool IsInside(const Pnt3Cart& pnt) const override
		{
			// Check if point is inside the pyramid
			// Pyramid has square base [-a/2, a/2] x [-a/2, a/2] at z=0 and apex at (0, 0, h)
			
			// Translate point relative to center
			Real x = pnt.X() - _center.X();
			Real y = pnt.Y() - _center.Y();
			Real z = pnt.Z() - _center.Z();
			
			// Check if point is below base or above apex
			if (z < 0.0 || z > _h)
				return false;
			
			// For a pyramid, at height z, the cross-section is a square
			// The side length decreases linearly from _a at z=0 to 0 at z=_h
			// At height z, half-side = (_a / 2) * (1 - z / _h)
			Real half_side = (_a / 2.0) * (1.0 - z / _h);
			
			// Check if point is within the square cross-section at height z
			return (std::abs(x) <= half_side && std::abs(y) <= half_side);
		}

		// IBody interface implementation
		Real Volume() const override
		{
			// Volume = (1/3) * base_area * height = (1/3) * a² * h
			return (_a * _a * _h) / 3.0;
		}

		Real SurfaceArea() const override
		{
			// Surface area = base + 4 triangular faces
			// Base area = a²
			// Each triangular face: slant height s = sqrt(h² + (a/2)²)
			// Area of one face = (1/2) * a * s
			// Total = a² + 4 * (1/2) * a * s = a² + 2as
			Real slant_height = std::sqrt(_h * _h + (_a / 2.0) * (_a / 2.0));
			return _a * _a + 2.0 * _a * slant_height;
		}

		Pnt3Cart GetCenter() const override
		{
			// Geometric center (centroid) is at 1/4 of height from base
			return Pnt3Cart(_center.X(), _center.Y(), _center.Z() + _h / 4.0);
		}

		BoundingBox3D GetBoundingBox() const override
		{
			// AABB from base corners to apex
			Real half_a = _a / 2.0;
			Pnt3Cart min(_center.X() - half_a, _center.Y() - half_a, _center.Z());
			Pnt3Cart max(_center.X() + half_a, _center.Y() + half_a, _center.Z() + _h);
			return BoundingBox3D(min, max);
		}

		BoundingSphere3D GetBoundingSphere() const override
		{
			// Bounding sphere centered at geometric center
			// Need to ensure it contains both base corners and apex
			Pnt3Cart geomCenter = GetCenter();
			Real half_a = _a / 2.0;
			
			// Distance from centroid (at h/4 from base) to base corner at (±a/2, ±a/2, 0)
			Real distToBaseCorner = std::sqrt(half_a * half_a + half_a * half_a + (_h / 4.0) * (_h / 4.0));
			
			// Distance from centroid (at h/4 from base) to apex at (0, 0, h)
			// dz = h - h/4 = 3h/4
			Real distToApex = (3.0 * _h) / 4.0;
			
			// Use the maximum distance as radius
			Real radius = std::max(distToBaseCorner, distToApex);
			return BoundingSphere3D(geomCenter, radius);
		}

		std::string ToString() const override
		{
			std::ostringstream oss;
			Pnt3Cart geomCenter = GetCenter();
			oss << "Pyramid3D{Center=(" << geomCenter.X() << ", " << geomCenter.Y() << ", " << geomCenter.Z()
				<< "), BaseSize=" << _a << ", Height=" << _h 
				<< ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}

		// Getters
		Real GetBaseSize() const { return _a; }
		Real GetHeight() const { return _h; }
		Vec3Cart GetBaseCenter() const { return _center; }
	};

	class PyramidEquilateral3D : public Pyramid3D
	{
	public:
		PyramidEquilateral3D(Real a) : Pyramid3D(a, a / sqrt(3))
		{ }

		PyramidEquilateral3D(Real a, const Vec3Cart& center) : Pyramid3D(a, a / sqrt(3), center)
		{ }

		std::string ToString() const override
		{
			std::ostringstream oss;
			Vec3Cart baseCenter = GetBaseCenter();
			oss << "PyramidEquilateral3D (Equilateral triangular base): Center=(" 
				<< baseCenter.X() << ", " << baseCenter.Y() << ", " << baseCenter.Z() << ")"
				<< ", BaseSize=" << GetBaseSize()
				<< ", Height=" << GetHeight()
				<< ", Volume=" << Volume()
				<< ", SurfaceArea=" << SurfaceArea();
			return oss.str();
		}
	};
}

#endif
