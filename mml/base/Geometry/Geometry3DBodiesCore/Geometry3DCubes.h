///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DCubes.h                                                   ///
///  Description: 3D cube classes (quad and triangle mesh versions)                   ///
///               Cube3D, CubeWithTriangles3D                                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry3DCubes.h
/// @brief Cube solid body representations.
/// - Cube3D: Axis-aligned cube with 6 rectangular (quad) faces
/// - CubeWithTriangles3D: Axis-aligned cube with 12 triangular faces
/// @see Geometry3DBodies.h for the aggregate header

#if !defined MML_GEOMETRY_3D_CUBES_H
#define MML_GEOMETRY_3D_CUBES_H

#include <cmath>
#include <sstream>

#include "mml/MMLBase.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry3D.h"
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBounding.h"
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBodyBase.h"

namespace MML {

	/// @brief Axis-aligned cube with rectangular (quad) surface mesh.
	/// Represents a cube with 6 rectangular faces (RectSurface3D).
	/// The cube is axis-aligned with faces parallel to coordinate planes.
	/// @note Vertex ordering produces outward-facing normals.
	/// @see CubeWithTriangles3D for triangulated version
	/// @see Box3D for AABB representation

	class Cube3D : public BodyWithRectSurfaces {
		Real _a;		  ///< Side length
		Pnt3Cart _center; ///< Center point
	public:
		/// @name Constructors
		/// @{

		/// @brief Construct cube centered at origin.
		/// @param a Side length

		Cube3D(Real a)
			: _a(a)
			, _center(0, 0, 0) {
			Pnt3Cart pnt1(a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt2(a / 2, a / 2, -a / 2);
			Pnt3Cart pnt3(-a / 2, a / 2, -a / 2);
			Pnt3Cart pnt4(-a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt5(a / 2, -a / 2, a / 2);
			Pnt3Cart pnt6(a / 2, a / 2, a / 2);
			Pnt3Cart pnt7(-a / 2, a / 2, a / 2);
			Pnt3Cart pnt8(-a / 2, -a / 2, a / 2);

			// Add all 6 faces
			_surfaces.push_back(RectSurface3D(pnt1, pnt4, pnt3, pnt2)); // lower side in xy plane
			_surfaces.push_back(RectSurface3D(pnt5, pnt6, pnt7, pnt8)); // upper side in xy plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt2, pnt6, pnt5)); // front side in yz plane
			_surfaces.push_back(RectSurface3D(pnt4, pnt8, pnt7, pnt3)); // back side in yz plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt5, pnt8, pnt4)); // left side in xz plane
			_surfaces.push_back(RectSurface3D(pnt2, pnt3, pnt7, pnt6)); // right side in xz plane
		}

		/// @brief Construct cube centered at specified point.
		/// @param a Side length
		/// @param center Center point of the cube

		Cube3D(Real a, const Pnt3Cart& center)
			: _a(a)
			, _center(center) {
			// Create 8 corner vertices centered at the given point
			Pnt3Cart pnt1(center.X() + a / 2, center.Y() - a / 2, center.Z() - a / 2);
			Pnt3Cart pnt2(center.X() + a / 2, center.Y() + a / 2, center.Z() - a / 2);
			Pnt3Cart pnt3(center.X() - a / 2, center.Y() + a / 2, center.Z() - a / 2);
			Pnt3Cart pnt4(center.X() - a / 2, center.Y() - a / 2, center.Z() - a / 2);
			Pnt3Cart pnt5(center.X() + a / 2, center.Y() - a / 2, center.Z() + a / 2);
			Pnt3Cart pnt6(center.X() + a / 2, center.Y() + a / 2, center.Z() + a / 2);
			Pnt3Cart pnt7(center.X() - a / 2, center.Y() + a / 2, center.Z() + a / 2);
			Pnt3Cart pnt8(center.X() - a / 2, center.Y() - a / 2, center.Z() + a / 2);

			// Add all 6 faces
			_surfaces.push_back(RectSurface3D(pnt1, pnt4, pnt3, pnt2)); // lower side in xy plane
			_surfaces.push_back(RectSurface3D(pnt5, pnt6, pnt7, pnt8)); // upper side in xy plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt2, pnt6, pnt5)); // front side in yz plane
			_surfaces.push_back(RectSurface3D(pnt4, pnt8, pnt7, pnt3)); // back side in yz plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt5, pnt8, pnt4)); // left side in xz plane
			_surfaces.push_back(RectSurface3D(pnt2, pnt3, pnt7, pnt6)); // right side in xz plane
		}
		/// @}

		/// @name Containment Test
		/// @{

		bool IsInside(const Pnt3Cart& pnt) const override {
			Real half_a = _a / 2.0;
			return (pnt.X() >= _center.X() - half_a && pnt.X() <= _center.X() + half_a && pnt.Y() >= _center.Y() - half_a &&
					pnt.Y() <= _center.Y() + half_a && pnt.Z() >= _center.Z() - half_a && pnt.Z() <= _center.Z() + half_a);
		}
		/// @}

		/// @name IBody Interface
		/// @{

		/// /** @brief Volume = a³ */

		Real Volume() const override { return _a * _a * _a; }

		/// /** @brief Surface area = 6a² */

		Real SurfaceArea() const override { return 6.0 * _a * _a; }

		Pnt3Cart GetCenter() const override { return _center; }

		Box3D GetBoundingBox() const override {
			Real half_a = _a / 2.0;
			Pnt3Cart min(_center.X() - half_a, _center.Y() - half_a, _center.Z() - half_a);
			Pnt3Cart max(_center.X() + half_a, _center.Y() + half_a, _center.Z() + half_a);
			return Box3D(min, max);
		}

		/// /** @brief Bounding sphere radius = (a√3)/2 (half space diagonal) */

		BoundingSphere3D GetBoundingSphere() const override {
			Real radius = _a * std::sqrt(3.0) / 2.0;
			return BoundingSphere3D(_center, radius);
		}

		std::string ToString() const override {
			std::ostringstream oss;
			oss << "Cube3D{Center=(" << _center.X() << ", " << _center.Y() << ", " << _center.Z() << "), Side=" << _a
				<< ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}
		/// @}

		/// @name Accessors
		/// @{

		/// /** @brief Get side length. */

		Real GetSide() const { return _a; }
		/// @}
	};

	/// @brief Axis-aligned cube with triangular surface mesh.
	/// Represents a cube with 12 triangular faces (2 per quad face).
	/// Triangle winding is CCW for outward-facing normals.
	/// @see Cube3D for quad-face version

	class CubeWithTriangles3D : public BodyWithTriangleSurfaces {
		Real _a;		  ///< Side length
		Pnt3Cart _center; ///< Center point
	public:
		/// @name Constructors
		/// @{

		/// @brief Construct cube centered at origin.
		/// @param a Side length

		CubeWithTriangles3D(Real a)
			: _a(a)
			, _center(0, 0, 0) {
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

		/// @brief Construct cube centered at specified point.
		/// @param a Side length
		/// @param center Center point of the cube

		CubeWithTriangles3D(Real a, const Pnt3Cart& center)
			: _a(a)
			, _center(center) {
			Real h = a / 2.0;
			// 8 vertices centered at _center
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
		/// @}

		/// @name Containment Test
		/// @{

		bool IsInside(const Pnt3Cart& pnt) const override {
			Real h = _a / 2.0;
			return (pnt.X() >= _center.X() - h && pnt.X() <= _center.X() + h && pnt.Y() >= _center.Y() - h && pnt.Y() <= _center.Y() + h &&
					pnt.Z() >= _center.Z() - h && pnt.Z() <= _center.Z() + h);
		}
		/// @}

		/// @name IBody Interface
		/// @{

		/// /** @brief Volume = a³ */

		Real Volume() const override { return _a * _a * _a; }

		/// /** @brief Surface area = 6a² */

		Real SurfaceArea() const override { return 6.0 * _a * _a; }

		Pnt3Cart GetCenter() const override { return _center; }

		Box3D GetBoundingBox() const override {
			Real h = _a / 2.0;
			Pnt3Cart min(_center.X() - h, _center.Y() - h, _center.Z() - h);
			Pnt3Cart max(_center.X() + h, _center.Y() + h, _center.Z() + h);
			return Box3D(min, max);
		}

		/// /** @brief Bounding sphere radius = (a/2)√3 (corner distance) */

		BoundingSphere3D GetBoundingSphere() const override {
			Real radius = (_a / 2.0) * std::sqrt(3.0);
			return BoundingSphere3D(_center, radius);
		}

		std::string ToString() const override {
			std::ostringstream oss;
			oss << "CubeWithTriangles3D: Center=(" << _center.X() << ", " << _center.Y() << ", " << _center.Z() << ")"
				<< ", Side=" << _a << ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << ", Triangles=" << _surfaces.size();
			return oss.str();
		}
		/// @}

		/// @name Accessors
		/// @{

		/// /** @brief Get side length. */

		Real GetSide() const { return _a; }
		/// @}
	};

} // namespace MML

#endif
