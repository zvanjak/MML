///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DSpherePyramid.h                                           ///
///  Description: Sphere and Pyramid 3D solid classes                                 ///
///               UV-sphere and square-based pyramid representations                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry3DSpherePyramid.h
/// @brief Sphere and Pyramid solid body representations.
/// - Sphere3D: UV-sphere with triangular surface mesh
/// - Pyramid3D: Square-based pyramid with triangular faces
/// - PyramidEquilateral3D: Equilateral pyramid (all edges equal)
/// @see Geometry3DBodies.h for the aggregate header

#if !defined MML_GEOMETRY_3D_SPHERE_PYRAMID_H
#define MML_GEOMETRY_3D_SPHERE_PYRAMID_H

#include <cmath>
#include <sstream>

#include "mml/MMLBase.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry3D.h"
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBounding.h"
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBodyBase.h"

namespace MML {

	/// @brief Sphere with triangular surface mesh.
	/// UV-sphere parameterization with latitude/longitude divisions.
	/// Special handling at poles to avoid degenerate triangles.
	/// @note Formulas:
	/// - Volume = (4/3)πR³
	/// - Surface Area = 4πR²

	class Sphere3D : public BodyWithTriangleSurfaces {
		Real _R;		   ///< Radius
		Pnt3Cart _center;  ///< Center point
		int _numLatitude;  ///< Latitude divisions (θ: 0 to π)
		int _numLongitude; ///< Longitude divisions (φ: 0 to 2π)
	public:
		/// @name Constructors
		/// @{

		/// @brief Construct sphere centered at origin.
		/// @param R Radius
		/// @param numLatitude Latitude divisions (default 16)
		/// @param numLongitude Longitude divisions (default 20)

		Sphere3D(Real R, int numLatitude = 16, int numLongitude = 20)
			: _R(R)
			, _center(0, 0, 0)
			, _numLatitude(numLatitude)
			, _numLongitude(numLongitude) {
			constructSurfaces();
		}

		/// @brief Construct sphere centered at specified point.
		/// @param R Radius
		/// @param center Center point
		/// @param numLatitude Latitude divisions (default 16)
		/// @param numLongitude Longitude divisions (default 20)

		Sphere3D(Real R, const Pnt3Cart& center, int numLatitude = 16, int numLongitude = 20)
			: _R(R)
			, _center(center)
			, _numLatitude(numLatitude)
			, _numLongitude(numLongitude) {
			constructSurfaces();
		}
		/// @}

		/// @name Surface Generation
		/// @{

		/// /** @brief Generate triangles for the sphere surface. */

		void constructSurfaces() {
			_surfaces.clear();
			Real dTheta = Constants::PI / _numLatitude;
			Real dPhi = 2.0 * Constants::PI / _numLongitude;

			for (int i = 0; i < _numLatitude; ++i) {
				Real theta1 = i * dTheta;
				Real theta2 = (i + 1) * dTheta;

				for (int j = 0; j < _numLongitude; ++j) {
					Real phi1 = j * dPhi;
					Real phi2 = (j + 1) * dPhi;

					Pnt3Cart p1 = spherePoint(theta1, phi1);
					Pnt3Cart p2 = spherePoint(theta1, phi2);
					Pnt3Cart p3 = spherePoint(theta2, phi2);
					Pnt3Cart p4 = spherePoint(theta2, phi1);

					// Special pole handling (avoid degenerate triangles)
					if (i == 0) {
						// Top pole: single triangle
						_surfaces.emplace_back(p1, p4, p3);
					} else if (i == _numLatitude - 1) {
						// Bottom pole: single triangle
						_surfaces.emplace_back(p1, p2, p3);
					} else {
						// Regular patch: two triangles
						_surfaces.emplace_back(p1, p2, p3);
						_surfaces.emplace_back(p1, p3, p4);
					}
				}
			}
		}

		/// @brief Compute point on sphere surface.
		/// @param theta Latitude angle [0, π] (0 = north pole)
		/// @param phi Longitude angle [0, 2π]

		Pnt3Cart spherePoint(Real theta, Real phi) const {
			Real x = _R * std::sin(theta) * std::cos(phi);
			Real y = _R * std::sin(theta) * std::sin(phi);
			Real z = _R * std::cos(theta);
			return Pnt3Cart(_center.X() + x, _center.Y() + y, _center.Z() + z);
		}
		/// @}

		/// @name Containment Test
		/// @{

		bool IsInside(const Pnt3Cart& pnt) const override {
			Real dx = pnt.X() - _center.X();
			Real dy = pnt.Y() - _center.Y();
			Real dz = pnt.Z() - _center.Z();
			Real distSq = dx * dx + dy * dy + dz * dz;
			return distSq <= _R * _R;
		}
		/// @}

		/// @name IBody Interface
		/// @{

		/// /** @brief Volume = (4/3)πR³ */

		Real Volume() const override { return (4.0 / 3.0) * Constants::PI * _R * _R * _R; }

		/// /** @brief Surface area = 4πR² */

		Real SurfaceArea() const override { return 4.0 * Constants::PI * _R * _R; }

		Pnt3Cart GetCenter() const override { return _center; }

		Box3D GetBoundingBox() const override {
			Pnt3Cart min(_center.X() - _R, _center.Y() - _R, _center.Z() - _R);
			Pnt3Cart max(_center.X() + _R, _center.Y() + _R, _center.Z() + _R);
			return Box3D(min, max);
		}

		/// /** @brief Bounding sphere = the sphere itself. */

		BoundingSphere3D GetBoundingSphere() const override { return BoundingSphere3D(_center, _R); }

		std::string ToString() const override {
			std::ostringstream oss;
			oss << "Sphere3D{Center=(" << _center.X() << ", " << _center.Y() << ", " << _center.Z() << "), Radius=" << _R
				<< ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}
		/// @}

		/// @name Accessors
		/// @{

		Real GetRadius() const { return _R; }
		int GetNumLatitude() const { return _numLatitude; }
		int GetNumLongitude() const { return _numLongitude; }
		/// @}
	};

	/// @brief Square-based pyramid with triangular surface mesh.
	/// Pyramid with square base in XY plane and apex at (0, 0, h).
	/// Base is centered at origin (or specified center), extending ±a/2 in X and Y.
	/// @note Formulas:
	/// - Volume = (1/3)a²h
	/// - Surface Area = a² + 2as (s = slant height = √(h² + (a/2)²))

	class Pyramid3D : public BodyWithTriangleSurfaces {
		Real _a;		  ///< Base side length
		Real _h;		  ///< Pyramid height
		Vec3Cart _center; ///< Base center
	public:
		/// @name Constructors
		/// @{

		/// @brief Construct pyramid with base centered at origin.
		/// @param a Base side length
		/// @param h Height (apex at z = h)

		Pyramid3D(Real a, Real h)
			: _a(a)
			, _h(h)
			, _center(0, 0, 0) {
			Pnt3Cart pnt1(a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt2(a / 2, a / 2, -a / 2);
			Pnt3Cart pnt3(-a / 2, a / 2, -a / 2);
			Pnt3Cart pnt4(-a / 2, -a / 2, -a / 2);
			Pnt3Cart pnt5(0, 0, h); // apex

			_surfaces.push_back(TriangleSurface3D(pnt1, pnt4, pnt3)); // base
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt1, pnt2)); // front face
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt2, pnt3)); // right face
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt3, pnt4)); // back face
			_surfaces.push_back(TriangleSurface3D(pnt5, pnt4, pnt1)); // left face
		}

		/// @brief Construct pyramid with base centered at specified point.
		/// @param a Base side length
		/// @param h Height
		/// @param center Base center point

		Pyramid3D(Real a, Real h, const Vec3Cart& center)
			: Pyramid3D(a, h) {
			_center = center;
		}
		/// @}

		/// @name Containment Test
		/// @{

		bool IsInside(const Pnt3Cart& pnt) const override {
			// Translate relative to center
			Real x = pnt.X() - _center.X();
			Real y = pnt.Y() - _center.Y();
			Real z = pnt.Z() - _center.Z();

			// Height bounds
			if (z < 0.0 || z > _h)
				return false;

			// Cross-section shrinks linearly from base to apex
			Real half_side = (_a / 2.0) * (1.0 - z / _h);

			return (std::abs(x) <= half_side && std::abs(y) <= half_side);
		}
		/// @}

		/// @name IBody Interface
		/// @{

		/// /** @brief Volume = (1/3)a²h */

		Real Volume() const override { return (_a * _a * _h) / 3.0; }

		/// /** @brief Surface area = a² + 2as where s = √(h² + (a/2)²) */

		Real SurfaceArea() const override {
			Real slant_height = std::sqrt(_h * _h + (_a / 2.0) * (_a / 2.0));
			return _a * _a + 2.0 * _a * slant_height;
		}

		/// /** @brief Centroid at 1/4 height from base. */

		Pnt3Cart GetCenter() const override { return Pnt3Cart(_center.X(), _center.Y(), _center.Z() + _h / 4.0); }

		Box3D GetBoundingBox() const override {
			Real half_a = _a / 2.0;
			Pnt3Cart min(_center.X() - half_a, _center.Y() - half_a, _center.Z());
			Pnt3Cart max(_center.X() + half_a, _center.Y() + half_a, _center.Z() + _h);
			return Box3D(min, max);
		}

		BoundingSphere3D GetBoundingSphere() const override {
			Pnt3Cart geomCenter = GetCenter();
			Real half_a = _a / 2.0;

			// Distance from centroid to base corner
			Real distToBaseCorner = std::sqrt(half_a * half_a + half_a * half_a + (_h / 4.0) * (_h / 4.0));

			// Distance from centroid to apex (dz = 3h/4)
			Real distToApex = (3.0 * _h) / 4.0;

			Real radius = std::max(distToBaseCorner, distToApex);
			return BoundingSphere3D(geomCenter, radius);
		}

		std::string ToString() const override {
			std::ostringstream oss;
			Pnt3Cart geomCenter = GetCenter();
			oss << "Pyramid3D{Center=(" << geomCenter.X() << ", " << geomCenter.Y() << ", " << geomCenter.Z() << "), BaseSize=" << _a
				<< ", Height=" << _h << ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}
		/// @}

		/// @name Accessors
		/// @{

		Real GetBaseSize() const { return _a; }
		Real GetHeight() const { return _h; }
		Vec3Cart GetBaseCenter() const { return _center; }
		/// @}
	};

	/// @brief Equilateral pyramid (height = a/√3).
	/// Special case of Pyramid3D with all edges equal length.

	class PyramidEquilateral3D : public Pyramid3D {
	public:
		/// @brief Construct equilateral pyramid centered at origin.
		/// @param a Base side length (height = a/√3)

		PyramidEquilateral3D(Real a)
			: Pyramid3D(a, a / sqrt(3)) {}

		/// @brief Construct equilateral pyramid centered at specified point.
		/// @param a Base side length
		/// @param center Base center point

		PyramidEquilateral3D(Real a, const Vec3Cart& center)
			: Pyramid3D(a, a / sqrt(3), center) {}

		std::string ToString() const override {
			std::ostringstream oss;
			Vec3Cart baseCenter = GetBaseCenter();
			oss << "PyramidEquilateral3D (Equilateral triangular base): Center=(" << baseCenter.X() << ", " << baseCenter.Y() << ", "
				<< baseCenter.Z() << ")"
				<< ", BaseSize=" << GetBaseSize() << ", Height=" << GetHeight() << ", Volume=" << Volume()
				<< ", SurfaceArea=" << SurfaceArea();
			return oss.str();
		}
	};

} // namespace MML

#endif
