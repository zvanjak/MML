///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GeometrySpherical.h                                                 ///
///  Description: Spherical coordinate geometry (SphericalTriangle, geodesics)        ///
///               Great circle calculations on sphere surface                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file GeometrySpherical.h
/// @brief Spherical geometry calculations and coordinate conversions.
/// Provides utilities for working with spherical coordinates including:
/// - Coordinate conversions (Cartesian ↔ Spherical ↔ Lat/Long)
/// - Great circle distance calculations (geodesics)
/// - Geographic coordinate support (latitude/longitude)
/// @section coord_conventions Coordinate Conventions
/// **Spherical coordinates (r, θ, φ):**
/// - r: radial distance from origin
/// - θ (theta): polar angle from +Z axis [0, π]
/// - φ (phi): azimuthal angle in XY plane from +X axis [0, 2π]
/// **Geographic coordinates:**
/// - Latitude: -90° (South Pole) to +90° (North Pole)
/// - Longitude: -180° to +180° (East positive)
/// @section formulas Key Formulas
/// **Cartesian from Spherical:**
/// @f[
/// x = r \sin\theta \cos\phi, \quad
/// y = r \sin\theta \sin\phi, \quad
/// z = r \cos\theta
/// @f]
/// **Spherical Law of Cosines (angular distance):**
/// @f[
/// \cos c = \cos\theta_1 \cos\theta_2 + \sin\theta_1 \sin\theta_2 \cos(\phi_1 - \phi_2)
/// @f]
/// @see Pnt3Sph for spherical point representation
/// @see Pnt3Cart for Cartesian point representation


#if !defined MML_GEOMETRY_SPHERICAL_H
#define MML_GEOMETRY_SPHERICAL_H

// Standard headers - include what we use
#include <algorithm>
#include <cmath>

#include "mml/MMLBase.h"

#include "mml/interfaces/IFunction.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry.h"



namespace MML
{
	/// @brief Static utility class for spherical geometry calculations.
/// Provides coordinate conversions and geodesic distance calculations
/// for points on a sphere. Supports both mathematical spherical coordinates
/// and geographic (latitude/longitude) coordinates.
/// @note All angle inputs/outputs are in radians unless explicitly noted.
/// @note Distance functions return angular distance; multiply by radius for arc length.

	class SphericalGeometryCalculator
	{
		public:
		/// @name Coordinate Conversions
		/// @{
		
		/// @brief Compute radial distance from Cartesian coordinates.
/// @param pnt Cartesian point
/// @return Distance from origin: √(x² + y² + z²)

		static Real RadiusFromCartesian(const Pnt3Cart& pnt)
		{
			return hypot(pnt.X(), hypot(pnt.Y(), pnt.Z()));
		}
		
		/// @brief Convert spherical coordinates to Cartesian.
/// @param pnt Spherical point (r, θ, φ)
/// @return Cartesian point (x, y, z)
/// Uses standard physics convention:
/// - θ (theta) = polar angle from +Z axis [0, π]
/// - φ (phi) = azimuthal angle in XY plane [0, 2π]

		static Pnt3Cart CartesianFromSpherical(const Pnt3Sph& pnt)
		{
			return Pnt3Cart(pnt.R() * sin(pnt.Theta()) * cos(pnt.Phi()),
											pnt.R() * sin(pnt.Theta()) * sin(pnt.Phi()),
											pnt.R() * cos(pnt.Theta()));
		}

		/// @brief Convert geographic coordinates to spherical (unit sphere).
/// @param latitudeDeg Latitude in degrees [-90 (South) to +90 (North)]
/// @param longitudeDeg Longitude in degrees [-180 to +180]
/// @return Spherical point on unit sphere (r=1)
/// Conversion: θ = 90° - latitude (so North Pole has θ=0)

		static Pnt3Sph SphericalFromLatLong(Real latitudeDeg, Real longitudeDeg)
		{
			Real latitudeRad = Utils::DegToRad(latitudeDeg);
			Real longitudeRad = Utils::DegToRad(longitudeDeg);
			// latitude 90° (North) -> theta = 0, latitude -90° (South) -> theta = PI
			Real theta = Constants::PI / 2.0 - latitudeRad;
			return Pnt3Sph(1.0, theta, longitudeRad);
		}
		/// @}

		/// @name Geodesic Distances
		/// @{
		
		/// @brief Compute great circle angular distance between two geographic points.
/// @param lat1Deg First point latitude (degrees)
/// @param long1Deg First point longitude (degrees)
/// @param lat2Deg Second point latitude (degrees)
/// @param long2Deg Second point longitude (degrees)
/// @return Angular distance in radians
/// Uses spherical law of cosines. The result is the central angle
/// subtended by the two points. Multiply by sphere radius to get arc length.
/// @note For Earth: multiply by 6371 km for approximate surface distance.

		static Real DistanceBetweenLatLong(Real lat1Deg, Real long1Deg, Real lat2Deg, Real long2Deg)
		{
			Pnt3Sph pnt1 = SphericalFromLatLong(lat1Deg, long1Deg);
			Pnt3Sph pnt2 = SphericalFromLatLong(lat2Deg, long2Deg);

			// Spherical law of cosines
			Real cosAngle = cos(pnt1.Theta()) * cos(pnt2.Theta()) +
											sin(pnt1.Theta()) * sin(pnt2.Theta()) *
											cos(pnt1.Phi() - pnt2.Phi());
			// Clamp to [-1, 1] for numerical stability
			cosAngle = std::max(Real(-1.0), std::min(Real(1.0), cosAngle));
			return acos(cosAngle);
		}
		
		/// @brief Compute great circle distance in physical units.
/// @param radius Sphere radius (in desired units, e.g., meters or km)
/// @param lat1Deg First point latitude (degrees)
/// @param long1Deg First point longitude (degrees)
/// @param lat2Deg Second point latitude (degrees)
/// @param long2Deg Second point longitude (degrees)
/// @return Arc length distance (same units as radius)
/// @note For Earth calculations, use radius ≈ 6371000 m or 6371 km.

		static Real DistanceBetweenLatLong(Real radius, Real lat1Deg, Real long1Deg, Real lat2Deg, Real long2Deg)
		{
      Real angleRad = DistanceBetweenLatLong(lat1Deg, long1Deg, lat2Deg, long2Deg);
      return radius * angleRad;
    }
		/// @}
	};
}

#endif