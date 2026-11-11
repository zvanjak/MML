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

#include <algorithm>
#include <cmath>

#include <mml/MMLBase.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/BaseUtils.h>
#include <mml/base/Geometry/Geometry.h>



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
		struct DistanceChords {
			Real latitude, longitude, weight, opposite;
		};

		static void ValidateGeographicCoordinates(Real latitude, Real longitude) {
			detail::ValidateGeometryCoordinates({latitude, longitude}, "Spherical geographic coordinates");
			if (latitude < -90 || latitude > 90)
				throw GeometryError("Spherical geographic latitude must be in [-90,90]");
		}

		static Real DegreeChord(Real degrees) {
			// Degree-scaled chords retain tiny input differences before radius multiplication.
			if (std::abs(degrees) <= Precision::RoundoffZeroThreshold) return degrees;
			const Real radiansPerDegree = Constants::PI / 180;
			return 2 * std::sin(degrees * radiansPerDegree / 2) / radiansPerDegree;
		}

		static DistanceChords GeographicChords(Real latitude1, Real longitude1, Real latitude2, Real longitude2) {
			ValidateGeographicCoordinates(latitude1, longitude1);
			ValidateGeographicCoordinates(latitude2, longitude2);
			Real longitudeDifference = longitude2 - longitude1;
			if (std::isfinite(longitudeDifference)) {
				const Real error = detail::GeometryAngleDifferenceError(longitude2, longitude1, longitudeDifference);
				longitudeDifference = std::remainder(longitudeDifference, Real(360)) + error;
			} else {
				longitudeDifference = std::remainder(longitude2, Real(360)) - std::remainder(longitude1, Real(360));
			}
			longitudeDifference = std::remainder(longitudeDifference, Real(360));
			const Real cosine1 = std::abs(latitude1) == 90 ? 0 : std::cos(latitude1 * (Constants::PI / 180));
			const Real cosine2 = std::abs(latitude2) == 90 ? 0 : std::cos(latitude2 * (Constants::PI / 180));
			const Real weight = std::sqrt(cosine1) * std::sqrt(cosine2);
			const Real opposite = std::hypot(DegreeChord(latitude1 + latitude2),
				weight * DegreeChord(180 - std::abs(longitudeDifference)));
			return {DegreeChord(latitude2 - latitude1), DegreeChord(longitudeDifference), weight, opposite};
		}

		public:
		/// @name Coordinate Conversions
		/// @{
		
		/// @brief Compute radial distance from Cartesian coordinates.
    /// @param pnt Cartesian point
    /// @return Distance from origin: √(x² + y² + z²)

		static Real RadiusFromCartesian(const Pnt3Cart& pnt)
		{
			return Pnt3Cart().Dist(pnt);
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
			ValidateGeographicCoordinates(latitudeDeg, longitudeDeg);
			Real latitudeRad = latitudeDeg * (Constants::PI / 180);
			Real longitudeRad = longitudeDeg * (Constants::PI / 180);
			// latitude 90° (North) -> theta = 0, latitude -90° (South) -> theta = PI
			Real theta = Constants::PI / MML::Real(2.0) - latitudeRad;
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
    /// Uses independently evaluated half-angle chords and atan2. The result is the central angle
    /// subtended by the two points. Multiply by sphere radius to get arc length.
    /// @note For Earth: multiply by 6371 km for approximate surface distance.

		static Real DistanceBetweenLatLong(Real lat1Deg, Real long1Deg, Real lat2Deg, Real long2Deg)
		{
			const auto chords = GeographicChords(lat1Deg, long1Deg, lat2Deg, long2Deg);
			const Real direct = std::hypot(chords.latitude, chords.weight * chords.longitude);
			return detail::CheckedGeometryDistance(2 * std::atan2(direct, chords.opposite), "Great-circle angle",
				chords.latitude != 0 || (chords.weight != 0 && chords.longitude != 0));
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
			if (!std::isfinite(radius) || radius < 0)
				throw GeometryError("Great-circle radius must be finite and nonnegative");
			const auto chords = GeographicChords(lat1Deg, long1Deg, lat2Deg, long2Deg);
			if (radius == 0) return 0;
			const Real direct = std::hypot(chords.latitude, chords.weight * chords.longitude);
			Real distance;
			if (direct <= Precision::RoundoffZeroThreshold * chords.opposite) {
				// atan(x)=x up to relative O(x^2); scale before tiny angular components underflow.
				const Real scale = radius / (chords.opposite / 2);
				distance = std::hypot(scale * chords.latitude, (scale * chords.weight) * chords.longitude);
			} else {
				distance = radius * (2 * std::atan2(direct, chords.opposite));
			}
			return detail::CheckedGeometryDistance(distance, "Great-circle arc length",
				chords.latitude != 0 || (chords.weight != 0 && chords.longitude != 0));
		}
		/// @}
	};
}

#endif