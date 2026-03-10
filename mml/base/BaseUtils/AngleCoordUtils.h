///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        AngleCoordUtils.h                                                   ///
///  Description: Angle conversion and coordinate transformation utilities            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_ANGLE_COORD_UTILS_H
#define MML_ANGLE_COORD_UTILS_H

#include <cmath>
#include "mml/MMLBase.h"

namespace MML
{
	namespace Utils
	{
		// ============================================================================
		// Angle Conversions
		// ============================================================================

		/// @brief Converts angle from degrees to radians.
		/// @param angleDeg Angle in degrees
		/// @return Angle in radians
		static constexpr Real DegToRad(Real angleDeg) noexcept { return angleDeg * Constants::PI / 180.0; }

		/// @brief Converts angle from radians to degrees.
		/// @param angleRad Angle in radians
		/// @return Angle in degrees
		static constexpr Real RadToDeg(Real angleRad) noexcept { return angleRad * 180.0 / Constants::PI; }

		/// @brief Converts decimal degrees to degrees, minutes, seconds.
		/// @param angle Angle in decimal degrees
		/// @param[out] deg Integer degrees
		/// @param[out] min Integer minutes
		/// @param[out] sec Seconds (fractional)
		static void AngleDegToExplicit(Real angle, Real& deg, Real& min, Real& sec)
		{
			deg = floor(angle);
			min = floor((angle - deg) * 60.0);
			sec = (angle - deg - min / 60.0) * 3600.0;
		}

		/// @brief Converts radians to degrees, minutes, seconds.
		/// @param angleRad Angle in radians
		/// @param[out] deg Integer degrees
		/// @param[out] min Integer minutes
		/// @param[out] sec Seconds (fractional)
		static void AngleRadToExplicit(Real angleRad, Real& deg, Real& min, Real& sec) 
		{ 
			AngleDegToExplicit(angleRad * 180.0 / Constants::PI, deg, min, sec); 
		}

		/// @brief Converts degrees, minutes, seconds to decimal degrees.
		/// @param deg Integer degrees
		/// @param min Integer minutes
		/// @param sec Seconds
		/// @return Angle in decimal degrees
		static constexpr Real ExplicitToAngleDeg(Real deg, Real min, Real sec) noexcept 
		{ 
			return deg + min / 60.0 + sec / 3600.0; 
		}

		/// @brief Converts degrees, minutes, seconds to radians.
		/// @param deg Integer degrees
		/// @param min Integer minutes
		/// @param sec Seconds
		/// @return Angle in radians
		static constexpr Real ExplicitToAngleRad(Real deg, Real min, Real sec) noexcept 
		{ 
			return ExplicitToAngleDeg(deg, min, sec) * Constants::PI / 180.0; 
		}

		/// @brief Normalizes angle to [0, 2π) range.
		/// @param rad Angle in radians
		/// @return Normalized angle in [0, 2π)
		static Real AngleTo2PiRange(Real rad)
		{
			while (rad < 0)
				rad += 2 * Constants::PI;
			while (rad >= 2 * Constants::PI)
				rad -= 2 * Constants::PI;
			return rad;
		}

		/// @brief Normalizes angle to [-π, π) range.
		/// @param rad Angle in radians
		/// @return Normalized angle in [-π, π)
		static Real AngleToPiPiRange(Real rad)
		{
			while (rad < -Constants::PI)
				rad += 2 * Constants::PI;
			while (rad >= Constants::PI)
				rad -= 2 * Constants::PI;
			return rad;
		}

		// ============================================================================
		// Coordinate Transformations
		// ============================================================================
		/// @details Basic coordinate transformation functions following Math/ISO conventions.
		/// These are the underlying formulas used by CoordTransf classes in Core layer.
		///
		/// 2D: Cartesian (x, y) <-> Polar (r, phi)
		///     phi = azimuthal angle from positive x-axis
		///
		/// 3D Spherical (r, theta, phi) - Math/ISO convention:
		///     r     = radial distance from origin
		///     theta = polar angle (inclination) from positive z-axis, θ ∈ [0, π]
		///     phi   = azimuthal angle in xy-plane from positive x-axis
		///
		/// 3D Cylindrical (r, phi, z):
		///     r   = distance from z-axis
		///     phi = azimuthal angle in xy-plane from positive x-axis
		///     z   = height along z-axis

		/// @brief Converts 2D Cartesian coordinates to polar coordinates.
		/// @param x X-coordinate
		/// @param y Y-coordinate
		/// @param[out] r Radial distance from origin
		/// @param[out] phi Azimuthal angle from positive x-axis (radians)
		static void CartesianToPolar(Real x, Real y, Real& r, Real& phi)
		{
			r   = std::sqrt(x * x + y * y);
			phi = std::atan2(y, x);
		}

		/// @brief Converts 2D polar coordinates to Cartesian coordinates.
		/// @param r Radial distance from origin
		/// @param phi Azimuthal angle from positive x-axis (radians)
		/// @param[out] x X-coordinate
		/// @param[out] y Y-coordinate
		static void PolarToCartesian(Real r, Real phi, Real& x, Real& y)
		{
			x = r * std::cos(phi);
			y = r * std::sin(phi);
		}

		/// @brief Converts 3D Cartesian coordinates to spherical coordinates (Math/ISO convention).
		/// @param x X-coordinate
		/// @param y Y-coordinate
		/// @param z Z-coordinate
		/// @param[out] r Radial distance from origin
		/// @param[out] theta Polar angle from positive z-axis (radians, [0, π])
		/// @param[out] phi Azimuthal angle in xy-plane from positive x-axis (radians)
		static void CartesianToSpherical(Real x, Real y, Real z, Real& r, Real& theta, Real& phi)
		{
			r = std::sqrt(x * x + y * y + z * z);
			if (r == 0.0) {
				theta = 0.0;
				phi   = 0.0;
			} else {
				theta = std::acos(z / r);
				phi   = std::atan2(y, x);
			}
		}

		/// @brief Converts 3D spherical coordinates to Cartesian coordinates (Math/ISO convention).
		/// @param r Radial distance from origin
		/// @param theta Polar angle from positive z-axis (radians, [0, π])
		/// @param phi Azimuthal angle in xy-plane from positive x-axis (radians)
		/// @param[out] x X-coordinate
		/// @param[out] y Y-coordinate
		/// @param[out] z Z-coordinate
		static void SphericalToCartesian(Real r, Real theta, Real phi, Real& x, Real& y, Real& z)
		{
			x = r * std::sin(theta) * std::cos(phi);
			y = r * std::sin(theta) * std::sin(phi);
			z = r * std::cos(theta);
		}

		/// @brief Converts 3D Cartesian coordinates to cylindrical coordinates.
		/// @param x X-coordinate
		/// @param y Y-coordinate
		/// @param z_in Z-coordinate (input)
		/// @param[out] r Radial distance from z-axis
		/// @param[out] phi Azimuthal angle in xy-plane from positive x-axis (radians)
		/// @param[out] z_out Z-coordinate (output, same as z_in)
		static void CartesianToCylindrical(Real x, Real y, Real z_in, Real& r, Real& phi, Real& z_out)
		{
			r     = std::sqrt(x * x + y * y);
			phi   = std::atan2(y, x);
			z_out = z_in;
		}

		/// @brief Converts 3D cylindrical coordinates to Cartesian coordinates.
		/// @param r Radial distance from z-axis
		/// @param phi Azimuthal angle in xy-plane from positive x-axis (radians)
		/// @param z_in Z-coordinate (input)
		/// @param[out] x X-coordinate
		/// @param[out] y Y-coordinate
		/// @param[out] z_out Z-coordinate (output, same as z_in)
		static void CylindricalToCartesian(Real r, Real phi, Real z_in, Real& x, Real& y, Real& z_out)
		{
			x     = r * std::cos(phi);
			y     = r * std::sin(phi);
			z_out = z_in;
		}

	} // namespace Utils
} // namespace MML

#endif // MML_ANGLE_COORD_UTILS_H
