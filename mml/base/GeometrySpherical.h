///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GeometrySpherical.h                                                 ///
///  Description: Spherical coordinate geometry (SphericalTriangle, geodesics)        ///
///               Great circle calculations on sphere surface                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_SPHERICAL_H
#define MML_GEOMETRY_SPHERICAL_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/BaseUtils.h"
#include "base/VectorTypes.h"
#include "base/Geometry.h"



namespace MML
{
	class SphericalGeometryCalculator
	{
		public:
		static Real RadiusFromCartesian(const Pnt3Cart& pnt)
		{
			return hypot(pnt.X(), hypot(pnt.Y(), pnt.Z()));
		}
		static Pnt3Cart CartesianFromSpherical(const Pnt3Sph& pnt)
		{
			// Standard spherical coordinates: (r, theta, phi)
			// theta = polar angle from z-axis [0, PI]
			// phi = azimuthal angle in xy-plane [0, 2*PI]
			return Pnt3Cart(pnt.R() * sin(pnt.Theta()) * cos(pnt.Phi()),
											pnt.R() * sin(pnt.Theta()) * sin(pnt.Phi()),
											pnt.R() * cos(pnt.Theta()));
		}

		// given latitude and longitude in degrees, return spherical point
		// latitude: -90 (South) to +90 (North), longitude: -180 to +180
		static Pnt3Sph SphericalFromLatLong(Real latitudeDeg, Real longitudeDeg)
		{
			Real latitudeRad = Utils::DegToRad(latitudeDeg);
			Real longitudeRad = Utils::DegToRad(longitudeDeg);
			// Convert latitude to theta (polar angle from z-axis)
			// latitude 90° (North) -> theta = 0, latitude -90° (South) -> theta = PI
			Real theta = Constants::PI / 2.0 - latitudeRad;
			return Pnt3Sph(1.0, theta, longitudeRad);
		}

		// given two points on sphere, with latitude and longitude in degrees,
		// return distance between them in radians
		static Real DistanceBetweenLatLong(Real lat1Deg, Real long1Deg, Real lat2Deg, Real long2Deg)
		{
			Pnt3Sph pnt1 = SphericalFromLatLong(lat1Deg, long1Deg);
			Pnt3Sph pnt2 = SphericalFromLatLong(lat2Deg, long2Deg);

			// using spherical law of cosines
			Real cosAngle = cos(pnt1.Theta()) * cos(pnt2.Theta()) +
											sin(pnt1.Theta()) * sin(pnt2.Theta()) *
											cos(pnt1.Phi() - pnt2.Phi());
			// Clamp to [-1, 1] to handle numerical precision issues
			cosAngle = std::max(Real(-1.0), std::min(Real(1.0), cosAngle));
			return acos(cosAngle);
		}
		// given two points on sphere of radius R in meters, with latitude and longitude in degrees,
		// return distance between them in meters
		static Real DistanceBetweenLatLong(Real radius, Real lat1Deg, Real long1Deg, Real lat2Deg, Real long2Deg)
		{
      Real angleRad = DistanceBetweenLatLong(lat1Deg, long1Deg, lat2Deg, long2Deg);
      return radius * angleRad;
    }    
	};
}

#endif