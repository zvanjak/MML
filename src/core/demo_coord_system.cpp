#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordSystem.h"
#endif

using namespace MML;

// TODO 0.9 - 2D verzije!!!
// carousel - rotirajuci sustav
// 1. simple- carousel promjera 10 km, i na njemu ljudi
// 2. carousel promjera 10 km, i na njemu ljudi, i na njemu je jos jedan carousel promjera 1 km
// na njemu se odvija streljačko natjecanje


void Demo_CoordSystem()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                         COORD SYSTEM                          ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	InertialFrame3D SolarSystemCMFrame;
	Real earthOrbitPeriod = 365 * 24 * 3600;
	CircleOrbitingCartesianFrame3DToCartesian EarthCMFrame(SolarSystemCMFrame, 150000, earthOrbitPeriod);    // radijus kruženja (jedinica 10^6 m, 1000km), period kruženja (s)

	std::cout << "Earth origin at T=0  : " << EarthCMFrame.GetOriginPositionAtTime(0) << std::endl;
	std::cout << "Earth origin at T/4  : " << EarthCMFrame.GetOriginPositionAtTime(earthOrbitPeriod / 4) << std::endl;
	std::cout << "Earth origin at T/2  : " << EarthCMFrame.GetOriginPositionAtTime(earthOrbitPeriod / 2) << std::endl;
	std::cout << "Earth origin at 3T/4T: " << EarthCMFrame.GetOriginPositionAtTime(earthOrbitPeriod * 3 / 4) << std::endl;

	// za bilo koji (x,y,z) u EarthCM frameu, mogu dobiti poziciju u SolarSystem frameu za vrijeme T
	Vector3Cartesian posAtEquator1({ 6.4, 0, 0 });
	Vector3Cartesian posAtEquator2({ -6.4, 0, 0 });
	std::cout << "Pos 1 at T=0: " << EarthCMFrame.GetLocalPositionAtTime(posAtEquator1, 0) << std::endl;
	std::cout << "Pos 1 at T/4: " << EarthCMFrame.GetLocalPositionAtTime(posAtEquator1, earthOrbitPeriod / 4) << std::endl;
	std::cout << "Pos 1 at T/2: " << EarthCMFrame.GetLocalPositionAtTime(posAtEquator1, earthOrbitPeriod / 2) << std::endl;
	std::cout << "Pos 1 at 3T/4T: " << EarthCMFrame.GetLocalPositionAtTime(posAtEquator1, earthOrbitPeriod * 3 / 4) << std::endl;

	// za bilo koji (r,theta,phi) u EarthCM frameu, mogu dobiti poziciju u SolarSystem frameu za vrijeme T
	// to me zanima, jer ce mi tako sustav Zemlje davati svoje lokalne koordinate
	// problem - Zemlja se rotira!

	// sad kreiramo rotirajuci sustav Zemlje, s osi rotacije pod kutom 23 stupnja prema z-osi sustava CM
	// u trenuttku t=0, os rotacije je nagnuta u y-z ravnini
	VectorN<Real, 3> axis({ 0, sin(23 * Constants::PI / 180), cos(23 * Constants::PI / 180) });
	HardSphereRotatingFrame EarthSurfaceFrame(EarthCMFrame, 6.4, 24 * 3600, axis);

	// pozicija na površini zemlje - izračunati u EarthCM i SolarSystemCM frameu
	VectorN<Real, 3> pos({ 0, 0, 0 });
	double latitude = 45;
	double longitude = 16;
	double height = 0;          // on the surface
	double T = 10;

	// pozicija tocke u lokalnom sustavu EarthCMFrame, ali je u SFERNIM koordinatama!!
	Vector3Spherical pos_EarthCM = EarthSurfaceFrame.GetLocalPositionAtTime(latitude, longitude, height, T);
	// koristeci PARENT referentni sustav, transformiramo to u kartezian kord. u EarthCM frameu
	Vector3Cartesian pos_EarthCM_cart = EarthCMFrame.GetLocalPositionAtTime(pos_EarthCM, T);

	std::cout << "Position on Earth surface in EarthCM frame: " << pos_EarthCM_cart << std::endl;

	// is sun visible from this position?

}