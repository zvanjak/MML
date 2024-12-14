#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "base/VectorN.h"

#include "core/CoordSystem.h"
#include "core/Derivation.h"
#endif


using namespace MML;

// Simulacija referencijalnih sistema, s ukljucenim efektima specijalne relativnosti
//   twin paradoks, kosi hitac, kretanje u orbiti, Coriolis, rotacija Zemlje

// TODO 0.9
/*
0. odrediti tocne koordinate!
1. lokalna ravnina za prva tri slucaja, i onda vec za treci Coriolis
2. 3 slucaja gdje je Coriolis presudan
3. 2 slucaja di smo u orbiti
4. 3 slucaja di je specijalna relativnost presudna
*/

void Example1_kosi_hitac()
{
  std::cout << "***********************************************************************" << std::endl;
  std::cout << "****                     EXAMPLE 1 - kosi hitac                    ****" << std::endl;
  std::cout << "***********************************************************************" << std::endl;

  InertialFrame3D SolarSystemCMFrame(nullptr, Vec3Cart(0, 0, 0), Vec3Cart(0, 0, 0));

  Real earthOrbitPeriod = 365 * 24 * 3600;
  CircleOrbitingFrame3DCartesian EarthCMFrame(&SolarSystemCMFrame, 150000, earthOrbitPeriod);    // radijus kruženja (jedinica 10^6 m, 1000km), period kruženja (s)

  // sad kreiramo rotirajuci sustav Zemlje, s osi rotacije pod kutom 23 stupnja prema z-osi sustava CM
  // u trenutku t=0, os rotacije je nagnuta u y-z ravnini
  VectorN<Real, 3> axis({ 0, sin(23 * Constants::PI / 180), cos(23 * Constants::PI / 180) });
  HardSphereRotatingFrame EarthSurfaceFrame(&EarthCMFrame, 6.4, 24 * 3600, axis);

  double h = 0.0;
  double lat = Utils::ExplicitToAngleDeg(45, 48, 47.4); // 45°48'47.4"N
  double lon = Utils::ExplicitToAngleDeg(15, 58, 38.3); // 15°58'38.3"E

  double T = 0;
  double V = 10; // m/s 50, 100, 300, 500, 1000, 5000

  // koristeci PARENT referentni sustav, transformiramo to u kartezian kord. u EarthCM frameu
  VectorN<Real, 3> init_pos({ lat, lon, h });

  Vector3Cartesian pos_EarthCM_cart = EarthSurfaceFrame.GetLocalPosInParentFrameAtTime(init_pos, T);
  std::cout << "Position on Earth surface in EarthCM frame: " << pos_EarthCM_cart << std::endl;

	Vector3Cartesian pos_SolarFrame = EarthCMFrame.GetLocalPosInParentFrameAtTime(pos_EarthCM_cart, 0);
  std::cout << "Position on Earth surface in Solar System frame at T=0: " << pos_SolarFrame << std::endl;

  RotatingSphereLocalCartesian local(EarthSurfaceFrame, lat, lon);

  double one_lat_deg_in_km = 2 * 6400 * cos(Utils::DegToRad(45)) * Constants::PI / 360;
  double one_km_in_long_deg = 1 / one_lat_deg_in_km;
  // transformirati lokalni (x, y, z) u parent sferni (lat, lon, h)

  // najprije kalkulacije u lokalnom frameu, sve dok nije preko 1. kozmičke brzine

  // PRECIZNE KALKULACIJE IDU U EarthSurfaceFrame

  // postaviti lokalni (sferni) koord sustav, fiksiran prema Suncu, u kojem se zemlja vrti, 
  // ali se ovisno o t, tocno zna transf. tocke iz Zemlj. koord u taj sferni sustav

  // imamo pocetnu tocku, za dani T tocno znamu njene sferne koordinate

  // Solarni sustav
  // 20.000, 1e5, 1e6, 1e7, 1e8, 2.5e8 m/s, where it will be in 1 hour
}