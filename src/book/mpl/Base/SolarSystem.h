#ifndef MPL_SOLAR_SYSTEM_H
#define MPL_SOLAR_SYSTEM_H

namespace MPL
{
	enum class SolarSystemBody
	{
		Sun = 0,
		Mercury = 1,
		Venus = 2,
		Earth = 3,
		Mars = 4,
		Jupiter = 5,
		Saturn = 6,
		Uranus = 7,
		Neptune = 8,
		Pluto = 9
	};

	class SolarSystemBodyInfo
	{
	public:
		const SolarSystemBody Id;
		const std::string Name;
		
		const double Mass_kg;
		const double Radius_km;
		const double OrbitalPeriod_days;
		const double RotationPeriod_hours;

		const double AxialTilt;
		const double OrbitalEccentricity;
		const double OrbitalInclination;

		SolarSystemBodyInfo(SolarSystemBody id, const std::string& name, double mass, double radius, double orbitalPeriod, double rotationPeriod,
			double axialTilt = 0.0, double orbitalEccentricity = 0.0, double orbitalInclination = 0.0)
			: Id(id), Name(name), Mass_kg(mass), Radius_km(radius), OrbitalPeriod_days(orbitalPeriod), RotationPeriod_hours(rotationPeriod),
			AxialTilt(axialTilt), OrbitalEccentricity(orbitalEccentricity), OrbitalInclination(orbitalInclination) 
		{	}
	};

	class SolarSystem
	{
	public:
		const static inline SolarSystemBodyInfo _bodies[] = {
			{	SolarSystemBody::Sun,			"Sun",			1.989e30, 695700.0,     0.0,     25.38,     7.25,      0.0,     0.0 },
			{ SolarSystemBody::Mercury, "Mercury",	3.285e23,		2439.7,    87.97,    58.646,    7.005,		 0.2056,  7.0049 },
			{ SolarSystemBody::Venus,		"Venus",		4.867e24,		6051.8,   224.701, -243,		  177.4,			 0.0068,  3.3946 },
			{ SolarSystemBody::Earth,		"Earth",		5.972e24,		6371.0,   365.256,   23.9345,  23.4392811, 0.0167, -7.155 },
			{ SolarSystemBody::Mars,		"Mars",			6.417e23,		3389.5,   687,      -24.6229,  25.19,			 0.0934, -1.850 },
			{ SolarSystemBody::Jupiter, "Jupiter",	1.898e27,  69911.0,  4332,       -9.925 ,  -3.13,			 0.0484, -1.303 },
			{ SolarSystemBody::Saturn,	"Saturn",		5.683e26,  58232.0, 10759,      -10.656 , -26.73,			 0.0542, -2.485 },
			{ SolarSystemBody::Uranus,	"Uranus",		8.681e25,  25362.0, 30689,      -17.24  ,  97.77,			 0.0472,  0.773 },
			{ SolarSystemBody::Neptune, "Neptune",	1.024e26,  24622.0, 60190,      -16.11  ,  28.32,			 0.0086,  1.769 },
			{ SolarSystemBody::Pluto,		"Pluto",		1.303e22,   1188.3, 90560,     -153.3   , 122.53,			 0.2488, 17.14 }
		};

		static const SolarSystemBodyInfo& GetBodyInfo(SolarSystemBody body)
		{
			return _bodies[static_cast<int>(body)];
		}
		static const SolarSystemBodyInfo& GetBodyInfo(std::string name)
		{
			for (const auto& body : _bodies)
			{
				if (body.Name == name)
					return body;
			}
			throw std::invalid_argument("Solar system body not found: " + name);
		}
	};
}

#endif