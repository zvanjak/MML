#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/GeometrySpherical.h"
#include "base/VectorTypes.h"
#endif

using namespace MML;
using namespace Catch::Matchers;

namespace MML::Tests::Base::GeometrySphericalTests
{
	/*********************************************************************/
	/*****     SphericalGeometryCalculator::RadiusFromCartesian     *****/
	/*********************************************************************/
	TEST_CASE("RadiusFromCartesian::Origin", "[GeometrySpherical]")
	{
		Pnt3Cart origin(0.0, 0.0, 0.0);
		Real radius = SphericalGeometryCalculator::RadiusFromCartesian(origin);
		REQUIRE_THAT(radius, WithinAbs(0.0, REAL(1e-10)));
	}

	TEST_CASE("RadiusFromCartesian::OnAxes", "[GeometrySpherical]")
	{
		Pnt3Cart xAxis(3.0, 0.0, 0.0);
		Pnt3Cart yAxis(0.0, 4.0, 0.0);
		Pnt3Cart zAxis(0.0, 0.0, 5.0);

		REQUIRE_THAT(SphericalGeometryCalculator::RadiusFromCartesian(xAxis), WithinAbs(3.0, REAL(1e-10)));
		REQUIRE_THAT(SphericalGeometryCalculator::RadiusFromCartesian(yAxis), WithinAbs(4.0, REAL(1e-10)));
		REQUIRE_THAT(SphericalGeometryCalculator::RadiusFromCartesian(zAxis), WithinAbs(5.0, REAL(1e-10)));
	}

	TEST_CASE("RadiusFromCartesian::ArbitraryPoints", "[GeometrySpherical]")
	{
		// 3-4-5 right triangle
		Pnt3Cart pnt1(3.0, 4.0, 0.0);
		REQUIRE_THAT(SphericalGeometryCalculator::RadiusFromCartesian(pnt1), WithinAbs(5.0, REAL(1e-10)));

		// 3D Pythagorean triple
		Pnt3Cart pnt2(1.0, 2.0, 2.0);
		REQUIRE_THAT(SphericalGeometryCalculator::RadiusFromCartesian(pnt2), WithinAbs(3.0, REAL(1e-10)));

		// Unit sphere point
		Real x = 1.0 / std::sqrt(3.0);
		Pnt3Cart unitPoint(x, x, x);
		REQUIRE_THAT(SphericalGeometryCalculator::RadiusFromCartesian(unitPoint), WithinAbs(1.0, REAL(1e-10)));
	}

	/*********************************************************************/
	/*****   SphericalGeometryCalculator::CartesianFromSpherical    *****/
	/*********************************************************************/
	TEST_CASE("CartesianFromSpherical::NorthPole", "[GeometrySpherical]")
	{
		// Standard spherical: Pnt3Sph(r, theta, phi)
		// theta = polar angle from z-axis, phi = azimuthal angle
		// North pole: theta = 0
		Pnt3Sph northPole(5.0, 0.0, 0.0);
		Pnt3Cart cart = SphericalGeometryCalculator::CartesianFromSpherical(northPole);
		
		// x = R*sin(theta)*cos(phi) = 5*0*1 = 0
		// y = R*sin(theta)*sin(phi) = 5*0*0 = 0
		// z = R*cos(theta) = 5*1 = 5
		REQUIRE_THAT(cart.X(), WithinAbs(0.0, REAL(1e-10)));
		REQUIRE_THAT(cart.Y(), WithinAbs(0.0, REAL(1e-10)));
		REQUIRE_THAT(cart.Z(), WithinAbs(5.0, REAL(1e-10)));
	}

	TEST_CASE("CartesianFromSpherical::SouthPole", "[GeometrySpherical]")
	{
		// South pole: theta = PI
		Pnt3Sph southPole(5.0, Constants::PI, 0.0);
		Pnt3Cart cart = SphericalGeometryCalculator::CartesianFromSpherical(southPole);
		
		// x = 5*sin(PI)*cos(0) ≈ 0, y = 5*sin(PI)*sin(0) ≈ 0, z = 5*cos(PI) = -5
		REQUIRE_THAT(cart.X(), WithinAbs(0.0, REAL(1e-10)));
		REQUIRE_THAT(cart.Y(), WithinAbs(0.0, REAL(1e-10)));
		REQUIRE_THAT(cart.Z(), WithinAbs(-5.0, REAL(1e-10)));
	}

	TEST_CASE("CartesianFromSpherical::Equator", "[GeometrySpherical]")
	{
		// Equator: theta = PI/2
		// Point on x-axis: phi = 0
		Pnt3Sph equatorX(3.0, Constants::PI / 2, 0.0);
		Pnt3Cart cartX = SphericalGeometryCalculator::CartesianFromSpherical(equatorX);
		// x = 3*sin(PI/2)*cos(0) = 3*1*1 = 3, y = 3*sin(PI/2)*sin(0) = 0, z = 3*cos(PI/2) = 0
		REQUIRE_THAT(cartX.X(), WithinAbs(3.0, REAL(1e-10)));
		REQUIRE_THAT(cartX.Y(), WithinAbs(0.0, REAL(1e-10)));
		REQUIRE_THAT(cartX.Z(), WithinAbs(0.0, REAL(1e-10)));

		// Point on y-axis: phi = PI/2
		Pnt3Sph equatorY(4.0, Constants::PI / 2, Constants::PI / 2);
		Pnt3Cart cartY = SphericalGeometryCalculator::CartesianFromSpherical(equatorY);
		// x = 4*sin(PI/2)*cos(PI/2) = 0, y = 4*sin(PI/2)*sin(PI/2) = 4, z = 4*cos(PI/2) = 0
		REQUIRE_THAT(cartY.X(), WithinAbs(0.0, REAL(1e-10)));
		REQUIRE_THAT(cartY.Y(), WithinAbs(4.0, REAL(1e-10)));
		REQUIRE_THAT(cartY.Z(), WithinAbs(0.0, REAL(1e-10)));
	}

	TEST_CASE("CartesianFromSpherical::RoundTripConversion", "[GeometrySpherical]")
	{
		// Test that converting spherical -> cartesian maintains radius
		Pnt3Sph sph(7.0, Constants::PI / 3, Constants::PI / 4);
		Pnt3Cart cart = SphericalGeometryCalculator::CartesianFromSpherical(sph);
		Real calculatedRadius = SphericalGeometryCalculator::RadiusFromCartesian(cart);
		
		REQUIRE_THAT(calculatedRadius, WithinAbs(7.0, REAL(1e-10)));
	}

	/*********************************************************************/
	/*****   SphericalGeometryCalculator::SphericalFromLatLong      *****/
	/*********************************************************************/
	TEST_CASE("SphericalFromLatLong::Equator", "[GeometrySpherical]")
	{
		// Latitude 0° (equator) -> theta = PI/2 (polar angle from z-axis)
		Pnt3Sph pnt = SphericalGeometryCalculator::SphericalFromLatLong(0.0, 0.0);
		
		REQUIRE_THAT(pnt.R(), WithinAbs(1.0, REAL(1e-10)));
		REQUIRE_THAT(pnt.Theta(), WithinAbs(Constants::PI / 2, REAL(1e-10)));  // equator
		REQUIRE_THAT(pnt.Phi(), WithinAbs(0.0, REAL(1e-10)));    // longitude
	}

	TEST_CASE("SphericalFromLatLong::PositiveLatitudes", "[GeometrySpherical]")
	{
		// 45° North latitude -> theta = PI/2 - PI/4 = PI/4
		Pnt3Sph pnt45N = SphericalGeometryCalculator::SphericalFromLatLong(45.0, 0.0);
		REQUIRE_THAT(pnt45N.Theta(), WithinAbs(Constants::PI / 4, REAL(1e-10)));

		// 90° North latitude (North Pole) -> theta = 0
		Pnt3Sph pntNorth = SphericalGeometryCalculator::SphericalFromLatLong(90.0, 0.0);
		REQUIRE_THAT(pntNorth.Theta(), WithinAbs(0.0, REAL(1e-10)));
	}

	TEST_CASE("SphericalFromLatLong::NegativeLatitudes", "[GeometrySpherical]")
	{
		// -45° latitude (45° South) -> theta = PI/2 + PI/4 = 3*PI/4
		Pnt3Sph pnt45S = SphericalGeometryCalculator::SphericalFromLatLong(-45.0, 0.0);
		REQUIRE_THAT(pnt45S.Theta(), WithinAbs(3 * Constants::PI / 4, REAL(1e-10)));

		// -90° latitude (South Pole) -> theta = PI
		Pnt3Sph pntSouth = SphericalGeometryCalculator::SphericalFromLatLong(-90.0, 0.0);
		REQUIRE_THAT(pntSouth.Theta(), WithinAbs(Constants::PI, REAL(1e-10)));
	}

	TEST_CASE("SphericalFromLatLong::Longitudes", "[GeometrySpherical]")
	{
		// Longitude is stored in Phi()
		// Various longitudes at equator
		Pnt3Sph pnt0 = SphericalGeometryCalculator::SphericalFromLatLong(0.0, 0.0);
		REQUIRE_THAT(pnt0.Phi(), WithinAbs(0.0, REAL(1e-10)));

		Pnt3Sph pnt90E = SphericalGeometryCalculator::SphericalFromLatLong(0.0, 90.0);
		REQUIRE_THAT(pnt90E.Phi(), WithinAbs(Constants::PI / 2, REAL(1e-10)));

		Pnt3Sph pnt180 = SphericalGeometryCalculator::SphericalFromLatLong(0.0, 180.0);
		REQUIRE_THAT(pnt180.Phi(), WithinAbs(Constants::PI, REAL(1e-10)));

		Pnt3Sph pnt270 = SphericalGeometryCalculator::SphericalFromLatLong(0.0, 270.0);
		REQUIRE_THAT(pnt270.Phi(), WithinAbs(3 * Constants::PI / 2, REAL(1e-10)));
	}

	/*********************************************************************/
	/*****  SphericalGeometryCalculator::DistanceBetweenLatLong     *****/
	/*********************************************************************/
	TEST_CASE("DistanceBetweenLatLong::SamePoint", "[GeometrySpherical]")
	{
		Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(45.0, 30.0, 45.0, 30.0);
		REQUIRE_THAT(distance, WithinAbs(0.0, REAL(1e-10)));
	}

	TEST_CASE("DistanceBetweenLatLong::EquatorPoints", "[GeometrySpherical]")
	{
		// Distance along equator: 0° to 90° longitude
		Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(0.0, 0.0, 0.0, 90.0);
		REQUIRE_THAT(distance, WithinAbs(Constants::PI / 2, REAL(1e-9)));
	}

	TEST_CASE("DistanceBetweenLatLong::PoleToPole", "[GeometrySpherical]")
	{
		// North pole to South pole
		Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(90.0, 0.0, -90.0, 0.0);
		REQUIRE_THAT(distance, WithinAbs(Constants::PI, REAL(1e-9)));
	}

	TEST_CASE("DistanceBetweenLatLong::QuarterSphere", "[GeometrySpherical]")
	{
		// From equator to North pole
		Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(0.0, 0.0, 90.0, 0.0);
		REQUIRE_THAT(distance, WithinAbs(Constants::PI / 2, REAL(1e-9)));
	}

	TEST_CASE("DistanceBetweenLatLong::RealWorldExample", "[GeometrySpherical]")
	{
		// Approximate coordinates: London (51.5°N, 0°E) to New York (40.7°N, 74°W)
		Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(51.5, 0.0, 40.7, -74.0);
		
		// Great circle distance London to NYC is about 51.4° or 0.897 radians
		// Allow reasonable range due to approximate coordinates
		REQUIRE(distance > 0.8);
		REQUIRE(distance < 1.4);
	}

	TEST_CASE("DistanceBetweenLatLong::Symmetry", "[GeometrySpherical]")
	{
		// Distance should be symmetric
		Real dist1 = SphericalGeometryCalculator::DistanceBetweenLatLong(30.0, 45.0, 60.0, 90.0);
		Real dist2 = SphericalGeometryCalculator::DistanceBetweenLatLong(60.0, 90.0, 30.0, 45.0);
		
		REQUIRE_THAT(dist1, WithinAbs(dist2, REAL(1e-10)));
	}

	/*********************************************************************/
	/*****              Integration Tests                           *****/
	/*********************************************************************/
	TEST_CASE("Integration::SphericalCartesianConsistency", "[GeometrySpherical][Integration]")
	{
		// Create point in spherical, convert to cartesian, verify radius
		Pnt3Sph sph(10.0, Constants::PI / 3, Constants::PI / 6);
		Pnt3Cart cart = SphericalGeometryCalculator::CartesianFromSpherical(sph);
		
		Real radiusFromCart = SphericalGeometryCalculator::RadiusFromCartesian(cart);
		
		REQUIRE_THAT(radiusFromCart, WithinAbs(10.0, REAL(1e-10)));
	}

	TEST_CASE("Integration::LatLongToCartesianConsistency", "[GeometrySpherical][Integration]")
	{
		// Convert lat/long to spherical, then to cartesian, verify unit radius
		Real lat = 37.7749;  // San Francisco latitude
		Real lon = -122.4194; // San Francisco longitude
		
		Pnt3Sph sph = SphericalGeometryCalculator::SphericalFromLatLong(lat, lon);
		Pnt3Cart cart = SphericalGeometryCalculator::CartesianFromSpherical(sph);
		Real radius = SphericalGeometryCalculator::RadiusFromCartesian(cart);
		
		REQUIRE_THAT(radius, WithinAbs(1.0, REAL(1e-10)));
	}

	TEST_CASE("Integration::DistanceTriangleInequality", "[GeometrySpherical][Integration]")
	{
		// Triangle inequality: d(A,C) <= d(A,B) + d(B,C)
		Real distAB = SphericalGeometryCalculator::DistanceBetweenLatLong(0.0, 0.0, 30.0, 30.0);
		Real distBC = SphericalGeometryCalculator::DistanceBetweenLatLong(30.0, 30.0, 60.0, 60.0);
		Real distAC = SphericalGeometryCalculator::DistanceBetweenLatLong(0.0, 0.0, 60.0, 60.0);
		
		REQUIRE(distAC <= distAB + distBC + 1e-9);
	}

	/*********************************************************************/
	/*****   DistanceBetweenLatLong with Radius (Earth distances)   *****/
	/*********************************************************************/
	TEST_CASE("DistanceBetweenLatLong::WithRadius_BasicTest", "[GeometrySpherical]")
	{
		// Test with unit sphere first
		Real distRadians = SphericalGeometryCalculator::DistanceBetweenLatLong(0.0, 0.0, 0.0, 90.0);
		Real distMeters = SphericalGeometryCalculator::DistanceBetweenLatLong(1.0, 0.0, 0.0, 0.0, 90.0);
		
		// For radius 1.0, distance in meters should equal distance in radians
		REQUIRE_THAT(distMeters, WithinAbs(distRadians, REAL(1e-10)));
		
		// Quarter of great circle on equator (90 degrees) = PI/2 radians
		REQUIRE_THAT(distRadians, WithinAbs(Constants::PI / 2, REAL(1e-10)));
	}

	TEST_CASE("DistanceBetweenLatLong::WithRadius_ScalingTest", "[GeometrySpherical]")
	{
		// Test that distance scales linearly with radius
		Real lat1 = 40.7128, lon1 = -74.0060;  // New York
		Real lat2 = 51.5074, lon2 = -0.1278;   // London
		
		Real distRadians = SphericalGeometryCalculator::DistanceBetweenLatLong(lat1, lon1, lat2, lon2);
		
		Real radius1 = 1000.0;
		Real radius2 = 5000.0;
		
		Real dist1 = SphericalGeometryCalculator::DistanceBetweenLatLong(radius1, lat1, lon1, lat2, lon2);
		Real dist2 = SphericalGeometryCalculator::DistanceBetweenLatLong(radius2, lat1, lon1, lat2, lon2);
		
		// Distance should scale linearly with radius
		REQUIRE_THAT(dist2 / dist1, WithinAbs(5.0, REAL(1e-10)));
		
		// Both should equal radius * angle
		REQUIRE_THAT(dist1, WithinAbs(radius1 * distRadians, REAL(1e-10)));
		REQUIRE_THAT(dist2, WithinAbs(radius2 * distRadians, REAL(1e-10)));
	}

	TEST_CASE("DistanceBetweenLatLong::RealWorldCities_MajorRoutes", "[GeometrySpherical][RealWorld]")
	{
		// Earth's mean radius in kilometers
		const Real EARTH_RADIUS_KM = 6371.0;
		
		// City coordinates (latitude, longitude) and known great circle distances
		struct City {
			std::string name;
			Real lat, lon;
		};
		
		City newYork    = {"New York",    40.7128,  -74.0060};
		City london     = {"London",      51.5074,   -0.1278};
		City moscow     = {"Moscow",      55.7558,   37.6173};
		City tokyo      = {"Tokyo",       35.6762,  139.6503};
		City sydney     = {"Sydney",     -33.8688,  151.2093};
		City saoPaulo   = {"Sao Paulo",  -23.5505,  -46.6333};
		
		SECTION("New York to London")
		{
			Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(
				EARTH_RADIUS_KM, newYork.lat, newYork.lon, london.lat, london.lon);
			
			// Known distance: approximately 5,585 km
			// Allow 1% tolerance for coordinate precision and Earth's non-spherical shape
			REQUIRE_THAT(distance, WithinAbs(5585.0, REAL(60.0)));
		}
		
		SECTION("London to Moscow")
		{
			Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(
				EARTH_RADIUS_KM, london.lat, london.lon, moscow.lat, moscow.lon);
			
			// Known distance: approximately 2,500 km
			REQUIRE_THAT(distance, WithinAbs(2500.0, REAL(50.0)));
		}
		
		SECTION("Tokyo to Sydney")
		{
			Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(
				EARTH_RADIUS_KM, tokyo.lat, tokyo.lon, sydney.lat, sydney.lon);
			
			// Known distance: approximately 7,820 km
			REQUIRE_THAT(distance, WithinAbs(7820.0, REAL(80.0)));
		}
		
		SECTION("New York to Sao Paulo")
		{
			Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(
				EARTH_RADIUS_KM, newYork.lat, newYork.lon, saoPaulo.lat, saoPaulo.lon);
			
			// Known distance: approximately 7,700 km
			REQUIRE_THAT(distance, WithinAbs(7700.0, REAL(80.0)));
		}
		
		SECTION("Moscow to Tokyo")
		{
			Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(
				EARTH_RADIUS_KM, moscow.lat, moscow.lon, tokyo.lat, tokyo.lon);
			
			// Known distance: approximately 7,500 km
			REQUIRE_THAT(distance, WithinAbs(7500.0, REAL(80.0)));
		}
		
		SECTION("Sydney to Sao Paulo")
		{
			Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(
				EARTH_RADIUS_KM, sydney.lat, sydney.lon, saoPaulo.lat, saoPaulo.lon);
			
			// Known distance: approximately 13,350 km (nearly half-way around Earth)
			REQUIRE_THAT(distance, WithinAbs(13350.0, REAL(200.0)));
		}
	}

	TEST_CASE("DistanceBetweenLatLong::RealWorldCities_ComprehensiveMatrix", "[GeometrySpherical][RealWorld]")
	{
		// Earth's mean radius in kilometers
		const Real EARTH_RADIUS_KM = 6371.0;
		
		// Test a comprehensive matrix of all city pairs
		struct CityData {
			std::string name;
			Real lat, lon;
		};
		
		std::vector<CityData> cities = {
			{"New York",   40.7128,  -74.0060},
			{"London",     51.5074,   -0.1278},
			{"Moscow",     55.7558,   37.6173},
			{"Tokyo",      35.6762,  139.6503},
			{"Sydney",    -33.8688,  151.2093},
			{"Sao Paulo", -23.5505,  -46.6333}
		};
		
		// Known distances in km (approximate, due to Earth's oblate shape)
		// Matrix is symmetric, storing upper triangle
		std::map<std::pair<int, int>, Real> knownDistances = {
			{{0, 1}, 5585.0},   // NYC - London
			{{0, 2}, 7510.0},   // NYC - Moscow
			{{0, 3}, 10850.0},  // NYC - Tokyo
			{{0, 4}, 16000.0},  // NYC - Sydney
			{{0, 5}, 7700.0},   // NYC - Sao Paulo
			{{1, 2}, 2500.0},   // London - Moscow
			{{1, 3}, 9560.0},   // London - Tokyo
			{{1, 4}, 17000.0},  // London - Sydney
			{{1, 5}, 9480.0},   // London - Sao Paulo
			{{2, 3}, 7500.0},   // Moscow - Tokyo
			{{2, 4}, 14450.0},  // Moscow - Sydney
			{{2, 5}, 11750.0},  // Moscow - Sao Paulo
			{{3, 4}, 7820.0},   // Tokyo - Sydney
			{{3, 5}, 18550.0},  // Tokyo - Sao Paulo
			{{4, 5}, 13350.0}   // Sydney - Sao Paulo
		};
		
		// Test all city pairs
		int testCount = 0;
		for (int i = 0; i < cities.size(); i++) {
			for (int j = i + 1; j < cities.size(); j++) {
				Real calculated = SphericalGeometryCalculator::DistanceBetweenLatLong(
					EARTH_RADIUS_KM, 
					cities[i].lat, cities[i].lon, 
					cities[j].lat, cities[j].lon);
				
				Real expected = knownDistances[{i, j}];
				Real tolerance = expected * 0.015;  // 1.5% tolerance
				
				INFO("Distance from " << cities[i].name << " to " << cities[j].name);
				REQUIRE_THAT(calculated, WithinAbs(expected, tolerance));
				testCount++;
			}
		}
		
		REQUIRE(testCount == 15);  // Should test all 15 pairs (6 choose 2)
	}

	TEST_CASE("DistanceBetweenLatLong::WithRadius_EdgeCases", "[GeometrySpherical]")
	{
		const Real EARTH_RADIUS_KM = 6371.0;
		
		SECTION("Same point - zero distance")
		{
			Real distRadians = SphericalGeometryCalculator::DistanceBetweenLatLong(40.0, 50.0, 40.0, 50.0);
			Real distMeters = SphericalGeometryCalculator::DistanceBetweenLatLong(EARTH_RADIUS_KM, 40.0, 50.0, 40.0, 50.0);
			
			REQUIRE_THAT(distRadians, WithinAbs(0.0, REAL(1e-10)));
			REQUIRE_THAT(distMeters, WithinAbs(0.0, REAL(1e-10)));
		}
		
		SECTION("Antipodal points - maximum distance")
		{
			// North and South poles are antipodal
			Real distRadians = SphericalGeometryCalculator::DistanceBetweenLatLong(90.0, 0.0, -90.0, 0.0);
			Real distMeters = SphericalGeometryCalculator::DistanceBetweenLatLong(EARTH_RADIUS_KM, 90.0, 0.0, -90.0, 0.0);
			
			// Should be PI radians (half circumference)
			REQUIRE_THAT(distRadians, WithinAbs(Constants::PI, REAL(1e-9)));
			REQUIRE_THAT(distMeters, WithinAbs(EARTH_RADIUS_KM * Constants::PI, REAL(1e-6)));
		}
		
		SECTION("Equator opposite sides")
		{
			// Points on equator, 180 degrees apart
			Real distRadians = SphericalGeometryCalculator::DistanceBetweenLatLong(0.0, 0.0, 0.0, 180.0);
			Real distMeters = SphericalGeometryCalculator::DistanceBetweenLatLong(EARTH_RADIUS_KM, 0.0, 0.0, 0.0, 180.0);
			
			REQUIRE_THAT(distRadians, WithinAbs(Constants::PI, REAL(1e-9)));
			REQUIRE_THAT(distMeters, WithinAbs(EARTH_RADIUS_KM * Constants::PI, REAL(1e-6)));
		}
		
		SECTION("Across International Date Line")
		{
			// Tokyo (Japan) to Anchorage (Alaska) crosses date line
			Real tokyo_lat = 35.6762, tokyo_lon = 139.6503;
			Real anchorage_lat = 61.2181, anchorage_lon = -149.9003;
			
			Real distance = SphericalGeometryCalculator::DistanceBetweenLatLong(
				EARTH_RADIUS_KM, tokyo_lat, tokyo_lon, anchorage_lat, anchorage_lon);
			
			// Known distance approximately 5,500 km
			REQUIRE_THAT(distance, WithinAbs(5500.0, REAL(100.0)));
		}
	}
}
