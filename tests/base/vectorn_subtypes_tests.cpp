#include "../catch/catch.hpp"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry.h"
#include "base/VectorTypes.h"
#endif

using namespace MML;

namespace MML::Tests::Base::VectorNTests
{
	/*********************************************************************/
	/*****                 Vector2Cartesian tests                    *****/
	/*********************************************************************/

	TEST_CASE("Test_Vector2Cartesian_GetUnitVector", "[simple]") {
		Vector2Cartesian a(4, 0);

		auto unit = a.GetUnitVector();

		REQUIRE(unit.X() == 1.0);
		REQUIRE(unit.Y() == 0.0);
	}

	TEST_CASE("Test_Point3Cartesian", "[simple]") {
		Point3Cartesian a(1, 1, 1);

		REQUIRE(a.X() == 1.0);
		REQUIRE(a.Y() == 1.0);
		REQUIRE(a.Z() == 1.0);
	}
	// Dist
	// oper + Vector
	// oper - Point

	/*********************************************************************/
	/*****                   Vector2Polar tests                      *****/
	/*********************************************************************/


	/*********************************************************************/
	/*****                 Vector3Cartesian tests                    *****/
	/*********************************************************************/
	TEST_CASE("Test_Vector3Cartesian_init", "[simple]") {
		Vector3Cartesian a(1, 1, 1);

		REQUIRE(a.X() == 1.0);
		REQUIRE(a.Y() == 1.0);
		REQUIRE(a.Z() == 1.0);

		Vector3Cartesian b{ 2,2,2 };

		REQUIRE(b.X() == 2.0);
		REQUIRE(b.Y() == 2.0);
		REQUIRE(b.Z() == 2.0);

		VectorN<Real, 3> c{ 3,3,3 };
		Vector3Cartesian d(c);

		REQUIRE(d.X() == 3.0);
		REQUIRE(d.Y() == 3.0);
		REQUIRE(d.Z() == 3.0);
	}

	TEST_CASE("Test_Vector3Cartesian_GetUnitVector", "[simple]") {
		Vector3Cartesian a(4, 0, 0);

		auto unit = a.GetAsUnitVector();

		REQUIRE(unit.X() == 1.0);
		REQUIRE(unit.Y() == 0.0);
		REQUIRE(unit.Z() == 0.0);
	}

	// IsParallelTo
	// IsPerpendicularTo
	// ScalarProd
	// VectorProd
	// VectorsAngle

	/*********************************************************************/
	/*****                 Vector3Spherical tests                    *****/
	/*********************************************************************/

	/*********************************************************************/
	/*****                Vector3Cylindrical tests                   *****/
	/*********************************************************************/


} // namespace MML::Tests::Base::GeometryTests