#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/base/Geometry/Geometry.h"
#include "base/Vector/VectorTypes.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::VectorNTests
{
	/*********************************************************************/
	/*****                 Vector2Cartesian tests                    *****/
	/*********************************************************************/

	TEST_CASE("Vector2Cartesian::GetUnitVector", "[simple]") {
		TEST_PRECISION_INFO();
		Vector2Cartesian a(4, 0);

		auto unit = a.GetAsUnitVector();

		REQUIRE(unit.X() == REAL(1.0));
		REQUIRE(unit.Y() == REAL(0.0));
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
	TEST_CASE("Vector3Cartesian::Init", "[simple]") {
		TEST_PRECISION_INFO();
		Vector3Cartesian a(1, 1, 1);

		REQUIRE(a.X() == REAL(1.0));
		REQUIRE(a.Y() == REAL(1.0));
		REQUIRE(a.Z() == REAL(1.0));

		Vector3Cartesian b{ 2,2,2 };

		REQUIRE(b.X() == REAL(2.0));
		REQUIRE(b.Y() == REAL(2.0));
		REQUIRE(b.Z() == REAL(2.0));

		VectorN<Real, 3> c{ 3,3,3 };
		Vector3Cartesian d(c);

		REQUIRE(d.X() == REAL(3.0));
		REQUIRE(d.Y() == REAL(3.0));
		REQUIRE(d.Z() == REAL(3.0));
	}

	TEST_CASE("Vector3Cartesian::GetUnitVector", "[simple]") {
		TEST_PRECISION_INFO();
		Vector3Cartesian a(4, 0, 0);

		auto unit = a.GetAsUnitVector();

		REQUIRE(unit.X() == REAL(1.0));
		REQUIRE(unit.Y() == REAL(0.0));
		REQUIRE(unit.Z() == REAL(0.0));
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
