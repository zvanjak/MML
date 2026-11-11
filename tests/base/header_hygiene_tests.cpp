#include <cmath>

#ifndef __STDCPP_WANT_MATH_SPEC_FUNCS__
#define MML_TEST_EXPECT_NO_MATH_OPT_IN
#endif

#include <mml/base/Vector/VectorN.h>

#if defined(MML_GEOMETRY_H) || defined(MML_GEOMETRY_POINTS_H) || \
	defined(MML_GEOMETRY_2D_SHAPES_H) || defined(MML_GEOMETRY_3D_SHAPES_H)
#error VectorN must not depend on geometry headers.
#endif

#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>

#if defined(MML_GEOMETRY_H) || defined(MML_GEOMETRY_2D_SHAPES_H) || defined(MML_GEOMETRY_3D_SHAPES_H)
#error Coordinate vectors need points, not the geometry aggregate or shapes.
#endif

#include <mml/base/StandardFunctions.h>

#if defined(MML_TEST_EXPECT_NO_MATH_OPT_IN) && defined(__STDCPP_WANT_MATH_SPEC_FUNCS__)
#error MML must not inject a special-math opt-in after cmath has already been included.
#endif

#include <mml/base/Matrix/Matrix.h>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Foundational headers work without geometry or a math opt-in", "[headers][VectorN]")
{
	const MML::VectorN<MML::Real, 3> fixed{1, -4, 2};
	const MML::Vector<MML::Real> dynamic{1, -4, 2};
	REQUIRE(fixed.NormLInf() == 4);
	REQUIRE(dynamic.NormLInf() == 4);

	const std::array<std::array<MML::Real, 2>, 2> values{{{{1, 2}}, {{3, 4}}}};
	const MML::Matrix<MML::Real> matrix(values);
	REQUIRE(matrix(0, 1) == 2);
	REQUIRE(matrix(1, 0) == 3);

	REQUIRE(MML::Functions::Legendre(0, MML::Real(0)) == 1);
	REQUIRE(MML::Functions::Hermite(1, MML::Real(2)) == 4);
}

TEST_CASE("Coordinate-vector headers supply point interoperability without shapes", "[headers][VectorN]")
{
	const MML::Point2Cartesian from2(1, 2);
	const MML::Point2Cartesian to2(4, 6);
	const MML::Vector2Cartesian displacement2(from2, to2);
	const auto translated2 = from2 + displacement2;
	REQUIRE(translated2.X() == to2.X());
	REQUIRE(translated2.Y() == to2.Y());

	const MML::Point3Cartesian from3(1, 2, 3);
	const MML::Point3Cartesian to3(4, 6, 8);
	const MML::Vector3Cartesian displacement3(from3, to3);
	const auto translated3 = from3 + displacement3;
	REQUIRE(translated3.X() == to3.X());
	REQUIRE(translated3.Y() == to3.Y());
	REQUIRE(translated3.Z() == to3.Z());
	REQUIRE(displacement3.ToPoint().X() == 3);
}
