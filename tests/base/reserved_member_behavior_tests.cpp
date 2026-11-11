#include <catch2/catch_all.hpp>
#include <mml/base/DiracDeltaFunction.h>
#include <mml/base/Geometry/Geometry3D/Geometry3DPlane.h>
#include <mml/base/Geometry/GeometryBase/Geometry3DShapes.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DSpherePyramid.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DTorusCylinder.h>
#include <mml/core/Surfaces.h>
#include <mml/core/LinAlgEqSolvers/LinAlgQR.h>
#include <mml/systems/DiscreteMaps.h>
#include <mml/systems/LinearSystem.h>

using namespace MML;

TEST_CASE("Plane equation coefficients retain mutable public accessors", "[reserved-members][geometry]")
{
	Plane3D plane(Pnt3Cart(0, 0, 2), Vec3Cart(0, 0, 1));
	REQUIRE(plane.A() == 0);
	REQUIRE(plane.B() == 0);
	REQUIRE(plane.C() == 1);
	REQUIRE(plane.D() == -2);
	plane.A() = 1;
	plane.B() = 2;
	plane.C() = 3;
	plane.D() = -6;
	const Plane3D& read = plane;
	REQUIRE(read.A() == 1);
	REQUIRE(read.B() == 2);
	REQUIRE(read.C() == 3);
	REQUIRE(read.D() == -6);
	const auto normal = plane.Normal();
	REQUIRE(normal[0] == 1);
	REQUIRE(normal[1] == 2);
	REQUIRE(normal[2] == 3);
}

TEST_CASE("Radius and height member renames preserve body and surface geometry", "[reserved-members][geometry][surfaces]")
{
	TorusGeom torus(3, 1);
	REQUIRE(torus.MajorRadius() == 3);
	REQUIRE(torus.MinorRadius() == 1);
	torus.MajorRadius() = 4;
	torus.MinorRadius() = 2;
	REQUIRE(torus.InnerRadius() == 2);
	REQUIRE(torus.OuterRadius() == 6);
	REQUIRE(IsNearlyEqual(torus.Volume(), Real(32) * Constants::PI * Constants::PI));
	Surfaces::CylinderSurface cylinder(3, 5);
	REQUIRE(cylinder.getMaxW() == 5);
	const auto point = cylinder(0, 2);
	REQUIRE(point[0] == 3);
	REQUIRE(point[1] == 0);
	REQUIRE(point[2] == 2);
	Surfaces::Torus surface(3, 1);
	const auto outer = surface(0, 0);
	REQUIRE(outer[0] == 4);
	REQUIRE(outer[1] == 0);
	REQUIRE(outer[2] == 0);
	Surfaces::Sphere sphere(2);
	REQUIRE(IsNearlyEqual(sphere(0, 0).NormL2(), Real(2)));
	Sphere3D sphereBody(2, 4, 6);
	REQUIRE(sphereBody.getRadius() == 2);
	REQUIRE(IsNearlyEqual(sphereBody.Volume(), Real(32) * Constants::PI / 3));
	Torus3D torusBody(3, 1, 4, 4);
	REQUIRE(torusBody.getMajorRadius() == 3);
	REQUIRE(IsNearlyEqual(torusBody.Volume(), Real(6) * Constants::PI * Constants::PI));
	Cylinder3D cylinderBody(2, 5, 6);
	REQUIRE(cylinderBody.getRadius() == 2);
	REQUIRE(cylinderBody.getHeight() == 5);
	REQUIRE(IsNearlyEqual(cylinderBody.Volume(), Real(20) * Constants::PI));
}

TEST_CASE("Protected Dirac sharpness parameter remains available to derived approximations", "[reserved-members][Dirac]")
{
	class CustomDelta : public DiracFunction {
	public:
		using DiracFunction::DiracFunction;
		Real operator()(Real) const override { return static_cast<Real>(_sharpness); }
	};
	const CustomDelta custom(7);
	REQUIRE(custom(0) == 7);
	const DiracStep step(4);
	REQUIRE(step(0) == 4);
	REQUIRE(step(1) == 0);
	REQUIRE(IsNearlyEqual(DiracExp(4)(0), Real(4) / std::sqrt(2 * Constants::PI)));
	REQUIRE(IsNearlyEqual(DiracSqr(4)(0), Real(4) / Constants::PI));
	REQUIRE(IsNearlyEqual(DiracSin(4)(0), Real(4) / Constants::PI));
	REQUIRE_THROWS_AS(CustomDelta(0), ArgumentError);
}

TEST_CASE("Renamed QR and linear-system storage preserves solve and multiple-RHS behavior", "[reserved-members][QR][LinearSystem]")
{
	Matrix<Real> matrix(2, 2, {3, 1, 1, 2});
	Vector<Real> rhs{9, 8};
	QRSolver<Real> qr(matrix);
	const auto solution = qr.Solve(rhs);
	REQUIRE(IsNearlyEqual(solution[0], Real(2)));
	REQUIRE(IsNearlyEqual(solution[1], Real(3)));
	Systems::LinearSystem<Real> system(matrix, rhs);
	const auto actual = system.Solve();
	REQUIRE(actual.IsEqualTo(solution, Precision::LooseTolerance));
	REQUIRE(system.getMatrix() == matrix);
	const Matrix<Real> multipleRHS(2, 2, {9, 18, 8, 16});
	const Systems::LinearSystem<Real> multiple(matrix, multipleRHS);
	const auto solutions = multiple.SolveMultiple();
	REQUIRE(IsNearlyEqual(solutions(0, 0), Real(2)));
	REQUIRE(IsNearlyEqual(solutions(1, 0), Real(3)));
	REQUIRE(IsNearlyEqual(solutions(0, 1), Real(4)));
	REQUIRE(IsNearlyEqual(solutions(1, 1), Real(6)));
}

TEST_CASE("Standard map retains K parameter naming and analytical Jacobian", "[reserved-members][DiscreteMap]")
{
	Systems::StandardMap map(2);
	REQUIRE(map.getParam(0) == 2);
	REQUIRE(map.getParamName(0) == "K");
	map.setParam(0, 3);
	REQUIRE(map.getParam(0) == 3);
	Vector<Real> next(2);
	map.map(Vector<Real>{0, 0}, next);
	REQUIRE(next[0] == 0);
	REQUIRE(next[1] == 0);
	Matrix<Real> jacobian;
	map.jacobian(Vector<Real>{0, 0}, jacobian);
	REQUIRE(jacobian(0, 0) == 1);
	REQUIRE(jacobian(0, 1) == 3);
	REQUIRE(jacobian(1, 0) == 1);
	REQUIRE(jacobian(1, 1) == 4);
}
