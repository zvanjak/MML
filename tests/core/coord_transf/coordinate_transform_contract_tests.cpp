#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/CoordTransf/CoordTransf2D.h>
#include <mml/core/CoordTransf/CoordTransf3D.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#include <mml/core/Fields/TensorFieldAdapters.h>
#include <mml/core/Fields/ScalarFieldOperations.h>
#include <mml/core/DifferentialGeometry/CoordinateMap.h>
#include <mml/core/DifferentialGeometry/Atlas.h>
#endif
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <optional>
#include <tuple>
#include <utility>

using namespace MML;
using Catch::Matchers::WithinAbs;
using namespace MML::DifferentialGeometry;

namespace {
	class Field4 : public ITensorField4<2> {
	public:
		Field4() : ITensorField4<2>(2, 2) {}
		Tensor4<2> operator()(const VectorN<Real, 2>& p) const override {
			Tensor4<2> result(2, 2);
			result.setContravar(0, false);
			result.setContravar(1, true);
			result.setContravar(2, false);
			result.setContravar(3, true);
			result(0, 1, 0, 1) = p[0] + 2 * p[1];
			return result;
		}
		Real component(int i, int j, int k, int l, const VectorN<Real, 2>& p) const override {
			return (*this)(p)(i, j, k, l);
		}
	};
	class Field5 : public ITensorField5<2> {
	public:
		Field5() : ITensorField5<2>(2, 3) {}
		Tensor5<2> operator()(const VectorN<Real, 2>& p) const override {
			Tensor5<2> result(3, 2);
			for (int slot = 0; slot < 5; ++slot) result.setContravar(slot, slot == 0 || slot == 3);
			result(0, 1, 0, 1, 0) = p[0] + 2 * p[1];
			return result;
		}
		Real component(int i, int j, int k, int l, int m, const VectorN<Real, 2>& p) const override {
			return (*this)(p)(i, j, k, l, m);
		}
	};

	template<class Transform, int N, class FirstArgs, class SecondArgs>
	void CheckValueCopies(const VectorN<Real, N>& point, const FirstArgs& first, const SecondArgs& second) {
		std::optional<Transform> original;
		std::apply([&](const auto&... args) { original.emplace(args...); }, first);
		const auto expected = original->transf(point);
		const auto expectedJac = original->jacobian(point);
		Transform copied(*original);
		Transform moved(std::move(*original));
		// Replace the live source before checking: no dangling callback is executed in the red baseline.
		original.reset();
		std::apply([&](const auto&... args) { original.emplace(args...); }, second);
		for (const auto* value : {&copied, &moved}) {
			const auto result = (*value)(point);
			const auto direct = value->transf(point);
			for (int axis = 0; axis < N; ++axis) {
				CHECK_THAT(result[axis], WithinAbs(expected[axis], Precision::LooseTolerance));
				CHECK_THAT(result[axis], WithinAbs(direct[axis], Precision::LooseTolerance));
				CHECK_THAT(value->inverseCoordTransfFunc(axis)(expected), WithinAbs(point[axis], Precision::LooseTolerance));
			}
		}
		original.reset();
		for (const auto* value : {&copied, &moved}) {
			const auto result = (*value)(point);
			const auto jac = value->jacobian(point);
			const auto back = value->transfInverse(expected);
			for (int i = 0; i < N; ++i) {
				CHECK_THAT(result[i], WithinAbs(expected[i], Precision::LooseTolerance));
				CHECK_THAT(back[i], WithinAbs(point[i], Precision::LooseTolerance));
				CHECK_THAT(value->inverseCoordTransfFunc(i)(expected), WithinAbs(point[i], Precision::LooseTolerance));
				for (int j = 0; j < N; ++j)
					CHECK_THAT(jac(i, j), WithinAbs(expectedJac(i, j), Precision::LooseTolerance));
			}
		}
	}
}

TEST_CASE("Coordinate contracts instantiate rank four and five transforms and field components",
	"[CoordTransf][coordinate-contract]") {
	CoordTransfCart2DRotation identity(0);
	const Vector2Cartesian p{2, 1};
	Field4 field4;
	Field5 field5;
	TransformedTensorField4<Vector2Cartesian, Vector2Cartesian, 2> adapted4(field4, identity);
	TransformedTensorField5<Vector2Cartesian, Vector2Cartesian, 2> adapted5(field5, identity);
	CHECK_THAT(identity.transfTensor4(field4(p), p)(0, 1, 0, 1), WithinAbs(Real(4), Precision::LooseTolerance));
	CHECK_THAT(identity.transfTensor5(field5(p), p)(0, 1, 0, 1, 0), WithinAbs(Real(4), Precision::LooseTolerance));
	CHECK_THAT(adapted4.component(0, 1, 0, 1, p), WithinAbs(Real(4), Precision::LooseTolerance));
	CHECK_THAT(adapted5.component(0, 1, 0, 1, 0, p), WithinAbs(Real(4), Precision::LooseTolerance));
}

TEST_CASE("Cartesian rotation Jacobians retain exact linear coefficients at every point",
	"[CoordTransf][coordinate-contract]") {
	for (Real angle : {Real(0), Constants::PI / Real(3), -Constants::PI / Real(4)}) {
		const CoordTransfCart2DRotation rotation(angle);
		const Real c = std::cos(angle), s = std::sin(angle);
		const MatrixNM<Real, 2, 2> expected{c, -s, s, c}, inverse{c, s, -s, c};
		for (const auto& point : {Vector2Cartesian{0, 0}, Vector2Cartesian{2, 1},
			Vector2Cartesian{std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::max()}}) {
			const auto forward = rotation.jacobian(point);
			const auto backward = rotation.inverseJacobian(point);
			for (int i = 0; i < 2; ++i)
				for (int j = 0; j < 2; ++j) {
					CHECK(forward(i, j) == expected(i, j));
					CHECK(backward(i, j) == inverse(i, j));
				}
		}
		CHECK_THROWS_AS(rotation.jacobian(Vector2Cartesian{std::numeric_limits<Real>::quiet_NaN(), 0}), DomainError);
		CHECK_THROWS_AS(rotation.inverseJacobian(Vector2Cartesian{0, std::numeric_limits<Real>::infinity()}), DomainError);
	}
}

TEST_CASE("Coordinate contracts expose identity one form pullback directly",
	"[CoordinateMap][coordinate-contract]") {
	auto identity = MakeIdentityCoordinateMap<Cartesian2, 2>();
	const Point<Real, 2, Cartesian2> p{1, 2};
	const auto dx = BasisOneForm<0, 2, Cartesian2>();
	const auto pulled = identity.pull_back(dx, p);
	CHECK(pulled.Component(0) == 1);
	CHECK(pulled.Component(1) == 0);
}

TEST_CASE("Coordinate contracts preserve scaled orthogonal axes and reject degenerate bases",
	"[CoordTransf][coordinate-contract]") {
	CoordTransf3DCartOrthogonal transform(Vec3Cart{2, 0, 0}, Vec3Cart{0, 3, 0}, Vec3Cart{0, 0, 4});
	const Vec3Cart point{1, 1, 1};
	const auto target = transform.transf(point);
	CHECK(target[0] == 2);
	CHECK(target[1] == 3);
	CHECK(target[2] == 4);
	const auto roundtrip = transform.transfInverse(target);
	for (int i = 0; i < 3; ++i) CHECK_THAT(roundtrip[i], WithinAbs(point[i], Precision::LooseTolerance));
	CHECK_THROWS_AS((CoordTransf3DCartOrthogonal(Vec3Cart{0, 0, 0}, Vec3Cart{0, 1, 0}, Vec3Cart{0, 0, 1})), GeometryError);
}

TEST_CASE("Coordinate contracts rebind copied and relocated transform components",
	"[CoordTransf][coordinate-contract]") {
	CheckValueCopies<CoordTransfCart2DRotation>(VectorN<Real, 2>{1, 2}, std::tuple{Real(0)}, std::tuple{Constants::PI / 2});
	const auto firstAngle = std::tuple{Constants::PI / 2};
	const auto secondAngle = std::tuple{Constants::PI / 3};
	CheckValueCopies<CoordTransfObliqueToCartesian2D>(VectorN<Real, 2>{1, 2}, firstAngle, secondAngle);
	CheckValueCopies<CoordTransfCartesianToOblique2D>(VectorN<Real, 2>{1, 2}, firstAngle, secondAngle);
	const VectorN<Real, 3> p{1, 2, 3};
	CheckValueCopies<CoordTransfCart3DRotationXAxis>(p, firstAngle, secondAngle);
	CheckValueCopies<CoordTransfCart3DRotationYAxis>(p, firstAngle, secondAngle);
	CheckValueCopies<CoordTransfCart3DRotationZAxis>(p, firstAngle, secondAngle);
	CheckValueCopies<CoordTransfCart3DRotationEuler>(p, std::tuple{Real(0), Real(0), Real(0)},
		std::tuple{Constants::PI / 2, Constants::PI / 3, Constants::PI / 4});
	const auto unitAxes = std::tuple{Vec3Cart{1, 0, 0}, Vec3Cart{0, 1, 0}, Vec3Cart{0, 0, 1}};
	CheckValueCopies<CoordTransf3DCartOrthogonal>(p, unitAxes,
		std::tuple{Vec3Cart{2, 0, 0}, Vec3Cart{0, 3, 0}, Vec3Cart{0, 0, 4}});
	CheckValueCopies<CoordTransf3DCartGeneral>(p, std::tuple{MatrixNM<Real, 3, 3>{1, 0, 0, 0, 1, 0, 0, 0, 1}},
		std::tuple{MatrixNM<Real, 3, 3>{2, 1, 0, 0, 3, 0, 0, 0, 4}});
	CheckValueCopies<CoordTransfCartesianToOblique3D>(p, unitAxes,
		std::tuple{Vec3Cart{1, 0, 0}, Vec3Cart{1, 2, 0}, Vec3Cart{0, 1, 3}});
	CheckValueCopies<CoordTransfCart3DRotationQuaternion>(p,
		std::tuple{Vec3Cart{0, 0, 1}, Real(0)}, std::tuple{Vec3Cart{0, 0, 1}, Constants::PI / 2});
}

TEST_CASE("Coordinate contracts polar inverse Jacobian is local across the atan2 branch cut",
	"[CoordTransf][coordinate-contract]") {
	CoordTransfPolarToCartesian2D transform;
	const auto jac = transform.inverseJacobian(VectorN<Real, 2>{-1, 0});
	CHECK_THAT(jac(0, 0), WithinAbs(Real(-1), Precision::LooseTolerance));
	CHECK_THAT(jac(0, 1), WithinAbs(Real(0), Precision::LooseTolerance));
	CHECK_THAT(jac(1, 0), WithinAbs(Real(0), Precision::LooseTolerance));
	CHECK_THAT(jac(1, 1), WithinAbs(Real(-1), Precision::LooseTolerance));
}

TEST_CASE("Coordinate contracts ordinary conversions and stereographic maps retain large finite inputs",
	"[CoordTransf][Atlas][coordinate-contract]") {
	const Real large = std::sqrt(std::numeric_limits<Real>::max()) * 16;
	const Vec3Cart p{large, 0, large};
	const auto spherical = CoordTransfCartesianToSpherical{}.transf(p);
	CHECK(std::isfinite(spherical[0]));
	CHECK_THAT(spherical[0] / large, WithinAbs(std::sqrt(Real(2)), Precision::LooseTolerance));
	CHECK_THAT(spherical[1], WithinAbs(Constants::PI / 4, Precision::LooseTolerance));
	const auto cylindrical = CoordTransfCartesianToCylindrical{}.transf(p);
	CHECK(std::isfinite(cylindrical[0]));
	CHECK_THAT(cylindrical[0] / large, WithinAbs(Real(1), Precision::LooseTolerance));
	CHECK(std::isfinite(CoordTransfPolarToCartesian2D{}.transfInverse(Vector2Cartesian{large, large})[0]));
	UnitSphereStereographicChart north(UnitSphereStereographicChart::Pole::North);
	const auto sphere = north.map_point(Point<Real, 2, Parameter2>{large, 0});
	for (int i = 0; i < 3; ++i) CHECK(std::isfinite(sphere[i]));
	CHECK_THAT(sphere[2], WithinAbs(Real(-1), Precision::LooseTolerance));
}

TEST_CASE("Coordinate contracts rank four and five preserve every variance slot independently",
	"[CoordTransf][TensorField][coordinate-contract]") {
	CoordTransfPolarToCartesian2D transform;
	const Vector2Polar point{2, Real(0.2)};
	const Real c = std::cos(point[1]), s = std::sin(point[1]), r = point[0];
	// Independent outer-product oracle: J for vectors, K^T for covectors.
	const Real j[2][2]{{c, -r * s}, {s, r * c}};
	const Real kt[2][2]{{c, -s / r}, {s, c / r}};
	for (int mask = 0; mask < 32; ++mask) {
		int count = 0;
		for (int slot = 0; slot < 5; ++slot) count += (mask >> slot) & 1;
		Tensor5<2> input(5 - count, count);
		for (int slot = 0; slot < 5; ++slot) input.setContravar(slot, ((mask >> slot) & 1) != 0);
		for (int flat = 0; flat < 32; ++flat) {
			const int a = flat / 16, b = flat / 8 % 2, d = flat / 4 % 2, e = flat / 2 % 2, f = flat % 2;
			input(a, b, d, e, f) = Real((1 + a) * (2 + b) * (3 + d) * (4 + e) * (5 + f));
		}
		const auto output = transform.transfTensor5(input, point);
		for (int slot = 0; slot < 5; ++slot) CHECK(output.isContravar(slot) == input.isContravar(slot));
		for (int flat = 0; flat < 32; ++flat) {
			const int axes[]{flat / 16, flat / 8 % 2, flat / 4 % 2, flat / 2 % 2, flat % 2};
			Real expected = 1, absoluteTerms = 1;
			for (int slot = 0; slot < 5; ++slot) {
				const auto& matrix = input.isContravar(slot) ? j : kt;
				expected *= matrix[axes[slot]][0] * (slot + 1) + matrix[axes[slot]][1] * (slot + 2);
				absoluteTerms *= std::abs(matrix[axes[slot]][0] * (slot + 1)) + std::abs(matrix[axes[slot]][1] * (slot + 2));
			}
			const Real roundoff = (32 + 2 * 5) * Precision::Eps * absoluteTerms;
			CHECK_THAT(output(axes[0], axes[1], axes[2], axes[3], axes[4]),
				WithinAbs(expected, roundoff + Precision::LooseTolerance * std::max(Real(1), std::abs(expected))));
		}
		if (mask >= 16) continue;
		count -= (mask >> 4) & 1;
		Tensor4<2> input4(4 - count, count);
		for (int slot = 0; slot < 4; ++slot) input4.setContravar(slot, ((mask >> slot) & 1) != 0);
		for (int flat = 0; flat < 16; ++flat) {
			const int a = flat / 8, b = flat / 4 % 2, d = flat / 2 % 2, e = flat % 2;
			input4(a, b, d, e) = Real((1 + a) * (2 + b) * (3 + d) * (4 + e));
		}
		const auto output4 = transform.transfTensor4(input4, point);
		for (int slot = 0; slot < 4; ++slot) CHECK(output4.isContravar(slot) == input4.isContravar(slot));
		for (int flat = 0; flat < 16; ++flat) {
			const int axes[]{flat / 8, flat / 4 % 2, flat / 2 % 2, flat % 2};
			Real expected = 1, absoluteTerms = 1;
			for (int slot = 0; slot < 4; ++slot) {
				const auto& matrix = input4.isContravar(slot) ? j : kt;
				expected *= matrix[axes[slot]][0] * (slot + 1) + matrix[axes[slot]][1] * (slot + 2);
				absoluteTerms *= std::abs(matrix[axes[slot]][0] * (slot + 1)) + std::abs(matrix[axes[slot]][1] * (slot + 2));
			}
			const Real roundoff = (16 + 2 * 4) * Precision::Eps * absoluteTerms;
			CHECK_THAT(output4(axes[0], axes[1], axes[2], axes[3]),
				WithinAbs(expected, roundoff + Precision::LooseTolerance * std::max(Real(1), std::abs(expected))));
		}
	}
	Field4 field4;
	Field5 field5;
	TransformedTensorField4<Vector2Polar, Vector2Cartesian, 2> adapted4(field4, transform);
	TransformedTensorField5<Vector2Polar, Vector2Cartesian, 2> adapted5(field5, transform);
	const auto target = transform.transf(point);
	const Real magnitude = point[0] + 2 * point[1];
	for (int flat = 0; flat < 32; ++flat) {
		const int a = flat / 16, b = flat / 8 % 2, d = flat / 4 % 2, e = flat / 2 % 2, f = flat % 2;
		const Real expected = magnitude * j[a][0] * kt[b][1] * kt[d][0] * j[e][1] * kt[f][0];
		CHECK_THAT(adapted5.component(a, b, d, e, f, target),
			WithinAbs(expected, Precision::LooseTolerance * std::max(Real(1), std::abs(expected))));
	}
	for (int flat = 0; flat < 16; ++flat) {
		const int a = flat / 8, b = flat / 4 % 2, d = flat / 2 % 2, e = flat % 2;
		const Real expected = magnitude * kt[a][0] * j[b][1] * kt[d][0] * j[e][1];
		CHECK_THAT(adapted4.component(a, b, d, e, target),
			WithinAbs(expected, Precision::LooseTolerance * std::max(Real(1), std::abs(expected))));
	}
}

TEST_CASE("Coordinate contracts preserve typed pairings and curve reparameterization",
	"[CoordinateMap][coordinate-contract]") {
	CoordTransfPolarToCartesian2D transform;
	const auto map = MakePolarToCartesian2DMap(transform);
	const Point<Real, 2, Polar2> p{2, Real(0.4)};
	const TangentVector<2, Polar2> velocity{1, 2}, rescaled{3, 6};
	const Covector<2, Cartesian2> covector{3, -4};
	const auto pushed = map.push_forward(velocity, p);
	const auto pushedRescaled = map.push_forward(rescaled, p);
	const auto pulled = map.pull_back(covector, p);
	CHECK_THAT(pushed[0], WithinAbs(std::cos(p[1]) - 4 * std::sin(p[1]), Precision::LooseTolerance));
	CHECK_THAT(pushed[1], WithinAbs(std::sin(p[1]) + 4 * std::cos(p[1]), Precision::LooseTolerance));
	const Real pairing = covector[0] * pushed[0] + covector[1] * pushed[1];
	CHECK_THAT(pulled[0] * velocity[0] + pulled[1] * velocity[1],
		WithinAbs(pairing, Precision::LooseTolerance * std::abs(pairing)));
	for (int i = 0; i < 2; ++i) CHECK_THAT(pushedRescaled[i], WithinAbs(3 * pushed[i], Precision::LooseTolerance));
}

TEST_CASE("Coordinate contracts finite errors are not converted into singularity sentinels",
	"[CoordTransf][Atlas][coordinate-contract]") {
	CoordTransfPolarToCartesian2D polar;
	CoordTransfSphericalToCartesian spherical;
	CoordTransfCylindricalToCartesian cylindrical;
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	for (const Real invalid : {nan, std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(polar.transfInverse(Vector2Cartesian{invalid, 0}), DomainError);
		CHECK_THROWS_AS(polar.inverseJacobian(VectorN<Real, 2>{invalid, 0}, SingularityPolicy::ReturnZero), DomainError);
		for (int axis = 0; axis < 3; ++axis) {
			VectorN<Real, 3> p{2, 1, 1};
			p[axis] = invalid;
			CHECK_THROWS_AS(spherical.transf(Vector3Spherical(p)), DomainError);
			CHECK_THROWS_AS(spherical.transfInverse(Vec3Cart(p)), DomainError);
			CHECK_THROWS_AS(cylindrical.transf(Vector3Cylindrical(p)), DomainError);
			CHECK_THROWS_AS(cylindrical.transfInverse(Vec3Cart(p)), DomainError);
			CHECK_THROWS_AS(cylindrical.coordTransfFunc(2)(p), DomainError);
			CHECK_THROWS_AS(spherical.jacobian(p), DomainError);
			CHECK_THROWS_AS(cylindrical.jacobian(p), DomainError);
		}
		const Point<Real, 2, Parameter2> p{invalid, 0};
		UnitSphereStereographicChart chart(UnitSphereStereographicChart::Pole::North);
		CHECK_THROWS_AS(chart.map_point(p), DomainError);
		CHECK_THROWS_AS(chart.jacobian(p), DomainError);
		CHECK_THROWS_AS(StereographicSphereTransition(p), DomainError);
	}
	const Real max = std::numeric_limits<Real>::max();
	CHECK_THROWS_AS(polar.transfInverse(Vector2Cartesian{max, max}), DomainError);
	CHECK_THROWS_AS(spherical.transfInverse(Vec3Cart{max, max, max}), DomainError);
	CHECK_THROWS_AS(cylindrical.transfInverse(Vec3Cart{max, max, 0}), DomainError);
	const VectorN<Real, 2> origin{0, 0};
	CHECK_THROWS_AS(polar.inverseJacobian(origin), SingularMatrixError);
	CHECK_THROWS_AS(polar.inverseJacobian(origin, SingularityPolicy::Clamp), DomainError);
	CHECK_THROWS_AS(polar.inverseJacobian(origin, static_cast<SingularityPolicy>(-1)), ArgumentError);
	const auto zeros = polar.inverseJacobian(origin, SingularityPolicy::ReturnZero);
	const auto nans = polar.inverseJacobian(origin, SingularityPolicy::ReturnNaN);
	const auto infs = polar.inverseJacobian(origin, SingularityPolicy::ReturnInf);
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j) {
			CHECK(zeros(i, j) == 0);
			CHECK(std::isnan(nans(i, j)));
			CHECK(std::isinf(infs(i, j)));
		}
}

TEST_CASE("Coordinate contracts small conversions preserve direction and scaled bases remain invertible",
	"[CoordTransf][coordinate-contract]") {
	const Real large = std::sqrt(std::numeric_limits<Real>::max()) * 16;
	const Real tiny = 1 / large;
	const Vec3Cart point{tiny, 0, tiny};
	const auto spherical = CoordTransfCartesianToSpherical{}.transf(point);
	CHECK_THAT(spherical[0] / tiny, WithinAbs(std::sqrt(Real(2)), Precision::LooseTolerance));
	CHECK_THAT(spherical[1], WithinAbs(Constants::PI / 4, Precision::LooseTolerance));
	CHECK(CoordTransfPolarToCartesian2D{}.transfInverse(Vector2Cartesian{tiny, tiny})[0] > 0);
	for (const Real scale : {tiny, large}) {
		CoordTransf3DCartOrthogonal transform(Vec3Cart{scale, 0, 0}, Vec3Cart{0, scale, 0}, Vec3Cart{0, 0, scale});
		const auto back = transform.transfInverse(transform.transf(Vec3Cart{1, 2, 3}));
		for (int i = 0; i < 3; ++i) CHECK_THAT(back[i], WithinAbs(Real(i + 1), Precision::LooseTolerance));
		CHECK(transform.IsRightHanded());
	}
	CHECK_THROWS_AS((CoordTransf3DCartOrthogonal(Vec3Cart{tiny, 0, 0}, Vec3Cart{tiny, tiny, 0}, Vec3Cart{0, 0, tiny})), GeometryError);
	CHECK_THROWS_AS((CoordTransf3DCartOrthogonal(Vec3Cart{std::numeric_limits<Real>::quiet_NaN(), 0, 0},
		Vec3Cart{0, 1, 0}, Vec3Cart{0, 0, 1})), DomainError);
	const Real max = std::numeric_limits<Real>::max();
	CoordTransf3DCartOrthogonal wide(Vec3Cart{max, max, 0}, Vec3Cart{max, -max, 0}, Vec3Cart{0, 0, max});
	const auto back = wide.transfInverse(wide.transf(Vec3Cart{Real(0.25), Real(0.25), Real(0.25)}));
	for (int i = 0; i < 3; ++i) CHECK_THAT(back[i], WithinAbs(Real(0.25), Precision::LooseTolerance));
	CHECK_FALSE(wide.IsRightHanded());
}

TEST_CASE("Coordinate contracts stereographic transitions preserve embeddings derivatives and metric",
	"[Atlas][coordinate-contract]") {
	auto atlas = MakeUnitSphereStereographicAtlas();
	const Point<Real, 2, Parameter2> p{Real(1.5), Real(0.5)};
	const auto q = atlas.transition(0, 1, p);
	const auto north = atlas.chart(0).map_point(p), south = atlas.chart(1).map_point(q);
	const auto jn = atlas.chart(0).jacobian(p), js = atlas.chart(1).jacobian(q);
	const Real r2 = p[0] * p[0] + p[1] * p[1];
	const MatrixNM<Real, 2, 2> dt{(p[1] * p[1] - p[0] * p[0]) / (r2 * r2),
		-2 * p[0] * p[1] / (r2 * r2), -2 * p[0] * p[1] / (r2 * r2),
		(p[0] * p[0] - p[1] * p[1]) / (r2 * r2)};
	const auto chained = js * dt;
	const auto metric = atlas.chart(0).induced_metric(p);
	for (int i = 0; i < 3; ++i) {
		CHECK_THAT(south[i], WithinAbs(north[i], Precision::LooseTolerance));
		for (int a = 0; a < 2; ++a) CHECK_THAT(chained(i, a), WithinAbs(jn(i, a), Precision::LooseTolerance));
	}
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			CHECK_THAT(metric(i, j), WithinAbs(i == j ? 4 / ((1 + r2) * (1 + r2)) : 0, Precision::LooseTolerance));
	for (const Real scale : {std::sqrt(std::numeric_limits<Real>::max()) * 16,
		1 / (std::sqrt(std::numeric_limits<Real>::max()) * 16)}) {
		const Point<Real, 2, Parameter2> extreme{scale, 0};
		const auto opposite = atlas.transition(0, 1, extreme);
		CHECK_THAT(opposite[0] * scale, WithinAbs(Real(1), Precision::LooseTolerance));
		const auto roundtrip = atlas.transition(1, 0, opposite);
		CHECK_THAT(roundtrip[0] / scale, WithinAbs(Real(1), Precision::LooseTolerance));
		const auto embedded = atlas.chart(0).map_point(extreme);
		CHECK_THAT(embedded[0] * (scale > 1 ? scale : 1 / scale), WithinAbs(Real(2), Precision::LooseTolerance));
		const auto jac = atlas.chart(0).jacobian(extreme);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 2; ++j) CHECK(std::isfinite(jac(i, j)));
		const Real expectedJac = scale > 1 ? (-2 / scale) / scale : Real(2);
		CHECK_THAT(jac(0, 0), WithinAbs(expectedJac, std::max(std::numeric_limits<Real>::denorm_min(),
			Precision::LooseTolerance * std::abs(expectedJac))));
	}
	CHECK_THROWS_AS(StereographicSphereTransition(Point<Real, 2, Parameter2>{0, 0}), ArgumentError);
	CHECK_THROWS_AS(StereographicSphereTransition(Point<Real, 2, Parameter2>{std::numeric_limits<Real>::denorm_min(), 0}), DomainError);
}

TEST_CASE("Coordinate contracts downstream detailed evaluations retain conversion failure evidence",
	"[CoordTransf][FieldOperations][coordinate-contract]") {
	CoordTransfPolarToCartesian2D transform;
	const Real scale = std::numeric_limits<Real>::max() * Real(0.75);
	ScalarFunctionFromStdFunc<2> radius([&](const VectorN<Real, 2>& p) {
		return transform.transfInverse(Vector2Cartesian{scale * p[0], scale * p[1]})[0];
	});
	FieldOperationConfig config;
	config.estimate_error = false;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const VectorN<Real, 2> p{1, 1};
	const auto failed = ScalarFieldOperations::GradientCartDetailed(radius, p, config);
	CHECK(failed.status == AlgorithmStatus::InvalidInput);
	CHECK_FALSE(failed.error_message.empty());
	CHECK(failed.function_evaluations == 1);
	for (int i = 0; i < 2; ++i) {
		CHECK(std::isnan(failed.value[i]));
		CHECK(std::isnan(failed.error[i]));
	}
	config.exception_policy = EvaluationExceptionPolicy::Propagate;
	CHECK_THROWS_AS(ScalarFieldOperations::GradientCartDetailed(radius, p, config), DomainError);
}
