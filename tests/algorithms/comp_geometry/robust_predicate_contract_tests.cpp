#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/CompGeometry/RobustPredicates.h>
#endif
#include <catch2/catch_all.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>

namespace {
using namespace MML;
using MML::CompGeometry::RobustPredicates;

template<std::size_t N>
std::int64_t IntegerDeterminant(const std::array<std::array<std::int64_t, N>, N>& matrix) {
	if constexpr (N == 1) return matrix[0][0];
	else {
		std::int64_t sum = 0;
		for (std::size_t column = 0; column < N; ++column) {
			std::array<std::array<std::int64_t, N - 1>, N - 1> minor{};
			for (std::size_t row = 1; row < N; ++row) {
				std::size_t next = 0;
				for (std::size_t c = 0; c < N; ++c)
					if (c != column) minor[row - 1][next++] = matrix[row][c];
			}
			sum += (column % 2 == 0 ? 1 : -1) * matrix[0][column] * IntegerDeterminant(minor);
		}
		return sum;
	}
}

int IntegerSign(std::int64_t value) { return (value > 0) - (value < 0); }

int PermutationSign(const std::array<int, 4>& permutation) {
	int sign = 1;
	for (int i = 0; i < 4; ++i)
		for (int j = i + 1; j < 4; ++j)
			if (permutation[i] > permutation[j]) sign = -sign;
	return sign;
}

template<typename T>
std::pair<T, T> OriginalMixedWitness() {
	if constexpr (std::numeric_limits<T>::max_exponent >= 1024)
		return {T(1e200), T(1e-300)};
	else return {std::ldexp(T(1), std::numeric_limits<T>::max_exponent / 2),
		std::ldexp(T(1), std::numeric_limits<T>::min_exponent + std::numeric_limits<T>::digits)};
}

constexpr std::array<int, 5> Exponents{
	0, std::numeric_limits<Real>::max_exponent / 2, std::numeric_limits<Real>::min_exponent / 2,
	std::numeric_limits<Real>::max_exponent - 4, std::numeric_limits<Real>::min_exponent - std::numeric_limits<Real>::digits + 2};
}

TEST_CASE("All predicate signs survive uniform scale beyond determinant range", "[predicate-contract][Predicates]") {
	for (int exponent : Exponents) {
		CAPTURE(exponent);
		const Real s = std::ldexp(Real(1), exponent);
		CHECK(RobustPredicates::Orientation2D({0, 0}, {s, 0}, {0, s}) == 1);
		CHECK(RobustPredicates::Orientation2D({0, 0}, {0, s}, {s, 0}) == -1);
		CHECK(RobustPredicates::InCircle2D({s, 0}, {0, s}, {-s, 0}, {0, 0}) == 1);
		CHECK(RobustPredicates::Orientation3D({0, 0, s}, {s, 0, 0}, {0, s, 0}, {0, 0, 0}) == 1);
		CHECK(RobustPredicates::Orientation3D({s, 0, 0}, {0, 0, s}, {0, s, 0}, {0, 0, 0}) == -1);
		CHECK(RobustPredicates::InSphere3D({0, 0, 0}, {0, s, 0}, {s, 0, 0}, {0, 0, s}, {s / 2, s / 2, s / 2}) == 1);
	}
}

TEST_CASE("Mixed-scale orientation retains small coordinates without uniform rescaling", "[predicate-contract][Predicates]") {
	const auto original = OriginalMixedWitness<Real>();
	for (const auto& [large, tiny] : {original,
		std::make_pair(std::numeric_limits<Real>::max() / 4, std::numeric_limits<Real>::denorm_min())}) {
		const std::array<Point3Cartesian, 4> points{{{0, 0, tiny}, {large, 0, 0}, {0, large, 0}, {0, 0, 0}}};
		std::array<int, 4> order{0, 1, 2, 3};
		do {
			CHECK(RobustPredicates::Orientation3D(points[order[0]], points[order[1]], points[order[2]], points[order[3]]) ==
				PermutationSign(order));
		} while (std::next_permutation(order.begin(), order.end()));
		CHECK(RobustPredicates::Orientation3D({0, 0, -tiny}, points[1], points[2], points[3]) == -1);
		CHECK(RobustPredicates::Orientation3D({0, 0, 0}, points[1], points[2], points[3]) == 0);
	}
}

TEST_CASE("In-sphere distinguishes cospherical and adjacent values at every scale", "[predicate-contract][Predicates]") {
	for (int exponent : Exponents) {
		CAPTURE(exponent);
		const Real s = std::ldexp(Real(1), exponent);
		const std::array<Point3Cartesian, 4> tetra{{{0, 0, 0}, {0, s, 0}, {s, 0, 0}, {0, 0, s}}};
		const Point3Cartesian on(s, s, s), inside(s, s, std::nextafter(s, Real(0)));
		const Point3Cartesian outside(s, s, std::nextafter(s, std::numeric_limits<Real>::infinity()));
		std::array<int, 4> order{0, 1, 2, 3};
		do {
			const auto& a = tetra[order[0]];
			const auto& b = tetra[order[1]];
			const auto& c = tetra[order[2]];
			const auto& d = tetra[order[3]];
			CHECK(RobustPredicates::InSphere3D(a, b, c, d, on) == 0);
			CHECK(RobustPredicates::InSphere3D(a, b, c, d, inside) == PermutationSign(order));
			CHECK(RobustPredicates::InSphere3D(a, b, c, d, outside) == -PermutationSign(order));
		} while (std::next_permutation(order.begin(), order.end()));
	}
}

TEST_CASE("Mixed-scale in-sphere retains thin tetrahedra and subnormal query offsets", "[predicate-contract][Predicates]") {
	const Real large = std::numeric_limits<Real>::max() / 4;
	const Real tiny = 4 * std::numeric_limits<Real>::denorm_min();
	const Point3Cartesian a(0, 0, 0), b(0, large, 0), c(large, 0, 0), d(0, 0, tiny);
	CHECK(RobustPredicates::InSphere3D(a, b, c, d, {large, large, tiny}) == 0);
	CHECK(RobustPredicates::InSphere3D(a, b, c, d, {large, large, tiny / 2}) == 1);
	CHECK(RobustPredicates::InSphere3D(a, b, c, d, {large, large, std::nextafter(tiny, Real(0))}) == 1);
	CHECK(RobustPredicates::InSphere3D(a, b, c, d, {large, large, std::nextafter(tiny, Real(1))}) == -1);
	CHECK(RobustPredicates::InSphere3D(a, c, b, d, {large, large, tiny / 2}) == -1);
}

TEST_CASE("Orientation cancellation remains exact with mixed and uniform exponents", "[predicate-contract][Predicates]") {
	const Real epsilon = std::numeric_limits<Real>::epsilon();
	for (int exponent : {0, std::numeric_limits<Real>::max_exponent / 2, std::numeric_limits<Real>::min_exponent / 2}) {
		const Real s = std::ldexp(Real(1), exponent);
		const Point3Cartesian a(0, 0, 0), b(s, s * (1 + epsilon), 0);
		const Point3Cartesian c(s * (1 + epsilon), s * (1 + 2 * epsilon), 0), d(0, 0, s);
		CHECK(RobustPredicates::Orientation3D(a, b, c, d) == 1);
		CHECK(RobustPredicates::Orientation3D(a, c, b, d) == -1);
		CHECK(RobustPredicates::Orientation3D(a, b, c, {0, 0, 0}) == 0);
	}
	const Real large = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 4);
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	CHECK(RobustPredicates::Orientation3D({0, 0, 0}, {large, large * (1 + epsilon), 0},
		{large * (1 + epsilon), large * (1 + 2 * epsilon), 0}, {0, 0, tiny}) == 1);
}

TEST_CASE("Adaptive predicates agree with independent integer determinant oracles", "[predicate-contract][Predicates]") {
	std::mt19937 random(4191);
	std::uniform_int_distribution<int> coordinate(-8, 8);
	for (int sample = 0; sample < 64; ++sample) {
		std::array<std::array<std::int64_t, 3>, 3> orientation2{};
		std::array<std::array<std::int64_t, 4>, 4> circle{}, orientation3{};
		std::array<std::array<std::int64_t, 5>, 5> sphere{};
		std::array<std::array<int, 3>, 5> input{};
		for (int i = 0; i < 5; ++i) {
			auto& p = input[i];
			for (auto& value : p) value = coordinate(random);
			sphere[i] = {p[0], p[1], p[2], p[0] * p[0] + p[1] * p[1] + p[2] * p[2], 1};
			if (i < 4) {
				circle[i] = {p[0], p[1], p[0] * p[0] + p[1] * p[1], 1};
				orientation3[i] = {p[0], p[1], p[2], 1};
			}
			if (i < 3) orientation2[i] = {p[0], p[1], 1};
		}
		const int o2 = IntegerSign(IntegerDeterminant(orientation2)), ic = IntegerSign(IntegerDeterminant(circle));
		const int o3 = IntegerSign(IntegerDeterminant(orientation3)), is = IntegerSign(IntegerDeterminant(sphere));
		for (int exponent : {0, 3 * std::numeric_limits<Real>::max_exponent / 4, 3 * std::numeric_limits<Real>::min_exponent / 4}) {
			for (Real shift : {Real(0), std::ldexp(Real(1), std::numeric_limits<Real>::digits / 4)}) {
				CAPTURE(sample, exponent, shift);
				std::array<Point2Cartesian, 4> xy{};
				std::array<Point3Cartesian, 5> xyz{};
				for (int i = 0; i < 5; ++i) {
					const auto& p = input[i];
					const Real x = std::ldexp(Real(p[0]) + shift, exponent), y = std::ldexp(Real(p[1]) + shift, exponent);
					const Real z = std::ldexp(Real(p[2]) + shift, exponent);
					xyz[i] = Point3Cartesian(x, y, z);
					if (i < 4) xy[i] = Point2Cartesian(x, y);
				}
				CHECK(RobustPredicates::Orientation2D(xy[0], xy[1], xy[2]) == o2);
				CHECK(RobustPredicates::InCircle2D(xy[0], xy[1], xy[2], xy[3]) == ic);
				CHECK(RobustPredicates::Orientation3D(xyz[0], xyz[1], xyz[2], xyz[3]) == o3);
				CHECK(RobustPredicates::InSphere3D(xyz[0], xyz[1], xyz[2], xyz[3], xyz[4]) == is);
			}
		}
	}
}

TEST_CASE("Every robust predicate rejects nonfinite coordinates in every argument", "[predicate-contract][Predicates]") {
	for (Real invalid : {std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity(),
		std::numeric_limits<Real>::quiet_NaN()}) {
		for (int point = 0; point < 5; ++point)
			for (int axis = 0; axis < 3; ++axis) {
				std::array<Point3Cartesian, 5> p{{{0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}, {1, 1, 1}}};
				p[point] = Point3Cartesian(axis == 0 ? invalid : p[point].X(),
					axis == 1 ? invalid : p[point].Y(), axis == 2 ? invalid : p[point].Z());
				CHECK_THROWS_AS(RobustPredicates::InSphere3D(p[0], p[1], p[2], p[3], p[4]), GeometryError);
				if (point < 4) CHECK_THROWS_AS(RobustPredicates::Orientation3D(p[0], p[1], p[2], p[3]), GeometryError);
				if (point < 4 && axis < 2)
					CHECK_THROWS_AS(RobustPredicates::InCircle2D({p[0].X(), p[0].Y()}, {p[1].X(), p[1].Y()},
						{p[2].X(), p[2].Y()}, {p[3].X(), p[3].Y()}), GeometryError);
				if (point < 3 && axis < 2)
					CHECK_THROWS_AS(RobustPredicates::Orientation2D({p[0].X(), p[0].Y()}, {p[1].X(), p[1].Y()}, {p[2].X(), p[2].Y()}),
						GeometryError);
			}
	}
}

TEST_CASE("Repeated zero signed-zero and lower-dimensional predicates stay degenerate", "[predicate-contract][Predicates]") {
	CHECK(RobustPredicates::Orientation2D({-Real(0), 0}, {0, -Real(0)}, {0, 0}) == 0);
	CHECK(RobustPredicates::InCircle2D({0, 0}, {0, 0}, {0, 0}, {0, 0}) == 0);
	CHECK(RobustPredicates::Orientation3D({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}) == 0);
	CHECK(RobustPredicates::InSphere3D({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 0}) == 0);
	CHECK(RobustPredicates::Orientation3D({-Real(0), 0, 0}, {0, -Real(0), 0}, {0, 0, -Real(0)}, {0, 0, 0}) == 0);
	CHECK(RobustPredicates::InSphere3D({0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-Real(0), 0, 0}) == 0);
}
