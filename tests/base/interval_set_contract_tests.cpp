#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Intervals.h>
#endif
#include <catch2/catch_all.hpp>
#include <climits>
#include <iterator>
#include <numeric>

using namespace MML;

namespace {
	std::shared_ptr<BaseInterval> Piece(int kind, Real lo, Real hi) {
		switch (kind) {
		case 0: return std::make_shared<OpenInterval>(lo, hi);
		case 1: return std::make_shared<ClosedInterval>(lo, hi);
		case 2: return std::make_shared<OpenClosedInterval>(lo, hi);
		default: return std::make_shared<ClosedOpenInterval>(lo, hi);
		}
	}
	void CheckSamples(const IInterval& interval, int budget, const std::vector<Real>& points) {
		CHECK(points.size() <= static_cast<size_t>(std::max(0, budget)));
		for (size_t i = 0; i < points.size(); ++i) {
			CHECK(std::isfinite(points[i]));
			CHECK(interval.contains(points[i]));
			if (i != 0) CHECK(points[i] > points[i - 1]);
		}
	}
	class ThrowingInterval : public BaseInterval {
	public:
		ThrowingInterval() : BaseInterval(0, EndpointType::CLOSED, 1, EndpointType::CLOSED) {}
		bool contains(Real) const override { throw std::logic_error("Membership callback witness"); }
	};
	class LegacyInterval : public IInterval {
	public:
		Real getLowerBound() const override { return 0; }
		Real getUpperBound() const override { return 1; }
		Real getLength() const override { return 1; }
		bool isContinuous() const override { return true; }
		bool contains(Real x) const override { return x >= 0 && x <= 1; }
		void getEquidistantCovering(int n, std::vector<Real>& points) const override {
			points.clear();
			if (n > 0) points.push_back(Real(0.5));
		}
	};
}

TEST_CASE("Interval contracts predicates respect endpoints empty sets and joint coverage",
	"[intervals][interval-contract]") {
	for (int kind = 0; kind < 4; ++kind) {
		const auto piece = Piece(kind, 0, 1);
		const Interval set(*piece);
		CHECK(set.contains(*piece));
		CHECK(set.contains(OpenInterval(0, 0)));
		CHECK_FALSE(set.intersects(OpenInterval(1, 2)));
	}
	Interval joined;
	joined.AddInterval(ClosedInterval(0, 1)).AddInterval(ClosedInterval(1, 2));
	CHECK(joined.contains(OpenInterval(0, 2)));
	Interval punctured;
	punctured.AddInterval(OpenInterval(0, 1)).AddInterval(OpenInterval(1, 2));
	CHECK_FALSE(punctured.contains(OpenInterval(0, 2)));
	CHECK_FALSE(punctured.intersects(OpenInterval(1, 1)));
	CHECK(Interval().contains(OpenInterval(2, 2)));
	CHECK_FALSE(Interval().contains(ClosedInterval(2, 2)));
}

TEST_CASE("Interval contracts plain algebra agrees with independent pointwise truth tables",
	"[intervals][interval-contract]") {
	for (int leftKind = 0; leftKind < 4; ++leftKind)
		for (int rightKind = 0; rightKind < 4; ++rightKind)
			for (const auto& limits : {std::pair{Real(0), Real(2)}, std::pair{Real(1), Real(3)},
				std::pair{Real(2), Real(3)}, std::pair{Real(1), Real(1)}}) {
				const auto left = Piece(leftKind, 0, 2), right = Piece(rightKind, limits.first, limits.second);
				const Interval a(*left), b(*right);
				const auto both = a & b, either = a | b, difference = a - b, complement = ~a;
				for (int i = -2; i <= 14; ++i) {
					const Real x = Real(i) / 4;
					const bool inA = left->contains(x), inB = right->contains(x);
					CHECK(both.contains(x) == (inA && inB));
					CHECK(either.contains(x) == (inA || inB));
					CHECK(difference.contains(x) == (inA && !inB));
					CHECK(complement.contains(x) == !inA);
				}
			}
}

TEST_CASE("Interval contracts hole bearing algebra fails explicitly instead of returning hull sets",
	"[intervals][interval-contract]") {
	const ClosedIntervalWithReccuringPointHoles source(0, 4, 1, 2);
	Interval holes;
	holes.AddInterval(source);
	CHECK(holes.contains(0));
	CHECK_FALSE(holes.contains(1));
	CHECK_THROWS_AS(holes & holes, NotImplementedError);
	CHECK_THROWS_AS(holes | ClosedInterval(0, 4), NotImplementedError);
	CHECK_THROWS_AS(holes - ClosedInterval(0, 1), NotImplementedError);
	CHECK_THROWS_AS(~holes, NotImplementedError);
	CHECK_THROWS_AS(holes.contains(ClosedInterval(0, 4)), NotImplementedError);
	CHECK_THROWS_AS(holes.intersects(ClosedInterval(0, 4)), NotImplementedError);
	CHECK_THROWS_AS(Interval::Intersection(source, source), NotImplementedError);
	CHECK_FALSE(holes.contains(1));
}

TEST_CASE("Interval contracts empty and singleton coverings respect distinct point budgets",
	"[intervals][interval-contract]") {
	for (int budget : {INT_MIN, -1, 0, 1, 2, 5}) {
		for (int kind = 0; kind < 4; ++kind) {
			const auto interval = Piece(kind, 3, 3);
			std::vector<Real> points{99};
			interval->getEquidistantCovering(budget, points);
			CheckSamples(*interval, budget, points);
			CHECK(points.size() == static_cast<size_t>(budget > 0 && kind == 1 ? 1 : 0));
		}
	}
	const Real next = std::nextafter(Real(1), Real(2));
	for (int kind = 0; kind < 4; ++kind) {
		const auto interval = Piece(kind, 1, next);
		std::vector<Real> points;
		interval->getEquidistantCovering(20, points);
		CheckSamples(*interval, 20, points);
		CHECK(points.size() == static_cast<size_t>(kind == 0 ? 0 : kind == 1 ? 2 : 1));
	}
}

TEST_CASE("Interval contracts composite allocation uses occupied lengths not hull gaps",
	"[intervals][interval-contract]") {
	Interval separated;
	separated.AddInterval(ClosedInterval(0, 1)).AddInterval(ClosedInterval(100, 101));
	std::vector<Real> points;
	separated.getEquidistantCovering(100, points);
	CheckSamples(separated, 100, points);
	REQUIRE(points.size() == 100);
	CHECK(std::count_if(points.begin(), points.end(), [](Real x) { return x <= 1; }) == 50);
	Interval many;
	for (int i = 0; i < 20; ++i) many.AddInterval(ClosedInterval(Real(2 * i), Real(2 * i + 1)));
	for (int budget : {0, 1, 2, 7, 20, 41}) {
		many.getEquidistantCovering(budget, points);
		CheckSamples(many, budget, points);
		CHECK(points.size() == static_cast<size_t>(budget));
	}
	Interval singleton;
	singleton.AddInterval(ClosedInterval(7, 7));
	singleton.getEquidistantCovering(10, points);
	CHECK(points == std::vector<Real>{7});
}

TEST_CASE("Interval contracts periodic sampling never emits excluded points",
	"[intervals][interval-contract]") {
	const ClosedIntervalWithReccuringPointHoles holes(0, 4, 0, 1);
	for (int budget : {1, 2, 5, 17}) {
		std::vector<Real> points;
		holes.getEquidistantCovering(budget, points);
		CheckSamples(holes, budget, points);
		CHECK(points.size() == static_cast<size_t>(budget));
		Interval composite;
		composite.AddInterval(holes);
		composite.getEquidistantCovering(budget, points);
		CheckSamples(composite, budget, points);
		CHECK(points.size() == static_cast<size_t>(budget));
	}
}

TEST_CASE("Interval contracts automatic unbounded window never reverses past a finite bound",
	"[intervals][interval-contract]") {
	const Real outside = 2 * Precision::UnboundedIntervalSamplingLimit;
	for (const auto& ray : {std::shared_ptr<BaseInterval>(std::make_shared<ClosedToInfInterval>(outside)),
		std::shared_ptr<BaseInterval>(std::make_shared<NegInfToClosedInterval>(-outside))}) {
		std::vector<Real> points;
		ray->getEquidistantCovering(5, points);
		CheckSamples(*ray, 5, points);
		CHECK(points.empty());
	}
}

TEST_CASE("Interval contracts finite extreme sampling avoids overflow and remains distinct",
	"[intervals][interval-contract]") {
	const Real max = std::numeric_limits<Real>::max();
	for (const auto& limits : {std::pair{-max, max}, std::pair{max / 2, max}}) {
		const ClosedInterval interval(limits.first, limits.second);
		std::vector<Real> points;
		interval.getEquidistantCovering(5, points);
		CheckSamples(interval, 5, points);
		REQUIRE(points.size() == 5);
		CHECK(points.front() == limits.first);
		CHECK(points.back() == limits.second);
		CHECK(points[2] == std::midpoint(limits.first, limits.second));
	}
}

TEST_CASE("Interval contracts explicit finite windows work through both concrete and abstract views",
	"[intervals][interval-contract]") {
	const Real outside = 2 * Precision::UnboundedIntervalSamplingLimit;
	const ClosedToInfInterval ray(outside);
	const IInterval& view = ray;
	std::vector<Real> points;
	view.getEquidistantCovering(5, points, outside, outside + 4);
	CheckSamples(view, 5, points);
	REQUIRE(points.size() == 5);
	CHECK(points.front() == outside);
	CHECK(points.back() == outside + 4);
	Interval composite;
	composite.AddInterval(OpenInterval(-2, 0)).AddInterval(ClosedInterval(1, 3));
	const IInterval& compositeView = composite;
	compositeView.getEquidistantCovering(20, points, -1, 2);
	CheckSamples(composite, 20, points);
	REQUIRE(points.size() == 20);
	CHECK(points.front() == -1);
	CHECK(points.back() == 2);
	compositeView.getEquidistantCovering(10, points, 0, 0);
	CHECK(points.empty());
	ray.getEquidistantCovering(3, points, outside + 1, outside + 1);
	CHECK(points == std::vector<Real>{outside + 1});
	const ClosedInterval bounded(-2, 4);
	bounded.getEquidistantCovering(1, points);
	CHECK(points == std::vector<Real>{1});
	OpenInterval(-2, 4).getEquidistantCovering(1, points);
	CHECK(points == std::vector<Real>{1});
}

TEST_CASE("Interval contracts unbounded composite sampling clips before integer budget allocation",
	"[intervals][interval-contract]") {
	// The old unbounded-composite NaN-to-int conversion was never executed in the red run.
	Interval whole;
	whole.AddInterval(CompleteRInterval());
	Interval rays;
	rays.AddInterval(NegInfToClosedInterval(-1)).AddInterval(ClosedToInfInterval(1));
	Interval periodic;
	periodic.AddInterval(CompleteRWithReccuringPointHoles(0, 1));
	for (const IInterval* domain : std::initializer_list<const IInterval*>{&whole, &rays, &periodic}) {
		for (int budget : {INT_MIN, 0, 1, 2, 5, 100}) {
			std::vector<Real> points{99};
			domain->getEquidistantCovering(budget, points);
			CheckSamples(*domain, budget, points);
			CHECK(points.size() == static_cast<size_t>(std::max(0, budget)));
			for (Real x : points) CHECK(std::abs(x) <= Precision::UnboundedIntervalSamplingLimit);
		}
	}
}

TEST_CASE("Interval contracts normalized sampling is independent of overlap insertion order and huge hulls",
	"[intervals][interval-contract]") {
	Interval overlap, permuted;
	overlap.AddInterval(ClosedInterval(0, 2)).AddInterval(ClosedInterval(1, 3)).AddInterval(ClosedInterval(0, 2));
	permuted.AddInterval(ClosedInterval(1, 3)).AddInterval(ClosedInterval(0, 2));
	std::vector<Real> first, second;
	overlap.getEquidistantCovering(31, first);
	permuted.getEquidistantCovering(31, second);
	CheckSamples(overlap, 31, first);
	CHECK(first.size() == 31);
	CHECK(first == second);
	const Real max = std::numeric_limits<Real>::max();
	Interval huge;
	huge.AddInterval(ClosedInterval(-max, -max / 2)).AddInterval(ClosedInterval(max / 2, max));
	huge.getEquidistantCovering(20, first);
	CheckSamples(huge, 20, first);
	CHECK(first.size() == 20);
	Interval enormous;
	enormous.AddInterval(ClosedInterval(-max, max));
	enormous.getEquidistantCovering(20, first);
	CheckSamples(enormous, 20, first);
	CHECK(first.size() == 20);
	Interval points;
	points.AddInterval(ClosedInterval(1, 1)).AddInterval(ClosedInterval(3, 3));
	points.getEquidistantCovering(INT_MAX, first);
	CHECK(first == std::vector<Real>{1, 3});
}

TEST_CASE("Interval contracts invalid definitions and windows fail without corrupting sample output",
	"[intervals][interval-contract]") {
	const Real nan = std::numeric_limits<Real>::quiet_NaN(), inf = std::numeric_limits<Real>::infinity();
	for (Real invalid : {nan, inf, -inf}) {
		CHECK_THROWS_AS(ClosedInterval(invalid, 1), NumericInputError);
		CHECK_THROWS_AS(OpenInterval(0, invalid), NumericInputError);
		CHECK_THROWS_AS(ClosedIntervalWithReccuringPointHoles(0, 1, invalid, 1), NumericInputError);
		CHECK_THROWS_AS(CompleteRWithReccuringPointHoles(0, invalid), NumericInputError);
	}
	CHECK_THROWS_AS(CompleteRWithReccuringPointHoles(0, 0), ArgumentError);
	CHECK_THROWS_AS(ClosedIntervalWithReccuringPointHoles(0, 1, 0, 0), ArgumentError);
	const ClosedInterval interval(0, 1);
	const Interval composite(interval);
	for (const IInterval* domain : std::initializer_list<const IInterval*>{&interval, &composite}) {
		std::vector<Real> points{7, 8};
		CHECK_THROWS_AS(domain->getEquidistantCovering(5, points, 1, 0), ArgumentError);
		CHECK(points == std::vector<Real>{7, 8});
		for (Real invalid : {nan, inf, -inf}) {
			CHECK_THROWS_AS(domain->getEquidistantCovering(5, points, invalid, 1), NumericInputError);
			CHECK_THROWS_AS(domain->getEquidistantCovering(5, points, 0, invalid), NumericInputError);
			CHECK(points == std::vector<Real>{7, 8});
		}
		CHECK_NOTHROW(domain->getEquidistantCovering(0, points, nan, inf));
		CHECK(points.empty());
	}
}

TEST_CASE("Interval contracts callback and unresolved hole failures propagate atomically",
	"[intervals][interval-contract]") {
	ThrowingInterval throwing;
	Interval composite;
	composite.AddInterval(throwing);
	for (const IInterval* domain : std::initializer_list<const IInterval*>{&throwing, &composite}) {
		std::vector<Real> points{7, 8};
		CHECK_THROWS_AS(domain->getEquidistantCovering(5, points), std::logic_error);
		CHECK(points == std::vector<Real>{7, 8});
	}
	const ClosedIntervalWithReccuringPointHoles unresolved(0, 1, 0, std::numeric_limits<Real>::min());
	std::vector<Real> points{7, 8};
	CHECK_THROWS_AS(unresolved.getEquidistantCovering(5, points), NumericalMethodError);
	CHECK(points == std::vector<Real>{7, 8});
	const ClosedIntervalWithReccuringPointHoles empty(1, 1, 1, 2);
	empty.getEquidistantCovering(INT_MAX, points);
	CHECK(points.empty());
	LegacyInterval legacy;
	const IInterval& view = legacy;
	CHECK_NOTHROW(view.getEquidistantCovering(2, points));
	CHECK_THROWS_AS(view.getEquidistantCovering(2, points, 0, 1), NotImplementedError);
	CHECK_NOTHROW(view.getEquidistantCovering(0, points, 0, 1));
	CHECK(points.empty());
}

TEST_CASE("Interval contracts infinite endpoint metadata does not erase finite extreme members",
	"[intervals][interval-contract]") {
	const Real max = std::numeric_limits<Real>::max();
	const Interval all{CompleteRInterval()};
	for (Real x : {-max, Real(-1), Real(0), Real(1), max}) {
		const ClosedInterval point(x, x);
		CHECK((all & point).contains(x));
		CHECK(all.contains(point));
		CHECK(all.intersects(point));
	}
}

TEST_CASE("Interval contracts narrow coverings agree with independently enumerated representable members",
	"[intervals][interval-contract]") {
	for (Real lower : {Real(0), Real(1)}) {
		Real upper = lower;
		std::vector<Real> representable{lower};
		for (int step = 0; step < 8; ++step) {
			upper = std::nextafter(upper, Real(2));
			representable.push_back(upper);
		}
		for (int kind = 0; kind < 4; ++kind) {
			const auto interval = Piece(kind, lower, upper);
			std::vector<Real> points, expected;
			interval->getEquidistantCovering(100, points);
			std::copy_if(representable.begin(), representable.end(), std::back_inserter(expected),
				[&interval](Real x) { return interval->contains(x); });
			CheckSamples(*interval, 100, points);
			CHECK(points == expected);
		}
	}
}
