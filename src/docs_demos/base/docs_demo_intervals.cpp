#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/Intervals.h>
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace MML;

namespace
{
void IntervalTypesExample()
{
	const ClosedInterval closed(0, 10);
	const OpenInterval open(0, 10);
	const OpenClosedInterval openClosed(0, 10);
	const ClosedOpenInterval closedOpen(0, 10);
	const ClosedToInfInterval right(0);
	const OpenToInfInterval openRight(0);
	const NegInfToClosedInterval left(0);
	const NegInfToOpenInterval openLeft(0);
	const CompleteRInterval line;
	const ClosedInterval singleton(3, 3);
	const OpenInterval empty(3, 3);

	if (!closed.contains(0) || !closed.contains(10) ||
		open.contains(0) || open.contains(10) || !open.contains(5) ||
		openClosed.contains(0) || !openClosed.contains(10) ||
		!closedOpen.contains(0) || closedOpen.contains(10) ||
		!right.contains(0) || !right.contains(100) || openRight.contains(0) ||
		!openRight.contains(100) || !left.contains(0) || !left.contains(-100) ||
		openLeft.contains(0) || !openLeft.contains(-100) ||
		!line.contains(-100) || !line.contains(100) ||
		!singleton.contains(3) || empty.contains(3))
		throw InvalidStateError("Interval type example failed");
	std::cout << "Finite endpoints, rays, real line, singleton and empty interval verified.\n";
}

void IntervalSetOperationsExample()
{
	const ClosedInterval sleeping(15, 20);
	const ClosedInterval working(18, 24);
	const Interval optimal = Interval::Intersection(sleeping, working);
	const Interval sleepOnly = Interval::Difference(sleeping, working);
	const Interval either = Interval::Union(sleeping, working);
	const Interval outside = Interval::Complement(sleeping);
	const Interval set(sleeping);

	if (!optimal.contains(ClosedInterval(18, 20)) || optimal.contains(21) ||
		!sleepOnly.contains(15) || sleepOnly.contains(18) ||
		!either.contains(ClosedInterval(15, 24)) ||
		!outside.contains(14) || outside.contains(15) || !outside.contains(21) ||
		!set.contains(sleeping) || !set.intersects(working) ||
		!(set & working).contains(18) || !(set - working).contains(15) ||
		!(set | working).contains(24) || !(~set).contains(14) ||
		!set.intersectWith(working).contains(18) ||
		!set.differenceWith(working).contains(15) ||
		!set.unionWith(working).contains(24) || !set.complement().contains(14))
		throw InvalidStateError("Interval set-operation example failed");

	const Interval touching = Interval::Intersection(OpenInterval(0, 1), OpenInterval(1, 2));
	std::vector<Real> points;
	touching.getEquidistantCovering(1, points);
	if (!points.empty() || Interval(OpenInterval(0, 1)).intersects(OpenInterval(1, 2)))
		throw InvalidStateError("Open intervals must not intersect at an excluded endpoint");
	std::cout << "Set operations: [18,20], [15,18), [15,24] and the complement of [15,20].\n";
}

void IntervalCompoundPropertiesExample()
{
	Interval compound;
	compound.AddInterval(ClosedInterval(0, 5))
		.AddInterval(ClosedInterval(10, 15))
		.AddInterval(ClosedInterval(20, 25));
	if (compound.getLowerBound() != 0 || compound.getUpperBound() != 25 ||
		compound.getLength() != 25 || compound.getMeasure() != 15 ||
		!compound.contains(3) || compound.contains(7) || !compound.contains(12))
		throw InvalidStateError("Compound interval properties failed");

	Interval overlapping;
	overlapping.AddInterval(ClosedInterval(0, 2)).AddInterval(ClosedInterval(1, 3));
	if (overlapping.getMeasure() != 4)
		throw InvalidStateError("Stored-piece measure must count an overlap before normalization");
	overlapping.Normalize();
	if (overlapping.getMeasure() != 3 || overlapping.getLength() != 3 ||
		!overlapping.contains(ClosedInterval(0, 3)) || overlapping.isContinuous())
		throw InvalidStateError("Normalized composite properties failed");
	overlapping = Interval{};
	std::vector<Real> points;
	overlapping.getEquidistantCovering(1, points);
	if (!points.empty()) throw InvalidStateError("Reset interval is not empty");
	std::cout << "Disjoint measure=15, hull length=25; overlapping stored measure=4, normalized=3.\n";
}

void IntervalDomainsExample()
{
	const Real pi = Constants::PI;
	const CompleteRWithReccuringPointHoles tangentDomain(pi / 2, pi);
	Interval threeBranches;
	threeBranches.AddInterval(OpenInterval(-pi / 2, pi / 2))
		.AddInterval(OpenInterval(pi / 2, 3 * pi / 2))
		.AddInterval(OpenInterval(3 * pi / 2, 5 * pi / 2));
	if (!tangentDomain.contains(0) || tangentDomain.contains(pi / 2) ||
		!tangentDomain.contains(-pi) || !threeBranches.contains(0) ||
		!threeBranches.contains(pi) || threeBranches.contains(pi / 2) ||
		threeBranches.contains(-pi))
		throw InvalidStateError("Tangent-domain example failed");
	std::cout << "Global periodic tangent domain and three finite open branches verified.\n";
}

void IntervalSetChainExample()
{
	const ClosedInterval universe(0, 100);
	const ClosedInterval exclude1(20, 30);
	const ClosedInterval exclude2(60, 70);
	const Interval step1 = Interval::Difference(universe, exclude1);
	const Interval remaining = Interval::Difference(step1, exclude2);
	if (!remaining.contains(0) || remaining.contains(20) || remaining.contains(30) ||
		!remaining.contains(40) || remaining.contains(60) || remaining.contains(70) ||
		!remaining.contains(100) || remaining.getMeasure() != 80)
		throw InvalidStateError("Chained interval difference failed");

	const ClosedToInfInterval ray(0);
	const Interval overlap = Interval::Intersection(ray, ClosedInterval(5, 10));
	const Interval cutRay = Interval::Difference(ray, ClosedInterval(5, 10));
	if (!overlap.contains(ClosedInterval(5, 10)) ||
		!cutRay.contains(0) || cutRay.contains(5) || cutRay.contains(10) || !cutRay.contains(11))
		throw InvalidStateError("Ray set operations failed");
	std::cout << "Chained difference: [0,20) U (30,60) U (70,100].\n";
}

void IntervalBracketingExample()
{
	const auto f = [](Real x) { return (x - 1) * (x - 3); };
	const ClosedInterval searchRange(0, 5);
	for (const int sampleCount : {6, 8}) {
		std::vector<Real> points;
		searchRange.getEquidistantCovering(sampleCount, points);
		std::vector<Real> exactRoots;
		std::vector<ClosedInterval> brackets;
		for (size_t i = 0; i < points.size(); ++i) {
			const Real value = f(points[i]);
			if (value == 0) exactRoots.push_back(points[i]);
			if (i == 0) continue;
			const Real previous = f(points[i - 1]);
			if (previous != 0 && value != 0 && std::signbit(previous) != std::signbit(value))
				brackets.emplace_back(points[i - 1], points[i]);
		}
		if (exactRoots.size() + brackets.size() != 2 ||
			(sampleCount == 6 && exactRoots != std::vector<Real>{1, 3}) ||
			(sampleCount == 8 && brackets.size() != 2))
			throw InvalidStateError("Root sampling missed the two roots in this example");
		for (const auto& bracket : brackets)
			if (bracket.contains(1) == bracket.contains(3))
				throw InvalidStateError("Expected exactly one known root in each bracket");
		std::cout << sampleCount << " samples: " << exactRoots.size()
			<< " exact roots, " << brackets.size() << " sign-change brackets.\n";
	}
}

void IntervalSamplingExample()
{
	std::vector<Real> points;
	ClosedInterval(0, 10).getEquidistantCovering(5, points);
	if (points != std::vector<Real>{0, Real(2.5), 5, Real(7.5), 10})
		throw InvalidStateError("Closed covering grid failed");
	OpenInterval(0, 10).getEquidistantCovering(5, points);
	if (points.size() != 5 || points.front() != std::nextafter(Real(0), Real(10)) ||
		points.back() != std::nextafter(Real(10), Real(0)))
		throw InvalidStateError("Open covering endpoints failed");
	ClosedInterval(3, 3).getEquidistantCovering(100, points);
	if (points != std::vector<Real>{3}) throw InvalidStateError("Singleton covering failed");

	const Real outside = 2 * Precision::UnboundedIntervalSamplingLimit;
	const ClosedToInfInterval ray(outside);
	const IInterval& view = ray;
	view.getEquidistantCovering(5, points);
	if (!points.empty()) throw InvalidStateError("Far ray must miss the default window");
	view.getEquidistantCovering(5, points, outside, outside + 4);
	if (points.size() != 5 || points.front() != outside || points.back() != outside + 4)
		throw InvalidStateError("Explicit window covering failed");
	std::cout << "Closed/open grids, distinct singleton and explicit ray window verified.\n";
}

void IntervalErrorsExample()
{
	bool reversedRejected = false;
	try {
		(void)ClosedInterval(10, 5);
	} catch (const ArgumentError& error) {
		reversedRejected = true;
		std::cout << error.what() << "\n";
	}
	bool nonfiniteRejected = false;
	try {
		(void)ClosedInterval(0, std::numeric_limits<Real>::infinity());
	} catch (const NumericInputError& error) {
		nonfiniteRejected = true;
		std::cout << error.what() << "\n";
	}
	if (!reversedRejected || !nonfiniteRejected)
		throw InvalidStateError("Invalid interval definitions were accepted");
}
}

void Docs_Demo_Intervals()
{
	IntervalTypesExample();
	IntervalSetOperationsExample();
	IntervalCompoundPropertiesExample();
	IntervalDomainsExample();
	IntervalSetChainExample();
	IntervalBracketingExample();
	IntervalSamplingExample();
	IntervalErrorsExample();
}

void Docs_Demo_IntervalContracts()
{
	const auto check = [](const IInterval& domain, int budget, const std::vector<Real>& samples) {
		if (samples.size() > static_cast<size_t>(budget))
			throw NumericalMethodError("Interval sample budget exceeded");
		for (size_t i = 0; i < samples.size(); ++i)
			if (!std::isfinite(samples[i]) || !domain.contains(samples[i]) ||
				(i > 0 && !(samples[i - 1] < samples[i])))
				throw NumericalMethodError("Interval sampling returned invalid or repeated points");
	};
	ClosedInterval interval(0, 10), sub(2, 8);
	const bool isContained = Interval(interval).contains(sub);
	const OpenInterval open(0, 1);
	if (!isContained || !Interval(open).contains(open) || Interval(open).intersects(OpenInterval(1, 2)))
		throw NumericalMethodError("Endpoint-aware set predicates failed");

	Interval separated;
	separated.AddInterval(ClosedInterval(0, 1)).AddInterval(ClosedInterval(100, 101));
	std::vector<Real> points;
	separated.getEquidistantCovering(100, points);
	check(separated, 100, points);
	if (points.size() != 100 ||
		std::count_if(points.begin(), points.end(), [](Real x) { return x <= 1; }) != 50)
		throw NumericalMethodError("Occupied-length allocation failed");
	ClosedInterval(3, 3).getEquidistantCovering(100, points);
	if (points != std::vector<Real>{3}) throw NumericalMethodError("Singleton samples were duplicated");

	const Real outside = 2 * Precision::UnboundedIntervalSamplingLimit;
	ClosedToInfInterval ray(outside);
	ray.getEquidistantCovering(5, points);
	if (!points.empty()) throw NumericalMethodError("Automatic window escaped the ray");
	const IInterval& view = ray;
	view.getEquidistantCovering(5, points, outside, outside + 4);
	check(view, 5, points);
	if (points.size() != 5 || points.front() != outside || points.back() != outside + 4)
		throw NumericalMethodError("Explicit sampling window failed");

	ClosedIntervalWithReccuringPointHoles source(0, 4, 1, 2);
	Interval punctured;
	punctured.AddInterval(source);
	punctured.getEquidistantCovering(9, points);
	check(punctured, 9, points);
	if (points.size() != 9 || punctured.contains(1) || punctured.contains(3))
		throw NumericalMethodError("Periodic holes were lost");
	bool rejected = false;
	try {
		(void)(punctured & punctured);
	} catch (const NotImplementedError& error) {
		rejected = true;
		std::cout << error.what() << "\n";
	}
	if (!rejected) throw InvalidStateError("Hole-bearing algebra was not rejected");
	std::cout << "Interval endpoint predicates, holes, distinct budgets and finite windows verified.\n";
}
