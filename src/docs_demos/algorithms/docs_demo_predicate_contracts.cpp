#include <mml/algorithms/CompGeometry/RobustPredicates.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

void Docs_Demo_PredicateContracts() {
	using namespace MML;
	using MML::CompGeometry::RobustPredicates;

	const Real large = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent / 2);
	const Real tiny = std::ldexp(Real(1), std::numeric_limits<Real>::min_exponent + std::numeric_limits<Real>::digits);
	const Point3Cartesian a(0, 0, tiny), b(large, 0, 0), c(0, large, 0), origin(0, 0, 0);
	if (RobustPredicates::Orientation3D(a, b, c, origin) != 1 ||
		RobustPredicates::Orientation3D(a, c, b, origin) != -1)
		throw std::runtime_error("Mixed-scale orientation must retain the tiny coordinate");

	const Point3Cartesian on(large, large, tiny);
	const Point3Cartesian inside(large, large, std::nextafter(tiny, Real(0)));
	const Point3Cartesian outside(large, large, std::nextafter(tiny, Real(1)));
	if (RobustPredicates::InSphere3D(origin, c, b, a, on) != 0 ||
		RobustPredicates::InSphere3D(origin, c, b, a, inside) != 1 ||
		RobustPredicates::InSphere3D(origin, c, b, a, outside) != -1)
		throw std::runtime_error("In-sphere signs must distinguish adjacent represented offsets");

	std::cout << "Exact predicate range, cancellation and mixed-scale contracts verified.\n";
}
