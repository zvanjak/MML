#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Matrix.h"
#include "mml/base/VectorTypes.h"
#include "mml/base/InterpolatedFunction.h"

#include "mml/core/Curves.h"
#include "mml/core/Integration/PathIntegration.h"

#include "mpl/Gravity/GravityBase.h"
#endif

#include <iostream>

using namespace MML;
using namespace MPL;

// Verify vector and scalar field operations for gravity fields
void Demo_Field_operations()
{
	// TODO: Implementation pending
}

// For two given points in space, calculates the potential at those points, 
// and then calculates the path integral of the force field along a line and a spline curve between those two points
void Verify_path_integrals()
{
	Real G = 1.0;
	Real mass = 100.0;

	GravityPotentialField		potentialField(G, mass, Vec3Cart(0, 0, 0));
	GravityForceField				forceField(G, mass, Vec3Cart(0, 0, 0));

	// our two points in space
	Vec3Cart point1(10.0, 10.0, -5.0);
	Vec3Cart point2(-5.0, 5.0, 8.0);

	Real potential1 = potentialField(point1);
	Real potential2 = potentialField(point2);

	std::cout << "Potential at point 1: " << potential1 << std::endl;
	std::cout << "Potential at point 2: " << potential2 << std::endl;
	std::cout << "Difference in potential: " << potential1 - potential2 << std::endl;

	// first, create a line between those two points
	Curves::LineCurve line(0.0, Pnt3Cart(point1.X(), point1.Y(), point1.Z()), 
												 1.0, Pnt3Cart(point2.X(), point2.Y(), point2.Z()));

	Real integral = PathIntegration::LineIntegral(forceField, line, 0.0, 1.0);

	std::cout << "Path integral between point 1 and point 2: " << integral << std::endl;

	// second, we'll create a spline parametric curve between those two points
	// with some points added in between to make path "curvaceous"
	Matrix<Real> curve_points1{ 5, 3,
														 { point1.X(), point1.Y(), point1.Z(),
															20.0, 5.0, -10.0,
															100.0, 150.5, -30.0,		// sending it far away
															5.0, 10.5, 1.0,
															point2.X(), point2.Y(), point2.Z() }
	};
	SplineInterpParametricCurve<3> curve(0.0, 1.0, curve_points1);

	Real integral2 = PathIntegration::LineIntegral(forceField, curve, 0.0, 1.0);

	std::cout << "Path integral along spline curve between point 1 and point 2: " << integral2 << std::endl;
}

// Simulate gravitational slingshot of Voyager 1 by Jupiter
void Demo_Voyager_Jupiter_Slingshot()
{
	// TODO: Implementation pending
}
