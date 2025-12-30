#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/core/Curves.h"
#endif

#include <iostream>

using namespace MML;
using namespace MML::Curves;

void Chapter20_Diff_geometry_curves()
{
	// testing curves
	Circle3DXY circle(1.0);

	auto tang = circle.getTangent(1.0);
	std::cout << "Tangent at t=1.0: " << tang << std::endl;
}
