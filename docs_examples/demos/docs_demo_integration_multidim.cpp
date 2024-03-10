#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/Integration.h"
#endif

using namespace MML;

void Docs_Demo_Integration_multidim()
{
	std::cout << "***********************************************************" << std::endl;
	std::cout << "*****               Multidim integration            *******" << std::endl;

	// 2D integration of constant scalar 2D function (ie. we'll get the area of the surface)
	ScalarFunction<2> f2([](const VectorN<Real, 2>& x) { return Real{ 1 }; });

	// we integrate over circle with radius 2
	Real val = IntegrateSurface(f2, IntegrationMethod::GAUSS10,
		-2, 2,              // x range
		[](Real x) { return -sqrt(4 - x * x); },   // y range lower limit
		[](Real x) { return sqrt(4 - x * x); });   // y range upper limit

	std::cout << "Calc. area = " << val << ", exact value: 4 * PI = " << 4 * Constants::PI << std::endl;

	// 3D integration of constant scalar 3D function (ie. we'll get the volume of the solid)
	ScalarFunction<3> f3([](const VectorN<Real, 3>& x) { return Real{ 1 }; });

	// integration over sphere of radius 1
	Real vol = IntegrateVolume(f3,
		-1, 1,
		[](Real x) { return -sqrt(1 - x * x); },
		[](Real x) { return sqrt(1 - x * x); },
		[](Real x, Real y) { return -sqrt(1 - x * x - y * y); },
		[](Real x, Real y) { return sqrt(1 - x * x - y * y); });

	std::cout << "Calc. vol. = " << vol << ", exact value: 4/3 * PI = " << 4.0 / 3.0 * Constants::PI << std::endl;

	/* OUTPUT
			Calc. area = 12.57211164, exact value: 4 * PI = 12.56637061
			Calc. vol. = 4.190703882, exact value: 4/3 * PI = 4.188790205
	*/
}