#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Function.h"
#include "mml/core/Integration.h"

#endif

using namespace MML;

void Demo_Integration_1D()
{
	RealFunction f1{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };
	RealFunction f1_integral{ [](Real x) { return (Real)(x * (-0.5 * x * cos(x) + sin(x))); } };

	double a = 0.0;
	double b = 10.0;
	double int_trap  = IntegrateTrap(f1, a, b);
	double int_simp  = IntegrateSimpson(f1, a, b);
	double int_gauss = IntegrateGauss10(f1, a, b);

	std::cout << "Integrating function f1 from " << a << " to " << b << std::endl;
	std::cout << "Exact integral   = " << f1_integral(b) - f1_integral(a) << std::endl;
	std::cout << "IntegrateTrap    = " << int_trap << std::endl;
	std::cout << "IntegrateSimpson = " << int_simp << std::endl;
	std::cout << "IntegrateGauss10 = " << int_gauss << std::endl;

	/* OUTPUT
	Integrating function f1 from 0 to 10
	Exact integral   = 36.51336534
	IntegrateTrap    = 36.51297408
	IntegrateSimpson = 36.51337665
	IntegrateGauss10 = 36.51336529
	*/
}

void Demo_Integration_2D()
{
	// 2D integration of constant 2D function over area (ie. we'll get the area of the surface)
	ScalarFunction<2> f2([](const VectorN<Real, 2>& x) { return Real{ 1 }; });

	// we integrate over circle with radius 2
	Real	val = Integrate2D(f2, IntegrationMethod::GAUSS10,
													-2, 2,              // x range
													[](Real x) { return -sqrt(4 - x * x); },   // y range lower limit
													[](Real x) { return sqrt(4 - x * x); });   // y range upper limit

	std::cout << "Circle area = " << val << ", exact value: 4 * PI = " << 4 * Constants::PI << std::endl;

	// calculate volume of half-sphere
	ScalarFunction<2> f3([](const VectorN<Real, 2>& x) { return Real{ sqrt(4 - POW2(x[0]) - POW2(x[1])) }; });
	Real val2 = Integrate2D(f3, IntegrationMethod::GAUSS10,
													-2, 2,              // x range
													[](Real x) { return -sqrt(4 - x * x); },   // y range lower limit
													[](Real x) { return sqrt(4 - x * x); });   // y range upper limit	

	std::cout << "Half-sphere volume  = " << val2 << ", exact value = " << 4.0 / 6 * 8 * Constants::PI << std::endl;
	
	/* OUTPUT
			Calc. area = 12.57211164, exact value: 4 * PI = 12.56637061
			Calc. half-sphere volume  = 16.76281553, exact value = 16.75516082
	*/
}

void Demo_Integration_3D()
{
	// 3D integration of constant scalar 3D function (ie. we'll get the volume of the solid)
	ScalarFunction<3> f3([](const VectorN<Real, 3>& x) { return Real{ 1 }; });

	// integration over cube with sides 2, 3, 4
	Real	cubeVol = Integrate3D(f3,
															-1, 1,
															[](Real x) -> Real { return Real(-1.5); },
															[](Real x) -> Real { return Real(1.5); },
															[](Real x, Real y) -> Real { return Real(-2.0); },
															[](Real x, Real y) -> Real { return Real(2.0); });

	std::cout << "Cube vol. = " << cubeVol << ", exact value: 2 * 3 * 4 = " << 2 * 3 * 4 << std::endl;

	// integration over sphere of radius 1
	Real	sphereVol = Integrate3D(f3,
																-1, 1,
																[](Real x) -> Real { return -sqrt(Real(1) - x * x); },
																[](Real x) -> Real { return sqrt(Real(1) - x * x); },
																[](Real x, Real y) -> Real { return -sqrt(Real(1) - x * x - y * y); },
																[](Real x, Real y) -> Real { return sqrt(Real(1) - x * x - y * y); });

	std::cout << "Sphere vol. = " << sphereVol << ", exact value: 4/3*PI = " << 4.0/3*Constants::PI << std::endl;

	/* OUTPUT
			Cube vol. = 24, exact value: 2 * 3 * 4 = 24
			Sphere vol. = 4.1907, exact value: 4/3 * PI = 4.18879
	*/
}

void Demo_Integration()
{
	Demo_Integration_1D();
	Demo_Integration_2D();
	Demo_Integration_3D();
}