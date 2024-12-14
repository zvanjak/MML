#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/FunctionHelpers.h"
#include "core/Derivation.h"
#include "core/Integration.h"

#include "tools/ConsolePrinter.h"
#endif


using namespace MML;

void Readme_integrating_functions()
{
	RealFunction f1{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };
	RealFunction f1_integral{ [](Real x) { return (Real)(x * (-0.5 * x * cos(x) + sin(x))); } };

	double a = 0.0;
	double b = 10.0;
	double int_trap = IntegrateTrap(f1, a, b);
	// we can use default Integrate routine (set to IntegrateSimpson), requires precision
	double int_def = Integrate(f1, a, b, 1e-04, nullptr);

	std::cout << "Integrating function f1 from " << a << " to " << b << std::endl;
	std::cout << "Exact integral   = " << f1_integral(b) - f1_integral(a) << std::endl;
	std::cout << "IntegrateTrap    = " << int_trap << std::endl;

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

void Readme_Integration_precision()
{
	// our test function
	RealFunction f{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };
	// and its exact integral
	RealFunction f_int{ [](Real x) { return (Real)(x * (-0.5 * x * cos(x) + sin(x))); } };

	double x1 = 0.0, x2 = 10.0;
	double err_sum1 = 0.0, err_sum2 = 0.0, err_sum3 = 0.0;
	const int numIntervals = 10;

	std::cout << "\nAVERAGE INTEGRATION ERROR FOR DIFFERENT INTEGRATORS" << std::endl;
	std::cout << "   Interval           Exact int.      Trap         Trap err.           Simpson      Simpson  err.       Romberg       Romberg err.          " << std::endl;
	std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

	for (int i = 1; i < numIntervals; i++) {
		double x = x1 + (x2 - x1) * i / (numIntervals - 1);

		double integral = f_int(x) - f_int(x1);

		double int_trap = MML::IntegrateTrap(f, x1, x, 1e-3, nullptr);

		double err1 = int_trap - integral;

		std::cout << "[" << std::fixed
			<< std::setw(6) << std::setprecision(3) << x1 << ", "
			<< std::setw(6) << std::setprecision(3) << x << "]  "
			<< std::setw(13) << std::setprecision(8) << integral << " "
			<< std::setw(13) << std::setprecision(8) << int_trap << "   "
			<< std::scientific << std::setw(15) << err1 << "   " << std::fixed
			<< std::endl;

		err_sum1 += std::abs(err1);
	}

	std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "Total abs error =                                "
		<< std::scientific << err_sum1 << "                    "
		<< std::endl;

	/* OUTPUT
	AVERAGE INTEGRATION ERROR FOR DIFFERENT INTEGRATORS
		 Interval           Exact int.      Trap         Trap err.           Simpson      Simpson  err.       Romberg       Romberg err.
	------------------------------------------------------------------------------------------------------------------------------------
	[ 0.000,  1.111]     0.72190120    0.72191911    1.79168330e-05      0.72190120   -2.20768637e-09      0.72190120    1.74171788e-12
	[ 0.000,  2.222]     3.26424438    3.26411014   -1.34238455e-04      3.26424432   -5.66896650e-08      3.26424438    5.37498046e-09
	[ 0.000,  3.333]     4.81851793    4.81806182   -4.56106873e-04      4.81851807    1.38713693e-07      4.81851833    4.00889971e-07
	[ 0.000,  4.444]    -1.67104024   -1.67124534   -2.05095264e-04     -1.67103857    1.67321678e-06     -1.67104024    4.24239088e-10
	[ 0.000,  5.556]   -15.21897376  -15.21421664    4.75711807e-03    -15.21897406   -2.96904277e-07    -15.21897377   -9.82410597e-09
	[ 0.000,  6.667]   -18.11382964  -18.10862331    5.20633435e-03    -18.11384814   -1.84979419e-05    -18.11382939    2.51947842e-07
	[ 0.000,  7.778]     5.45250452    5.45320716    7.02639753e-04      5.45247119   -3.33281723e-05      5.45250452    1.91555927e-09
	[ 0.000,  8.889]    38.50671730   38.49414238   -1.25749240e-02     38.50675166    3.43570339e-05     38.50673638    1.90788434e-05
	[ 0.000, 10.000]    36.51336534   36.50710489   -6.26045832e-03     36.51354664    1.81295066e-04     36.51339597    3.06248120e-05
	------------------------------------------------------------------------------------------------------------------------------------
	Total abs error =                                3.03148319e-02                    2.69645946e-04                    5.03740338e-05
	*/
}

void Readme_integration()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                      README - integration                     ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Readme_integrating_functions();
	Readme_Integration_precision();
}