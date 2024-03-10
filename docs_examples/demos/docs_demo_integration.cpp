#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/Integration.h"
#endif

using namespace MML;

void Docs_Demo_Integration()
{
	RealFunction f1{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };
	RealFunction f1_integral{ [](Real x) { return (Real)(x * (-0.5 * x * cos(x) + sin(x))); } };

	double a = 0.0;
	double b = 10.0;
	double int_trap  = IntegrateTrap(f1, a, b);
	double int_simp  = IntegrateSimpson(f1, a, b);
	double int_romb  = IntegrateRomberg(f1, a, b);
	double int_gauss = IntegrateGauss10(f1, a, b);
	// we can use default Integrate routine (set to IntegrateSimpson), requires precision
	double int_def = Integrate(f1, a, b, 1e-04);

	std::cout << "Integrating function f1 from " << a << " to " << b << std::endl;
	std::cout << "Exact integral   = " << f1_integral(b) - f1_integral(a) << std::endl;
	std::cout << "IntegrateTrap    = " << int_trap << std::endl;
	std::cout << "IntegrateSimpson = " << int_simp << std::endl;
	std::cout << "IntegrateRomberg = " << int_romb << std::endl;
	std::cout << "IntegrateGauss10 = " << int_gauss << std::endl;

/* OUTPUT
Integrating function f1 from 0 to 10
Exact integral   = 36.5134
IntegrateTrap    = 36.5133
IntegrateSimpson = 36.5134
IntegrateRomberg = 36.5134
IntegrateGauss10 = 36.5134
*/
}