#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "algorithms/PathIntegration.h"
#endif

using namespace std;
using namespace MML;

void Calc_curve_length()
{
	ParametricCurve<3> circle([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), 0}; });
	ParametricCurve<3> helix([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), t}; });

	Real len = PathIntegration::ParametricCurveLength(circle, 0, 2 * 3.14);

	std::cout << "Length of circle is " << PathIntegration::ParametricCurveLength(circle, 0, 2 * 3.14) << endl;
	std::cout << "Length of helix is " << PathIntegration::ParametricCurveLength(helix, 0, 2 * 3.14) << endl;
}

void Calc_work_integral()
{
	ScalarFunction<3>  potential([](const VectorN<Real, 3>& x) { return Real{ 10.0 } / x.NormL2(); });
	ParametricCurve<3> circle([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), 1}; });

	for (auto phi = 0.1; phi < 2 * 3.14159; phi += 0.25)
	{
		std::cout << "Work integral for phi : " << phi << " is : " << PathIntegration::WorkIntegral(potential, circle, 0, phi, 1e-03) << endl;
	}
}

void Demo_Path_Integration()
{
	std::cout << endl;
	std::cout << "***********************************************************************" << endl;
	std::cout << "****                     PATH INTEGRATION                          ****" << endl;
	std::cout << "***********************************************************************" << endl;

	Calc_curve_length();
	Calc_work_integral();
}