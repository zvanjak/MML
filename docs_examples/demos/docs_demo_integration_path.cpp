#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "algorithms/PathIntegration.h"
#endif

using namespace MML;


void Docs_Demo_Calc_curve_length()
{
  std::cout << "***********************************************************" << std::endl;
  std::cout << "*****            Calculating curve length           *******" << std::endl;

  ParametricCurve<3> circle([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), 0}; });
  ParametricCurve<3> helix([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), t}; });

  Real len = PathIntegration::ParametricCurveLength(circle, 0, 2 * 3.14);

  std::cout << "Length of circle is " << PathIntegration::ParametricCurveLength(circle, 0, 2 * 3.14) << std::endl;
  std::cout << "Length of helix is  " << PathIntegration::ParametricCurveLength(helix, 0, 2 * 3.14) << std::endl;
}

void Docs_Demo_Calc_work_integral()
{
  std::cout << "***********************************************************" << std::endl;
  std::cout << "*****            Calculating work integral          *******" << std::endl;

  ScalarFunction<3>  potential([](const VectorN<Real, 3>& x) { return Real{ 10.0 } / x.NormL2(); });
  ParametricCurve<3> circle([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), 1}; });

  Real phi = 6;
  //for (auto phi = 0.1; phi < 2 * 3.14159; phi += 0.25)
  {
    std::cout << "Work integral for phi : " << phi << " is : " << PathIntegration::WorkIntegral(potential, circle, 0, phi, 1e-03) << std::endl;
  }
}

void Docs_Demo_Integration_path()
{
  Docs_Demo_Calc_curve_length();
	Docs_Demo_Calc_work_integral();
}
