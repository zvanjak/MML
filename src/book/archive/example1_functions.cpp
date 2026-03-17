#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#endif

using namespace MML;


/////////////////////////////////////////////////////////////////////////////////////
// CASE 1 - we already have standalone function providing needed calculation
// 
// examples of already existing functions

Real Example1_functions_RealFunc(Real x) {
  return sin(x) * (1 + x * x / 2);
}

// simple scalar function, returning 1 / r, for given point
Real Example1_functions_ScalarFunc(const VectorN<Real, 3>& x) {
  return 1 / x.NormL2();
}

// rotating vector field in 2D
VectorN<Real, 2> Example1_functions_VectorFunc(const VectorN<Real, 2>& x) {
  return VectorN<Real, 2>{-x[1], x[0]};
}

// helix curve in 3D
VectorN<Real, 3> Example1_functions_ParamCurve(Real t) {
  return  VectorN<Real, 3>{ cos(t), sin(t), t};
}

// simple surface
VectorN<Real, 3> Example1_functions_ParamSurface(Real x, Real y) {
  return VectorN<Real, 3>{x* y, 2 * x * y, 3 * x};
}

// little bit more realistic scalar function :)
Real Example1_functions_two_masses_gravity_potential(const VectorN<Real, 3>& x)
{
  const VectorN<Real, 3> x1{ 10.0, 0.0, 0.0 };
  const VectorN<Real, 3> x2{ -10.0, 0.0, 0.0 };
  const Real m1 = 1000.0;
  const Real m2 = 1000.0;
  const Real G = 1.0;
  return -G * m1 / (x - x1).NormL2() - G * m2 / (x - x2).NormL2();
}

// EXAMPLE WITH MAGNETIC FILED

// something even more complex?

void Example1_functions_case_1_usage()
{
  // creating a function objects from an already existing (standalone) function
  RealFunction            fReal(Example1_functions_RealFunc);
  ScalarFunction<3>       fScalar(Example1_functions_ScalarFunc);
  ScalarFunction<3>       fScalar2(Example1_functions_two_masses_gravity_potential);
  VectorFunction<2>       fVector(Example1_functions_VectorFunc);
  ParametricCurve<3>      paramCurve(Example1_functions_ParamCurve);
  ParametricSurfaceRect<3>    paramSurface(Example1_functions_ParamSurface);
}

/////////////////////////////////////////////////////////////////////////////////////
// CASE 2 - creating function directly from lambda

void Example1_functions_case_2_usage()
{
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "***  CASE 2 - creating function directly from lambda" << std::endl;

  RealFunction      f2{ [](Real x) { return (Real)sin(x) * (1 + x * x / 2); } };
  ScalarFunction<3> fScalar([](const VectorN<Real, 3>& x) { return 1 / x.NormL2(); });
  VectorFunction<3> fVector([](const VectorN<Real, 3>& x)
    {
      return VectorN<Real, 3>{-x[2], x[0] * x[1], x[0]};
    });
  ParametricCurve<3>    paramCurve([](Real t) { return VectorN<Real, 3>{t, 2 * t, 3 * t}; });
  ParametricSurfaceRect<3>  paramSurface([](Real x, Real y) 
    { 
      return VectorN<Real, 3>{x* y, 2 * x * y, 3 * x}; 
    });

  ScalarFunction<3> two_masses_gravity_field_potential{ [](const VectorN<Real, 3>& x)
              {
                  const VectorN<Real, 3> x1{ 10.0, 0.0, 0.0 };
                  const VectorN<Real, 3> x2{ -10.0, 0.0, 0.0 };
                  const Real m1 = 1000.0;
                  const Real m2 = 1000.0;
                  const Real G = 1.0;
                  return -G * m1 / (x - x1).NormL2() - G * m2 / (x - x2).NormL2();
              }
  };

  // evaluating the functions
  std::cout << "f2(0.0) = " << f2(0.0) << std::endl;
  std::cout << "fScalar({1, 1, 1})   = " << fScalar(VectorN<Real, 3>{1, 1, 1}) << std::endl;
  std::cout << "fScalar({-1, -1, 1}) = " << fScalar(VectorN<Real, 3>{-1, -1, 1}) << std::endl;
  std::cout << "fScalar({2, 2, 2})   = " << fScalar(VectorN<Real, 3>{2, 2, 2}) << std::endl;
  std::cout << "fVector({1, 1, 1})   = " << fVector(VectorN<Real, 3>{1, 1, 1}) << std::endl;
  std::cout << "fVector({-1, -1, 1}) = " << fVector(VectorN<Real, 3>{-1, -1, 1}) << std::endl;
  std::cout << "fVector({2, 2, 2})   = " << fVector(VectorN<Real, 3>{2, 2, 2}) << std::endl;
  std::cout << "paramCurve(1.0)        = " << paramCurve(1.0) << std::endl;
  std::cout << "paramSurface(1.0, 1.0) = " << paramSurface(1.0, 1.0) << std::endl;

  /* OUTPUT
  f2(0.0) = 0
  fScalar({1, 1, 1})   = 1
  fScalar({-1, -1, 1}) = -1
  fScalar({2, 2, 2})   = 2
  fVector({1, 1, 1})   = [              0,               1,               0]
  fVector({-1, -1, 1}) = [              0,               1,               0]
  fVector({2, 2, 2})   = [              0,               4,               0]
  paramCurve(1.0)        = [              1,               2,               3]
  paramSurface(1.0, 1.0) = [              1,               2,               3]
  */
}

/////////////////////////////////////////////////////////////////////////////////////
// CASE 3 - we have function that has some additional parameters that it depends on 
//           (and that we might wanna vary easily)
double funcWithParams(double x, double p1, double p2)
{
  return p1 * cos(x) + p2 * sin(x);
}

double twoMassesGravityPotential(const VectorN<Real, 3>& x, 
                              Real m1, Real m2, Real G,
                              const VectorN<Real, 3>& x1, const VectorN<Real, 3>& x2)
{
  return -G * (m1 / (x - x1).NormL2() + m2 / (x - x2).NormL2());
}

class TwoMassesGravityPotential : public IScalarFunction<3>
{
  Real _m1, _m2, _G;
  VectorN<Real, 3> _x1, _x2;

public:
  TwoMassesGravityPotential(Real m1, Real m2, Real G,
                            const VectorN<Real, 3>& x1, const VectorN<Real, 3>& x2)
    : _m1(m1), _m2(m2), _G(G), _x1(x1), _x2(x2) {
  }

  void SetM1(Real m1) { _m1 = m1; }
  void SetM2(Real m2) { _m2 = m2; }
  void SetX1(const VectorN<Real, 3>& x1) { _x1 = x1; }
  void SetX2(const VectorN<Real, 3>& x2) { _x2 = x2; }

  Real operator()(const VectorN<Real, 3>& x) const
  {
    return -_G * (_m1 / (x - _x1).NormL2() + _m2 / (x - _x2).NormL2());
  }
};

// primjer param krivulje tako definirane
// infinite line - magnetsko polje, primjer vector field

void Example1_functions_case_3_usage()
{
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "***  CASE 3 - we have simple function, but that has some additional parameters that it depend on (and that we might wanna vary easily)" << std::endl;

	// creating a function object from a standalone function with additional parameters
  RealFunctionFromStdFunc f1{ [](Real x) { return funcWithParams(x, Real(1.0), Real(2.0)); } };
	
  std::cout << "f1(0.0) = " << f1(0.0) << std::endl;
	std::cout << "f1(1.0) = " << f1(1.0) << std::endl;
	std::cout << "f1(2.0) = " << f1(2.0) << std::endl;


  TwoMassesGravityPotential f3(1000.0, 1000.0, 1.0, VectorN<Real, 3>{10.0, 0.0, 0.0}, VectorN<Real, 3>{-10.0, 0.0, 0.0});

  std::cout << "f3({1, 1, 1})   = " << f3(VectorN<Real, 3>{1, 1, 1}) << std::endl;
  std::cout << "f3({-1, -1, 1}) = " << f3(VectorN<Real, 3>{-1, -1, 1}) << std::endl;
  std::cout << "f3({2, 2, 2})   = " << f3(VectorN<Real, 3>{2, 2, 2}) << std::endl;

  /* OUTPUT
  f3({1, 1, 1})   = -199.9312235
  f3({-1, -1, 1}) = -199.9312235
  f3({2, 2, 2})   = -198.9618408
  */
}

////////////////////////////////////////////////////////////////////////////////////
// CASE 4 - there is an external class that has EXACT FUNCTION you want to use (but you can't change the class)
class BigComplexClassYouCantChange1 {
public:
  double _param1, _param2;
  double Weight(double x) const { return 58 * _param1; }

  double UsefulRealFunc(double x) const { return (_param1 * cos(x) + _param2 * sin(x)) / Weight(x); }
};

void Example1_functions_case_4_usage()
{
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "***  CASE 4 - there is an external class that has EXACT FUNCTION you want to use (but you can't change the class)" << std::endl;

  BigComplexClassYouCantChange1 bigObj;

  bigObj._param1 = 1.0;
  bigObj._param2 = 2.0;

  RealFunctionFromStdFunc f4(std::function<Real(Real)>{ [&bigObj](Real x)
    {
      return bigObj.UsefulRealFunc(x);
    }
  });

  std::cout << "f4(0.0) = " << f4(0.0) << std::endl;
  std::cout << "f4(1.0) = " << f4(1.0) << std::endl;
  std::cout << "f4(2.0) = " << f4(2.0) << std::endl;

  bigObj._param1 = -2.0;
  bigObj._param2 = 5.0;

  std::cout << "After change of object parameters:" << std::endl;

  std::cout << "f4(0.0) = " << f4(0.0) << std::endl;
  std::cout << "f4(1.0) = " << f4(1.0) << std::endl;
  std::cout << "f4(2.0) = " << f4(2.0) << std::endl;

  /* OUTPUT
  f4(0.0) = 0.01724137931
  f4(1.0) = 0.03833179785
  f4(2.0) = 0.02418013823
  After change of object parameters:
  f4(0.0) = 0.01724137931
  f4(1.0) = -0.02695474407
  f4(2.0) = -0.04636880006
  */
}

////////////////////////////////////////////////////////////////////////////////////
// CASE 5 - there is an external class that has data relevant to calculation (but you can't change the class)
class BigComplexClassYouCantChange2 {
public:
  double _param1, _param2;
  double Weight(double x) const { return 58 * _param1; }
};

// Create a helper wrapper class, inherit it from IRealFunction and use it as RealFunction
class BigComplexRealFunc2 : public IRealFunction
{
  const BigComplexClassYouCantChange2& _ref;
public:
  BigComplexRealFunc2(const BigComplexClassYouCantChange2& bigClass)
    : _ref(bigClass) {
  }

  Real operator()(Real x) const {
    return (_ref._param1 * cos(x) + _ref._param2 * sin(x)) / _ref.Weight(x);
  }
};

void Example1_functions_case_5_usage()
{
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "***  CASE 5 - there is an external class that has data relevant to calculation (but you can't change the class)" << std::endl;
  BigComplexClassYouCantChange2 bigObj;

  bigObj._param1 = 1.0;
  bigObj._param2 = 2.0;

  BigComplexRealFunc2    f5(bigObj);     // usable RealFunction object

  std::cout << "f5(0.0) = " << f5(0.0) << std::endl;
  std::cout << "f5(1.0) = " << f5(1.0) << std::endl;
  std::cout << "f5(2.0) = " << f5(2.0) << std::endl;

  bigObj._param1 = -2.0;
  bigObj._param2 = 5.0;

  // our f5 function is now CHANGED!!!
  std::cout << "After change of object parameters:" << std::endl;

  std::cout << "f5(0.0) = " << f5(0.0) << std::endl;
  std::cout << "f5(1.0) = " << f5(1.0) << std::endl;
  std::cout << "f5(2.0) = " << f5(2.0) << std::endl;

  /* OUTPUT
  f5(0.0) = 0.01724137931
  f5(1.0) = 0.03833179785
  f5(2.0) = 0.02418013823
  After change of object parameters:
  f5(0.0) = 0.01724137931
  f5(1.0) = -0.02695474407
  f5(2.0) = -0.04636880006
  */
}

void Example1_Functions()
{
  Example1_functions_case_1_usage();
  Example1_functions_case_2_usage();
  Example1_functions_case_3_usage();
  Example1_functions_case_4_usage();
  Example1_functions_case_5_usage();
}