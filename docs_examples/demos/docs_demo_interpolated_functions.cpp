#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/InterpolatedFunction.h"
#include "core/Serializer.h"
#include "core/Visualizer.h"
#endif

using namespace MML;

void Docs_Demo_Interpolating_Real_function1_equally_spaced()
{
  Vector<Real> x_val({ 0.0,  1.0, 2.0, 3.0, 4.0,  5.0,  6.0,  7.0, 8.0, 9.0, 10.0 });
  Vector<Real> y_val({ 5.0, -2.0, 2.0, 3.0, 7.0, -6.0, -3.0, -2.0, 3.0, 0.0, -5.0 });

  int   NumInterpPnt = x_val.size();
  Real  x1 = x_val[0];
  Real  x2 = x_val[x_val.size() - 1];

  // these are the ways we can interpolate Real function
  LinearInterpRealFunc    f_linear(x_val, y_val);
  PolynomInterpRealFunc   f_polynom(x_val, y_val, 3);
  SplineInterpRealFunc    f_spline(x_val, y_val);
  BaryRatInterpRealFunc   f_baryrat(x_val, y_val, 3);

  Visualizer::VisualizeMultiRealFunction({ &f_linear, &f_polynom, &f_spline, &f_baryrat }, "docs_dmos_multi_real_func1", x1, x2, 500, "docs_dmos_multi_real_func1.txt");
}

void Docs_Demo_Interpolating_Real_function2()
{
  Vector<Real> x_val({ 0.0, 1.0, 5.0, 6.0, 7.0, 10.0, 15.0, 17.0, 19.0, 20.0 });
  Vector<Real> y_val({ 5.0, -2.0, 2.0, 3.0, 7.0, -6.0, -3.0, -2.0, 3.0, 0.0 });

  int   NumInterpPnt = x_val.size();
  Real  x1 = x_val[0];
  Real  x2 = x_val[x_val.size() - 1];

  // these are the ways we can interpolate Real function
  LinearInterpRealFunc    f_linear(x_val, y_val);
  PolynomInterpRealFunc   f_polynom(x_val, y_val, 3);
  SplineInterpRealFunc    f_spline(x_val, y_val);
  BaryRatInterpRealFunc   f_baryrat(x_val, y_val, 3);

  Visualizer::VisualizeMultiRealFunction({ &f_linear, &f_polynom, &f_spline, &f_baryrat }, "docs_dmos_multi_real_func2", x1, x2, 500, "docs_dmos_multi_real_func2.txt");
}

void Docs_Demo_Interpolating_Real_function3()
{
  // we have function, but we want to approximate it!!!
  const int NumInterpPnt = 12;
  const Real x1 = 0, x2 = 10.0;

  // we will use this as test func
  RealFunction test_func{ [](Real x) { return (Real)sin(x) * (1 + x * x / 2); } };

  // and using our helper, available for all real functions, create data for interpolation
  Vector<Real> x_val(NumInterpPnt), y_val(NumInterpPnt);
  test_func.GetValues(x1, x2, NumInterpPnt, x_val, y_val);

  // these are the ways we can interpolate Real function
  LinearInterpRealFunc    f_linear(x_val, y_val);
  PolynomInterpRealFunc   f_polynom(x_val, y_val, 3);
  SplineInterpRealFunc    f_spline(x_val, y_val);
  BaryRatInterpRealFunc   f_baryrat(x_val, y_val, 3);

  // situation - we need different number of points for different functions
  Serializer::SaveRealFuncEquallySpacedDetailed(test_func, "docs_demo_interp_test_func", x1, x2, 100, "..\\..\\results\\docs_demo_interp_test_func.txt");
  Serializer::SaveRealFuncEquallySpacedDetailed(f_linear, "docs_demo_interp_linear_5_pnt", x1, x2, 500, "..\\..\\results\\docs_demo_interp_linear_5_pnt.txt");
  Serializer::SaveRealFuncEquallySpacedDetailed(f_polynom, "docs_demo_interp_polynom_5_pnt", x1, x2, 100, "..\\..\\results\\docs_demo_interp_polynom_5_pnt.txt");
  Serializer::SaveRealFuncEquallySpacedDetailed, (x1, x2, 100, "..\\..\\results\\docs_demo_interp_spline_5_pnt.txt");
  Serializer::SaveRealFuncEquallySpacedDetailed(f_baryrat, "docs_demo_interp_baryrat_5_pnt", x1, x2, 100, "..\\..\\results\\docs_demo_interp_baryrat_5_pnt.txt");

  const char* cmd = "..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe"
    " ..\\..\\results\\docs_demo_interp_test_func.txt"
    " ..\\..\\results\\docs_demo_interp_linear_5_pnt.txt"
    " ..\\..\\results\\docs_demo_interp_polynom_5_pnt.txt"
    " ..\\..\\results\\docs_demo_interp_spline_5_pnt.txt"
    " ..\\..\\results\\docs_demo_interp_baryrat_5_pnt.txt";
  std::system(cmd);
}

void Docs_Demo_Interpolating_Parametric_curve()
{
  Matrix<Real> values(10, 3, { 0.0,  1.0,  5.0,
                              -1.0,  2.0,  3.0, 
                              -4.0,  3.0, -2.0, 
                              20.0,  5.0, -2.0, 
                               2.0,  3.0,  7.0, 
                              -6.0, -3.0, -2.0, 
                               3.0,  0.0,  5.0, 
                               5.0, -2.0,  1.0, 
                               6.0, -1.0, -2.0, 
                              -2.0,  3.0,  0.0 });

  SplineInterpParametricCurve<3>  f_spline(values);

  Visualizer::VisualizeParamCurve3D(f_spline, "docs_dmos_param_curve", 0.0, 1.0, 100, "docs_dmos_param_curve.txt");
}

void Docs_Demo_Interpolating_Scalar_function2()
{
  Vector<Real> x1v{ -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
  Vector<Real> x2v{ -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
  Matrix<Real> ym(11, 11, {  10, 10, 20, 10, 40, 30, 20, -10, -20, -60, -40,
                             10, 50, 30, 20, 10, 10, 20, -10, -30, -50,  30,
                             10, 40, 10, 00, 20, 40, 20,  10, -20, -30, -10,
                             10, 30, 20, 10, 20, 60, 00,  10,  00, -10,   0,
                             20, 20, 30, 40, 50, 70, 10,  20,  10,  10,  10,
                             30, 10, 30, 40, 60, 90, 40,  40,  20,  10,  10,
                             40, 20, 30, 40, 50, 60, 50,  10,  30,  20,  10,
                             20, 30, 30, 40, 40, 50, 30,  70,  50,  30,  10,
                             20, 40, 30, 30, 20, 10, 20,  50,  80,  40,  10,
                             10, 20, 30, 20, 10, 10, 10,  40,  30,  20,   0,
                             10, 10, 20, 10, 00, 10, 00,  10,   0,  10,   0 });

  BilinInterpScalarFunction2D	f_bilin(x1v, x2v, ym);

	//Visualizer::VisualizeScalarFunc2DCartesian(f_bilin, "docs_dmos_surface1", -5.0, 5.0, 11, -5.0, 5.0, 11, "docs_dmos_surface1.txt");

  PolynomInterpScalarFunction2D	f_poly(x1v, x2v, ym, 3, 3);

  Visualizer::VisualizeScalarFunc2DCartesian(f_poly, "docs_dmos_surface2", -5.0, 5.0, 11, -5.0, 5.0, 11, "docs_dmos_surface2.txt");

  SplineInterpScalarFunction2D	f_spline(x1v, x2v, ym);

  //Visualizer::VisualizeScalarFunc2DCartesian(f_spline, "docs_dmos_surface3", -5.0, 5.0, 11, -5.0, 5.0, 11, "docs_dmos_surface3.txt");
}

void Docs_Demo_Interpolated_functions()
{
  //Docs_Demo_Interpolating_Real_function1_equally_spaced();
  //Docs_Demo_Interpolating_Real_function2();
  //Docs_Demo_Interpolating_Real_function2();
  // Docs_Demo_Interpolating_Parametric_curve();
  Docs_Demo_Interpolating_Scalar_function2();
}