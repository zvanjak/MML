# Interpolated functions

Types of functions that can be interpolated:
- IRealFunction
- IScalarFunction<2>
- IParametricCurve

## Real function interpolation

Types of interpolation for Real functions:
- Linear (LinearInterpRealFunc)
- Polynomial (PolynomInterpRealFunc)
- Rational polynomial (RationalInterpRealFunc)
- Splines (SplineInterpRealFunc)
- Barycentric rational (BaryRatInterpRealFunc)

## Scalar function interpolation

Types of interpolation for Scalar functions:
- Bilinear (BilinInterpScalarFunction2D)
- Polynom(PolynomInterpScalarFunction2D)
- Splines (SplineInterpScalarFunction2D)

## Parametric curve interpolation

Types of interpolation for Parametric curves:
- Splines (SplineInterpParametricCurve)

## Example usage

### Real function interpolation

~~~C++
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
~~~

### Parametric curve interpolation

~~~C++
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
~~~

### Surface (ScalarFunction<2>) interpolation

~~~C++
void Docs_Demo_Interpolating_Scalar_function2()
{
  Vector<Real> x1v{ -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
  Vector<Real> x2v{ -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
  Matrix<Real> ym(11, 11, {  10, 10, 20, 10, 40, 30, 20, -1, -2, -6, -4,
                             10, 50, 30, 20, 10, 10, 20, -1, -3, -5,  3,
                             10, 40, 10, 00, 20, 40, 20,  1, -2, -3, -1,
                             10, 30, 20, 10, 20, 60, 00,  1,  0, -1,  0,
                             20, 20, 30, 40, 50, 70, 10,  2,  1,  1,  1,
                             30, 10, 30, 40, 60, 90, 40,  4,  2,  1,  1,
                             40, 20, 30, 40, 50, 60, 50,  1,  3,  2,  1,
                             20, 30, 30, 40, 40, 50, 30,  7,  5,  3,  1,
                             20, 40, 30, 30, 20, 10, 20,  5,  8,  4,  1,
                             10, 20, 30, 20, 10, 10, 10,  4,  3,  2,  0,
                             10, 10, 20, 10, 00, 10, 00,  1,  0,  1,  0 });

  BilinInterpScalarFunction2D	f_bilin(x1v, x2v, ym);

	Visualizer::VisualizeScalarFunc2DCartesian(f_bilin, "docs_dmos_surface1", -5.0, 5.0, 11, -5.0, 5.0, 11, "docs_dmos_surface1.txt");

  PolynomInterpScalarFunction2D	f_poly(x1v, x2v, ym, 3, 3);

  Visualizer::VisualizeScalarFunc2DCartesian(f_poly, "docs_dmos_surface2", -5.0, 5.0, 11, -5.0, 5.0, 11, "docs_dmos_surface2.txt");

  SplineInterpScalarFunction2D	f_spline(x1v, x2v, ym);

  Visualizer::VisualizeScalarFunc2DCartesian(f_spline, "docs_dmos_surface3", -5.0, 5.0, 11, -5.0, 5.0, 11, "docs_dmos_surface3.txt");
}
~~~
