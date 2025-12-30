#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/VectorN.h"

#include "mml/base/Function.h"
#include "mml/base/InterpolatedFunction.h"

#include "mml/core/FunctionHelpers.h"
#include "mml/core/Curves.h"

#include "mml/algorithms/ODESystemSolver.h"
#include "mml/algorithms/ODESystemStepCalculators.h"
#include "mml/algorithms/ODESystemSteppers.h"

#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"
#endif

//#include "../test_data/parametric_curves_test_bed.h"
//#include "../test_data/Fields.h"


using namespace MML;
using namespace MML::Curves;

void Example3_WPF_Real_function_visualization()
{
	RealFunctionFromStdFunc f1{ [](Real x) { return sin(x) * (x - Real(3)) * (x + Real(5)) / exp(Real(0.2) / std::abs(Real(2) - x)); } };

	Visualizer::VisualizeRealFunction(f1, "Simple function", -10.0, 10.0, 501, 
																				"example3_wpf_real_func1.txt");

	RealFuncDerived4 f1_der(f1);      // 4th order derivation of f1 function

	Visualizer::VisualizeRealFunction(f1_der, "Derivation", -10.0, 10.0, 501, 
																						"example3_wpf_real_func2.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction({ &f1, &f1_der }, 
																					"Function and derivation together", { "Func", "Derivation" },
																					-10.0, 10.0, 501, "example3_wpf_multi_real_func.txt");

	// Solving Lorentz system as MULTI_REAL_FUNCTION example
	ODESystem lorentz_system(3, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		{
			const Real sigma = 10.0;
			const Real rho = 28.0;
			const Real beta = 8.0 / 3.0;
			dxdt[0] = sigma * (x[1] - x[0]);
			dxdt[1] = x[0] * (rho - x[2]) - x[1];
			dxdt[2] = x[0] * x[1] - beta * x[2];
		} 
	);

	Vector<Real> ystart({ 2.0, 1.0, 1.0 });
	ODESystemSolver<RK5_CashKarp_Stepper> solver(lorentz_system);
	ODESystemSolution sol = solver.integrate(ystart, 0.0, 50.0, 0.1, 1e-08, 0.01);

	Visualizer::VisualizeODESysSolAsMultiFunc(sol, "Lorentz system", std::vector<std::string>{"x", "y", "z"}, "example3_wpf_lorentz.txt");

	PolynomInterpRealFunc 	comp1 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	comp2 = sol.getSolAsPolyInterp(1, 3);
	PolynomInterpRealFunc 	comp3 = sol.getSolAsPolyInterp(2, 3);

	Visualizer::VisualizeRealFunction(comp1, "Lorentz system - component 1",
																		0.0, 50.0, 500, "example3_wpf_lorentz1.txt");
	Visualizer::VisualizeRealFunction(comp2, "Lorentz system - component 2",
																		0.0, 50.0, 500, "example3_wpf_lorentz2.txt");
	Visualizer::VisualizeRealFunction(comp3, "Lorentz system - component 3",
																		0.0, 50.0, 500, "example3_wpf_lorentz3.txt");

}

void Example3_WPF_Scalar_function_2D_visualization()
{
	// Monkey saddle surface
	ScalarFunction<2> testFunc1{ [](const VectorN<Real, 2>& x) 
		{ 
			return 3 * x[0] / 5 * (x[0] * x[0] / 25 - 3 * x[1] * x[1] / 25); 
		} 
	};

	Visualizer::VisualizeScalarFunc2DCartesian(testFunc1, "Monkey saddle", 
																						 -10.0, 10.0, 20, -10.0, 10.0, 20, 
																						 "example3_wpf_surface1.txt");

	ScalarFunction<2> testFunc2{ [](const VectorN<Real, 2>& x) 
		{ 
			return (std::abs(x[0]) - 10) * (std::abs(x[1]) - 10) * sin(x[0]) * cos(x[1]); 
		} 
	};

	Visualizer::VisualizeScalarFunc2DCartesian(testFunc2, "Example surface", 
																						 -10.0, 10.0, 50, -10.0, 10.0, 50, 
																						 "example3_wpf_surface2.txt");
}

void Example3_WPF_Parametric_curve_2D_visualization()
{
	// using predefine 2D curves for visualization example
	Circle2DCurve   circle(5.0);

	LemniscateCurve lemniscate;
	Visualizer::VisualizeParamCurve2D(lemniscate, "Lemniscate", 0.0, 2 * Constants::PI, 101, 
																		"example3_wpf_curve_lemniscate.txt");
	
	LogSpiralCurve  log_spiral2(-0.2), log_spiral3(-0.3), log_spiral4(-0.4), log_spiral5(-0.5), log_spiral6(-0.6);
	LogSpiralCurve  log_spiral7(-0.7), log_spiral8(-0.8), log_spiral9(-0.9);

	Visualizer::VisualizeParamCurve2D(log_spiral2, "Log spiral", 0.0, 10.0, 101, 
																		"example3_wpf_curve_log_spiral.txt");

	Visualizer::VisualizeMultiParamCurve2D({ &log_spiral2, &log_spiral3, &log_spiral4, &log_spiral5, 
																					 &log_spiral6, &log_spiral7, &log_spiral8, &log_spiral9 },
																					 "Log spiral", 0.0, 20.0, 201, "example3_wpf_curve_log_spiral_multi");

	DeltoidCurve deltoid1(1), deltoid2(2), deltoid3(5), deltoid4(8), deltoid5(10), deltoid6(15);

	Visualizer::VisualizeParamCurve2D(deltoid1, "Deltoid", 0.0, 2 * Constants::PI, 101, 
																		"example3_wpf_curve_deltoid.txt");

	Visualizer::VisualizeMultiParamCurve2D({ &deltoid1, &deltoid2, &deltoid3, &deltoid4 },
																					"Deltoid", 0.0, 2 * Constants::PI, 101, 
																					"example3_wpf_curve_deltoid_multi");

	AstroidCurve astroid1(1), astroid2(2.0), astroid3(3.0), astroid4(4.0);
	AstroidCurve astroid5(5.0), astroid6(6.0), astroid7(7.0), astroid8(8.0), astroid9(9.0);

	Visualizer::VisualizeParamCurve2D(astroid1, "Astroid", 0.0, 2 * Constants::PI, 101, 
																		"example3_wpf_curve_astroid.txt");

	Visualizer::VisualizeMultiParamCurve2D({ &astroid1, &astroid2, &astroid3, &astroid4, &astroid5, 
																					 &astroid6, &astroid7, &astroid8, &astroid9 },
																					"Astroid", 0.0, 2 * Constants::PI, 101, 
																					"example3_wpf_curve_astroid_multi");
	
	ArchimedeanSpiralCurve arch_spiral1(1.0), arch_spiral2(2.0), arch_spiral3(3.0), arch_spiral4(4.0), arch_spiral5(5.0);
	ArchimedeanSpiralCurve arch_spiral6(6.0), arch_spiral7(7.0), arch_spiral8(8.0), arch_spiral9(9.0);

	Visualizer::VisualizeParamCurve2D(arch_spiral1, "Archimedean spiral", 0.0, 11.0, 101, "example3_wpf_curve_archimedean_spiral.txt");

	Visualizer::VisualizeMultiParamCurve2D({ &arch_spiral1, &arch_spiral2, &arch_spiral3, &arch_spiral4,
																					 & arch_spiral5,& arch_spiral6,& arch_spiral7,& arch_spiral8,& arch_spiral9 },
																					"Archimedean spiral", 0.0, 3.5 * Constants::PI, 101, 
																					"example3_wpf_curve_archimedean_spiral_multi");
}

void Example3_WPF_Parametric_curve_3D_visualization()
{
	// using predefined 3D curves for visualization example
	HelixCurve          helix(20.0, 2.0);
	ToroidalSpiralCurve toroid1(5, 20.0), toroid2(5, 10.0), toroid3(5, 5.0);

	Visualizer::VisualizeParamCurve3D(helix, "helix", -50.0, 50.0, 1000, 
																		"example3_wpf_curve_helix.txt");
	
	Visualizer::VisualizeParamCurve3D(toroid1, "toroid", 0.0, 5 * Constants::PI, 1000, 
																		"example3_wpf_curve_toroid.txt");

	// visualize all three curves
	Visualizer::VisualizeMultiParamCurve3D({ &toroid1, &toroid2, &toroid3 }, "toroidal spiral", 
																					0.0, 5 * Constants::PI, 1000, 
																					"example3_wpf_curve_toroidal_spiral.txt");

	// Constructing curve in 3d from 10 points
	// First case - we already have points in Matrix
	Matrix<Real> curve_points1{ 10, 3,
														 {0.0, 50.0, -100.0,
															-100.0,-100.0, 50.0,
															-150.0,-180.0, 150.0,
															-80.0, -80.0, 200.0,
															130.0,-150.0, 100.0,
															150.0,-100.0, -50.0,
															100.0, 180.0, 90.0,
															 50.0, 140.0, 100.0,
															-70.0, 100.5, 120.0,
															 20.0, -100.5, -50.0} };

	// we can easily create spline parametric curve
	SplineInterpParametricCurve<3> curve(0.0, 1.0, curve_points1);

	Visualizer::VisualizeParamCurve3D(curve, "Curve from points", 0.0, 1.0, 1000, 
																		"example3_wpf_curve_spline.txt");

	// second case - we have points in a vector
	std::vector<VectorN<Real, 3>> vec_points{ { 50.0, -50.0, 20.0},
																						{ 10, -30.0, 40.0},
																						{120,  50, -30},
																						{ 50, 100,  0},
																						{ 20, -50,  10},
																						{-50, -100, 80},
																						{-90, -120, 50},
																						{ 50, -20, -70},
																						{-20, -25,  10},
																						{ 80,  50,  60} };

	// create Matrix from vector of points
	Matrix<Real> curve_points2(vec_points.size(), 3);
	for (int i = 0; i < vec_points.size(); i++)
	{
		curve_points2(i, 0) = vec_points[i][0];
		curve_points2(i, 1) = vec_points[i][1];
		curve_points2(i, 2) = vec_points[i][2];
	}

	SplineInterpParametricCurve<3> curve2(0.0, 1.0, curve_points2);

	Visualizer::VisualizeParamCurve3D(curve2, "Curve from points", 0.0, 1.0, 1000, 
																		"example3_wpf_curve_spline2.txt");

	// shown together
	Visualizer::VisualizeMultiParamCurve3D({ &curve, &curve2 }, "Curve from points", 0.0, 1.0, 1000, "example3_wpf_curve_spline_multi.txt");

	// Lorentz system as 3D parametric curve
	ODESystem lorentz_system(3, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		{
			const Real sigma = 10.0;
			const Real rho = 28.0;
			const Real beta = 8.0 / 3.0;
			dxdt[0] = sigma * (x[1] - x[0]);
			dxdt[1] = x[0] * (rho - x[2]) - x[1];
			dxdt[2] = x[0] * x[1] - beta * x[2];
		}
	);

	Vector<Real> ystart({ -2.0, 3.0, 1.0 });
	ODESystemSolver<RK5_CashKarp_Stepper> solver(lorentz_system);
	
	ODESystemSolution sol = solver.integrate(ystart, 0.0, 50.0, 0.001, 1e-10, 0.001);
	
	SplineInterpParametricCurve<3> solCurve1 = sol.getSolAsParamCurve3D(0, 1, 2);
	
	Visualizer::VisualizeODESysSolAsParamCurve3(sol, 0, 1, 2, "Lorentz system", 
																							"example3_wpf_lorentz_system1.txt");

	Vector<Real> ystart1({15.0,-35.0,-10.0}), ystart2({-10.0,25.0,-30.0}), ystart3({25.0,-10.0,40.0});
	
	sol = solver.integrate(ystart1, 0.0, 50.0, 0.001, 1e-10, 0.001);
	SplineInterpParametricCurve<3> solCurve2 = sol.getSolAsParamCurve3D(0, 1, 2);
	sol = solver.integrate(ystart2, 0.0, 50.0, 0.001, 1e-10, 0.001);
	SplineInterpParametricCurve<3> solCurve3 = sol.getSolAsParamCurve3D(0, 1, 2);
	sol = solver.integrate(ystart3, 0.0, 50.0, 0.001, 1e-10, 0.001);
	SplineInterpParametricCurve<3> solCurve4 = sol.getSolAsParamCurve3D(0, 1, 2);


	// visualize together
	Visualizer::VisualizeMultiParamCurve3D({ &solCurve1, &solCurve2, &solCurve3, &solCurve4 }, 
																					"Lorentz attractor", 
																					0.0, 50.0, 10000, "example3_wpf_lorentz_system_multi.txt");

	Visualizer::VisualizeODESysSolAsParamCurve3(sol, 0, 1, 2, "Lorentz system", "example3_wpf_lorentz_system2.txt");
}

void Example3_WPF_Vector_field_2D_visualization()
{
	VectorFunction<2> rotational_field{ [](const VectorN<Real, 2>& x)
		{
			Real norm = x.NormL2();
			return VectorN<Real, 2>{-x[1] / norm, x[0] / norm};
		} 
	};

	Visualizer::VisualizeVectorField2DCartesian(rotational_field, "rotational_field", 
																							-10.0, 10.0, 20, -10.0, 10.0, 20, 
																							"example3_wpf_vector_field_2d.txt");
}

void Example3_WPF_Vector_field_3D_visualization()
{
	// let's define potential gravity field of two masses as ScalarFunction
	ScalarFunction<3> gravity_field{ [](const VectorN<Real, 3>& x)
	{
			const VectorN<Real, 3> x1{ 100.0, 0.0, 0.0 };
			const VectorN<Real, 3> x2{ -100.0, 50.0, 0.0 };
			const Real m1 = 100.0;
			const Real m2 = -20.0;
			const Real G = 10.0;
			return -G * m1 / (x - x1).NormL2() - G * m2 / (x - x2).NormL2();
	} };
	
	// defining force field of two masses as VectorFunction
	VectorFunction<3> gravity_force_field{ [](const VectorN<Real, 3>& x)
	{
			const VectorN<Real, 3> x1{ 100, 0, 0}, x2{ -100, 0, 0};
			const Real m1 = 100.0, m2 = 100.0;
			const Real G = 10;
			return -G * m1 * (x - x1) / std::pow((x - x1).NormL2(), 3) -
						  G * m2 * (x - x2) / std::pow((x - x2).NormL2(), 3);
	} };

	Visualizer::VisualizeVectorField3DCartesian(gravity_force_field, "Gravity field", 
																							-200, 200, 30, -200, 200, 30, -200, 200, 30, 
																							"example3_wpf_vector_field_3d.txt");
}

void Example3_WPF_visualizators()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                  EXAMPLE - visualizators                      ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Example3_WPF_Real_function_visualization();
	Example3_WPF_Scalar_function_2D_visualization();
	Example3_WPF_Parametric_curve_2D_visualization();
	Example3_WPF_Parametric_curve_3D_visualization();
	Example3_WPF_Vector_field_2D_visualization();
	Example3_WPF_Vector_field_3D_visualization();
}
