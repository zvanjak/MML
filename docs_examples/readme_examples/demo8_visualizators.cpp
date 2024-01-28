#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#include "core/Function.h"
#include "core/FunctionHelpers.h"

#include "core/Fields.h"
#include "core/Curves.h"
#include "core/Surfaces.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"


using namespace MML;


void Readme_Real_function_visualization()
{
    RealFunction f1{[](double x) { return sin(x) * (x-3)*(x+5) / sqrt(std::abs(2 - x)); } };

    f1.SerializeEquallySpacedDetailed(-10.0, 10.0, 500, "..\\..\\results\\readme_real_func1.txt");
    std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe ..\\..\\results\\readme_real_func1.txt");

    RealFuncDerived4 f2(f1);        // derivation of f1 function

    f2.SerializeEquallySpacedDetailed(-10.0, 10.0, 100, "..\\..\\results\\readme_real_func2.txt");
    std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe ..\\..\\results\\readme_real_func2.txt");
    
    // shown together
    std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe ..\\..\\results\\readme_real_func1.txt ..\\..\\results\\readme_real_func2.txt");

}

void Readme_Scalar_function_2D_visualization()
{
    // Monkey saddle surface
    ScalarFunction<2> testFunc1{[](const VectorN<Real, 2> &x) { return 3*x[0]/5 * (x[0]*x[0]/25 - 3 * x[1]*x[1]/25); } };

    testFunc1.Serialize2DCartesian(-10.0, 10.0, 20, -10.0, 10.0, 20, "..\\..\\results\\readme_surface1.txt");
    std::system("..\\..\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe ..\\..\\results\\readme_surface1.txt");

    ScalarFunction<2> testFunc2{[](const VectorN<Real, 2> &x) { return (std::abs(x[0])-10) * (std::abs(x[1])-10) * sin(x[0]) * cos(x[1]); } };

    testFunc2.Serialize2DCartesian(-10.0, 10.0, 50, -10.0, 10.0, 50, "..\\..\\results\\readme_surface2.txt");
    std::system("..\\..\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe ..\\..\\results\\readme_surface2.txt");    
}

// void Readme_Scalar_function_3D_visualization()
// {
//     ScalarFunction<3> testFunc{[](const VectorN<Real, 3> &x) { return x[0] * x[1] * x[2]; } };

//     testFunc.Serialize3DCartesian(-10.0, 10.0, 10, -10.0, 10.0, 10, -10.0, 10.0, 10, "func_cart_3d.txt");
// }

void Readme_Parametric_curve_visualization()
{
    // using predefine 3D curves for visualization example
    Curves::HelixCurve              helix(20.0, 2.0);
    Curves::ToroidalSpiralCurve     toroid(20.0);
    
    helix.SerializeCartesian3D(-50.0, 50.0, 1000, "..\\..\\results\\readme_curve_helix.txt");
    std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe ..\\..\\results\\readme_curve_helix.txt");

    toroid.SerializeCartesian3D(0.0, 2 * Constants::PI, 5000, "..\\..\\results\\readme_curve_toroid.txt");
    std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe ..\\..\\results\\readme_curve_toroid.txt");
}

void Readme_Vector_field_visualization()
{
    // let's define potential gravity field of two masses as ScalarFunction
    ScalarFunction<3> gravity_field{ [](const VectorN<Real, 3> &x) 
    {
        const VectorN<Real, 3> x1{ 10.0, 0.0, 0.0 };
        const VectorN<Real, 3> x2{ -10.0, 0.0, 0.0 };
        const Real m1 = 1000.0;
        const Real m2 = 1000.0;
        const Real G = 1.0;
        return -G * m1 / (x - x1).NormL2() - G * m2 / (x - x2).NormL2();
    } };
    // defining force field of two masses as VectorFunction
    VectorFunction<3> gravity_force_field{ [](const VectorN<Real, 3> &x) 
    {
        const VectorN<Real, 3> x1{ 100.0, 0.0, 0.0 }, x2{ -100.0, 0.0, 0.0 };
        const Real m1 = 1000.0, m2 = 1000.0;
        const Real G = 10;
        return -G * m1 * (x - x1) / std::pow((x - x1).NormL2(), 3) - G * m2 * (x - x2) / std::pow((x - x2).NormL2(), 3);
    } };

    gravity_force_field.Serialize3DCartesian(-200.0, 200.0, 15, -200.0, 200.0, 15, -200.0, 200.0, 15, "..\\..\\results\\readme_vector_field.txt", 5);
    auto ret2 = std::system("..\\..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe ..\\..\\results\\readme_vector_field.txt");
}

void Readme_visualizators()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                   README - visualizators                      ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    //Readme_Real_function_visualization();
    //Readme_Scalar_function_2D_visualization();
    Readme_Parametric_curve_visualization();
    //Readme_Vector_field_visualization();
}