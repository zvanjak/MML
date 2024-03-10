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
#include "core/Serializer.h"
#include "core/Visualizer.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"


using namespace MML;


void Readme_Real_function_visualization()
{
    RealFunction f1{[](Real x) { return sin(x) * (x-3)*(x+5) / sqrt(std::abs(2 - x)); } };

    Visualizer::VisualizeRealFunction(f1, "readme_real_func1" , -10.0, 10.0, 500, "readme_real_func1.txt");

    RealFuncDerived4 f2(f1);        // derivation of f1 function

    Visualizer::VisualizeRealFunction(f2, "readme_real_func2", -10.0, 10.0, 100, "readme_real_func2.txt");
    
    // shown together
    Visualizer::VisualizeMultiRealFunction({&f1, &f2}, "readme_multi_real_func", -10.0, 10.0, 500, "readme_multi_real_func.txt" );
}

void Readme_Scalar_function_2D_visualization()
{
    // Monkey saddle surface
    ScalarFunction<2> testFunc1{[](const VectorN<Real, 2> &x) { return 3*x[0]/5 * (x[0]*x[0]/25 - 3 * x[1]*x[1]/25); } };

    Visualizer::VisualizeScalarFunc2DCartesian(testFunc1, "readme_surface1" , -10.0, 10.0, 20, -10.0, 10.0, 20, "readme_surface1.txt");

    ScalarFunction<2> testFunc2{[](const VectorN<Real, 2> &x) { return (std::abs(x[0])-10) * (std::abs(x[1])-10) * sin(x[0]) * cos(x[1]); } };

    Visualizer::VisualizeScalarFunc2DCartesian(testFunc2, "readme_surface2" , -10.0, 10.0, 50, -10.0, 10.0, 50, "readme_surface2.txt");
}

// void Readme_Scalar_function_3D_visualization()
// {
//     ScalarFunction<3> testFunc{[](const VectorN<Real, 3> &x) { return x[0] * x[1] * x[2]; } };

//     testFunc.Serialize3DCartesian(-10.0, 10.0, 10, -10.0, 10.0, 10, -10.0, 10.0, 10, "func_cart_3d.txt");
// }

void Readme_Parametric_curve_visualization()
{
    // using predefine 3D curves for visualization example
    Curves3D::HelixCurve              helix(20.0, 2.0);
    Curves3D::ToroidalSpiralCurve     toroid(Real{20.0});
    
    Visualizer::VisualizeParamCurve3D(helix, "helix" , -50.0, 50.0, 1000, "readme_curve_helix.txt");
    Visualizer::VisualizeParamCurve3D(toroid, "toroid", 0.0, 2 * Constants::PI, 5000, "readme_curve_toroid.txt");
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

    Visualizer::VisualizeVectorField3DCartesian(gravity_force_field, "gravity_force_field" , -200.0, 200.0, 15, -200.0, 200.0, 15, -200.0, 200.0, 15, "readme_vector_field.txt");
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