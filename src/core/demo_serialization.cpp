#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "utilities/StdFunctions.h"

#include "base/VectorN.h"
#include "core/Function.h"
#include "core/FieldOperations.h"

#include "basic_types/Fields.h"
#include "basic_types/Curves.h"
#include "basic_types/Surfaces.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

void Demo_Real_function_serialization()
{
    RealFunction f1{[](double x) { return sin(x) * x; } };

    f1.SerializeEquallySpacedDetailed(-10.0, 10.0, 100, "func_sin_x_x.txt");
    auto ret1 = std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe func_sin_x_x.txt");
}

void Demo_Scalar_function_2D_serialization()
{
    ScalarFunction<2> testFunc{[](const VectorN<Real, 2> &x) { return x[0] * x[1]; } };

    testFunc.Serialize2DCartesian(-10.0, 10.0, 20, -10.0, 10.0, 20, "func_cart_2d.txt");
    auto ret2 = std::system("..\\..\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe func_cart_2d.txt");
}

void Demo_Scalar_function_3D_serialization()
{
    ScalarFunction<3> testFunc{[](const VectorN<Real, 3> &x) { return x[0] * x[1] * x[2]; } };

    testFunc.Serialize3DCartesian(-10.0, 10.0, 10, -10.0, 10.0, 10, -10.0, 10.0, 10, "func_cart_3d.txt");
}

void Demo_Parametric_curve_serialization()
{
    Curves::HelixCurve              helix(20.0, 20.0);
    Curves::ToroidalSpiralCurve     toroid(20.0);

    // helix.SerializeCartesian3D(0.0, 2.0 * Constants::PI, 100, "helix.txt");
    // std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe helix.txt");

    toroid.SerializeCartesian3D(0.0, 5.0 * Constants::PI, 500, "toroid.txt");
    auto ret = std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe toroid.txt");
}

void Demo_Vector_field_serialization()
{
    InverseRadialForceFieldCart field_cart(1000.0);

    field_cart.Serialize3DCartesian(-30.0, 30.0, 4, -30.0, 30.0, 4, -30.0, 30.0, 4, "vector_field.txt");
    auto ret2 = std::system("..\\..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe vector_field.txt");
    std::cout << "Return code = " << ret2 << std::endl;    

    InverseRadialForceFieldSpher field_spher(1000.0);

    field_spher.SerializeSpherical(0.1, 50.1, 5, -30.0, 30.0, 4, -30.0, 30.0, 4, "vector_field_spherical.txt");
    auto ret3 = std::system("..\\..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe vector_field_spherical.txt");
    std::cout << "Return code = " << ret3 << std::endl; 
}

void Demo_Serialization()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         SERIALIZATION                         ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Demo_Real_function_serialization();
    // Demo_Scalar_function_2D_serialization();
    Demo_Scalar_function_3D_serialization();
    Demo_Parametric_curve_serialization();
    // Demo_Vector_field_serialization();
}