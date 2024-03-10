#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/StdFunctions.h"

#include "base/VectorN.h"

#include "core/Function.h"
#include "core/FieldOperations.h"
#include "core/Fields.h"
#include "core/Curves.h"
#include "core/Surfaces.h"
#include "core/Serializer.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

void Demo_Real_function_serialization()
{
    RealFunction f1{[](Real x) { return sin(x) * x; } };

    //f1.SerializeEquallySpacedDetailed(-10.0, 10.0, 100, "..\\..\\results\\func_sin_x_x.txt");
    Serializer::SaveRealFuncEquallySpacedDetailed(f1, "Func y = sin(x) * x", - 10.0, 10.0, 100, "..\\..\\results\\func_sin_x_x.txt");
    auto ret1 = std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe ..\\..\\results\\func_sin_x_x.txt");
}

void Demo_Scalar_function_2D_serialization()
{
    ScalarFunction<2> testFunc{[](const VectorN<Real, 2> &x) { return x[0] * x[1]; } };

    Serializer::SaveScalarFunc2DCartesian(testFunc , "Simple func z = x * y",  - 10.0, 10.0, 20, -10.0, 10.0, 20, "..\\..\\results\\func_cart_2d.txt");
    auto ret2 = std::system("..\\..\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe ..\\..\\results\\func_cart_2d.txt");
}

void Demo_Scalar_function_3D_serialization()
{
    ScalarFunction<3> testFunc{[](const VectorN<Real, 3> &x) { return x[0] * x[1] * x[2]; } };

    Serializer::SaveScalarFunc3DCartesian(testFunc, "Simple func f = x * y * z", - 10.0, 10.0, 10, -10.0, 10.0, 10, -10.0, 10.0, 10, "..\\..\\results\\func_cart_3d.txt");
}

void Demo_Parametric_curve_serialization()
{
    Curves3D::HelixCurve              helix(20.0, 20.0);
    Curves3D::ToroidalSpiralCurve     toroid(Real{20.0});

    // helix.SerializeCartesian3D(0.0, 2.0 * Constants::PI, 100, "helix.txt");
    // std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe helix.txt");

    Serializer::SaveParamCurveCartesian3D(toroid, "Toroid", 0.0, 5.0 * Constants::PI, 500, "..\\..\\results\\toroid.txt");
    auto ret = std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe ..\\..\\results\\toroid.txt");
}

void Demo_Vector_field_serialization()
{
    Fields::InverseRadialForceFieldCart field_cart(1000.0);

    Serializer::SaveVectorFunc3DCartesian(field_cart, "Simple vector field",  - 30.0, 30.0, 4, -30.0, 30.0, 4, -30.0, 30.0, 4, "..\\..\\results\\vector_field.txt");
    auto ret2 = std::system("..\\..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe ..\\..\\results\\vector_field.txt");
    std::cout << "Return code = " << ret2 << std::endl;    

    Fields::InverseRadialForceFieldSpher field_spher(1000.0);

    Serializer::SaveVectorFuncSpherical(field_spher, "Simple vector field - spherical", 0.1, 50.1, 5, -30.0, 30.0, 4, -30.0, 30.0, 4, "..\\..\\results\\vector_field_spherical.txt");
    auto ret3 = std::system("..\\..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe ..\\..\\results\\vector_field_spherical.txt");
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