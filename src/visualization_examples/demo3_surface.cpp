#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "utilities/StdFunctions.h"

#include "base/VectorN.h"
#include "core/Function.h"
#include "core/FieldOperations.h"

#include "core/Surfaces.h"

#endif

using namespace MML;

void Demo3_surface()
{

    ScalarFunction<2> testFunc{[](const VectorN<Real, 2> &x) { return (10 - std::abs(x[0])) * (10 - std::abs(x[1])) * sin(x[0]) * cos(x[1]); } };

    testFunc.Serialize2DCartesian(-10.0, 10.0, 50, -10.0, 10.0, 50, "..\\..\\results\\demo3_surface.txt");
    auto ret2 = std::system("..\\..\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe ..\\..\\results\\demo3_surface.txt");


}

