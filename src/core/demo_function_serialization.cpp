#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "utilities/StdFunctions.h"

#include "core/VectorN.h"
#include "core/Function.h"

#include "basic_types/Curves.h"
#include "basic_types/Surfaces.h"
#endif

using namespace MML;


void Demo_Function_serialization()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                    FUNCTION SERIALIZATION                     ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    RealFunction f2{[](double x) { return sin(x); } };

    f2.SerializeVariableSpaced(0.0, 10.0, 100, "test.txt");
    auto ret = std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe");
}