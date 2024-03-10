#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/StdFunctions.h"

#include "base/VectorN.h"

#include "core/Function.h"
#include "core/FieldOperations.h"
#include "core/Surfaces.h"
#include "core/Serializer.h"
#include "core/Visualizer.h"

#endif

using namespace MML;

void Demo3_surface()
{

    ScalarFunction<2> testFunc{[](const VectorN<Real, 2> &x) { return (10 - std::abs(x[0])) * (10 - std::abs(x[1])) * sin(x[0]) * cos(x[1]); } };

    Visualizer::VisualizeScalarFunc2DCartesian(testFunc, "demo3_surface" , -10.0, 10.0, 50, -10.0, 10.0, 50, "demo3_surface.txt");
}

