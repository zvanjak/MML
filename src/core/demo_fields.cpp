#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#include "core/Function.h"
#include "core/FieldOperations.h"

#include "core/Fields.h"
#endif

using namespace MML;

Real MyPotentialCart(const VectorN<Real, 3> &x )   { return x[0]*x[1]; }

VectorN<Real, 3> MyPotentialForceField(const VectorN<Real, 3> &x)
{
    ScalarFunction<3> gravPot([](const VectorN<Real, 3> &x) { return MyPotentialCart(x); });
    return -1000.0 * ScalarFieldOperations::GradientCart(gravPot,x);
}

void Demo_Fields()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          FIELDS                               ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // defining fields directly from lambda
    ScalarFunction<3> scalarField([](const VectorN<Real, 3> &x) { return x[0]; });
    VectorFunction<3> vectorField([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });

    // using user defined function (function pointer)
    VectorFunction<3> user_fun(MyPotentialForceField);

    // predefined scalar fields
    Fields::InverseRadialFieldCart  inverseRadialPotentialField(10);
    Fields::InverseRadialFieldSpher inverseRadialPotentialFieldSpher(10);

    // predefined vector fields 
    Fields::InverseRadialForceFieldCart  inverseRadialPotentialForceField(10);
    Fields::InverseRadialForceFieldSpher inverseRadialPotentialForceFieldSpher(10);

    // no need for this, but you can do it
    Fields::InverseRadialForceFieldCart test_fun(1000.0);
    VectorFunction<3> test_fun1(Fields::InverseRadialPotentialForceFieldCart);
}