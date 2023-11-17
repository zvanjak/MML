#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/VectorN.h"
#include "core/Function.h"

#include "basic_types/Fields.h"

#include "core/FieldOperations.h"
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

/*
    - TODO - definirati različite potencijalne funkcije (ScalarFunction)
        - gravitacija
    - definirati vektorske funcije polja
        - polje EM naboja u pokretu
        - polje vodiča kroz koji teče struja
*/
    
    ScalarFunction<3> funcScalar([](const VectorN<Real, 3> &x) { return x[0]; });
    VectorFunction<3> funcVector([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });

    // potencijali - ScalarFunction
    InverseRadialFieldCart inverseRadialPotentialField(10);
    InverseRadialFieldSpher inverseRadialPotentialFieldSpher(10);
    InverseRadialFieldCyl inverseRadialPotentialFieldCyl(10);
    
    GravityPotentialFieldCart gravityPotentialFieldCart(100);
    GravityPotentialFieldSpher gravityPotentialFieldSpher(100);
    GravityPotentialFieldCyl gravityPotentialFieldCyl(100);

    // force fields - VectorFunction
    GravityForceFieldCart gravityCart(100);
    GravityForceFieldSpher gravitySpher(100);
    GravityForceFieldCyl gravityCyl(100);

    // no need for this, but you can do it
    InverseRadialForceFieldCart test_fun(1000.0);
    VectorFunction<3> test_fun1(InverseRadialPotentialForceFieldCart);

    // using user defined function
    VectorFunction<3> user_fun(MyPotentialForceField);

    test_fun.SerializeCartesian(-30.0, 30.0, 4, -30.0, 30.0, 4, -30.0, 30.0, 4, "vector_field.txt");
    auto ret = std::system("..\\..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe vector_field.txt");
    std::cout << "Return code = " << ret << std::endl;
}