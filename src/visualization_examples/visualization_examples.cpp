#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#endif

void Demo1_Lorenz_multi_func();
void Demo1_Lorenz_parametric_curve();
void Demo3_surface();
void Demo4_vector_field_gravity();
void Demo5_vector_field_EM();

void Demo_Visualization_Examples()
{
    // Demo1_Lorenz_multi_func();
    // Demo1_Lorenz_parametric_curve();
    Demo3_surface();
    Demo4_vector_field_gravity();
    Demo5_vector_field_EM();
}