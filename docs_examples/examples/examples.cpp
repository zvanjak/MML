#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/Function.h"
#include "core/Curves.h"
#include "base/Geometry3D.h"
#endif

using namespace MML;

void Main_example();

void Example1_kosi_hitac();
void Example2_collision_calculator();
void Example3_tensor_of_inertia();
void Example4_Gravity_field_visualization();
void Example5_Voyager_travels();
void Example6_electric_charge_distribution();


void Docs_Examples()
{
    //Main_example();

    //Example1_kosi_hitac();
    //Example2_collision_calculator();
    Example3_tensor_of_inertia();
    Example4_Gravity_field_visualization();
    //Example5_Voyager_travels();
    //Example6_electric_charge_distribution();
}