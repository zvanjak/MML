#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "core/Function.h"
#include "basic_types/Curves.h"
#include "basic_types/Geometry3D.h"
#endif

using namespace MML;

void Readme_vectors_matrices();
void Readme_functions();
void Readme_linear_system_solvers();
void Readme_ode_solvers();
void Readme_vector_field_operations();
void Readme_parametric_curves();
void Readme_interpolating_functions();
void Readme_visualizators();
void Readme_function_analyzer();

void Readme_coordinate_transformations()
{
}
void Readme_tensors()
{
}
void Readme_path_integration()
{
}
void Readme_geometry()
{
}

void Demo_Readme_Examples()
{
    Readme_vectors_matrices();
    Readme_linear_system_solvers();
    Readme_functions();
    Readme_ode_solvers();
    Readme_vector_field_operations();
    Readme_interpolating_functions();
    Readme_visualizators();
    Readme_function_analyzer();

    Readme_coordinate_transformations();
    Readme_tensors();
    Readme_path_integration();
    Readme_geometry();
}