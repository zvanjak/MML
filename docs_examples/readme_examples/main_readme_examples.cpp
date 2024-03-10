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

void Readme_vectors_matrices();
void Readme_linear_system_solvers();
void Readme_defining_functions();
void Readme_working_with_functions();
void Readme_ode_solvers();
void Readme_vector_field_operations();
void Readme_parametric_curves();
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

void Main_Readme_Examples()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                        README EXAMPLES                        ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Readme_vectors_matrices();
	Readme_linear_system_solvers();
	Readme_defining_functions();
	Readme_working_with_functions();
	Readme_ode_solvers();
	// Readme_vector_field_operations();
	// Readme_visualizators();
	Readme_function_analyzer();

	//Readme_coordinate_transformations();
	//Readme_tensors();
	//Readme_path_integration();
	//Readme_geometry();
}