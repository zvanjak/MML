#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "base/Function.h"
#include "core/Surfaces.h"
#include "base/Geometry3D.h"
#endif

using namespace MML;

void Readme_demonstration();
void Readme_vectors();
void Readme_matrices();
void Readme_tensors();
void Readme_geometry();
void Readme_polynoms();
void Readme_functions();
void Readme_derivation();
void Readme_integration();
void Readme_vector_field_operations();
void Readme_curves_surfaces();
void Readme_coord_system_transf();
void Readme_linear_system_solvers();
void Readme_ode_solvers();
void Readme_function_analyzer();
void Readme_visualizators();


int main()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                        README EXAMPLES                        ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Readme_demonstration();
	// Readme_vectors();
	// Readme_matrices();
	// Readme_tensors();
	// Readme_geometry();
	// Readme_polynoms();
	// Readme_functions();
	// Readme_derivation();
	// Readme_integration();
	// Readme_vector_field_operations();
	// Readme_curves_surfaces();
	// Readme_coord_system_transf();

	// Readme_linear_system_solvers();
	// Readme_ode_solvers();
	// Readme_visualizators();
	// Readme_function_analyzer();

	// Readme_tensors();
	// Readme_geometry();
}