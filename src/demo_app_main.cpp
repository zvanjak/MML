#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif

// Utilities
void Demo_Intervals();

// Core
void Demo_Algebra();
void Demo_CoordTransf();
void Demo_CoreUtils();
void Demo_Derivation();
void Demo_Function();
void Demo_Geometry();
void Demo_Interpolated_Function();
void Demo_Integration();
void Demo_LinearAlgEqSolvers();
void Demo_LinearFunctionals();
void Demo_LinearOperators();
void Demo_Matrix();
void Demo_Matrix_Other();
void Demo_MatrixNM();
void Demo_Polynom();
void Demo_QuadraticForms();
void Demo_Serialization();
void Demo_Tensors();
void Demo_VectorSpaces();
void Demo_Vector();
void Demo_VectorN();

// Basic types
void Demo_CoordSystem();
void Demo_Covar_Contravar_transformations();
void Demo_Dirac_function();
void Demo_Fields();
void Demo_Function_Space();
void Demo_Geometry_2D();
void Demo_Geometry_3D();
void Demo_Interpolators();
void Demo_Metric_Tensors();
void Demo_Surfaces();

// Algorithms
void Demo_Diff_geometry();
void Demo_EigenSolvers();
void Demo_Field_operations();
void Demo_ODESystemSolvers();
void Demo_Path_Integration();
void Demo_Function_analyzer();
void Demo_Statistics();
void Demo_Root_finding();
void Demo_Surface_integration();
void Demo_Volume_integration();

// Speed
void Test_Speed_Functions();
void Test_Speed_Derivation();
void Test_Speed_Linear_alg_eq_solvers();

// Precision
void Test_Precision_Derivation();
void Test_Precision_Integration();

void Demo_Visualization_Examples();


int main(int, char**)
{
	// double max = std::numeric_limits<double>::max();
	// double min = -std::numeric_limits<double>::max();
	// double inf = std::numeric_limits<double>::infinity();
	// double min_inf = -std::numeric_limits<double>::infinity();

	// if (inf > max)
	//     std::cout << inf << " is greater than " << max << '\n';

	// if (min_inf < min)
	//     std::cout << min_inf << " is less than " << min << '\n';

	// Demo_Intervals();

	// Demo_Algebra();
	// Demo_CoreUtils();
	// Demo_Matrix();
	// Demo_Matrix_Other();
	// Demo_MatrixNM();
	// Demo_Tensors();
	// Demo_Vector();
	// Demo_VectorN();

	// Demo_CoordTransf();
	Demo_CoordSystem();
	// Demo_Covar_Contravar_transformations();
	// Demo_Dirac_function();
	// Demo_Fields();
	// Demo_Function_Space();
	// Demo_Function();
	// Demo_Serialization();

	// Demo_LinearFunctionals();
	// Demo_LinearOperators();
	// Demo_VectorSpaces();
	// Demo_Geometry();
	// void Demo_Geometry_2D();
	// void Demo_Geometry_3D();    
	// Demo_Interpolated_Function();
	// Demo_Interpolators();
	// Demo_Metric_Tensors();
	// Demo_Polynom();
	// Demo_QuadraticForms();
	// Demo_Surfaces();

	// Demo_Derivation();
	// Demo_Integration();
	// Demo_LinearAlgEqSolvers();

	// Demo_Diff_geometry();
	// Demo_EigenSolvers();
	// Demo_Field_operations();
	// Demo_ODESystemSolvers();
	// Demo_Path_Integration();
	// Demo_Function_analyzer();
	// Demo_Statistics();
	// Demo_Root_finding();
	Demo_Surface_integration();
	// Demo_Volume_integration();

	// Test_Speed_Functions();
	// Test_Speed_Derivation();
	// Test_Speed_Linear_alg_eq_solvers();

	// Test_Precision_Derivation();
	// Test_Precision_Integration();

	// Demo_Visualization_Examples();
}