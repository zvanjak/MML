#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif


// Base
void Demo_BaseUtils();
void Demo_Intervals();
void Demo_Vector();
void Demo_VectorN();
void Demo_Matrix();
void Demo_Matrix_Other();
void Demo_MatrixNM();
void Demo_Polynom();
void Demo_Tensors();
void Demo_Geometry();
void Demo_Geometry_2D();
void Demo_Geometry_3D();

// Core
void Demo_CoordSystem();
void Demo_Covar_Contravar_transformations();
void Demo_CoordTransf();
void Demo_Function();
void Demo_Interpolated_Function();
void Example2_Derivation();
void Demo_Integration();
void Demo_LinearAlgEqSolvers();
void Demo_Serialization();
void Demo_Dirac_function();
void Demo_Fields();
void Demo_Metric_Tensors();
void Demo_Curves();
void Demo_Surfaces();

// Algorithms
void Demo_Diff_geometry();
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

// Note: Visualization examples are now in their own app: MML_VisualizationApp
// See src/visualization_examples/


int main(int, char**)
{
	// Demo_Intervals();

	// Demo_Algebra();
	// Demo_BaseUtils();
	// Demo_Vector();
	// Demo_VectorN();
	// Demo_Matrix();
	// Demo_MatrixNM();
	// Demo_Matrix_Other();
	// Demo_Tensors();
	// Demo_Polynom();
	// Demo_Geometry();
	// Demo_Geometry_2D();
	// Demo_Geometry_3D();  

	// Demo_CoordTransf();
	
	// Demo_CoordSystem();
	
	// Demo_Covar_Contravar_transformations();
	// Demo_Dirac_function();
	// Demo_Fields();
	// Demo_Function_Space();
	// Demo_Function();
	// Demo_Serialization();

	// Demo_LinearFunctionals();
	// Demo_LinearOperators();
  
	// Demo_Interpolated_Function();
	// Demo_Interpolators();
	// Demo_Metric_Tensors();

	// Demo_QuadraticForms();
	// Demo_Curves();
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
	// Demo_Surface_integration();
	// Demo_Volume_integration();

	// Test_Speed_Functions();
	// Test_Speed_Derivation();
	// Test_Speed_Linear_alg_eq_solvers();

	// Test_Precision_Derivation();
	// Test_Precision_Integration();

	// Note: To run visualization examples, use MML_VisualizationApp instead
	// Demo_Visualization_Examples();  // Moved to separate app
}