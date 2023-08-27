#include "MMLBase.h"

void Demo_Algebra();

void Demo_Vector();
void Demo_VectorN();
void Demo_Matrix();
void Demo_MatrixNM();
void Demo_Tensors();
void Demo_Geometry();
void Demo_CoordTransf();
void Demo_CoordSystem();

void Demo_Polynom();
void Demo_Function();

void Demo_Derivation();
void Demo_Field_operations();
void Demo_Fields();

void Demo_Integration();
void Demo_Path_Integration();
void Demo_Interpolators();
void Demo_Interpolated_Function();

void Demo_LinearAlgEqSolvers();
void Demo_EigenSolvers();
void Demo_DiffEqSolvers();

void Demo_Diff_geometry();

void Test_Speed_Functions();
void Test_Speed_Derivation();
void Test_Speed_Linear_alg_eq_solvers();

void Test_Precision_Derivation();
void Test_Precision_Integration();

void Demo_Readme_Examples();

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
    
    // Demo_Algebra();

    // Demo_Vector();
    // Demo_VectorN();
    // Demo_Matrix();
    // Demo_MatrixNM();
    // Demo_Polynom();
    // Demo_Tensors();
    // Demo_Function();
    // Demo_Interpolated_Function();

    // Demo_Geometry();
    Demo_CoordSystem();
    Demo_CoordTransf();

    // Demo_Derivation();
    // Demo_Field_operations();
    // Demo_Fields();
    // Demo_Integration();
    // Demo_Path_Integration();
    // Demo_Interpolators();
    
    // Demo_LinearAlgEqSolvers();
    // Demo_EigenSolvers();
    // Demo_DiffEqSolvers();
    // Demo_Diff_geometry();

    // Test_Speed_Functions();
    // Test_Speed_Derivation();
    // Test_Speed_Linear_alg_eq_solvers();

    // Test_Precision_Derivation();
    // Test_Precision_Integration();

    // Demo_Readme_Examples();
}