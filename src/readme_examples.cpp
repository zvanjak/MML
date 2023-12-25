#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "core/Function.h"
#include "basic_types/Curves.h"
#include "basic_types/Geometry3D.h"

#include "algorithms/ODESystemSolvers.h"
#include "algorithms/ODESystemSteppers.h"
#include "algorithms/ODESystemSteppers_Stiff.h"
#include "algorithms/DiffGeometryAlgorithms.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"
#include "../test_data/diff_eq_systems_test_bed.h"
#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

void Readme_vectors_matrices()
{
    Vector<double>  vec_dbl_3({ 1.0, 2.0, 3.0 }); 
    VectorComplex   vec_cmplx_2({ Complex(1,1), Complex(-1,2), Complex(2, -0.5) });

    Matrix<Real>   mat_3x3{ 3, 3, { 1.0, 2.0, -1.0, 
                                   -1.0, 5.0, 6.0, 
                                    3.0, 1.0, 1.0 }};  
    MatrixComplex  mat_cmplx(2,2, { Complex(1,1),  Complex(-1,2), 
                                    Complex(2, -0.5), Complex(1,1) });
    Matrix<Real>   unit_mat3 = MML::Matrix<Real>::GetUnitMatrix(3);

}

double Readme_functions_TestFunc(double x) 
{ 
    return sin(x)*(1.0 + 0.5*x*x); 
}
void Readme_functions()
{
    // creating a function object from a already existing (standalone) function
    RealFunction f1(Readme_functions_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
}

void Readme_parametric_curves()
{
    // creating curve directly with lambda
    ParametricCurve<3>        test_curve1( [](double t) -> VectorN<Real, 3> { return VectorN<Real, 3>{t, t*t, t*t*t}; } );
    
    // using predefined curve
    Curves::HelixCurve        helix(2.0, 2.0);
    
    // using curve from TestData
    const ParametricCurve<3> &test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve(0)._curve;

    double t = 0.5;
    auto tangent   = Vector3Cartesian( DiffGeometry::getTangent(test_curve, t) );
    auto unit_tang = Vector3Cartesian( DiffGeometry::getTangentUnit(test_curve, t) );
    auto normal    = Vector3Cartesian( DiffGeometry::getNormal(test_curve, t) );
    auto unit_norm = Vector3Cartesian( DiffGeometry::getNormalUnit(test_curve, t) );
    auto binormal  = VectorProd(unit_tang, unit_norm);
    
    auto curv_vec   = DiffGeometry::getCurvatureVector(test_curve, t);
    auto curvature  = DiffGeometry::getCurvature(test_curve, t);
    auto curvature3 = DiffGeometry::getCurvature3(test_curve, t);
}

void Readme_linear_system_solvers()
{
 	Matrix<Real>    mat{5, 5, { 1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7 } };
	Vector<Real> 	rhs{1.1, 4.7, 0.1, 9.3, 0.4};
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);

	Vector<Real>    res_rhs = mat * vecSol;
}

void Readme_ode_systems()
{
    auto sys0 = TestBeds::ODESystemTestBed::getODESystem(1);

    const double atol=1.0e-3, rtol=atol;
    const double h1=0.01, hmin=0.0;
    const double x1=0.0, x2=2.0;
    
    Vector<Real> ystart01{2.0, 0.0};
    Output out(20);            

    ODESystemSolver<StepperDopr853> ode_solver01(sys0, atol, rtol, out);
    ODESystemSolution               sol01 = ode_solver01.integrate(ystart01, 0.0, 2.0, h1, hmin);

    std::cout << "x values:\n";  sol01._xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - "; sol01._yval.Print(std::cout, 6, 3);
}

void Readme_coordinate_transformations()
{

}

void Readme_tensors()
{

}

void Readme_field_operations()
{

}

void Readme_path_integration()
{

}

void Readme_geometry()
{

}

void Readme_()
{

}

void Demo_Readme_Examples()
{
    Readme_vectors_matrices();
    Readme_functions();
    Readme_parametric_curves();
    Readme_linear_system_solvers();
    Readme_ode_systems();
    Readme_coordinate_transformations();
    Readme_tensors();
    Readme_field_operations();
    Readme_path_integration();
    Readme_geometry();
}