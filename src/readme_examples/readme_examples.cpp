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

void Readme_vectors_matrices();
void Readme_functions();
void Readme_derivation();
void Readme_integration();
void Readme_linear_system_solvers();
void Readme_ode_solvers();
void Readme_vector_field_operations();

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
    Readme_derivation();
    Readme_integration();
    Readme_ode_solvers();
    Readme_vector_field_operations();

    Readme_parametric_curves();
    Readme_coordinate_transformations();
    Readme_tensors();
    Readme_path_integration();
    Readme_geometry();
}