#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/Curves.h"

#include "algorithms/ParametricCurveAnalyzer.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

void Docs_Demo_Diff_geometry_curves()
{
  std::cout << "***********************************************************" << std::endl;
  std::cout << "*****               Curves diff. geometry           *******" << std::endl;

  // using Helix curve from TestData
  const ParametricCurve<3> &test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
     
  std::cout << "t                          tangent                             unit_tangent                            normal                             unit_normal                              binormal" << std::endl;
  for( int i=-10; i<=10; ++i)
  {
    Real t = i / 10.0;

    auto tangent   = ParametricCurveAnalyzer::getTangent(test_curve, t);
    auto unit_tang = ParametricCurveAnalyzer::getTangentUnit(test_curve, t);
    auto normal    = ParametricCurveAnalyzer::getNormal(test_curve, t);
    auto unit_norm = ParametricCurveAnalyzer::getNormalUnit(test_curve, t);
    auto binormal  = ParametricCurveAnalyzer::getBinormal(test_curve, t);
    //auto binormal = VectorProd(Vector3Cartesian(unit_tang), Vector3Cartesian(unit_norm));
    //auto princNorm = ParametricCurveAnalyzer::getPrincipalNormal(test_curve, t);

    std::cout << "t = " << std::setw(4) << t << " : ";

    tangent.Print(std::cout, 10, 6, 1e-10); std::cout << "  ";
    unit_tang.Print(std::cout, 10, 6, 1e-10); std::cout << "  ";
    normal.Print(std::cout, 10, 6, 1e-10); std::cout << "  ";
    unit_norm.Print(std::cout, 10, 6, 1e-10); std::cout << "  ";
    binormal.Print(std::cout, 10, 6, 1e-10);

    std::cout << std::endl;

    auto curv_vec   = ParametricCurveAnalyzer::getCurvatureVector(test_curve, t);
    auto curvature  = ParametricCurveAnalyzer::getCurvature(test_curve, t);
    auto curvature3 = ParametricCurveAnalyzer::getCurvature3(test_curve, t);
  }    

/* OUTPUT
t                          tangent                             unit_tangent                            normal                             unit_normal                              binormal
t =   -1 : [  0.841471,   0.540302,          1]  [   0.59501,   0.382051,   0.707107]  [ -0.540302,   0.841471,          0]  [ -0.540302,   0.841471,          0]  [   0.59501,   0.382051,  -0.707107]
t = -0.9 : [  0.783327,    0.62161,          1]  [  0.553896,   0.439545,   0.707107]  [  -0.62161,   0.783327,          0]  [  -0.62161,   0.783327,          0]  [  0.553896,   0.439545,  -0.707107]
t = -0.8 : [  0.717356,   0.696707,          1]  [  0.507247,   0.492646,   0.707107]  [ -0.696707,   0.717356,          0]  [ -0.696707,   0.717356,          0]  [  0.507247,   0.492646,  -0.707107]
t = -0.7 : [  0.644218,   0.764842,          1]  [  0.455531,   0.540825,   0.707107]  [ -0.764842,   0.644218,          0]  [ -0.764842,   0.644218,          0]  [  0.455531,   0.540825,  -0.707107]
t = -0.6 : [  0.564642,   0.825336,          1]  [  0.399263,     0.5836,   0.707107]  [ -0.825336,   0.564642,          0]  [ -0.825336,   0.564642,          0]  [  0.399263,     0.5836,  -0.707107]
t = -0.5 : [  0.479426,   0.877583,          1]  [  0.339005,   0.620545,   0.707107]  [ -0.877583,   0.479426,          0]  [ -0.877583,   0.479426,          0]  [  0.339005,   0.620545,  -0.707107]
t = -0.4 : [  0.389418,   0.921061,          1]  [   0.27536,   0.651288,   0.707107]  [ -0.921061,   0.389418,          0]  [ -0.921061,   0.389418,          0]  [   0.27536,   0.651288,  -0.707107]
t = -0.3 : [   0.29552,   0.955336,          1]  [  0.208964,   0.675525,   0.707107]  [ -0.955336,    0.29552,          0]  [ -0.955336,    0.29552,          0]  [  0.208964,   0.675525,  -0.707107]
t = -0.2 : [  0.198669,   0.980067,          1]  [   0.14048,   0.693012,   0.707107]  [ -0.980067,   0.198669,          0]  [ -0.980067,   0.198669,          0]  [   0.14048,   0.693012,  -0.707107]
t = -0.1 : [ 0.0998334,   0.995004,          1]  [ 0.0705929,   0.703574,   0.707107]  [ -0.995004,  0.0998334,          0]  [ -0.995004,  0.0998334,          0]  [ 0.0705929,   0.703574,  -0.707107]
t =    0 : [         0,          1,          1]  [         0,   0.707107,   0.707107]  [        -1,          0,          0]  [        -1,          0,          0]  [         0,   0.707107,  -0.707107]
t =  0.1 : [-0.0998334,   0.995004,          1]  [-0.0705929,   0.703574,   0.707107]  [ -0.995004, -0.0998334,          0]  [ -0.995004, -0.0998334,          0]  [-0.0705929,   0.703574,  -0.707107]
t =  0.2 : [ -0.198669,   0.980067,          1]  [  -0.14048,   0.693012,   0.707107]  [ -0.980067,  -0.198669,          0]  [ -0.980067,  -0.198669,          0]  [  -0.14048,   0.693012,  -0.707107]
t =  0.3 : [  -0.29552,   0.955336,          1]  [ -0.208964,   0.675525,   0.707107]  [ -0.955336,   -0.29552,          0]  [ -0.955336,   -0.29552,          0]  [ -0.208964,   0.675525,  -0.707107]
t =  0.4 : [ -0.389418,   0.921061,          1]  [  -0.27536,   0.651288,   0.707107]  [ -0.921061,  -0.389418,          0]  [ -0.921061,  -0.389418,          0]  [  -0.27536,   0.651288,  -0.707107]
t =  0.5 : [ -0.479426,   0.877583,          1]  [ -0.339005,   0.620545,   0.707107]  [ -0.877583,  -0.479426,          0]  [ -0.877583,  -0.479426,          0]  [ -0.339005,   0.620545,  -0.707107]
t =  0.6 : [ -0.564642,   0.825336,          1]  [ -0.399263,     0.5836,   0.707107]  [ -0.825336,  -0.564642,          0]  [ -0.825336,  -0.564642,          0]  [ -0.399263,     0.5836,  -0.707107]
t =  0.7 : [ -0.644218,   0.764842,          1]  [ -0.455531,   0.540825,   0.707107]  [ -0.764842,  -0.644218,          0]  [ -0.764842,  -0.644218,          0]  [ -0.455531,   0.540825,  -0.707107]
t =  0.8 : [ -0.717356,   0.696707,          1]  [ -0.507247,   0.492646,   0.707107]  [ -0.696707,  -0.717356,          0]  [ -0.696707,  -0.717356,          0]  [ -0.507247,   0.492646,  -0.707107]
t =  0.9 : [ -0.783327,    0.62161,          1]  [ -0.553896,   0.439545,   0.707107]  [  -0.62161,  -0.783327,          0]  [  -0.62161,  -0.783327,          0]  [ -0.553896,   0.439545,  -0.707107]
t =    1 : [ -0.841471,   0.540302,          1]  [  -0.59501,   0.382051,   0.707107]  [ -0.540302,  -0.841471,          0]  [ -0.540302,  -0.841471,          0]  [  -0.59501,   0.382051,  -0.707107]
*/
}

void Docs_Demo_Diff_geometry_surfaces()
{
}

void Docs_Demo_Diff_geometry()
{
	Docs_Demo_Diff_geometry_curves();
	Docs_Demo_Diff_geometry_surfaces();
}
