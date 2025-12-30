#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Curves.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"
#include "../test_data/parametric_surfaces_test_bed.h"

using namespace MML;

void Demo_curves()
{
	std::cout << std::endl;
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                    DIFFERENTIAL GEOMETRY                      ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Curves::CurveCartesian3D test_curve1([](Real t) -> VectorN<Real, 3> { return VectorN<Real, 3>{Real(t), Real(t* t), Real(t* t* t)}; });
	Curves::HelixCurve      helix(Real(2.0), Real(2.0));
	const Curves::CurveCartesian3D& test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve(0)._curve;

	std::cout << "          Tangent                   Tangent unit                   Normal                  Normal unit                    Binormal                   Curv.vec.                Curv.vec.norm.        Curvature\n";

	std::cout << std::fixed;

	for (double t = 0.0; t < 2 * Constants::PI; t += 0.4)
	{
		auto tangent = Vector3Cartesian(test_curve1.getTangent(t));
		auto unit_tang = Vector3Cartesian(test_curve1.getTangentUnit(t));
		auto normal = Vector3Cartesian(test_curve1.getNormal(t));
		auto unit_norm = Vector3Cartesian(test_curve1.getNormalUnit(t));
		auto binormal = test_curve1.getBinormal(t);

		auto curv_vec = test_curve1.getCurvatureVector(t);
		auto curvature = test_curve1.getCurvature(t);

		tangent.Print(std::cout, 7, 3); std::cout << " ";
		unit_tang.Print(std::cout, 7, 3); std::cout << " ";
		normal.Print(std::cout, 7, 3); std::cout << " ";
		unit_norm.Print(std::cout, 7, 3); std::cout << " ";
		binormal.Print(std::cout, 7, 3); std::cout << " ";

		curv_vec.Print(std::cout, 7, 3); std::cout << " ";
		auto curv_vec_norm = curv_vec / curv_vec.NormL2();
		curv_vec_norm.Print(std::cout, 7, 3);
		std::cout << "   " << curvature << "   " << std::endl;
	}

	// TODO: isArcLengthParametrized method not available on IParametricCurve<3> interface
	// bool b = test_curve.isArcLengthParametrized(0.0, 2 * Constants::PI);
	// std::cout << "Is arc length parametrized : " << b << std::endl;
}

void Demo_surfaces()
{
}

void Demo_Diff_geometry()
{
	Demo_curves();
	Demo_surfaces();
}