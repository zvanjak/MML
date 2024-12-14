#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/CurvesSurfaces.h"

#include "algorithms/PathIntegration.h"
#endif

using namespace std;
using namespace MML;

void Calc_curve_length()
{
	ParametricCurve<3> circle([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), 0}; });
	ParametricCurve<3> helix([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), t}; });

	Real len = PathIntegration::ParametricCurveLength(circle, 0, 2 * 3.14);

	std::cout << "Length of circle is " << PathIntegration::ParametricCurveLength(circle, 0, 2 * 3.14) << endl;
	std::cout << "Length of helix is " << PathIntegration::ParametricCurveLength(helix, 0, 2 * 3.14) << endl;
}

void Calc_line_integral()
{
	ScalarFunction<3>  pot_simple([](const VectorN<Real, 3>& x) { return x[0] + x[1]; });
	ParametricCurve<3> halfcircle([](Real t) { return VectorN<Real, 3>{3*cos(t), 3*sin(t), 0}; });
		
	std::cout << "Line integral is : " << PathIntegration::LineIntegral(pot_simple, halfcircle, 0, Constants::PI/2, 1e-03) << endl;
	std::cout << "Exact values is  : " << 18.0000 << endl;

	// simple inverse radial potential (ie. conservative field)
	ScalarFunction<3>  potential_field([](const VectorN<Real, 3>& x) { return Real{ -10.0 } / x.NormL2(); });
	VectorFunction<3>  force_field_exact([](const VectorN<Real, 3>& x) { return x * (10.0 / (x.NormL2() * x.NormL2() * x.NormL2())); });
	VectorFunction<3>  force_field_calc([](const VectorN<Real, 3>& x) 
		{ 
			ScalarFunction<3>  potential_field([](const VectorN<Real, 3>& x) { return Real{ -10.0 } / x.NormL2(); });
			return ScalarFieldOperations::GradientCart(potential_field, x);
		}
	);

	ParametricCurve<3> circle([](Real t) { return VectorN<Real, 3>{cos(t), sin(t), 1}; });

	//std::cout << "Work integral for circle is : " << PathIntegration::WorkIntegral(potential, circle, 0, 2 * Constants::PI, 1e-03) << endl;
	//for (auto phi = 0.1; phi < 2 * 3.14159; phi += 0.25)
	//{
	//	std::cout << "Work integral for phi : " << phi << " is : " << PathIntegration::WorkIntegral(potential, circle, 0, phi, 1e-03) << endl;
	//}

	Point3Cartesian pnt1{3,-5, 2};
	Point3Cartesian pnt2{5, 5, 2};
	Point3Cartesian pntmid((pnt1 + pnt2)/2.0);
	Vector3Cartesian vec1{pnt1.X(), pnt1.Y(), pnt1.Z() };
	Vector3Cartesian vec2{pnt2.X(), pnt2.Y(), pnt2.Z() };
	double t1 = 0.0;
	double t2 = 1.0;

	Curves3D::LineCurve  line(t1, pnt1, t2, pnt2);
	int numPnt = 10;
	double dt = (t2 - t1) / numPnt;
	double sum = 0.0;
	for(int i=0; i<numPnt; i++ )
	{
		auto pnt = line(t1 + i * (t2 - t1) / numPnt);
		auto pntNext = line(t1 + (i+1) * (t2 - t1) / numPnt);
		auto dr = pntNext - pnt;
		double ds = dr.NormL2();
		double part_sum = potential_field(pnt) * ds;
		sum += part_sum;

		auto t = t1 + i * (t2 - t1) / numPnt;
		std::cout << "T = " << t << "  ";
		std::cout << "Point on line is : " << pnt << "  ";
		std::cout << "Potential at point is : " << potential_field(pnt) << "  ";
		std::cout << "Part. sum = " << part_sum << endl;

		auto p1 = line(t1);
		auto p2 = line(t2);
	}
	std::cout << "Total sum = " << sum << endl;
	std::cout << "Line integral between points is : " << PathIntegration::LineIntegral(potential_field, line, t1, t2, 1e-03) << endl;
	std::cout << "Work integral between points is : " << PathIntegration::LineIntegral(force_field_exact, line, t1, t2, 1e-03) << endl;
	std::cout << "Work integral 2 between points is : " << PathIntegration::LineIntegral(force_field_calc, line, t1, t2, 1e-03) << endl;

	Curves3D::LineCurve  line1(t1, pnt1, t2, pntmid);
	Curves3D::LineCurve  line2(t1, pntmid, t2, pnt2);

	std::cout << "Line 1 : " << PathIntegration::LineIntegral(potential_field, line1, t1, t2, 1e-03) << endl;
	std::cout << "Line 2 : " << PathIntegration::LineIntegral(potential_field, line2, t1, t2, 1e-03) << endl;
	
	double diff = potential_field(vec2) - potential_field(vec1);
	std::cout << "Difference in potential is : " << diff << endl;
}

void Demo_Path_Integration()
{
	std::cout << endl;
	std::cout << "***********************************************************************" << endl;
	std::cout << "****                     PATH INTEGRATION                          ****" << endl;
	std::cout << "***********************************************************************" << endl;

	Calc_curve_length();
	Calc_line_integral();
}