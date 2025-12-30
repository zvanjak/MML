#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Matrix.h"
#include "mml/base/VectorN.h"
#include "mml/base/Function.h"
#include "mml/base/InterpolatedFunction.h"

#include "mml/core/Derivation.h"
#include "mml/core/MetricTensor.h"

#include "mml/core/Integration/PathIntegration.h"

#include "mml/tools/Visualizer.h"

#include "mpl/SpecialRelativity/LorentzTransformation.h"
#include "mpl/SpecialRelativity/RelativisticMechanicsSimulator.h"
#endif

#include <iostream>
#include <stdexcept>

using namespace MML;
using namespace MPL;

// Helper class for calculating proper time along curves in Minkowski space-time
class HelperCurveProperTimeMinkowski : public IRealFunction
{
	const IParametricCurve<4>& _curve;
public:
	HelperCurveProperTimeMinkowski(const IParametricCurve<4>& curve) : _curve(curve){}

	Real operator()(Real t) const
	{
		auto tangent_vec = Derivation::DeriveCurve<4>(_curve, t, nullptr);

		MetricTensorMinkowski metricMinkowski;
		Tensor2<4> metricTensor = metricMinkowski(VectorN<Real, 4>{0.0, 0.0, 0.0, 0.0}); // get Minkowski metric tensor

		Real intVal = -metricTensor(tangent_vec, tangent_vec); // this is the Minkowski metric applied to the tangent vector

		if(intVal < 0.0)
			throw std::runtime_error("HelperCurveProperTime: negative value under square root, spacelike point/vector");

		return sqrt(intVal);
	}
};

// Investigating Minkowski spacetime:
// - calculating point properties (timelike/spacelike/lightlike)
// - calculating distance between two points in Minkowski space-time
// - calculating integral of proper time along curves in Minkowski space-time
void Demo1_Vector4Minkowski()
{
	// lets create a point in Minkowski space-time
	Vector4Minkowski point1{ 10.0, 5.0, 3.0, 0.0 }; // (t, x, y, z)

	// check if this point is timelike, spacelike or lightlike
	if (point1.isTimelike())
	{
		std::cout << "Point is timelike." << std::endl;
	}
	else if (point1.isSpacelike())
	{
		std::cout << "Point is spacelike." << std::endl;
	}
	else if (point1.isLightlike())
	{
		std::cout << "Point is lightlike." << std::endl;
	}

	Vector4Minkowski point2{ 20.0, 10.0, 6.0, 0.0 }; // (t, x, y, z)
	// calculate distance between two points in Minkowski space-time
	Real distance = point1.Distance(point2);
	std::cout << "Distance between points: " << distance << std::endl;


	// let's compare integral of proper time along two curves
	// first one is a straight line in Minkowski space-time
	ParametricCurve<4> line([](Real t) { return VectorN<Real, 4>{ t, Real(0.8) * t, Real(0.0), Real(0.0) }; });
	ParametricCurve<3> line3D([](Real t) { return VectorN<Real,3>{ Real(0.8) * t, Real(0.0), Real(0.0) }; } );


	// second one is a curve in Minkowski space-time, that we'll form from a 4D spline
	Vector<Vector4Minkowski> vec_curve_pnts{ {0.0, 0.0, 0.0, 0.0},
																	{20.0, 10.0, 0.0, 0.0}, 
																	{40.0, 26.0, 0.0, 0.0}, 
																	{60.0, 43.0, 0.0, 0.0}, 
																	{80.0, 61.0, 0.0, 0.0}, 
																	{100.0, 80.0, 0.0, 0.0} };

	// create Matrix from vector of points
	Matrix<Real> curve_points(vec_curve_pnts.size(), 4);
	for (int i = 0; i < vec_curve_pnts.size(); i++)
	{
		curve_points(i, 0) = vec_curve_pnts[i][0];
		curve_points(i, 1) = vec_curve_pnts[i][1];
		curve_points(i, 2) = vec_curve_pnts[i][2];			// z coordinate is actually time T
		curve_points(i, 3) = vec_curve_pnts[i][3];
	}

	SplineInterpParametricCurve<4> curve(0.0, 100.0, curve_points);

	// create Matrix from vector of points
	Matrix<Real> curve_points_viz(vec_curve_pnts.size(), 3);
	for (int i = 0; i < vec_curve_pnts.size(); i++)
	{
		curve_points_viz(i, 0) = vec_curve_pnts[i][1];
		curve_points_viz(i, 1) = vec_curve_pnts[i][2];
		curve_points_viz(i, 2) = vec_curve_pnts[i][0];			// z coordinate is actually time T
	}

	SplineInterpParametricCurve<3> curve3D(0.0, 100.0, curve_points_viz);

	// visualize curve3D
	//Visualizer::VisualizeParamCurve3D(curve3D, "Curve in Minkowski space-time",	0.0, 100.0, 100, "Chapter15_curve3D.txt");

	// let's visualize these two curves, but with z coordinate changed with T 
	Visualizer::VisualizeMultiParamCurve3D({ &line3D, &curve3D }, "Two worldlines",
																					0.0, 100.0, 100,
																					"Chapter15_two_worldlines");

	// let's calculate path integral for both these two curves, for proper time
	MetricTensorMinkowski metricMinkowski;
	Tensor2<4> metricTensor = metricMinkowski(VectorN<Real, 4>{0.0, 0.0, 0.0, 0.0}); // get Minkowski metric tensor

	for(int i = 0; i < 100; ++i)
	{
		Real t = i * 1.0;
		
		auto tangent_vec_line	= Derivation::DeriveCurve<4>(line, t, nullptr);
		auto tangent_vec_curve = Derivation::DeriveCurve<4>(curve, t, nullptr);
		
		Vec4Mink v1{ tangent_vec_line[0], tangent_vec_line[1], tangent_vec_line[2], tangent_vec_line[3] };

		Real intValLine = metricTensor(tangent_vec_line, tangent_vec_line); 
		Real intValCurve = metricTensor(tangent_vec_curve, tangent_vec_curve);
		
		std::cout << "Proper time along line at t=" << t << ": " << intValLine;
		std::cout << "   Proper time along curve at t=" << t << ": " << intValCurve << std::endl;
	}
	// now, let's calculate integral of proper time along the curve
	Real lineInt = IntegrateTrap(HelperCurveProperTimeMinkowski(line), 0.0, 100.0);
	Real curveInt = IntegrateTrap(HelperCurveProperTimeMinkowski(curve), 0.0, 100.0);

	std::cout << "Integral of proper time along line: " << lineInt << std::endl;
	std::cout << "Integral of proper time along curve: " << curveInt << std::endl;
}
