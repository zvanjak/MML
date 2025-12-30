#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Matrix.h"
#include "base/VectorN.h"
#include "base/Function.h"
#include "base/InterpolatedFunction.h"

#include "core/Derivation.h"
#include "core/MetricTensor.h"

#include "core/Integration/PathIntegration.h"

#include "tools/Visualizer.h"

#include "mpl/SpecialRelativity/LorentzTransformation.h"
#include "mpl/SpecialRelativity/RelativisticMechanicsSimulator.h"
#endif


using namespace MML;
using namespace MPL;

// modeling point trajectory in 4D space-time, with everything happening in x-y plane
// effectively supressing z-coordinate and showing t as 3rd coordinate

// visualization of trajectories
// 1. observer at rest, with clocks synchronized in the rest frame
// 2. passanger moving with speed v along x-axis, with clocks synchronized in the moving frame

// kako parametrizirati tu krivulju?

// TWO perspectives
// one is based on a set of clocks throughout the space, which are synchronized in the rest frame
// and afterwards we analyze how these clocks see events in the moving frame
// second is based on a passanger moving with speed v along x-axis, and observer at rest

// what sees a "global" observer, who is at rest, and has clocks synchronized in the rest frame


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

// investigating Minkowski spacetime
// calculating points properties
// calculating distance between two points in Minkowski space-time
// calculating integral of proper time along a curve in Minkowski space-time
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
	// frist one is a straight line in Minkowski space-time
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
	//Visualizer::VisualizeParamCurve3D(curve3D, "Curve in Minkowski space-time",	0.0, 100.0, 100, "example15_curve3D.txt");

	// let's visualize these two curves, but with z coordinate changed with T 
	Visualizer::VisualizeMultiParamCurve3D({ &line3D, &curve3D }, "Two wordlines",
																					0.0, 100.0, 100,
																					"example15_two_worldlines");

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

// check that Lorentz transformation is equal to Galilean transformation for low velocities
void Example15_Lorentz_as_Gallilean_transf()
{
	Real speed = 0.8; // speed in units of c, so 0.8 means 80% of the speed of light

	// moving frame actually!
	CoordTransfLorentzXAxis ctLorentzX(speed);

	// we have event happening at T=10 s at the origin of rest (observer) frame
	Vector4Minkowski eventLocalRestFrame1{ 10.0, 0.0, 0.0, 0.0 };

	// how passanger at origin of moving frame sees this event?
	Vector4Minkowski eventMovingFrame1 = ctLorentzX.transf(eventLocalRestFrame1);

	std::cout << "Event in the local rest frame of reference: " << eventLocalRestFrame1 << std::endl;
	std::cout << "Event in the moving frame of reference    : " << eventMovingFrame1 << std::endl;

	// now, inverse transformation
	// we have an event at T'=10 s at the origin of moving frame
	Vector4Minkowski eventMovingFrame2{ 10.0, 0.0, 0.0, 0.0 };

	// how observer at rest sees this event?
	Vector4Minkowski eventLocalRestFrame2 = ctLorentzX.transfInverse(eventMovingFrame2);

	std::cout << "Event in the moving frame of reference    : " << eventMovingFrame2 << std::endl;
	std::cout << "Event in the moving frame of reference    : " << eventLocalRestFrame2 << std::endl;

	// event at x = 8.0 m, T = 10.0 s in the local rest frame
	// at that time, our passanger origin passes that point and we are interested in how observer at rest sees this event
	// so, we have event at T=10.0 s, x=8.0 m, y=0.0 m, z=0.0 m in the local rest frame
	Vector4Minkowski eventLocalRestFrame3{ Real(10.0), Real(10.0) * speed, Real(0.0), Real(0.0) };

	// how passanger at origin of moving frame sees this event?
	Vector4Minkowski eventMovingFrame3 = ctLorentzX.transf(eventLocalRestFrame3);
	std::cout << "Event in the local rest frame of reference: " << eventLocalRestFrame3 << std::endl;
	std::cout << "Event in the moving frame of reference    : " << eventMovingFrame3 << std::endl;


	// for 10 point along x-axis, following movement of the passanger
	std::cout << "Moving to point B with speed 0.8c:" << std::endl;
	for (int i = 0; i <= 10; ++i)
	{
		// event at T=10.0 s, x=i m, y=0.0 m, z=0.0 m in the local rest frame
		Vector4Minkowski xPos{ Real(i * 1.0), Real(i * speed), Real(0.0), Real(0.0) };

		// how passanger at origin of moving frame sees this event?
		Vector4Minkowski eventMovingFrame = ctLorentzX.transf(xPos);
		std::cout << "Local rest frame point    : " << xPos.T() << " " << xPos.X() << "  ";
		std::cout << "Moving frame of reference : " << eventMovingFrame.T() << " " << eventMovingFrame.X() << std::endl;
	}

	// now, turning back to the observer at rest
	CoordTransfLorentzXAxis ctLorentzXBack(-speed);
	std::cout << "Moving back to point A:" << std::endl;
	for (int i = 10; i <= 20; ++i)
	{
		// event at T=10.0 s, x=i m, y=0.0 m, z=0.0 m in the observer's frame
		Vector4Minkowski xPos{ Real(i * 1.0), Real(10 * speed - (i - 10) * speed), Real(0.0), Real(0.0) };

		// how passanger at origin of moving frame sees this event?
		Vector4Minkowski eventMovingFrame = ctLorentzXBack.transf(xPos);
		std::cout << "Local rest frame point   : " << xPos.T() << " " << xPos.X() << "  ";
		std::cout << "Moving frame of reference: " << eventMovingFrame.T() << " " << eventMovingFrame.X() << std::endl;
	}
}

// demo with passanger moving along x-axis with speed v, and observer at rest
// calculate Lorentz transformation of coordinates and time for this observer
void Example15_Moving_x_axis()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "***  EXAMPLE 15 - Moving frame in x-axis drection transformation   ***" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Real speed = 0.9; // speed in units of c, so 0.8 means 80% of the speed of light

	std::cout << std::fixed << std::setprecision(2);

	// moving frame actually!
	CoordTransfLorentzXAxis ctLorentzX(speed);		

	// we have event happening at T=10 s at the origin of rest (observer) frame
	Vector4Minkowski eventLocalRestFrame1{ Real(10.0), Real(0.0), Real(0.0), Real(0.0) };	

	// how passanger at origin of moving frame sees this event?
	Vector4Minkowski eventMovingFrame1 = ctLorentzX.transf(eventLocalRestFrame1);

	std::cout << "Event in the local rest frame of reference: " << eventLocalRestFrame1 << std::endl;
	std::cout << "Event in the moving frame of reference    : " << eventMovingFrame1 << std::endl;

	// now, inverse transformation
	// we have an event at T'=10 s at the origin of moving frame
	Vector4Minkowski eventMovingFrame2{ Real(10.0), Real(0.0), Real(0.0), Real(0.0) };	

	// how observer at rest sees this event?
	Vector4Minkowski eventLocalRestFrame2 = ctLorentzX.transfInverse(eventMovingFrame2);
	
	std::cout << "Event in the moving frame of reference    : " << eventMovingFrame2 << std::endl;
	std::cout << "Event in the moving frame of reference    : " << eventLocalRestFrame2 << std::endl;

	// event at x = 8.0 m, T = 10.0 s in the local rest frame
	// at that time, our passanger origin passes that point and we are interested in how observer at rest sees this event
	// so, we have event at T=10.0 s, x=8.0 m, y=0.0 m, z=0.0 m in the local rest frame
	Vector4Minkowski eventLocalRestFrame3{ Real(10.0), Real(10.0) * speed, Real(0.0), Real(0.0) };

	// how passanger at origin of moving frame sees this event?
	Vector4Minkowski eventMovingFrame3 = ctLorentzX.transf(eventLocalRestFrame3);
	std::cout << "Event in the local rest frame of reference: " << eventLocalRestFrame3 << std::endl;
	std::cout << "Event in the moving frame of reference    : " << eventMovingFrame3 << std::endl;


	std::cout << std::fixed << std::setprecision(2);

	// for 10 point along x-axis, following movement of the passanger
	std::cout << "Moving to point B with speed 0.8c:" << std::endl;
	for (int i = 0; i <= 10; ++i)
	{
		// event at T=10.0 s, x=i m, y=0.0 m, z=0.0 m in the local rest frame
		Vector4Minkowski xPos{ Real(i * 1.0), Real(i * speed), Real(0.0), Real(0.0) };
		
		// how passanger at origin of moving frame sees this event?
		Vector4Minkowski eventMovingFrame = ctLorentzX.transf(xPos);
		std::cout << "Local rest frame : " << "T = " << std::setw(5) << xPos.T() << "  X = " << xPos.X() << "  ";
		std::cout << "Moving frame of reference : " << "T = " << eventMovingFrame.T() << "  X = " << eventMovingFrame.X() << std::endl;
	}

	// now, turning back to the observer at rest
	CoordTransfLorentzXAxis ctLorentzXBack(speed*2*10, -speed);
	
	std::cout << "Moving back to point A:" << std::endl;
	for (int i = 10; i <= 20; ++i)
	{
		// event at T=10.0 s, x=i m, y=0.0 m, z=0.0 m in the observer's frame
		Vector4Minkowski xPos{ Real(i * 1.0), Real(10 * speed - (i - 10) * speed), Real(0.0), Real(0.0) };
		
		// how passanger at origin of moving frame sees this event?
		Vector4Minkowski eventMovingFrame = ctLorentzXBack.transf(xPos);
		std::cout << "Local rest frame : " << "T = " << std::setw(5) << xPos.T() << "  X = " << xPos.X() << "  ";
		std::cout << "Moving frame of reference : " << "T = " << eventMovingFrame.T() << "  X = " << eventMovingFrame.X() << std::endl;
	}
}

// transforming between local frame at rest and inertial frame moving in general direction
// given with vector direction (but both system have the same orientation, so no rotation)
void Example15_Moving_frame_general_vector_direction()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "***  EXAMPLE 15 - Moving frame in general direction transformation  ***" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Real speed = 0.8;
	Vec3Cart direction{ Real(1.0), Real(0.0), Real(0.0) }; // moving along x-axis
	Vec3Cart velocity = direction * speed; // velocity vector in units of c
	
	Pnt3Cart pntA(0, 0, 0);
	Pnt3Cart pntB = pntA + direction * 10.0;	// point B is 10 seconds away from point A in the direction of movement

	CoordTransfLorentzGeneral ctGeneral(0.8, direction);
	
	std::cout << std::fixed << std::setprecision(2);

	// for 10 point along x-axis, following movement of the passanger
	std::cout << "Moving to point B with speed 0.8c:" << std::endl;
	// simulating ofr 10 time-units
	for (int i = 0; i <= 10; ++i)
	{
		Pnt3Cart posAtT = pntA + velocity * i; // position at T=i seconds
		Vector4Minkowski xPos{ Real(i * 1.0), Real(i * speed), Real(0.0), Real(0.0) };

		// how passanger at origin of moving frame sees this event?
		Vector4Minkowski eventMovingFrame = ctGeneral.transf(xPos);
		std::cout << "Local rest frame : " << "T = " << std::setw(5) << xPos.T() << "  X = " << xPos.X() << "  ";
		std::cout << "Moving frame of reference : " << "T = " << eventMovingFrame.T() << "  X = " << eventMovingFrame.X() << std::endl;

		// at the origin of local frame there is a clock, and on the passenger ship
		// there isa a telescope pointed at that clock at all times
		// when passenger ship passes midway point between A and B,
		// what time does he see on the clock at point A?
		if (i == 5)
		{
			// at T=5 seconds, the ship is at the point B
			Vector4Minkowski eventAtMidway = ctGeneral.transfInverse(Vector4Minkowski{ Real(5.0), Real(5.0) * speed, Real(0.0), Real(0.0) });
			std::cout << "At T=5 seconds, the ship sees clock at point A shows T = " << eventAtMidway.T() << " seconds." << std::endl;
		}
	}

	// now, turning back to the observer at rest
	CoordTransfLorentzGeneral ctLorentzXBack(0.8, -direction);
	std::cout << "Moving back to point A:" << std::endl;


	for (int i = 10; i <= 20; ++i)
	{
		// event at T=10.0 s, x=i m, y=0.0 m, z=0.0 m in the observer's frame
		Vector4Minkowski xPos{ Real(i * 1.0), Real(10 * speed - (i - 10) * speed), Real(0.0), Real(0.0) };

		// how passanger at origin of moving frame sees this event?
		Vector4Minkowski eventMovingFrame = ctLorentzXBack.transf(xPos);
		std::cout << "Local rest frame : " << "T = " << std::setw(5) << xPos.T() << "  X = " << xPos.X() << "  ";
		std::cout << "Moving frame of reference : " << "T = " << eventMovingFrame.T() << "  X = " << eventMovingFrame.X() << std::endl;
	}
}


void Example15_Spherical_ball()
{

}

void Example15_Lorentz_transformation()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 15 - Lorentz transformation         ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo1_Vector4Minkowski();
	//Example15_Lorentz_as_Gallilean_transf();
	//Example15_Moving_x_axis();
	//Example15_Moving_frame_general_vector_direction();
}