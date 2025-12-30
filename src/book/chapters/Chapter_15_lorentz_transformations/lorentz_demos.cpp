#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/VectorN.h"

#include "mpl/SpecialRelativity/LorentzTransformation.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;
using namespace MPL;

// Check that Lorentz transformation is equal to Galilean transformation for low velocities
void Chapter15_Lorentz_as_Gallilean_transf()
{
	Real speed = 0.8; // speed in units of c, so 0.8 means 80% of the speed of light

	// moving frame actually!
	CoordTransfLorentzXAxis ctLorentzX(speed);

	// we have event happening at T=10 s at the origin of rest (observer) frame
	Vector4Minkowski eventLocalRestFrame1{ REAL(10.0), REAL(0.0), REAL(0.0), REAL(0.0) };

	// how passanger at origin of moving frame sees this event?
	Vector4Minkowski eventMovingFrame1 = ctLorentzX.transf(eventLocalRestFrame1);

	std::cout << "Event in the local rest frame of reference: " << eventLocalRestFrame1 << std::endl;
	std::cout << "Event in the moving frame of reference    : " << eventMovingFrame1 << std::endl;

	// now, inverse transformation
	// we have an event at T'=10 s at the origin of moving frame
	Vector4Minkowski eventMovingFrame2{ REAL(10.0), REAL(0.0), REAL(0.0), REAL(0.0) };

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

// Demo with passenger moving along x-axis with speed v, and observer at rest
// Calculate Lorentz transformation of coordinates and time for this observer
void Chapter15_Moving_x_axis()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "***  EXAMPLE 15 - Moving frame in x-axis direction transformation   ***" << std::endl;
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
	Vector4Minkowski eventLocalRestFrame3{ REAL(10.0), REAL(10.0) * speed, REAL(0.0), REAL(0.0) };

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

// Transforming between local frame at rest and inertial frame moving in general direction
// (Both systems have the same orientation, so no rotation)
void Chapter15_Moving_frame_general_vector_direction()
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
	// simulating for 10 time-units
	for (int i = 0; i <= 10; ++i)
	{
		Pnt3Cart posAtT = pntA + velocity * i; // position at T=i seconds
		Vector4Minkowski xPos{ Real(i * 1.0), Real(i * speed), Real(0.0), Real(0.0) };

		// how passanger at origin of moving frame sees this event?
		Vector4Minkowski eventMovingFrame = ctGeneral.transf(xPos);
		std::cout << "Local rest frame : " << "T = " << std::setw(5) << xPos.T() << "  X = " << xPos.X() << "  ";
		std::cout << "Moving frame of reference : " << "T = " << eventMovingFrame.T() << "  X = " << eventMovingFrame.X() << std::endl;

		// at the origin of local frame there is a clock, and on the passenger ship
		// there is a telescope pointed at that clock at all times
		// when passenger ship passes midway point between A and B,
		// what time does he see on the clock at point A?
		if (i == 5)
		{
			// at T=5 seconds, the ship is at the point B
			Vector4Minkowski eventAtMidway = ctGeneral.transfInverse(Vector4Minkowski{ REAL(5.0), REAL(5.0) * speed, REAL(0.0), REAL(0.0) });
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

// Placeholder for spherical ball demonstration
void Chapter15_Spherical_ball()
{
	// TODO: Implementation pending
}
