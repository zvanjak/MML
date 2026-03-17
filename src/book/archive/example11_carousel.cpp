#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "core/Derivation.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"

#include "mpl/Mechanics/Carousel.h"
#include "mpl/Mechanics/MovingFrames.h"
#include "mpl/Mechanics/MovingFramesTransf.h"
#include "mpl/Mechanics/Projectiles3D.h"
#endif


using namespace MML;
using namespace MPL;


// We have carousel with R=1000 m radius, rotating counterclockwise in x-y plane with angular velocity w=0.1 rad/s
// At the distance of L=500 m from the center of carousel, and from the height of H = 2 m 
// we throw a ball with initial velocity v0 = 30 m/s at vertical angle alpha = 30 degrees
// and at angle to eastward direction (i.e. in the direction of carousel rotation) of beta = 30 degrees to the north

// We will form local Cartesian coordinate system ON THE CAROUSEL with the origin at the point of launch,
// with x-axis pointing outwards of the center of the carousel, y-axis pointing in the direction of carousel rotation, and z-axis pointing upwards
// This local coordinate system evidently moves with the carousel, and rotates with it, so the x-axis is always pointing outwards of the center of the carousel,

// We will also create global Cartesian coordinate system, where the origin is at the center of the carousel,
// and at the moment of launch, we will assume the center of the local coordinate system is at the point (L, 0, 0) in global Cartesian coordinates.
//  meaning that at the moment of launch, the local coordinate system axis are aligned with the global coordinate system axis

// We want to find correct landing point in both coordinate systems, and the time of flight.
void Demo_GiantCarousel()
{
	Real rotationPeriodInMinutes = 30;
	Real w = 2 * Constants::PI / (rotationPeriodInMinutes * 60); // angular velocity in rad/s
	
	Real R = 1000.0;			// radius in meters
	Real L = 500.0;				// distance from the center of carousel in meters
	Real H = 2.0;					// height in meters
	Real v0 = 30.0;				// magnitude of initial velocity in m/s
	
	Real alpha = Utils::DegToRad(45.0);			// vertical launch angle in radians
	Real beta = Utils::DegToRad(0);			// angle to eastward direction in radians

	// calculating launch velocity in local Cartesian system
	//	our angle is given with respect to eastward direction (ie local y-axis), with positive angle looking towards north
	//	so components of launch velocity in local system are:
	Real v0_x_local = -v0 * sin(beta) * cos(alpha);	
	Real v0_y_local = v0 * cos(beta) * cos(alpha);	
	Real v0_z_local = v0 * sin(alpha);

	// calculating time of flight (the same in both coordinate systems)
	Real g = 9.81;
	Real t_flight = (v0_z_local + std::sqrt(v0_z_local * v0_z_local + 2 * g * H)) / g; 

	std::cout << "Giant carousel parameters:" << std::endl;
	std::cout << "Radius of carousel                  : " << R << " meters" << std::endl;
	std::cout << "Distance from the center of carousel: " << L << " meters" << std::endl;
	std::cout << "Height of launch point              : " << H << " meters" << std::endl;
	std::cout << "Period of carousel rotation         : " << rotationPeriodInMinutes << " minutes" << std::endl;
	std::cout << "Angular velocity of carousel        : " << w << " rad/s" << std::endl;
	std::cout << "Additional velocity at pnt.of launch: " << L * w << " m/s" << std::endl;

	std::cout << "Initial velocity of launch:  " << v0 << " m/s" << std::endl;
	std::cout << "Vertical launch angle:       " << Utils::RadToDeg(alpha) << " degrees" << std::endl;
	std::cout << "Angle to eastward direction: " << Utils::RadToDeg(beta) << " degrees" << std::endl;

	std::cout << "\nTime of flight: " << t_flight << " seconds" << std::endl;

	std::cout << "\nLaunch point in local Cartesian system        : (" << L << ", 0.0, " << H << ")" << std::endl;
	std::cout << "Launch point in global Cartesian system       : (" << L << ", 0.0, " << H << ")" << std::endl;
	std::cout << "Launch point in local cylindrical coord.sys.  : ("
						<< L << ", " << Utils::RadToDeg(0.0) << " deg, " << H << ")" << std::endl;

	ProjectileInVacuumODE projectileSys;		// we'll use this for calculations in Cartesian coordinate system

	std::cout << "\n*********    SOLVING NAIVELY IN LOCAL COORDINATE SYSTEM    *********" << std::endl;

	Pnt3Cart landPosLocalCart1 = projectileSys.calculateLandingPos(L, 0.0, H, Vec3Cart(v0_x_local, v0_y_local, v0_z_local));

	std::cout << "\nLanding point in local Cartesian system       : (" << landPosLocalCart1.X() << ", " << landPosLocalCart1.Y() << ")" << std::endl;
	
	Real dx = landPosLocalCart1.X() - L;
	Real dy = landPosLocalCart1.Y();

	std::cout << "Landing point relative to starting point      : (" << dx << ", " << dy << ")" << std::endl;
	
	Vec3Cyl landPosLocalCyl1 = CoordTransfCartToCyl.transf(Vec3Cart(landPosLocalCart1.X(), landPosLocalCart1.Y(), landPosLocalCart1.Z()));

	std::cout << "Landing point in local cylindrical coord.sys. : ("
						<< landPosLocalCyl1.R() << ", " << Utils::RadToDeg(landPosLocalCyl1.Phi()) << " deg, " << landPosLocalCyl1.Z() << ")" << std::endl;
	
	std::cout << "\n*********    SOLVING EXACTLY IN GLOBAL COORDINATE SYSTEM    *********" << std::endl;
	RotatingCarouselFrame frame(w);
	RotatingCarouselLocalCartToFixedTransf transfLocalCartToFixed(frame);
	RotatingCarouselLocalCylToFixedTransf  transfLocalCylToFixed(frame);

	// additional velocity due to carousel rotation (which is in the plus y-direction of global Cartesian system)
	Real v_add_y = L * w;

	Pnt3Cart landPosGlobalCart = projectileSys.calculateLandingPos(L, 0.0, H, Vec3Cart(v0_x_local, v0_y_local + v_add_y, v0_z_local));
	
	std::cout << "\nLanding point in global Cartesian system      : (" << landPosGlobalCart.X() << ", " << landPosGlobalCart.Y() << ")" << std::endl;

	dx = landPosGlobalCart.X() - L;	
	dy = landPosGlobalCart.Y();	
	std::cout << "Landing point relative to starting point      : (" << dx << ", " << dy << ")" << std::endl;
	std::cout << "transforming to local system at t = t_flight" << std::endl;

	// transforming from global Cartesian coordinates to local cartesian coordinates
	Vec3Cart landPosLocal = transfLocalCartToFixed.transfInverse(t_flight, Vec3Cart(landPosGlobalCart.X(), landPosGlobalCart.Y(), landPosGlobalCart.Z()));

	std::cout << "Landing point in local Cart.system from global: ("
						<< landPosLocal.X() << ", " << landPosLocal.Y() << ", " << landPosLocal.Z() << ")" << std::endl;
	
	// transforming from global Cartesian coordinates to local cylindrical coordinates
	Vec3Cyl landPosLocalCyl2 = transfLocalCylToFixed.transfInverse(t_flight, Vec3Cart(landPosGlobalCart.X(), landPosGlobalCart.Y(), landPosGlobalCart.Z()));	
	
	std::cout << "Landing point in local cyl. system from global: ("
						<< landPosLocalCyl2.R() << ", " << Utils::RadToDeg(landPosLocalCyl2.Phi()) << " deg, " << landPosLocalCyl2.Z() << ")" << std::endl;
	
	std::cout << "\n*********    SOLVING IN ROTATING LOCAL CARTESIAN COORDINATE SYSTEM       *********" << std::endl;
	ProjectileInVacuumLaunchFromCarouselODE odeSys(w);

	Real dt = 0.01;
	int numSteps = 5000;			// expected number of steps in the solution
	ODESystemFixedStepSolver	odeSolver(odeSys, StepCalculators::EulerStepCalc);

	Vector<Real> initCond = odeSys.getInitCond(L, 0.0, H, Vec3Cart(v0_x_local, v0_y_local, v0_z_local));
	ODESystemSolution sol = odeSolver.integrate(initCond, 0.0, t_flight*1.00088, numSteps);

	// final coordinates in local Cartesian system
	// TODO - get exact coordinates in local system at the moment of landing
	Pnt3Cart landPosLocalCart3 = Pnt3Cart(sol.getXValues(0).back(), sol.getXValues(1).back(), sol.getXValues(2).back());

	std::cout << "\nFinal coordinates in local Cart.system        : (" 
						<< landPosLocalCart3.X() << ", " << landPosLocalCart3.Y() << ", " << landPosLocalCart3.Z() << ")" << std::endl;

	// convert to local cylindrical coordinates
	Vec3Cyl landPosLocalCyl3 = CoordTransfCartToCyl.transf(Vec3Cart(landPosLocalCart3.X(), landPosLocalCart3.Y(), landPosLocalCart3.Z()));

	std::cout << "Final coord in local cyl. system              : ("
						<< landPosLocalCyl3.R() << ", " << Utils::RadToDeg(landPosLocalCyl3.Phi()) << " deg, " << landPosLocalCyl3.Z() << ")" << std::endl;
	
	// calculate difference between simple inertial result and correct rotating frame result
	Pnt3Cyl pos1(landPosLocalCyl1.R(), landPosLocalCyl1.Phi(), landPosLocalCyl1.Z());
	Pnt3Cyl pos2(landPosLocalCyl2.R(), landPosLocalCyl2.Phi(), landPosLocalCyl2.Z());

	Real dist12 = pos1.Dist(pos2);

	std::cout << "\nDistance between local cylindrical coordinates from naive and correct solution: "
						<< dist12 << " meters" << std::endl;
	
	// calculate hom much did the launch point move during the flight?
	// we have to calculate the position of the launch point in global system at the moment of landing
	// transforming from local Cartesian coordinates to global cartesian coordinates
	Vec3Cart launchPosEnd = transfLocalCartToFixed.transf(t_flight, Vec3Cart(L, 0, 0));

	std::cout << "\nLaunch point position in global Cartesian system at the moment of landing: ("
						<< launchPosEnd.X() << ", " << launchPosEnd.Y() << ")" << std::endl;

	Pnt3Cart launchPosStart(L, 0.0, H); // position of launch point in global system at the moment of launch
	Real dx_launch = launchPosEnd.X() - launchPosStart.X();
	Real dy_launch = launchPosEnd.Y() - launchPosStart.Y();

	std::cout << "Launch point move during flight   : " << std::sqrt(dx_launch * dx_launch + dy_launch * dy_launch) << " meters" << std::endl;
	std::cout << "Movement in global x and y coords : (" << dx_launch << ", " << dy_launch << ")" << std::endl;



	/*
	// position of launch point in global system at the moment of landing (how much it moved during flight in global system)
	// we have to calculate roation for the time of flight
	Real angle_rotated = w * t_flight; // angle rotated by the carousel during the flight
	Real x_launch_point_global = xstart_global * cos(angle_rotated) - ystart_global * sin(angle_rotated);
	Real y_launch_point_global = xstart_global * sin(angle_rotated) + ystart_global * cos(angle_rotated);

	std::cout << "\nLaunch point position in global Cartesian system at the moment of landing: (" 
						<< x_launch_point_global << ", " << y_launch_point_global << ")" << std::endl;

	// we need to setup local rotating coordinate system (which is essential cylindrical coordinate system)
	// coordinates are R, phi, z, and at the moment of launch, we have:
	Real r_launch = L;
	Real phi_launch = 0;			// aligning local y-axis with global y-axis at the moment of launch
	Real z_launch = H;				// height of launch point in global system

	// and we need transformation from (static) global coordinate system to (moving/rotating) local coordinate system
	// where additional parameter is t, time of transformation
	Real x = landing_pos_global.X();
	Real y = landing_pos_global.Y();

	// TRANSFORM this to local rotating coordinate system
	Real t = t_flight; // time of transformation

	Real r_landing = std::sqrt(x * x + y * y); // radial coordinate in local system
	Real phi_landing = std::atan2(y, x) - w * t; // azimuthal angle in local system (subtract the angle rotated by the carousel during the flight)
	Real phi_landing_no_rotation = std::atan2(y, x); // azimuthal angle in local system without rotation

	std::cout << "\nLanding point in local cylindrical coordinate system: ("
		<< r_landing << ", " << Utils::RadToDeg(phi_landing) << " deg, " << 0.0 << ")" << std::endl;

	std::cout << "Landing point in local cylindrical coordinate system without rotation: ("
		<< r_landing << ", " << Utils::RadToDeg(phi_landing_no_rotation) << " deg, " << 0.0 << ")" << std::endl;

	// now we have to see how our local coordinate system is related to global coordinate system at the moment of landing
	// we have to calculate the position of the local coordinate system at the moment of landing
	Real angle_rotated_landing = w * t; // angle rotated by the carousel during the flight
	Real x_local_center = r_launch * cos(phi_launch + angle_rotated_landing); // x-coordinate of local center at the moment of landing
	Real y_local_center = r_launch * sin(phi_launch + angle_rotated_landing); // y-coordinate of local center at the moment of landing

	std::cout << "\nCenter of local coordinate system at the moment of launch in global Cartesian system: ("
		<< xstart_global << ", " << ystart_global << ")" << std::endl;

	std::cout << "Center of local coordinate system at the moment of landing in global Cartesian system: ("
		<< x_local_center << ", " << y_local_center << ")" << std::endl;
		*/
}

void Demo_CoriolisForceOnEarth()
{

}

void Demo_CarouselOnEarth()
{
}

// Demo - Earth is orbiting Sun, and we have a carousel on Earth, at given latitude and longitude.
// We will launch a projectile from the carousel, and calculate its landing point in both local and global coordinate systems.
// calculate landing point in carousel local coordinate system (both Cart and Cyl), in Earth local cartesian system with origin 
// at the center of the carousel, and in global Earth rotating frame with origin at the center of the Earth.
void Demo_CarouselOnEarthOrbitingSun()
{

}

void Example11_Carousel()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 11 - Carousel                        ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo_GiantCarousel();
	Demo_CoriolisForceOnEarth();
	Demo_CarouselOnEarth();
	Demo_CarouselOnEarthOrbitingSun();
}