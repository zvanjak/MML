#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/BaseUtils.h"

#include "mpl/Mechanics/Projectiles2D.h"
#endif

#include <iostream>

using namespace MML;
using namespace MPL;

void Projectile_trajectory(IProjectile2DODE &projSys, Real angle, Real initVel, Real initHeight, std::string title);
void Projectile_trajectory_multiple_angles(IProjectile2DODE& projSys, Vector<Real> angles, Real initVel, Real initHeight, std::string title);
void Compare_solutions_for_set_of_angles(IProjectile2DODE& sys1, std::string label1, IProjectile2DODE& sys2, std::string label2,
																				 Real initVel, Real initHeight, Real dt, Vector<Real> angles);

void Chapter5_Throwing_things_in_the_air()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****               CHAPTER 5 - Throwing things in the air          ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Projectile2DInVacuumODE						projSys1;
	Projectile2DWithAirResistanceODE	projSys2(4e-5);
	
	Projectile2DChangingAirDensityODE projSys3isothermal(4e-5);
	Projectile2DChangingAirDensityODE projSys3adiabatic(4e-5, AirDensityModel::Adiabatic);

	BaseballWithDragCoeffDependentOnSpeedODE	projSys4smooth(BaseballType::Smooth); 
	BaseballWithDragCoeffDependentOnSpeedODE	projSys4rough(BaseballType::Rough);
	BaseballWithDragCoeffDependentOnSpeedODE	projSys4normal(BaseballType::Normal);

	Real angle = Utils::DegToRad(45);			// launching angle for our projectile
	Real initVel		= 700;								// initial velocity in m/s	
	Real initHeight = 0.0;								// initial height in m

	//Projectile_trajectory(projSys1, angle, initVel, intHeight, "Projectile in vacuum - 700 m_s");

	Vector<Real> angles{ Utils::DegToRad(20), Utils::DegToRad(30), Utils::DegToRad(35), Utils::DegToRad(40), Utils::DegToRad(45), 
											 Utils::DegToRad(50), Utils::DegToRad(55), Utils::DegToRad(60), Utils::DegToRad(65), Utils::DegToRad(70) };

	//Projectile_trajectory_multiple_angles(projSys1, angles, 700, 0, "Projectile in vacuum");


	Vector<Real> angles2{ Utils::DegToRad(30), Utils::DegToRad(40),Utils::DegToRad(45), Utils::DegToRad(50), Utils::DegToRad(55) };
	
	//Compare_solutions_for_set_of_angles(projSys1, "Vacuum", projSys2, "Air resist.", initVel, initHeight, 0.1, angles2);

	//Compare_solutions_for_set_of_angles(projSys2, "Air resist.", projSys3adiabatic, "Chang.air dens.", initVel, initHeight, 0.1, angles2);

	//Compare_solutions_for_set_of_angles(projSys2, "Air resist.", projSys3isothermal, "Chang.air dens. isothermal", initVel, initHeight, 0.1, angles2);
	
	initVel = 35;		// initial velocity in m/s

	Compare_solutions_for_set_of_angles(projSys4normal, "Normal", projSys4smooth, "Smooth", initVel, initHeight, 0.1, angles2);
	Compare_solutions_for_set_of_angles(projSys4normal, "Normal", projSys4rough, "Rough", initVel, initHeight, 0.1, angles2);
	Compare_solutions_for_set_of_angles(projSys4smooth, "Smooth", projSys4rough, "Rough", initVel, initHeight, 0.1, angles2);

	// COMPARE EULER I RK5 SOLUTION, VERLET VELOCITY

	// PRECISION

	// REQUIRED NUMBER OF STEPS FOR GIVEN PRECISION
}
