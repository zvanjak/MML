#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/VectorTypes.h"

#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"

#include "mpl/Base/PhysicalConstants.h"
#include "mpl/Base/GasConstants.h"

#include "mpl/CollisionSimulator3D/ContainerFactory3D.h"
#include "mpl/CollisionSimulator3D/CollisionSimulator3D.h"

#include "mpl/MolecularDynamics/IdealGasCalculator.h"
#endif

using namespace MML;
using namespace MPL;

void Collision_simulator_3D_N_random_balls()
{
	int N = 1000; // Number of random balls

	const std::string& color1 = "Red";
	const std::string& color2 = "Blue";
	double minMass1 = 1.0, maxMass1 = 2.0, minRad1 = 2.0, maxRad1 = 3.0, velocityRange1 = 5.0;
	double minMass2 = 1.0, maxMass2 = 2.0, minRad2 = 2.0, maxRad2 = 3.0, velocityRange2 = 5.0;

	auto box = ContainerFactory3D::CreateConfig1(1000, 1000, 1000, N,
								color1, minMass1, maxMass1, minRad1, maxRad1, velocityRange1,
								color2, minMass2, maxMass2, minRad2, maxRad2, velocityRange2);

	CollisionSimulator3D simulator(box);

	int	 numSteps = 100;
	double dT = 1.0;
	auto simResult = simulator.Simulate(numSteps, dT);

	Serializer::SaveParticleSimulation3D(GetResultFilesPath() + "collision_sim_3d_example1.txt", box._balls.size(), box._width, box._height, box._depth, 
		simResult.BallPosList, box.getBallColors(), box.getBallRadii(), dT);

	Visualizer::VisualizeParticleSimulation3D("collision_sim_3d_example1.txt");
}

void Collision_simulator_3D_two_types_initially_separated()
{
	int		 nBall1 = 5000;
	double mass1 = 5, rad1 = 4, vel1 = 10;
	int		 nBall2 = 5000;
	double mass2 = 5, rad2 = 4, vel2 = 10;

	auto box = ContainerFactory3D::CreateConfig2(1000, 1000, 1000, 
																							nBall1, mass1, rad1, vel1,
																							nBall2, mass2, rad2, vel2);

	CollisionSimulator3D simulator(box, 3, 3, 3);

	int	 numSteps = 500;
	double dT = 0.5;
	auto simResult = simulator.Simulate(numSteps, dT);

	Serializer::SaveParticleSimulation3D(GetResultFilesPath() + "collision_sim_3d_example2.txt", 
																			 box._balls.size(), box._width, box._height, box._depth,
																			 simResult.BallPosList, box.getBallColors(), box.getBallRadii(), dT);

	Visualizer::VisualizeParticleSimulation3D("collision_sim_3d_example2.txt");
}
 
void Collision_Simulator_3D_Test_fast()
{

}

void Collision_Simulator_3D_Comparing_pressure_for_diff_N()
{

}

void Collision_simulator_3D_ideal_gas_in_100_nm_cube()
{
	// for 100 nm at STP, we get 27.000 balls of ideal gas!!!
	double boxSide_nm = 100.0;	
	double boxSide		= boxSide_nm * 1e-9; 
	double px_per_nm	= 10;									// 10 pixels per nm, so 1 pixel = 0.1 nm

	// STP constants
	double T = 273.0;							// temperature in K
	double P = 101325.0;					// pressure in Pa (1 atm)

	// oxygen (O2) gas at 273 K
	double radius_O2_nm		= 0.3;						// radius of O2 molecule in nm
	double molar_mass_O2	= 0.032;					// molar mass in kg/mol
	double mass_molecule_O2_kg = molar_mass_O2 / PhyConst::AvogardoConstant; 
	double avg_speed_O2_m_s    = sqrt(3 * PhyConst::BoltzmannConstant * T / mass_molecule_O2_kg); // in m/s

	// calculate the number of molecules in a 100 nm cube
	double volume_m3	 = boxSide * boxSide * boxSide; 
	double n_molecules = P * volume_m3 / (PhyConst::BoltzmannConstant * T); // number of molecules using ideal gas law

	// simulating real (ideal) gas in 100 nm cube, as if it contained molecules of Oxygen (O2) at 273 K
	int numBalls = (int) n_molecules;
	int radius = radius_O2_nm * px_per_nm;			

	double avg_speed = avg_speed_O2_m_s;						// in m/s, but also in nm/ns, 
	double avg_speed_px = avg_speed * px_per_nm;		// in pixels/s

	// this means that our time scale is in nanoseconds, so 1 timestep is 1 ns
	// but, as that would mean that particles would move too fast (1460 pixels per ns)
	// we will use 100 timesteps per ns, so 1 timestep is 0.01 ns

	int boxSide_px  = boxSide_nm * px_per_nm;					// in pixels, so 100 nm in real world
	double mass			= molar_mass_O2;
	double rad			= radius;
	double velocity = avg_speed_px;

	auto box = ContainerFactory3D::CreateConfigNSameBalls(boxSide_px, boxSide_px, boxSide_px, numBalls,
																												mass, rad, velocity, "Red");

	int numSteps = 100;				// number of timesteps
	double dT = 0.001;					// in ns, so 0.001 ns per timestep

	CollisionSimulator3D simulator(box, 5, 5, 5);

	auto simResult = simulator.Simulate(numSteps, dT);

	Serializer::SaveParticleSimulation3D(GetResultFilesPath() + "ideal_gas.txt", box._balls.size(), box._width, box._height, box._depth,
																				simResult.BallPosList, box.getBallColors(), box.getBallRadii(), dT);

	Visualizer::VisualizeParticleSimulation3D("ideal_gas.txt");
}

void Collision_simulator_3D_simulating_air_as_ideal_gas_in_100_nm_cube()
{
	// for 100 nm at STP, we get 27.000 molecules of ideal gas!!!
	double boxSide_nm = 100.0;
	double boxSide = boxSide_nm * 1e-9;
	double px_per_nm = 10;									// 10 pixels per nm, so 1 pixel = 0.1 nm

	// STP constants
	double k_B = PhyConst::BoltzmannConstant;		// Boltzmann constant in J/K
	double T = 273.0;							// temperature in K
	double P = 101325.0;					// pressure in Pa (1 atm)

	const GasData& O2 = GasConstants::gasByFormula("O2");
	const GasData& N2 = GasConstants::gasByFormula("N2");
	const GasData& Ar = GasConstants::gasByFormula("Ar");
	const GasData& CO2 = GasConstants::gasByFormula("CO2");

	// oxygen (O2) gas at 273 K
	double radius_O2_nm				 = O2.RadiusKinetic_nm();		
	double molar_mass_O2			 = O2.molarMass_g_mol / 1000;			// molar mass in kg/mol
	double mass_molecule_O2_kg = O2.MoleculeMass_kg();
	double avg_speed_O2_m_s		 = IdealGasCalculator::AvgSpeed(molar_mass_O2, T); 

	// nitrogen (N2) gas at 273 K
	double radius_N2_nm = N2.RadiusKinetic_nm();						// radius of N2 molecule in nm
	double molar_mass_N2 = 0.028;					// molar mass in kg/mol
	double mass_molecule_N2_kg = molar_mass_N2 / 6.02214076e23; // mass in kg (using Avogadro's number)
	double avg_speed_N2_m_s = sqrt(3 * k_B * T / mass_molecule_N2_kg); // in m/s

	// argon (Ar) gas at 273 K
	double radius_Ar_nm = 0.2;						// radius of Ar atom in nm
	double molar_mass_Ar = 0.040;					// molar mass in kg/mol
	double mass_molecule_Ar_kg = molar_mass_Ar / 6.02214076e23; // mass in kg (using Avogadro's number)
	double avg_speed_Ar_m_s = sqrt(3 * k_B * T / mass_molecule_Ar_kg); // in m/s

	// carbon dioxide (CO2) gas at 273 K
	double radius_CO2_nm = 0.33;					// radius of CO2 molecule in nm
	double molar_mass_CO2 = 0.044;				// molar mass in kg/mol
	double mass_molecule_CO2_kg = molar_mass_CO2 / 6.02214076e23; // mass in kg (using Avogadro's number)
	double avg_speed_CO2_m_s = sqrt(3 * k_B * T / mass_molecule_CO2_kg); // in m/s

	// calculate TOTAL number of molecules in a 100 nm cube
	double volume_m3 = boxSide * boxSide * boxSide;
	double numBalls = (int) P * volume_m3 / (k_B * T); // number of molecules using ideal gas law

	double rad_O2 = radius_O2_nm * px_per_nm;
	double rad_N2 = radius_N2_nm * px_per_nm;
	double rad_Ar = radius_Ar_nm * px_per_nm;
	double rad_CO2 = radius_CO2_nm * px_per_nm;

	double avg_speed_O2_px = avg_speed_O2_m_s * px_per_nm;		// in pixels/s
	double avg_speed_N2_px = avg_speed_N2_m_s * px_per_nm;		// in pixels/s
	double avg_speed_Ar_px = avg_speed_Ar_m_s * px_per_nm;		// in pixels/s
	double avg_speed_CO2_px = avg_speed_CO2_m_s * px_per_nm;	// in pixels/s

	int boxSide_px = boxSide_nm * px_per_nm;					// in pixels, so 100 nm in real world

	auto box = ContainerFactory3D::CreateConfigGases( boxSide_px, boxSide_px, boxSide_px, numBalls,
																										rad_N2, rad_O2, rad_Ar, rad_CO2,
																										avg_speed_N2_px, avg_speed_O2_px, avg_speed_Ar_px, avg_speed_CO2_px,
																										"Blue", "Red", "Yellow", "Green");
	int numSteps = 10;				// number of timesteps
	double dT = 0.001;				// in ns, so 0.001 ns per timestep

	CollisionSimulator3D simulator(box, 5, 5, 5);

	auto simResult = simulator.Simulate(numSteps, dT);

	Serializer::SaveParticleSimulation3D(GetResultFilesPath() + "air_as_ideal_gas.txt", box._balls.size(), 
																			 box._width, box._height, box._depth,
																			 simResult.BallPosList, box.getBallColors(), box.getBallRadii(), dT);

	Visualizer::VisualizeParticleSimulation3D("air_as_ideal_gas.txt");
}

void Collision_simulator_3D()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                 EXAMPLE 4 - collision calculator 3D           ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Collision_simulator_3D_N_random_balls();
	//Collision_simulator_3D_two_types_initially_separated();
	//Collision_Simulator_3D_Test_fast();
	//Collision_Simulator_3D_Comparing_pressure_for_diff_N();
	//Collision_simulator_3D_ideal_gas_in_100_nm_cube();
  Collision_simulator_3D_simulating_air_as_ideal_gas_in_100_nm_cube();
}