#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "algorithms/Statistics.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"
#include "tools/Timer.h"

#include "mpl/CollisionSimulator2D/CollisionSimulator2D.h"
#include "mpl/CollisionSimulator2D/ContainerFactory2D.h"
#endif

using namespace MML;
using namespace MPL;


// simulation with two types of balls, with randomized masses, radii, positions and velocities
void Collision_Simulator_2D_N_random_balls()
{
	int N = 10000; // Number of random balls

	const std::string& color1 = "Red";
	const std::string& color2 = "Blue";
	double minMass1 = 1.0, maxMass1 = 2.0, minRad1 = 1.0, maxRad1 = 2.0, velocityRange1 = 5.0;
	double minMass2 = 1.0, maxMass2 = 2.0, minRad2 = 1.0, maxRad2 = 2.0, velocityRange2 = 5.0;

	auto box = ContainerFactory2D::CreateConfig1(1000, 800, N, 
																						color1, minMass1, maxMass1, minRad1, maxRad1, velocityRange1,
																						color2, minMass2, maxMass2, minRad2, maxRad2, velocityRange2);

	CollisionSimulator2D simulator(box, 4, 4);

	int numSteps = 100;
	double dT = 1.0; // time step for simulation
	auto simResults = simulator.Simulate(numSteps, dT);

	std::string fileName = "collision_sim_" + std::to_string(N) + "_random_balls.txt";
	simulator.Serialize(GetResultFilesPath() + fileName, simResults, dT);

	Visualizer::VisualizeParticleSimulation2D(fileName);
}

// simulation with two types of balls, initially each type occupies one half of the container
// random balls position within each half
void Collision_Simulator_2D_two_types_initially_separated()
{
	int		 nBall1 = 5000;
	double mass1 = 5, rad1 = 2.5, vel1 = 10;
	int		 nBall2 = 5000;
	double mass2 = 5, rad2 = 2.5, vel2 = 50;		// five times faster than red balls

	auto box = ContainerFactory2D::CreateConfig2(1000, 800, "Red", nBall1, mass1, rad1, vel1, "Blue", nBall2, mass2, rad2, vel2);

	CollisionSimulator2D simulator(box, 10, 10);

	int numSteps = 20000;
	double dT = 0.05;										// time step for simulation	
	double totalTime = numSteps * dT;		// maximum time for simulation
	
	std::cout << "Starting simulation with " << numSteps << " steps" << std::endl;
	auto simResults = simulator.Simulate(numSteps, dT);

	std::cout << "Simulation finished." << std::endl;

	simulator.Serialize(GetResultFilesPath() + "collision_sim_2d_1000_x_800_1000_steps.txt", simResults, dT, 100);

	Visualizer::VisualizeParticleSimulation2D("collision_sim_2d_1000_x_800_1000_steps.txt");

/*

	//for (int i = 0; i < numSteps; i += 10)
	//{
	//	double minSpeed, maxSpeed, avgSpeed, speedDev;
	//	simResults.CalcAllBallsStatistic(i, minSpeed, maxSpeed, avgSpeed, speedDev);

	//	// print step number and statistics
	//	std::cout << "Step " << i << ": Min Speed = " << minSpeed << ", Max Speed = " << maxSpeed
	//		<< ", Avg Speed = " << avgSpeed << ", Speed Dev = " << speedDev << std::endl;
	//}

	// Average speed of balls by color visualization
	Vector<Real> vecTime, avgSpeedRed, avgSpeedBlue;
	for (int i = 0; i < numSteps; i += 10)
	{
		auto avgs = simResults.AvgSpeedByColor(i);

		if( avgs[0].first == "Red" ) {
			avgSpeedRed.push_back(avgs[0].second);
			avgSpeedBlue.push_back(avgs[1].second);
		}
		else {
			avgSpeedRed.push_back(avgs[1].second);
			avgSpeedBlue.push_back(avgs[0].second);
		}
		vecTime.push_back(i * dT);
	}

	// create linear interpolation functions for average speeds
	LinearInterpRealFunc avgSpeedRedFunc(vecTime, avgSpeedRed);
	LinearInterpRealFunc avgSpeedBlueFunc(vecTime, avgSpeedBlue);

	Visualizer::VisualizeMultiRealFunction({ avgSpeedRedFunc, avgSpeedBlueFunc },
																					"Average speeds of balls by color", { "Red balls", "Blue balls" },
																					vecTime[0], vecTime[vecTime.size() - 1], 100,
																					"avg_speed_by_color.txt");

	// Let's visualize evolution of x-coordinate of center of mass evolution for both types of balls
	Vector<Real> comRed, comBlue;
	for (int i = 0; i < numSteps; i += 10)
	{
		auto avgs = simResults.CentreOfMassByColor(i);

		if (avgs[0].first == "Red") {
			comRed.push_back(avgs[0].second.X());
			comBlue.push_back(avgs[1].second.X());
		}
		else {
			comRed.push_back(avgs[1].second.X());
			comBlue.push_back(avgs[0].second.X());
		}
	}

	// create linear interpolation functions for center of mass x-coordinates
	LinearInterpRealFunc comRedFunc(vecTime, comRed);
	LinearInterpRealFunc comBlueFunc(vecTime, comBlue);

	Visualizer::VisualizeMultiRealFunction({ comRedFunc, comBlueFunc },
																					"Centre of mass x-coordinate by color", { "Red balls", "Blue balls" },
																					vecTime[0], vecTime[vecTime.size() - 1], 100,
																					"com_x_by_color.txt");
*/
}

// energetic core of N particles, with M particles around it, which are moving randomly
void Collision_Simulator_2D_3()
{
	int		 numBalls = 50000;
	double mass1 = 5, rad1 = 1, vel1 = 2;

	int		 numEnergetic = 100;
	double mass2 = 5, rad2 = 1, vel2 = 500;

	double posRadius = 25;

	auto box = ContainerFactory2D::CreateConfig3(1000, 800, numBalls, mass1, rad1, vel1, numEnergetic, posRadius, mass2, rad2, vel2);
	
	CollisionSimulator2D simulator(box, 10, 10);

	int numSteps = 1000;
	double dT = 0.005; // time step for simulation
	
	auto simResults = simulator.Simulate(numSteps, dT);

	simulator.Serialize(GetResultFilesPath() + "collision_sim_2d_3.txt", simResults, dT);

	Visualizer::VisualizeParticleSimulation2D("collision_sim_2d_3.txt");
}

void Collision_Simulator_2D_Brown_motion()
{
	auto box = ContainerFactory2D::CreateBrownMotionConfig(1000, 800, 300);

	CollisionSimulator2D simulator(box);

	int numSteps = 100;
	double dT = 1.0;
	auto simResults = simulator.Simulate(numSteps, dT);

	simulator.Serialize(GetResultFilesPath() + "collision_sim_brown_motion.txt", simResults, dT);

	//Visualizer::VisualizeParticleSimulation2D("collision_sim_brown_motion.txt");

	auto path = simResults.getPathForBall(0);

	std::vector<VectorN<Real, 2>> res;
	for (int i = 0; i < path.size(); i++)
		res.push_back(VectorN<Real, 2>{path[i].X(), path[i].Y()});

	Serializer::SaveAsParamCurve<2>(res, "PARAMETRIC_CURVE_CARTESIAN_2D", "brown_motion_main_body_path.txt",
																	0, numSteps, path.size(),
																	GetResultFilesPath() + "brown_motion_main_body_path.txt");

	Visualizer::VisualizeMultiParamCurve2D({ "brown_motion_main_body_path.txt", "brown_motion_main_body_path1.txt" , "brown_motion_main_body_path2.txt", "brown_motion_main_body_path3.txt", 
																					 "brown_motion_main_body_path4.txt" , "brown_motion_main_body_path5.txt", "brown_motion_main_body_path6.txt" });

	//ParametricCurve2 curve(path, 1.0);
}

void Collision_Simulator_2D_Brown_motion_analysis()
{
	int NumSmallBalls = 500;
	std::vector<int> listStepNum{ 50, 100, 200, 500, 1000 };

	std::cout << "NUM SMALL BALLS = " << NumSmallBalls << std::endl;

	for (int numSteps : listStepNum)
	{
		std::cout << "NUM STEPS = " << numSteps << std::endl;

		double sumDist = 0.0;
		int NumRepeats = 500;
		for (int i = 0; i < NumRepeats; i++)
		{
			auto box = ContainerFactory2D::CreateBrownMotionConfig(1000, 800, NumSmallBalls);

			CollisionSimulator2D simulator(box);

			auto simResults = simulator.Simulate(numSteps, 1.0);

			auto path = simResults.getPathForBall(0);

			Pnt2Cart lastPos = path[path.size() - 1];

			double dist = lastPos.Dist(Pnt2Cart(500, 400));
			sumDist += dist;

			//std::cout << lastPos.X() << " " << lastPos.Y() << "  -  " << dist << std::endl;
		}
		double avgDist = sumDist / NumRepeats;

		std::cout << "AVG DIST = " << avgDist << std::endl;
	}

	// for given number of steps in simulation, investigate how average distance behaves 
	// with large number of trials
	//int numSteps = 50;
	//std::cout << "NUM STEPS = " << numSteps << std::endl;

	//double sumDist = 0.0;
	//std::vector<int> listNumRepeats{ 10 , 25, 50, 100, 200, 500, 1000 };
	//for (int numRepeats : listNumRepeats)
	//{
	//	for (int i = 0; i < numRepeats; i++)
	//	{
	//		auto box = ConfigFactory2D::CreateBrownMotionConfig(1000, 800, NumSmallBalls);

	//		CollisionSimulator2D simulator(box);

	//		auto simResults = simulator.Simulate(numSteps, 1.0);

	//		auto path = simResults.getPathForBall(0);

	//		Pnt2Cart lastPos = path[path.size() - 1];

	//		double dist = lastPos.Dist(Pnt2Cart(500, 400));
	//		sumDist += dist;

	//		// keep max and min distances
	//		// std::cout << lastPos.X() << " " << lastPos.Y() << "  -  " << dist << std::endl;
	//	}
	//	double avgDist = sumDist / numRepeats;
	//	std::cout << "NUM REPEATS = " << numRepeats << "  AVG DIST = " << avgDist << std::endl;
	//}
}

// for a given number of balls, and number of steps, compare simulation times for exact and fast implementation
void Collision_Simulator_2D_Test_fast()
{ 
	int N = 1000;						// Number of random balls
	int numSteps = 100;
	double dT = 1.0;				// time step for simulation

	const std::string& color1 = "Red";
	double minMass1 = 1.0, maxMass1 = 10.0, minRad1 = 1.0, maxRad1 = 4.0, velocityRange1 = 5.0;
	const std::string& color2 = "Blue";
	double minMass2 = 1.0, maxMass2 = 10.0, minRad2 = 1.0, maxRad2 = 4.0, velocityRange2 = 5.0;

	auto origConfig = ContainerFactory2D::CreateConfig1(1000, 800, N,
																									color1, minMass1, maxMass1, minRad1, maxRad1, velocityRange1,
																									color2, minMass2, maxMass2, minRad2, maxRad2, velocityRange2);

	Timer timer;

	timer.Start();

	auto box1 = origConfig; 
	CollisionSimulator2D simulator1(box1);
	auto simResults1 = simulator1.Simulate(numSteps, dT, CollisionSimulatorRunType::RunTypeExact);

	timer.MarkTime("Simulation time exact");

	auto box2 = origConfig;
	CollisionSimulator2D simulator2(box2, 4, 4);
	auto simResults2 = simulator2.Simulate(numSteps, dT, CollisionSimulatorRunType::RunTypeFastMultithread);

	timer.MarkTime("Simulation time fast");

	std::cout << "Serializing results..." << std::endl;
	std::string fileName = "collision_sim_" + std::to_string(N) + "_random_balls.txt";
	simulator2.Serialize(GetResultFilesPath() + fileName, simResults2, dT);

	timer.MarkTime("Serialization time");

	Visualizer::VisualizeParticleSimulation2D(fileName);

	timerprint();

	// usporediti rješenja simulacija
}

// TODO - investigate how different number of rows/cols influences simulation time
void Collision_Simulator_2D_Test_scaling_with_rows_cols()
{
	std::vector<std::string> strNumSteps;
	std::vector<LinearInterpRealFunc> simTimeFuncs;

	// std::map containing as key the number of balls, and as value vector containing
	// nmber of subdivision that we wanna test for that number of balls
	std::map<int, Vector<int>> mapNumBallsToNumSubdiv{  {500, Vector<int>{1, 2, 3, 4, 5} },
																										 {1000, Vector<int>{1, 2, 3, 4, 5, 6} },
																										 {2000, Vector<int>{1, 2, 3, 4, 5, 6, 7, 8} },
																										 {5000, Vector<int>{2, 3, 4, 5, 6, 7, 8, 9} },
																										{10000, Vector<int>{3, 4, 5, 6, 7, 8, 9, 10} }
																										//{20000, Vector<int>{4, 5, 6, 7, 8, 9, 10} }
	};

	for(auto & it : mapNumBallsToNumSubdiv)
	{
		std::cout << "Testing for number of balls: " << it.first << std::endl;
		
		int N = it.first;
		int numSteps = 100;
		double dT = 1.0;

		const std::string& color1 = "Red";
		double minMass1 = 1.0, maxMass1 = 10.0, minRad1 = 1.0, maxRad1 = 4.0, velocityRange1 = 5.0;
		const std::string& color2 = "Blue";
		double minMass2 = 1.0, maxMass2 = 10.0, minRad2 = 1.0, maxRad2 = 4.0, velocityRange2 = 5.0;

		auto box = ContainerFactory2D::CreateConfig1(1000, 800, N,
			color1, minMass1, maxMass1, minRad1, maxRad1, velocityRange1,
			color2, minMass2, maxMass2, minRad2, maxRad2, velocityRange2);

		Timer timer;

		Vector<Real> xCoordNumSubdivs;
		Vector<Real> yCoordSimTimes;
		for (int numSubdivs : it.second)
		{
			std::cout << "Num subdivisions = " << numSubdivs << " ";

			if( numSubdivs == 1 )
			{
				CollisionSimulator2D simulator1(box);				// no rows/cols

				timer.Start();
				auto simResults1 = simulator1.Simulate(numSteps, dT, CollisionSimulatorRunType::RunTypeExact);
				timer.MarkTime("Simulation time with " + std::to_string(numSubdivs) + " subdivs");
				std::cout << ", Simulation time: " << timer.GetTotalTime() << " seconds." << std::endl;

				xCoordNumSubdivs.push_back(numSubdivs);
				yCoordSimTimes.push_back(timer.GetTotalTime());
			}
			else
			{
				CollisionSimulator2D simulator1(box, numSubdivs, numSubdivs);

				timer.Start();

				auto simResults2 = simulator1.Simulate(numSteps, dT, CollisionSimulatorRunType::RunTypeFast);
				timer.MarkTime("Simulation time with " + std::to_string(numSubdivs) + " subdivs");
				std::cout << ", Simulation time: " << timer.GetTotalTime() << " seconds." << std::endl;
				
				xCoordNumSubdivs.push_back(numSubdivs);
				yCoordSimTimes.push_back(timer.GetTotalTime());
			}
		}

		// create linear interpolation function for simulation time
		LinearInterpRealFunc simTimeFunc(xCoordNumSubdivs, yCoordSimTimes);

		simTimeFuncs.push_back(simTimeFunc);
		strNumSteps.push_back("N - " + std::to_string(N));
	}

	// visualize all simulation time functions on one plot
	Visualizer::VisualizeMultiRealFunctionSeparately(simTimeFuncs,
																				"Simulation time vs number of subdivs for different number of balls", strNumSteps,
																				0, 10, 100, "sim_time_multi_N.txt");

}

// visualizing pressure on walls during simulation
void Collision_Simulator_2D_Pressure_on_walls()
{
	int numBalls = 10000;

	const std::string& color1 = "Red";
	double minMass1 = 1.0, maxMass1 = 10.0, minRad1 = 1.0, maxRad1 = 2.0, velocityRange1 = 5.0;
	const std::string& color2 = "Blue";
	double minMass2 = 1.0, maxMass2 = 10.0, minRad2 = 1.0, maxRad2 = 2.0, velocityRange2 = 5.0;

	auto box = ContainerFactory2D::CreateConfig1(1000, 1000, numBalls,
																							color1, minMass1, maxMass1, minRad1, maxRad1, velocityRange1,
																							color2, minMass2, maxMass2, minRad2, maxRad2, velocityRange2);

	BoxContainerWallPressureRecorder2D recorder(numBalls);
	box.setWallPressureRecorder(&recorder);

	CollisionSimulator2D simulator(box, 8, 8, &recorder);

	int numSteps = 100;
	double dT = 1;
	auto simResults = simulator.Simulate(numSteps, dT, CollisionSimulatorRunType::RunTypeFastMultithread);

	auto leftWallPressure = recorder.getTranferedMomentumPerStepForWall(0);
	auto rightWallPressure = recorder.getTranferedMomentumPerStepForWall(1);
	auto bottomWallPressure = recorder.getTranferedMomentumPerStepForWall(2);
	auto topWallPressure = recorder.getTranferedMomentumPerStepForWall(3);
	
	std::cout << "Step      Left Wall      Right Wall       Bottom Wall       Top Wall      Total" << std::endl;
	std::cout << "----------------------------------------------------------------------------------------" << std::endl;

	// print pressure on walls for each step
	for (int i = 0; i < leftWallPressure.size(); i++)
	{
		std::cout.precision(5);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout.setf(std::ios::showpoint);
		std::cout.setf(std::ios::right);
		std::cout.width(3);
		std::cout << i << " ";
		std::cout.width(15);
		std::cout << leftWallPressure[i] << " ";
		std::cout.width(15);
		std::cout << rightWallPressure[i] << " ";
		std::cout.width(15);
		std::cout << bottomWallPressure[i] << " " ;
		std::cout.width(15);
		std::cout << topWallPressure[i] << " " ;
		// print total pressure on walls
		std::cout.width(15);
		std::cout << leftWallPressure[i] + rightWallPressure[i] + bottomWallPressure[i] + topWallPressure[i] << std::endl;
	}
	std::cout << "----------------------------------------------------------------------------------------" << std::endl;


	auto totalWallPressure = recorder.getTotalTransferedMomentumPerStep();

	Vector<Real> x, y;
	for (int i = 0; i < totalWallPressure.size(); i++)
	{
		x.push_back(i);
		y.push_back(totalWallPressure[i]);
	}
	LinearInterpRealFunc funcTotalPressure(x, y);

	// visualize function of total pressure on walls
	Visualizer::VisualizeRealFunction(funcTotalPressure, "Total pressure on walls during simulation", 
																		0, numSteps - 1, numSteps-1, "total_pressure_on_walls.txt");

	// visualize pressure on all four walls separately, in same graph
	
	
	//simulator.Serialize(GetResultFilesPath() + "collision_sim_pressure_on_walls.txt", simResults.BallPosList);

	//Visualizer::VisualizeParticleSimulation2D("collision_sim_pressure_on_walls.txt");
}

// for a given number of steps, compare average pressure on walls for different number of balls
void Collision_Simulator_2D_Comparing_pressure_for_diff_N()
{
	int numSteps = 500;

	Vector<int>  numBallsList{ 1000, 2000, 4000, 8000, 16000 };
	Vector<Real> avgPressure;

	std::vector<std::string> strNumBalls;
	std::vector<LinearInterpRealFunc> totalPressureFuncs;

	std::cout << numSteps << " steps for each simulation." << std::endl;
	for (int i=0; i<numBallsList.size(); i++)
	{
		int numBalls = numBallsList[i];

		const std::string color = "Red";
		double mass = 1.0, rad = 1.0, velocityRange = 5.0;

		auto box = ContainerFactory2D::CreateConfigNSameBalls(1000, 1000, numBalls,
																													mass, rad, velocityRange, color);

		BoxContainerWallPressureRecorder2D recorder(numBalls);
		box.setWallPressureRecorder(&recorder);

		CollisionSimulator2D *pSimulator;
		if( numBalls < 10000 )
			pSimulator = new CollisionSimulator2D(box, 4, 4, &recorder);
		else
			pSimulator = new CollisionSimulator2D(box, 8, 8, &recorder);

		auto simResults = pSimulator->Simulate(numSteps, 1.0, CollisionSimulatorRunType::RunTypeFastMultithread);

		auto totalWallPressure = recorder.getTotalTransferedMomentumPerStep();

		Vector<Real> x, y;
		for (int i = 0; i < totalWallPressure.size(); i++)
		{
			x.push_back(i);
			y.push_back(totalWallPressure[i]);
		}

		LinearInterpRealFunc funcTotalPressure(x, y);

		totalPressureFuncs.push_back(funcTotalPressure);
		strNumBalls.push_back("N - " + std::to_string(numBallsList[i]));

		// print average pressure on walls for each N
		Real avg, stdDev;
		Statistics::AvgStdDev(y, avg, stdDev);
		avgPressure.push_back(avg);
		std::cout << "N = " << numBalls << "   Avg = " << avg << "   " << "StdDev = " << stdDev << "StdDev in perc.of avg = " << stdDev / avg * 100.0 << "%" << std::endl;


		// visualize function of total pressure on walls
		//Visualizer::VisualizeRealFunction(totalPressure, "Total pressure on walls", 0, numSteps - 1, numSteps - 1, "total_pressure_on_walls.txt");

		delete pSimulator;		// delete simulator to free memory
	}

	// visualize all total pressure functions on one plot
	Visualizer::VisualizeMultiRealFunction(totalPressureFuncs,
																				"Total pressure on walls for different number of balls", strNumBalls,
																				0, numSteps, numSteps, "total_pressure_multi_N.txt");

	// now visualize graph showing dependence of average pressure on N
	Vector<Real> x;
	for(int i = 0; i < numBallsList.size(); i++)
		x.push_back(numBallsList[i]);
	
	LinearInterpRealFunc avgPressureFunc(x, avgPressure);

	Visualizer::VisualizeRealFunction(avgPressureFunc, "Average pressure on walls for different number of balls", 
																		x[0], numBallsList[numBallsList.size() - 1], numSteps-1, "avg_pressure_on_walls_diff_N.txt");
}

// for a given number of steps, compare average pressure on walls for different size of the container
void Collision_Simulator_2D_Comparing_pressure_for_diff_V()
{
	int numSteps = 300;
	int numBalls = 5000;		// number of balls in each simulation

	Vector<int>  containerSizeList{ 500, 700, 1000, 1300, 1500, 2000 };
	Vector<Real> avgPressure;

	std::vector<std::string> strContSize;
	std::vector<LinearInterpRealFunc> totalPressureFuncs;

	std::cout << "Number of steps = " << numSteps << std::endl;
	std::cout << "Number of balls = " << numBalls << std::endl;
	for (int i = 0; i < containerSizeList.size(); i++)
	{
		int contSize = containerSizeList[i];

		const std::string color = "Red";
		double mass = 1.0, rad = 1.0, velocityRange = 5.0;

		auto box = ContainerFactory2D::CreateConfigNSameBalls(contSize, contSize, numBalls,
																													mass, rad, velocityRange, color);

		BoxContainerWallPressureRecorder2D recorder(numBalls);
		box.setWallPressureRecorder(&recorder);

		CollisionSimulator2D simulator(box, 6, 6, &recorder);;

		auto simResults = simulator.Simulate(numSteps, 1.0, CollisionSimulatorRunType::RunTypeFastMultithread);

		auto totalWallPressure = recorder.getTotalTransferedMomentumPerStep();

		Vector<Real> x, y;
		for (int i = 0; i < totalWallPressure.size(); i++)
		{
			x.push_back(i);
			y.push_back(totalWallPressure[i]);
		}

		LinearInterpRealFunc funcTotalPressure(x, y);

		totalPressureFuncs.push_back(funcTotalPressure);
		strContSize.push_back("Cont.size - " + std::to_string(containerSizeList[i]));
		
		Real avg, stdDev;
		Statistics::AvgStdDev(y, avg, stdDev);
		avgPressure.push_back(avg);
		std::cout << "Container size = " << contSize << "   Avg = " << avg << "   " << "StdDev = " << stdDev 
							<< " StdDev in perc.of avg = " << stdDev / avg * 100.0 << "%" << std::endl;
	}

	// visualize all total pressure functions on one plot
	Visualizer::VisualizeMultiRealFunction(totalPressureFuncs,
																				"Total pressure on walls for different size of container", strContSize,
																				0, numSteps-1, numSteps-1, "total_pressure_multi_V.txt");

	// now visualize graph showing dependence of average pressure on V
	Vector<Real> x_squared;
	for (int i = 0; i < containerSizeList.size(); i++)
		x_squared.push_back(containerSizeList[i]*containerSizeList[i]);

	LinearInterpRealFunc avgPressureFunc(x_squared, avgPressure);

	double maxX = POW2(containerSizeList[containerSizeList.size() - 1]);
	Visualizer::VisualizeRealFunction(avgPressureFunc, "Average pressure on walls for different area of container",
																		x_squared[0], maxX, 100, "avg_pressure_on_walls_diff_V.txt");
}

// for a given number of steps, compare average pressure on walls for different speed of balls
void Collision_Simulator_2D_Comparing_pressure_for_diff_speed_of_balls()
{
	int numSteps = 300;
	int numBalls = 5000;		// number of balls in each simulation
	
	Vector<Real>  speedList{ 1.0, 2.0, 5.0, 8.0, 10.0, 12.0, 15.0 };
	Vector<Real> avgPressure;
	
	std::vector<std::string> strSpeed;
	std::vector<LinearInterpRealFunc> totalPressureFuncs;
	
	std::cout << "Number of steps = " << numSteps << std::endl;
	std::cout << "Number of balls = " << numBalls << std::endl;
	for (int i = 0; i < speedList.size(); i++)
	{
		double speed = speedList[i];

		const std::string color = "Red";
		double mass = 1.0, rad = 1.0, velocityRange = speed;
		auto box = ContainerFactory2D::CreateConfigNSameBalls(1000, 1000, numBalls,
																													mass, rad, velocityRange, color);
		BoxContainerWallPressureRecorder2D recorder(numBalls);
		box.setWallPressureRecorder(&recorder);

		CollisionSimulator2D simulator(box, 6, 6, &recorder);;
		
		auto simResults = simulator.Simulate(numSteps, 0.5, CollisionSimulatorRunType::RunTypeFastMultithread);
		
		auto totalWallPressure = recorder.getTotalTransferedMomentumPerStep();
		Vector<Real> x, y;
		for (int j = 0; j < totalWallPressure.size(); j++)
		{
			x.push_back(j);
			y.push_back(totalWallPressure[j]);
		}
		LinearInterpRealFunc funcTotalPressure(x, y);
	
		totalPressureFuncs.push_back(funcTotalPressure);
		strSpeed.push_back("Speed - " + std::to_string(speedList[i]));
		Real avg, stdDev;
		Statistics::AvgStdDev(y, avg, stdDev);
		avgPressure.push_back(avg);
		std::cout << "Speed = " << speed << "   Avg = " << avg << "   " << "StdDev = " << stdDev 
							<< " StdDev in perc.of avg = " << stdDev / avg * 100.0 << "%" << std::endl;
	}
	// visualize all total pressure functions on one plot
	Visualizer::VisualizeMultiRealFunction(totalPressureFuncs,
																				"Total pressure on walls for different speed of balls", strSpeed,
																				0, numSteps - 1, numSteps - 1, "total_pressure_multi_speed.txt");

	// now visualize graph showing dependence of average pressure on speed
	Vector<Real> x_speed;
	for (int i = 0; i < speedList.size(); i++)
		x_speed.push_back(speedList[i]);

	LinearInterpRealFunc avgPressureFunc(x_speed, avgPressure);

	double maxX = speedList[speedList.size() - 1];
	Visualizer::VisualizeRealFunction(avgPressureFunc, "Average pressure on walls for different speed of balls",
																		x_speed[0], maxX, 100, "avg_pressure_on_walls_diff_speed.txt");

	// and visualization of square root of average pressure on walls depending on speed of molecules
	Vector<Real> x_speed_sqrt, y_speed_sqrt;
	for (int i = 0; i < speedList.size(); i++)
	{
		x_speed_sqrt.push_back(speedList[i]);
		y_speed_sqrt.push_back(sqrt(avgPressure[i]));
	}

	LinearInterpRealFunc avgPressureSqrtFunc(x_speed_sqrt, y_speed_sqrt);

	Visualizer::VisualizeRealFunction(avgPressureSqrtFunc, "Square root of average pressure on walls for different speed of balls",
		x_speed_sqrt[0], maxX, 100, "avg_pressure_sqrt_on_walls_diff_speed.txt");

}

// TODO - investigate how number of steps influences calculated average pressure on walls

// Investigate distibutions of velocities and positions of balls
void Collision_Simulator_2D_Test_statistics()
{
	int N = 200; // Number of random balls

	const std::string& color1 = "Red";
	const std::string& color2 = "Blue";
	double minMass1 = 1.0, maxMass1 =10.0, minRad1 = 3.0, maxRad1 = 5.0, velocityRange1 = 5.0;
	double minMass2 = 1.0, maxMass2 = 5.0, minRad2 = 3.0, maxRad2 = 5.0, velocityRange2 = 15.0;

	auto box = ContainerFactory2D::CreateConfig1(1000, 800, N,
		color1, minMass1, maxMass1, minRad1, maxRad1, velocityRange1,
		color2, minMass2, maxMass2, minRad2, maxRad2, velocityRange2);

	CollisionSimulator2D simulator(box);

	int numSteps = 500;
	SimResultsCollSim2D simResults = simulator.Simulate(numSteps, 1.0, CollisionSimulatorRunType::RunTypeExact);

	// print average and standard deviation of velocities at each step
	for(int i=0; i < numSteps; i+=20)
	{
		double minSpeed, maxSpeed, avgSpeed, speedDev;
		simResults.CalcAllBallsStatistic(i, minSpeed, maxSpeed, avgSpeed, speedDev);

		// print step number and statistics
		std::cout << "Step " << i << ": Min Speed = " << minSpeed << ", Max Speed = " << maxSpeed 
			<< ", Avg Speed = " << avgSpeed << ", Speed Dev = " << speedDev << std::endl;
	}

	// print average speed by color
	for (int i = 0; i < numSteps; i += 10)
	{
		auto avgs = simResults.AvgSpeedByColor(i);
		std::cout << "Step " << i << ": Avg Speed " << avgs[0].first << " = " << avgs[0].second 
			<< ", Avg Speed " << avgs[1].first << " = " << avgs[1].second << std::endl;
	}

	//std::string fileName = "collision_sim_" + std::to_string(N) + "_test_statistics.txt";
	//simulator.Serialize(GetResultFilesPath() + fileName, simResults.BallPosList);

	//Visualizer::VisualizeParticleSimulation2D(fileName);
}

// TODO - setup framework for evaluating and analyzing results of collision simulations with different dT
//   naci onu tocku nakon koje smanjivanje dT ne vodi drugacijim rezultatima

void Example4_collision_calculator_2D()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                 EXAMPLE 4 - collision calculator              ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Collision_Simulator_2D_N_random_balls();
	//Collision_Simulator_2D_two_types_initially_separated();
	//Collision_Simulator_2D_3();
	//Collision_Simulator_2D_Brown_motion();
	//Collision_Simulator_2D_Brown_motion_analysis();
	//Collision_Simulator_2D_Test_fast();
	Collision_Simulator_2D_Test_scaling_with_rows_cols();
	//Collision_Simulator_2D_Pressure_on_walls();
	//Collision_Simulator_2D_Comparing_pressure_for_diff_N();
	//Collision_Simulator_2D_Comparing_pressure_for_diff_V();
	//Collision_Simulator_2D_Comparing_pressure_for_diff_speed_of_balls();
	//Collision_Simulator_2D_Test_statistics();
}