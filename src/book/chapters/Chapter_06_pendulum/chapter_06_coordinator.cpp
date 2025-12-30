#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif

using namespace MML;

// Simple pendulum demonstrations
void Demo_SimplePendulum();
void Demo_SimplePendulum_Euler();

// Damped pendulum demonstrations  
void Demo_DampedPendulum();
void Demo_DampedPendulum_compare_three_cases();

// Forced pendulum demonstrations
void Demo_ForcedPendulum();
void Demo_Damped_forced_pendulum();
void Demo_Damped_forced_pendulum_compare_close_init_cond();

// Harmonic oscillator demonstrations
void Demo_HarmonicOscilator();
void Demo_DampedHarmonicOscilator();
void Demo_ForcedHarmonicOscilator();
void Demo_DampedForcedHarmonicOscilator();

void Chapter06_Pendulum()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 5 - Pendulum                        ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Demo_SimplePendulum();
	//Demo_SimplePendulumEuler();
	//Demo_DampedPendulum();
	Demo_DampedPendulum_compare_three_cases();
	//Demo_ForcedPendulum();
	//Demo_Damped_forced_pendulum();
	//Demo_Damped_forced_pendulum_compare_close_init_cond();

	//Demo_HarmonicOscilator();
	//Demo_DampedHarmonicOscilator();
	//Demo_ForcedHarmonicOscilator();
	//Demo_DampedForcedHarmonicOscilator();
}
