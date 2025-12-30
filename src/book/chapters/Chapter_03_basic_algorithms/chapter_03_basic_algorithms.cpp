#include "MMLBase.h"


using namespace MML;

void Demo_Derivation();
void Demo_Integration();
void Demo_Interpolation();
void Demo_Lin_alg_sys_solvers();
void Demo_Root_finding();


int Chapter_03_basic_algorithms()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                    EXAMPLE 2 - Basic algorithms               ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo_Derivation();
	Demo_Integration();
	Demo_Interpolation();
	Demo_Lin_alg_sys_solvers();
	Demo_Root_finding();

	return 0;
}

