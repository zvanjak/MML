#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#endif


using namespace MML;

void Example2_Derivation();
void Example2_Integration();
void Example2_Interpolation();
void Example2_Lin_alg_sys_solvers();
void Example2_Root_finding();


void Example2_Basic_algorithms()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                    EXAMPLE 2 - Basic algorithms               ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Example2_Derivation();
	//Example2_Integration();
	//Example2_Interpolation();
	Example2_Lin_alg_sys_solvers();
	//Example2_Root_finding();
}