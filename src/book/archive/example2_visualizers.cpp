#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#endif


using namespace MML;

void Example3_WPF_visualizators();

void Example3_Visualizers()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                      EXAMPLE 3 - Visualizers                  ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Example3_WPF_visualizators();
}