#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#endif

using namespace MML;

void Example1_Vector();
void Example1_Matrix();
void Example1_Tensor();
void Example1_Functions();
void Example1_Geometry();


void Example1_Basic_objects()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                      EXAMPLE 1 - Basic objects                ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Example1_Vector();
	Example1_Matrix();
	//Example1_Tensor();
	//Example1_Functions();
	//Example1_Geometry();
}