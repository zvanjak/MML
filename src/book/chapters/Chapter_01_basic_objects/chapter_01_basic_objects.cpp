#include "MMLBase.h"

#include <iostream>

void Demo_BaseUtils();
void Demo_Functions();
void Demo_Geometry();
void Demo_Geometry_2D();
void Demo_Geometry_3D();
void Demo_Matrix();
void Demo_MatrixNM();
void Demo_Vector();
void Demo_Tensor();
void Demo_Quaternion();

int Chapter01_basic_objects()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                      CHAPTER 1 - Basic Objects                ****" << std::endl;
	std::cout << "****                                                               ****" << std::endl;
	std::cout << "****    Vectors, Matrices, Geometry, Tensors, and Quaternions      ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo_BaseUtils();

	Demo_Functions();

	Demo_Geometry();
	Demo_Geometry_2D();
	Demo_Geometry_3D();

	Demo_Matrix();
	Demo_MatrixNM();

	Demo_Vector();

	Demo_Tensor();

	Demo_Quaternion();

	return 0;
}

