#include <iostream>

// Forward declarations of demo functions
void Demo_TwoMasses2D();
void Demo_Check_Kepler_laws();
void Demo_TwoMasses3D();
void Demo_Field_operations();
void Verify_path_integrals();
void Demo_Voyager_Jupiter_Slingshot();

void Chapter08_Gravity()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                       EXAMPLE 8 - Gravity                     ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo_TwoMasses2D();
	//Demo_Check_Kepler_laws();
	//Demo_TwoMasses3D();
	//Demo_Field_operations();
	//Verify_path_integrals();
	//Demo_Voyager_Jupiter_Slingshot();
}
