// Chapter 14 main.cpp

#include <iostream>
#include <exception>

void Chapter14_rotations_quaternions();

int main()
{
	try
	{
		Chapter14_rotations_quaternions();
		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
