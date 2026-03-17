// Chapter 22 main.cpp

#include <iostream>
#include <exception>

void Chapter22_General_relativity();

int main()
{
	try
	{
		Chapter22_General_relativity();
		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
