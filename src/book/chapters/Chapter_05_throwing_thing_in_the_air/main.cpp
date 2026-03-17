#include <iostream>
#include <exception>

void Chapter5_Throwing_things_in_the_air();

int main()
{
	try
	{
		Chapter5_Throwing_things_in_the_air();
		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
