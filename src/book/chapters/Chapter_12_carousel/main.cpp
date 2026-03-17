#include <iostream>
#include <exception>

void Chapter12_Carousel();

int main()
{
	try
	{
		Chapter12_Carousel();
		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
