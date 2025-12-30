#include <iostream>
#include <exception>

void Chapter11_Carousel();

int main()
{
	try
	{
		Chapter11_Carousel();
		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
