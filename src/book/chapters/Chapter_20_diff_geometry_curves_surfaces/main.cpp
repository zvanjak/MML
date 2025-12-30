#include <iostream>
#include <exception>

void Chapter20_Diff_geometry_curves_surfaces();

int main()
{
	try
	{
		Chapter20_Diff_geometry_curves_surfaces();
		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
