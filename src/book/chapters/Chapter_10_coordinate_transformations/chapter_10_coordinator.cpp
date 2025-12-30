#include <iostream>

// Forward declarations of demo functions
void ObliqueTransf_transforming_points_2D();
void ObliqueTransf_transforming_vectors_2D();
void ObliqueTransf_transforming_points_3D();
void ObliqueTransf_transforming_vectors_3D();

void SphericalTransf_transforming_points();
void SphericalTransf_transforming_vectors();

void CylindricalTransf_transforming_points();
void CylindricalTransf_transforming_vectors();

void Chapter10_testing_basis_vectors();

void Chapter10_Coordinate_transformations()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 10 - Coordinate transformations      ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	ObliqueTransf_transforming_points_2D();
	ObliqueTransf_transforming_vectors_2D();
	ObliqueTransf_transforming_points_3D();
	ObliqueTransf_transforming_vectors_3D();

	SphericalTransf_transforming_points();
	SphericalTransf_transforming_vectors();
	//CylindricalTransf_transforming_points();
	//CylindricalTransf_transforming_vectors();

	//Chapter10_testing_basis_vectors();
}
