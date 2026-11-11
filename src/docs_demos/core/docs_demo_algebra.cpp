#include <mml/base/Algebra_base.h>
#include <mml/algorithms/MatrixAlg/FieldMatrixAlgorithms.h>

#include <iostream>

using namespace MML;
using namespace MML::Algebra;

void Docs_Demo_Algebra()
{
	const SO2 planar(Constants::PI / Real(2));
	const auto rotation = SO3::FromAxisAngle({0, 0, 1}, Constants::PI / Real(2));
	const auto rotated = rotation.apply({1, 0, 0});
	const SE2 planarMotion(planar, {2, 3});
	const SE3 transform(rotation, {2, 3, 4});
	const auto movedPoint = transform.apply_point(VectorN<Real, 3>{1, 0, 0});
	const auto movedVector = transform.apply_vector(VectorN<Real, 3>{1, 0, 0});
	std::cout << "SO3 quarter-turn of e1: " << rotated << '\n';
	std::cout << "SE2 homogeneous matrix:\n" << planarMotion.homogeneous_matrix();
	std::cout << "SE3 point/vector: " << movedPoint << " / " << movedVector << '\n';

	const MatrixNM<Real, 2, 2> matrix{{0, 2}, {1, 3}};
	std::cout << "Generic field-matrix determinant: " << FieldMatrixDeterminant(matrix)
	          << "\nInverse:\n" << FieldMatrixInverse(matrix);
}
