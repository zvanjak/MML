/******************************************************************************
 * MML Example: Vector Spaces
 * ==========================
 *
 * Demonstrates the finite-dimensional vector-space layer:
 *
 *   1. basis changes keep VectorN as coordinate storage;
 *   2. subspaces expose semantic projection and membership;
 *   3. linear maps carry domain/codomain dimensions;
 *   4. dual spaces evaluate vectors through covectors;
 *   5. inner-product spaces own Gram-dependent orthogonality;
 *   6. affine spaces distinguish points from vectors.
 *
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/core/VectorSpaces.h>
#endif

#include <iomanip>
#include <iostream>

using namespace MML;
using namespace MML::VectorSpaces;

template<class VectorLike>
void PrintVector(const char* label, const VectorLike& vector)
{
	std::cout << "  " << std::left << std::setw(28) << label << " = (";
	for (int i = 0; i < vector.size(); i++) {
		if (i > 0)
			std::cout << ", ";
		std::cout << std::right << std::setw(9) << vector[i];
	}
	std::cout << ")\n";
}

int main()
{
	std::cout << std::fixed << std::setprecision(6);
	std::cout << "Vector spaces: bases, maps, duals, metrics, and affine frames\n\n";

	Basis<Real, 2> scaled({
		VectorN<Real, 2>{ REAL(2.0), REAL(0.0) },
		VectorN<Real, 2>{ REAL(0.0), REAL(3.0) }
	});
	VectorN<Real, 2> scaledCoords{ REAL(4.0), REAL(5.0) };
	VectorN<Real, 2> standardCoords = scaled.coordinatesInStandard(scaledCoords);
	PrintVector("scaled coordinates", scaledCoords);
	PrintVector("standard coordinates", standardCoords);

	Matrix<Real> xAxisBasis(2, 1);
	xAxisBasis(0, 0) = REAL(1.0);
	Subspace<Real> xAxis(2, xAxisBasis);
	Vector<Real> sample{ REAL(3.0), REAL(4.0) };
	Vector<Real> projection = xAxis.project(sample);
	PrintVector("projection onto x-axis", projection);

	MatrixNM<Real, 2, 2> rotation{
		REAL(0.0), -REAL(1.0),
		REAL(1.0), REAL(0.0)
	};
	LinearMap<Real, 2, 2> rotate90(rotation);
	PrintVector("rotated vector", rotate90(VectorN<Real, 2>{ REAL(1.0), REAL(0.0) }));

	DualSpace<Real, 2> dual;
	CovectorInSpace<Real, 2> alpha = dual.covector(VectorN<Real, 2>{ REAL(2.0), REAL(3.0) });
	VectorInSpace<Real, 2> v(VectorN<Real, 2>{ REAL(4.0), REAL(5.0) });
	std::cout << "  " << std::left << std::setw(28) << "alpha(v)" << " = " << alpha(v) << "\n";

	MatrixNM<Real, 2, 2> gram{
		REAL(2.0), REAL(0.0),
		REAL(0.0), REAL(1.0)
	};
	InnerProductSpace<Real, 2> weightedPlane(gram);
	std::cout << "  " << std::left << std::setw(28) << "weighted norm" << " = "
		<< weightedPlane.norm(VectorN<Real, 2>{ REAL(3.0), REAL(4.0) }) << "\n";

	AffineFrame<Real, 2> frame(
		PointInSpace<Real, 2>(VectorN<Real, 2>{ REAL(10.0), REAL(20.0) }),
		scaled);
	PointInSpace<Real, 2> point = frame.pointFromCoordinates(VectorN<Real, 2>{ REAL(1.0), REAL(1.0) });
	AffineMap<Real, 2, 2> translate = AffineMap<Real, 2, 2>::Translation(VectorN<Real, 2>{ REAL(5.0), REAL(6.0) });
	PointInSpace<Real, 2> moved = translate(point);
	PrintVector("point in standard coords", point.coordinates());
	PrintVector("translated point", moved.coordinates());

	return 0;
}