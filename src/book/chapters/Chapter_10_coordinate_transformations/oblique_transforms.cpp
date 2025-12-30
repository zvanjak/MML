#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Vector.h"
#include "mml/base/VectorTypes.h"

#include "mml/core/CoordTransf/CoordTransf2D.h"
#include "mml/core/CoordTransf/CoordTransf3D.h"
#endif

#include <iostream>

using namespace MML;

// Oblique coordinate transformations in 2D - transforming points
void ObliqueTransf_transforming_points_2D()
{
	std::cout << "\n-----------------------------------------------------------------" << std::endl;
	std::cout << "---------       Transforming points in 2D - OBLIQUE       -------" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;

	Real angle = 85;

	// testing oblique coordinate transformation in 2D
	std::cout << "Testing transformation from Cartesian to oblique coordinates:" << std::endl;
	std::cout << "Angle of oblique coordinates: " << angle << " degrees" << std::endl;
	Vector2Cartesian p1{ 5.0, 2.0 };
	CoordTransfCartesianToOblique2D ctCartToObl(Utils::DegToRad(angle)); 

	auto p1Obl = ctCartToObl.transf(p1); // transform to oblique coordinates
	std::cout << "Cartesian point:   " << p1 << std::endl;
	std::cout << "Oblique point:     " << p1Obl << std::endl;

	// and convert back to Cartesian coordinates
	auto p1Cart = ctCartToObl.transfInverse(p1Obl); // transform back to Cartesian coordinates
	std::cout << "Back to Cartesian: " << p1Cart << std::endl;

	std::cout << "\nTesting inverse transformation:" << std::endl;
	Vector2Cartesian p2Obl{ 3.0, 4.0 };
	CoordTransfObliqueToCartesian2D ctOblToCart(Utils::DegToRad(89.0));

	auto p2Cart = ctOblToCart.transf(p2Obl); // transform to Cartesian coordinates
	std::cout << "Oblique point:     " << p2Obl << std::endl;
	std::cout << "Cartesian point:   " << p2Cart << std::endl;

	// and convert back to oblique coordinates
	auto p2Ob = ctOblToCart.transfInverse(p2Cart); // transform back to oblique coordinates
	std::cout << "Back to oblique:   " << p2Ob << std::endl;
}

// Oblique coordinate transformations in 2D - transforming vectors
void ObliqueTransf_transforming_vectors_2D()
{
	std::cout << "\n-----------------------------------------------------------------" << std::endl;
	std::cout << "--------        Transforming vectors in 2D - OBLIQUE       -------" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;

	Real angle = 85;

	// testing transformation of vectors in 2D
	Vector2Cartesian v1{ 3.0, 4.0 };
	
	CoordTransfCartesianToOblique2D ctCartToObl(Utils::DegToRad(angle));
	
	// contravariant transformation of vector
	auto v1Obl = ctCartToObl.transfVecContravariant(v1, Vector2Cartesian{ 0.0, 0.0 }); // transform to oblique coordinates
	
	v1.PrintLine(std::cout, "Cartesian vector:         ", 12, 7);
	v1Obl.PrintLine(std::cout, "Oblique vector:           ", 12, 7);

	// covariant transformation of vector
	auto v1OblCov = ctCartToObl.transfVecCovariant(v1, Vector2Cartesian{ 0.0, 0.0 }); // transform to oblique coordinates
	v1OblCov.PrintLine(std::cout, "Oblique covariant vector: ", 12, 7);

	std::cout << "\nTransforming vector as point in oblique coordinates:" << std::endl;
	auto p1Obl = ctCartToObl.transf(v1); 
	v1.PrintLine(std::cout, "Cartesian point:          ", 12, 7);
	p1Obl.PrintLine(std::cout, "Oblique point:            ", 12, 7);
}

// Oblique coordinate transformations in 3D - transforming points
void ObliqueTransf_transforming_points_3D()
{
	std::cout << "\n-----------------------------------------------------------------" << std::endl;
	std::cout << "----------     Transforming points in 3D - OBLIQUE       --------" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;

	Vector3Cartesian p1{ 5.0, 2.0, -3.0 };

	Vec3Cart xAxisNew{ 1.0, 0.0, 0.0 }; // new x-axis direction
	Vec3Cart yAxisNew{ 0.0, 1.0, 0.0 }; // new y-axis direction
	Vec3Cart zAxisNew{ 0.1, 0.1, 1.0 }; // new z-axis direction

	std::cout	<< "Testing transformation from Cartesian to oblique coordinates:" << std::endl;
	std::cout	<< "New x-axis direction: " << xAxisNew << std::endl;
	std::cout << "New y-axis direction: " << yAxisNew << std::endl;
	std::cout << "New z-axis direction: " << zAxisNew << std::endl;
	std::cout << std::endl;

	CoordTransfCartesianToOblique3D ctCartToObl(xAxisNew, yAxisNew, zAxisNew);

	auto p1Obl = ctCartToObl.transf(p1); // transform to oblique coordinates
	p1.PrintLine(std::cout, "Cartesian point:          ", 12, 7);
	p1Obl.PrintLine(std::cout, "Oblique point:            ", 12, 7);

	// and convert back to Cartesian coordinates
	auto p1Cart = ctCartToObl.transfInverse(p1Obl); // transform back to Cartesian coordinates
	p1Cart.PrintLine(std::cout, "Back to Cartesian:        ", 12, 7);
}

// Oblique coordinate transformations in 3D - transforming vectors
void ObliqueTransf_transforming_vectors_3D()
{
	std::cout << "\n-----------------------------------------------------------------" << std::endl;
	std::cout << "----------     Transforming vectors in 3D - OBLIQUE      --------" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;

	Vector3Cartesian v1{ 5.0, 2.0, -3.0 };

	Vec3Cart xAxisNew{ 1.0, 0.0, 0.0 }; // new x-axis direction
	Vec3Cart yAxisNew{ 0.0, 1.0, 0.0 }; // new y-axis direction
	Vec3Cart zAxisNew{ 0.1, 0.1, 1.0 }; // new z-axis direction

	CoordTransfCartesianToOblique3D ctCartToObl(xAxisNew, yAxisNew, zAxisNew);

	// contravariant transformation of vector
	auto v1Obl = ctCartToObl.transfVecContravariant(v1, Vec3Cart(0,0,0)); // transform to oblique coordinates
	v1.PrintLine(std::cout, "Cartesian vector:         ", 12, 7);
	v1Obl.PrintLine(std::cout, "Oblique vector:           ", 12, 7);

	// covariant transformation of vector
	auto v1OblCov = ctCartToObl.transfVecCovariant(v1, Vec3Cart(0, 0, 0)); // transform to oblique coordinates
	v1OblCov.PrintLine(std::cout, "Oblique covariant vector: ", 12, 7);
}
