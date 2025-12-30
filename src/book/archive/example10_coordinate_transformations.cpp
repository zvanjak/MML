#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "base/Vector.h"
#include "base/VectorTypes.h"

#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/CoordTransf/CoordTransf2D.h"
#include "core/CoordTransf/CoordTransf3D.h"

#include "core/Derivation.h"
#include "core/FieldOperations.h"
#include "core/Fields.h"
#endif

using namespace MML;

// active vs passive transformation?

// transformacija tocke
// covar i contravar bazni vektori u tocki
// kako dobiti unit vektore
// metricki tenzor u tocki
// transformacija vektora

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

void ObliqueTransf_transforming_vectors_2D()
{
	std::cout << "\n-----------------------------------------------------------------" << std::endl;
	std::cout << "--------        Transforming vectors in 2D - OBLIUE       -------" << std::endl;
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

void SphericalTransf_transforming_points()
{
	std::cout << "\n***********************************************************************\n";
	std::cout << "****              SPHERICAL COORDINATE TRANSFORMATION              ****\n";
	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Point (position) transformation - direct \n";

	//	Vector3Cartesian p1{ 5.0, 2.0, -3.0 };
	Vector3Cartesian p1{ 2.0, 2.0, 1.0 };
	auto p1Spher = CoordTransfCartToSpher.transf(p1);
	auto p1BackTransf = CoordTransfSpherToCart.transf(p1Spher);

	std::cout << "Cartesian   : " << p1 << std::endl;
	std::cout << "Spherical   : " << p1Spher << std::endl;
	std::cout << "Back transf : " << p1BackTransf << std::endl;

	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Point (position) transformation - inverse \n";

	Vector3Cartesian p2{ 5.0, 2.0, -3.0 };
	auto p2Spher = CoordTransfSpherToCart.transfInverse(p2);
	auto p2BackTransf = CoordTransfCartToSpher.transfInverse(p2Spher);

	std::cout << "Cartesian   : " << p2 << std::endl;
	std::cout << "Spherical   : " << p2Spher << std::endl;
	std::cout << "Back transf : " << p2BackTransf << std::endl;
}

void SphericalTransf_transforming_vectors()
{
	std::cout << "\n-----------------------------------------------------------------------\n";
	std::cout << "Vector transformation - contravariant\n";

	Vector3Cartesian p1{ 2.0, 2.0, 1.0 };
	auto p1Spher = CoordTransfCartToSpher.transf(p1);

	Vec3Cart v1{ 2.0, 2.0, 0.0 };
	Vec3Cart p2_1{ 2.0, 2.0, 1.0 }, p2_2{ -2.0, -2.0, -1.0 };

	Vector<Vec3Cart> pntList{ Vec3Cart{2.0, 2.0, 1.0}, Vec3Cart{-2.0, -2.0, -1.0},
														Vec3Cart{2.0, -2.0, 1.0}, Vec3Cart{-2.0, 2.0, -1.0},
														Vec3Cart{2.0, 4.0, 2.0}, Vec3Cart{-2.0, -4.0, -2.0}
	};
	std::cout << "Vector (contravariant)       : " << v1 << std::endl;
	for (const auto& pnt : pntList)
	{
		auto v1Spher = CoordTransfCartToSpher.transfVecContravariant(v1, pnt);
		std::cout << "contravar transf. to spher at cart.pnt : " << pnt << " = " << v1Spher << std::endl;
	}
	//auto v2_1Spher = CoordTransfCartToSpher.transfVecContravariant(v1, p2_1);
	//std::cout << "contravar transf. to spher at cart.pnt : " << p2_1 << " = " << v2_1Spher << std::endl;
	//auto v2_2Spher = CoordTransfCartToSpher.transfVecContravariant(v1, p2_2);
	//std::cout << "contravar transf. to spher at cart.pnt : " << p2_2 << " = " << v2_2Spher << std::endl;

	std::cout << "\n-----------------------------------------------------------------------\n";
	std::cout << "Vector transformation - covariant\n";

	auto v2Spher = CoordTransfCartToSpher.transfVecCovariant(v1, p1Spher);
	std::cout << "covariant transf. to spher at cart.pnt : " << p1 << " = " << v2Spher << std::endl;

	std::cout << "Point (cartesian): " << p1 << std::endl;
	std::cout << "Point (spherical): " << p1Spher << std::endl;

	// calculate gradient in Cartesian system at given point
	ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);
	Vector3Cartesian grad_cart = ScalarFieldOperations::GradientCart<3>(fPotCart, p1);
	std::cout << "Gradient in Cartesian coords:" << grad_cart << std::endl;

	// calculate gradient in Spherical system at given point
	ScalarFunction<3> fPotSpher(Fields::InverseRadialPotentialFieldSpher);
	Vector3Spherical grad_spher = ScalarFieldOperations::GradientSpher(fPotSpher, p1Spher);
	std::cout << "Gradient in Spherical coords:" << grad_spher << std::endl;

	auto grad_cart_transf = CoordTransfSpherToCart.transfVecCovariant(grad_spher, p1);
	std::cout << "Cartesian gradient from spherical (covariant transf.):" << grad_cart_transf << std::endl;

	auto grad_spher_transf = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p1Spher);
	std::cout << "Spherical gradient from Cartesian (covariant transf.):" << grad_spher_transf << std::endl;

	auto grad_cart_transf_inv = CoordTransfCartToSpher.transfInverseVecCovariant(grad_spher, p1);
	std::cout << "Cartesian gradient from spherical (inverse covariant transf.):" << grad_cart_transf_inv << std::endl;

	auto grad_spher_transf_inv = CoordTransfSpherToCart.transfInverseVecCovariant(grad_cart, p1Spher);
	std::cout << "Spherical gradient from spherical (inverse covariant transf.):" << grad_spher_transf_inv << std::endl;
}

void CylindricalTransf_transforming_points()
{
}

void CylindricalTransf_transforming_vectors()
{
}

void Example10_testing_basis_vectors()
{
	// calculate basis vector for spherical coordinates

	VectorN<Real, 3> pos{ 1.0, Utils::DegToRad(45), Utils::DegToRad(30) }; // r, theta, phi

	VectorN<Real, 3> basisVec1x = CoordTransfSpherToCart.getBasisVec(0, pos); // basis vector for r

	//std::cout << "Basis vector for r: " << basisVec1x << std::endl;
}

void Example10_transforming_points_2D()
{
	// all available 2D coordinate transformations
	CoordTransfPolarToCartesian2D ctPolarToCart2D;
	CoordTransfCart2DRotation ctCart2DRot(Utils::DegToRad(30.0)); // 30 degrees rotation
	
	// oblique coordinate transformation in 2D
	CoordTransfObliqueToCartesian2D ctOblToCart(Utils::DegToRad(89.0));
}

void Example10_transforming_points_3D()
{
	// all available 3D coordinate transformations
	CoordTransfCart3DRotationXAxis ctCart3DRotX(Utils::DegToRad(30.0)); // 30 degrees rotation around X-axis
	CoordTransfCart3DRotationYAxis ctCart3DRotY(Utils::DegToRad(30.0)); // 30 degrees rotation around Y-axis
	CoordTransfCart3DRotationZAxis ctCart3DRotZ(Utils::DegToRad(30.0)); // 30 degrees rotation around Z-axis

	CoordTransf3DCartOrthogonal ctCartOrthogonal(
		Vec3Cart{ 1.0, 0.0, 0.0 }, // new X-axis
		Vec3Cart{ 0.0, 1.0, 0.0 }, // new Y-axis
		Vec3Cart{ 0.0, 0.0, 1.0 }  // new Z-axis
	);

	CoordTransf3DCartGeneral ctCartGeneral(MatrixNM<Real, 3, 3> {
																						{ 1.0, 0.0, 0.0 }, // new X-axis
																						{ 0.0, 1.0, 0.0 }, // new Y-axis
																						{ 0.0, 0.0, 1.0 }  // new Z-axis
																					}
																				);

	// oblique coordinate transformation in 3D
	CoordTransfCartesianToOblique3D ctCartToObl(
		Vec3Cart{ 1.0, 0.0, 0.0 }, // new X-axis
		Vec3Cart{ 0.0, 1.0, 0.0 }, // new Y-axis
		Vec3Cart{ 0.1, 0.1, 1.0 }  // new Z-axis, slightly tilted
	);

	CoordTransfSphericalToCartesian      ctSpherToCart;
	CoordTransfCylindricalToCartesian    ctCylToCart;

	// start with Cartesian coordinates
	Vector3Cartesian p1{ 2.0, 2.0, 1.0 };

	// and use INVERSE transformations to convert from Cart. to given coord system


	// and convert Cartesian point in all of them
}


void Example10_Coordinate_transformations()
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
	CylindricalTransf_transforming_points();
	CylindricalTransf_transforming_vectors();

	Example10_testing_basis_vectors();

	Example10_transforming_points_2D();
	Example10_transforming_points_3D();
}