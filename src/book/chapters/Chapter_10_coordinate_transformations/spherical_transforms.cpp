#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector.h"
#include "mml/base/VectorTypes.h"

#include "mml/core/CoordTransf.h"
#include "mml/core/CoordTransf/CoordTransfSpherical.h"

#include "mml/core/Derivation.h"
#include "mml/core/FieldOperations.h"
#include "mml/core/Fields.h"
#endif

#include <iostream>

using namespace MML;

// Spherical coordinate transformations - transforming points
void SphericalTransf_transforming_points()
{
	std::cout << "\n***********************************************************************\n";
	std::cout << "****              SPHERICAL COORDINATE TRANSFORMATION              ****\n";
	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Point (position) transformation - direct \n";

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

// Spherical coordinate transformations - transforming vectors
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
