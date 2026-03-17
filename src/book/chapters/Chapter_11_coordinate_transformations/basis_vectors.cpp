#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Vector/VectorN.h"

#include "mml/core/CoordTransf.h"
#include "mml/core/CoordTransf/CoordTransfSpherical.h"
#endif

using namespace MML;

// Calculate basis vectors for different coordinate systems
void Chapter11_testing_basis_vectors()
{
	// calculate basis vector for spherical coordinates
	VectorN<Real, 3> pos{ 1.0, Utils::DegToRad(45), Utils::DegToRad(30) }; // r, theta, phi

	VectorN<Real, 3> basisVec1x = CoordTransfSpherToCart.getBasisVec(0, pos); // basis vector for r

	//std::cout << "Basis vector for r: " << basisVec1x << std::endl;
}
