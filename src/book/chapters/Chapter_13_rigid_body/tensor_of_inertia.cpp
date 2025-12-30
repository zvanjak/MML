#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/VectorN.h"
#include "mml/base/Geometry3DBodies.h"

#include "mml/core/CoordTransf.h"
#include "mml/core/CoordTransf/CoordTransf3D.h"

#include "mpl/RigidBody/MomentOfInertiaCalculator.h"
#endif

using namespace MML;
using namespace MPL;

void Example13_tensor_of_inertia_discrete_masses()
{
	// define set of discrete masses
	double a = 1;
	Vector3Cartesian pos1(a, a, 0);
	Vector3Cartesian pos2(-a, a, 0);
	Vector3Cartesian pos3(-a, -a, 0);
	Vector3Cartesian pos4(a, -a, 0);
	Vector3Cartesian pos5(0, 0, 4 * a);
	double m1 = 1;
	double m2 = 1;
	double m3 = 1;
	double m4 = 1;
	double m5 = 1;

	DiscreteMass mass1(pos1, m1);
	DiscreteMass mass2(pos2, m2);
	DiscreteMass mass3(pos3, m3);
	DiscreteMass mass4(pos4, m4);
	DiscreteMass mass5(pos5, m5);

	std::vector<DiscreteMass> masses = { mass1, mass2, mass3, mass4, mass5 };
	DiscreteMassesConfig massesConfig(masses);

	DiscreteMassMomentOfInertiaTensorCalculator calculator(massesConfig);
	Tensor2<3> tensor_orig = calculator.calculate();

	std::cout << "Tensor of inertia: " << std::endl;
	std::cout << tensor_orig << std::endl;

	// investigating what happens if we change coord.system, in two cases:
	// 1. using coord.system transform we calculate TRANSFORMED (original) tensor
	// 2. using coord.system transform we calculate NEW set of masses and then calculate tensor

	// new coord.system is rotated around x axis for 30 degrees
	CoordTransfCart3DRotationXAxis coord_transf(Utils::DegToRad(30.0));

	// 1) - calculated using tensor transformation
	Tensor2<3> tensor_transf = coord_transf.transfTensor2(tensor_orig, Vector3Cartesian(1, 1, 1));

	std::cout << "Tensor of inertia transformed: " << std::endl;
	std::cout << tensor_transf << std::endl;

	// 2) - change masses position and calculate new tensor directly
	DiscreteMassesConfig massesTransConfig(masses);
	for (auto& mass : massesTransConfig._masses)
		mass._position = coord_transf.transf(mass._position);

	DiscreteMassMomentOfInertiaTensorCalculator calculator2(massesTransConfig);
	Tensor2<3> tensor_changed = calculator2.calculate();

	std::cout << "Tensor of inertia rotated masses: " << std::endl;
	std::cout << tensor_changed << std::endl;
}

void Example13_tensor_of_inertia_continuous_mass()
{
	// create ContinuousMass representation for cube of constant density
	SolidBodyWithBoundaryConstDensity cube(Real(-0.5), Real(0.5),
		[](Real x) { return Real(-0.5); },
		[](Real x) { return Real(0.5); },
		[](Real x, Real y) { return Real(-0.5); },
		[](Real x, Real y) { return Real(0.5); },
		1.0);

	ContinuousMassMomentOfInertiaTensorCalculator calculator(cube);

	Tensor2<3> tensor = calculator.calculate();

	std::cout << "Tensor of inertia for cube: " << std::endl;
	std::cout << tensor << std::endl;

	// let's do the same for sphere
	SolidBodyWithBoundaryConstDensity sphere(-1.0, 1.0,
		[](Real x) { return -sqrt(1 - x * x); },
		[](Real x) { return sqrt(1 - x * x); },
		[](Real x, Real y) { return -sqrt(1 - x * x - y * y); },
		[](Real x, Real y) { return sqrt(1 - x * x - y * y); },
		1.0);

	ContinuousMassMomentOfInertiaTensorCalculator calculator2(sphere);

	Tensor2<3> tensor2 = calculator2.calculate();

	std::cout << "Tensor of inertia for sphere: " << std::endl;
	std::cout << tensor2 << std::endl;

	std::cout << "Theoretical value for sphere: " << std::endl;
	std::cout << 2.0 / 5.0 * (4.0 / 3.0 * Constants::PI) << std::endl;
}

void Example13_tensor_of_inertia()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 13 - tensor of inertia               ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Example13_tensor_of_inertia_discrete_masses();
	Example13_tensor_of_inertia_continuous_mass();
}