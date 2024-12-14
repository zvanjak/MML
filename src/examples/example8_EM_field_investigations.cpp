#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/VectorTypes.h"
#include "base/Tensor.h"

#include "core/Derivation.h"
#include "core/CoordTransf/CoordTransfLorentz.h"
#endif


using namespace MML;


// Verify Gauss' law for a point charge
// integracija po liniji, calculate work done by electric field

// Biot-Savart law - calculate exact magnetic field for a wire loop

// Given EM field, calculate test charge trajectory, field lines

Tensor2<4> GetEMTensorContravariant(Vector3Cartesian E_field, Vector3Cartesian B_field)
{
	double c = 1.0;

	Tensor2<4> EM_tensor(2, 0);

	EM_tensor(0, 0) = 0.0;
	EM_tensor(0, 1) = -E_field.X() / c;
	EM_tensor(0, 2) = -E_field.Y() / c;
	EM_tensor(0, 3) = -E_field.Z() / c;

	EM_tensor(1, 0) = E_field.X() / c;
	EM_tensor(1, 1) = 0.0;
	EM_tensor(1, 2) = -B_field.Z();
	EM_tensor(1, 3) = B_field.Y();

	EM_tensor(2, 0) = E_field.Y() / c;
	EM_tensor(2, 1) = B_field.Z();
	EM_tensor(2, 2) = 0.0;
	EM_tensor(2, 3) = -B_field.X();

	EM_tensor(3, 0) = E_field.Z() / c;
	EM_tensor(3, 1) = -B_field.Y();
	EM_tensor(3, 2) = B_field.X();
	EM_tensor(3, 3) = 0.0;

	return EM_tensor;
}

Tensor2<4> GetEMTensorCovariant(Vector3Cartesian E_field, Vector3Cartesian B_field)
{
	double c = 1.0;

	Tensor2<4> EM_tensor(2, 0);

	EM_tensor(0, 0) = 0.0;
	EM_tensor(0, 1) = E_field.X() / c;
	EM_tensor(0, 2) = E_field.Y() / c;
	EM_tensor(0, 3) = E_field.Z() / c;

	EM_tensor(1, 0) = -E_field.X() / c;
	EM_tensor(1, 1) = 0.0;
	EM_tensor(1, 2) = -B_field.Z();   // check this!!! (this is what Wikipedia and Student Guide to Vectors and Tensors say - contrary to Github Copilot!)
	EM_tensor(1, 3) = B_field.Y();    // check this!!!

	EM_tensor(2, 0) = -E_field.Y() / c;
	EM_tensor(2, 1) = B_field.Z();    // check this!!!
	EM_tensor(2, 2) = 0.0;
	EM_tensor(2, 3) = -B_field.X();   // check this!!!

	EM_tensor(3, 0) = -E_field.Z() / c;
	EM_tensor(3, 1) = -B_field.Y();   // check this!!!
	EM_tensor(3, 2) = B_field.X();    // check this!!!
	EM_tensor(3, 3) = 0.0;

	return EM_tensor;
}

// Given t and r in lab frame of the observer, calculate Lienard-Wiechert potentials for charge moving along x-axis with given velocity
Real calcLienardWiechertScalarPotential(Vector3Cartesian r_at_point, Real t, Vector3Cartesian rs_charge_pos, Real q, Real charge_velocity)
{
	double c = 3e8;
	Vec3Cart charge_v(charge_velocity, 0, 0);
	
	// calculate retarded time
	Real r = (r_at_point - rs_charge_pos).NormL2();
	Real tr = t - (r_at_point - rs_charge_pos).NormL2() / c;
	
	// calculate retarded position
	Vec3Cart r_ret_charge_pos = r_at_point - charge_v * (t - tr);

	Real beta = charge_v.NormL2() / c;
	Vec3Cart ns = (r_at_point - r_ret_charge_pos) / (r_at_point - r_ret_charge_pos).NormL2();

	Real scalar_potential = q / ((r_at_point - r_ret_charge_pos).NormL2() * (1 - beta * ns.ScalarProductCartesian(charge_v) / c));

	return scalar_potential;
}
void Example8_EM_field_investigations()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 8 - EM field investigations         ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// 1) stationary charge at origin, and second observer moving along x-axis with given velocity
	Vector3Cartesian E_field;
	Vector3Cartesian B_field(0.0, 0.0, 0.0);

	Tensor2<4> EM_tensor(2,0);

	// define point in space
	Vector3Cartesian point(5.0, 3.0, -2.0);

	// calculate E field at that point of stationary charge
	Vector3Cartesian E_field_stationary = point / (point.NormL2() * point.NormL2() * point.NormL2());

	// form EM tensor at that point
	EM_tensor = GetEMTensorContravariant(E_field_stationary, B_field);

	EM_tensor.Print(std::cout, 12, 5);

	// transform EM tensor to moving observer
	CoordTransfLorentzXAxis lorentzTransf(0.5);
	
	Vector4Lorentz point4({point.X(), point.Y(), point.Z(), 0.0});

	auto transfEMtensor = lorentzTransf.transfTensor2(EM_tensor, point4);
	
	transfEMtensor.Print(std::cout, 12, 5);


	// calculate E and B fields for moving observer using retarded Lienard-Wiechert potentials



}

/*
  VERIFY EM TENSOR!!!!
  1) stationary charge at origin, and second observer moving along x-axis with given velocity
		- calculate E and B fields, then calculate EM tensor at couple of points
		- transform EM tensor to moving observer
		- calculate E and B fields for moving observer
		- compare with direct calculation for E & M fields of moving charge
	2) infinite line current along x-axis, and second observer moving along x-axis with given velocity
		- calculate E and B fields, then calculate EM tensor at couple of points
		- transform EM tensor to moving observer
		- calculate E and B fields for moving observer
		- compare with direct calculation? (is it possible?)
*/

// ispalis elektron u kompleksno EM polje

// egzaktno magnetsko polje zavojnice, s N navoja, promjerom D, i strujom I
// Biot-Savartov zakon ... za svaku točku integrirati po cijeloj liniji zavojnice

// TREBA MI VIZUALIZACIJA SILNICA b POLJA!!!
// dobro je što su zatvorene petlje