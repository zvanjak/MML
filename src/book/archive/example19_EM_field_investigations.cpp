#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/VectorTypes.h"
#include "base/Tensor.h"

#include "core/Derivation.h"

#include "mpl/Electromagnetism/EMTensor.h"
#include "mpl/Electromagnetism//LienardWiechertPotential.h"

#include "mpl/SpecialRelativity/LorentzTransformation.h"
#endif


using namespace MML;
using namespace MPL;


// Verify Gauss' law for a point charge
// integracija po liniji, calculate work done by electric field

// Biot-Savart law - calculate exact magnetic field for a wire loop

// Given EM field, calculate test charge trajectory, field lines

void Example19_EM_field_investigations()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 19 - EM field investigations        ****" << std::endl;
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
	
	Vector4Minkowski point4({point.X(), point.Y(), point.Z(), 0.0});

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
// Biot-Savartov zakon ... za svaku tocku integrirati po cijeloj liniji zavojnice

// TREBA MI VIZUALIZACIJA SILNICA b POLJA!!!
// dobro je što su zatvorene petlje