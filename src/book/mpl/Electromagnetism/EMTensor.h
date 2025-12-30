#if !defined MPL_EM_TENSOR_H
#define MPL_EM_TENSOR_H

#include "MMLBase.h"


using namespace MML;

namespace MPL
{
	Tensor2<4> GetEMTensorContravariant(const Vector3Cartesian &E_field, const Vector3Cartesian &B_field)
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

	Tensor2<4> GetEMTensorCovariant(const Vector3Cartesian &E_field, const Vector3Cartesian &B_field)
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
}

#endif // MPL_EM_TENSOR_H