#if !defined MPL_EM_TENSOR_H
#define MPL_EM_TENSOR_H

#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#include <cmath>

using namespace MML;

namespace MPL
{
	inline Tensor2<4> GetEMTensorContravariant(const Vector3Cartesian &E_field, const Vector3Cartesian &B_field)
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

	inline Tensor2<4> GetEMTensorCovariant(const Vector3Cartesian &E_field, const Vector3Cartesian &B_field)
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

	/// @brief Extract E and B field vectors from the contravariant EM tensor F^μν
	/// @param F_tensor The contravariant electromagnetic field tensor
	/// @param E_field Output: Electric field vector (in units where c=1)
	/// @param B_field Output: Magnetic field vector
	/// 
	/// The contravariant EM tensor has the form:
	///        0      -Ex/c   -Ey/c   -Ez/c
	///       Ex/c     0      -Bz      By
	///       Ey/c     Bz      0      -Bx
	///       Ez/c    -By      Bx      0
	///
	/// This is the inverse operation of GetEMTensorContravariant()
	inline void GetEandBFromEMTensorContravariant(const Tensor2<4>& F_tensor, 
	                                       Vector3Cartesian& E_field, 
	                                       Vector3Cartesian& B_field)
	{
		double c = 1.0;
		
		// Extract E field from F^0i components (with sign flip)
		// F^01 = -Ex/c, F^02 = -Ey/c, F^03 = -Ez/c
		E_field = Vector3Cartesian(
			-F_tensor(0, 1) * c,
			-F_tensor(0, 2) * c,
			-F_tensor(0, 3) * c
		);
		
		// Extract B field from spatial components
		// F^12 = -Bz, F^13 = By, F^23 = -Bx
		B_field = Vector3Cartesian(
			-F_tensor(2, 3),   // -Bx at (2,3)
			 F_tensor(1, 3),   //  By at (1,3)
			-F_tensor(1, 2)    // -Bz at (1,2)
		);
	}

	/// @brief Extract E and B field vectors from the covariant EM tensor F_μν
	/// @param F_tensor The covariant electromagnetic field tensor
	/// @param E_field Output: Electric field vector (in units where c=1)
	/// @param B_field Output: Magnetic field vector
	///
	/// The covariant EM tensor has the form:
	///        0       Ex/c    Ey/c    Ez/c
	///      -Ex/c     0      -Bz      By
	///      -Ey/c     Bz      0      -Bx
	///      -Ez/c    -By      Bx      0
	inline void GetEandBFromEMTensorCovariant(const Tensor2<4>& F_tensor, 
	                                   Vector3Cartesian& E_field, 
	                                   Vector3Cartesian& B_field)
	{
		double c = 1.0;
		
		// Extract E field from F_0i components (no sign flip for covariant)
		// F_01 = Ex/c, F_02 = Ey/c, F_03 = Ez/c
		E_field = Vector3Cartesian(
			F_tensor(0, 1) * c,
			F_tensor(0, 2) * c,
			F_tensor(0, 3) * c
		);
		
		// Extract B field from spatial components (same as contravariant)
		// F_12 = -Bz, F_13 = By, F_23 = -Bx
		B_field = Vector3Cartesian(
			-F_tensor(2, 3),
			 F_tensor(1, 3),
			-F_tensor(1, 2)
		);
	}

	//////////////////////////////////////////////////////////////////////////
	/// MOVING CHARGE FIELDS - Heaviside-Feynman Formula
	//////////////////////////////////////////////////////////////////////////

	/// @brief Calculate E and B fields of a uniformly moving point charge
	/// 
	/// This is the Heaviside-Feynman formula for the fields of a charge
	/// moving with constant velocity. The field lines are compressed
	/// ("pancaked") in the direction of motion due to Lorentz contraction.
	///
	/// @param chargePos Current position of the charge (at time t)
	/// @param fieldPoint Point where we want to calculate the field
	/// @param velocity Velocity vector of the charge (in units of c, |v| < 1)
	/// @param charge The charge q (in appropriate units, use q/(4πε₀) for SI)
	/// @param E_field Output: Electric field at fieldPoint
	/// @param B_field Output: Magnetic field at fieldPoint
	///
	/// Formulas (with c=1):
	///   E = q * (1 - β²) / (1 - β²sin²θ)^(3/2) * r̂ / r²
	///   B = v × E  (in units where c=1)
	///
	/// where θ is the angle between v and r (from charge to field point)
	inline void GetFieldsOfMovingCharge(
		const Vector3Cartesian& chargePos,
		const Vector3Cartesian& fieldPoint,
		const Vector3Cartesian& velocity,
		Real charge,
		Vector3Cartesian& E_field,
		Vector3Cartesian& B_field)
	{
		// Vector from charge to field point
		Vector3Cartesian r_vec(
			fieldPoint.X() - chargePos.X(),
			fieldPoint.Y() - chargePos.Y(),
			fieldPoint.Z() - chargePos.Z()
		);
		
		Real r = r_vec.NormL2();
		if (r < 1e-15) {
			E_field = Vector3Cartesian(0, 0, 0);
			B_field = Vector3Cartesian(0, 0, 0);
			return;
		}
		
		// Unit vector r̂
		Vector3Cartesian r_hat = r_vec * (1.0 / r);
		
		// β = v/c (already in units of c)
		Real beta = velocity.NormL2();
		Real beta2 = beta * beta;
		
		if (beta < 1e-15) {
			// Static case: Coulomb field
			E_field = r_hat * (charge / (r * r));
			B_field = Vector3Cartesian(0, 0, 0);
			return;
		}
		
		// cos(θ) = v̂ · r̂
		Vector3Cartesian v_hat = velocity * (1.0 / beta);
		Real cosTheta = ScalarProduct(v_hat, r_hat);
		Real sin2Theta = 1.0 - cosTheta * cosTheta;
		
		// Relativistic factor: (1 - β²) / (1 - β²sin²θ)^(3/2)
		Real denominator = 1.0 - beta2 * sin2Theta;
		Real factor = (1.0 - beta2) / (denominator * std::sqrt(denominator));
		
		// Electric field: E = q * factor * r̂ / r²
		E_field = r_hat * (charge * factor / (r * r));
		
		// Magnetic field: B = v × E (with c=1)
		B_field = VectorProduct(velocity, E_field);
	}

	/// @brief Simplified version for charge moving along x-axis
	/// @param fieldPoint Point where we want the field (charge at origin)
	/// @param beta Velocity as fraction of c (v/c), positive = +x direction
	/// @param charge The charge q
	/// @param E_field Output: Electric field
	/// @param B_field Output: Magnetic field
	inline void GetFieldsOfMovingChargeAlongX(
		const Vector3Cartesian& fieldPoint,
		Real beta,
		Real charge,
		Vector3Cartesian& E_field,
		Vector3Cartesian& B_field)
	{
		Vector3Cartesian chargePos(0, 0, 0);
		Vector3Cartesian velocity(beta, 0, 0);
		GetFieldsOfMovingCharge(chargePos, fieldPoint, velocity, charge, E_field, B_field);
	}
}

#endif // MPL_EM_TENSOR_H