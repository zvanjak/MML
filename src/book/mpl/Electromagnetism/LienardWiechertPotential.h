#if !defined MPL_LIENARDWIECHERT_POTENTIAL_H
#define MPL_LIENARDWIECHERT_POTENTIAL_H

#include "MMLBase.h"

#include "base/Vector/VectorN.h"
#include "base/Vector/VectorTypes.h"

#include "base/BaseUtils.h"

using namespace MML;
using namespace MML::Utils;

namespace MPL
{
	///////////////////////////////////////////////////////////////////////////
	///                   LIÉNARD-WIECHERT POTENTIALS & FIELDS              ///
	///////////////////////////////////////////////////////////////////////////
	///
	/// The Liénard-Wiechert formulas give the exact EM fields of a point charge
	/// in arbitrary motion. For UNIFORM velocity, they reduce to the same result
	/// as the Heaviside-Feynman formula.
	///
	/// Key insight: The fields depend on the RETARDED position and velocity -
	/// where the charge WAS when the "information" left to reach the field point.
	///
	/// For uniform velocity, the retarded position can be computed exactly,
	/// and the formula simplifies considerably.
	///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Calculate Liénard-Wiechert E and B fields for a uniformly moving charge
	/// 
	/// For uniform velocity, the general L-W formula simplifies to the same result
	/// as the Heaviside-Feynman formula:
	///   E = q(1-β²) / [r²(1 - β²sin²θ)^(3/2)] * r̂
	///   B = β × E
	/// 
	/// where r̂ points from the PRESENT position of charge to field point
	/// and θ is the angle between v and r.
	///
	/// This implementation uses the Heaviside-Feynman form to ensure
	/// consistency, but represents the same physics as L-W for uniform motion.
	///
	/// @param chargePos Current position of the charge
	/// @param fieldPoint Point where we want the field
	/// @param velocity Velocity vector (in units where c=1)
	/// @param charge The charge q
	/// @param E_field Output: Electric field
	/// @param B_field Output: Magnetic field
	void GetLienardWiechertFieldsUniformVelocity(
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
		
		// β = v (in units of c=1)
		Real beta = velocity.NormL2();
		Real beta2 = beta * beta;
		
		if (beta < 1e-15) {
			// Static case: Coulomb field
			E_field = r_hat * (charge / (r * r));
			B_field = Vector3Cartesian(0, 0, 0);
			return;
		}
		
		// For uniform velocity, the L-W formula reduces to Heaviside-Feynman:
		// E = q(1-β²) / [r²(1-β²sin²θ)^(3/2)] * r̂
		
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
		// Equivalently: B = β × E where β = v/c
		B_field = VectorProduct(velocity, E_field);
	}

	/// @brief Simplified version for charge at origin moving along x-axis
	void GetLienardWiechertFieldsAlongX(
		const Vector3Cartesian& fieldPoint,
		Real beta,
		Real charge,
		Vector3Cartesian& E_field,
		Vector3Cartesian& B_field)
	{
		Vector3Cartesian chargePos(0, 0, 0);
		Vector3Cartesian velocity(beta, 0, 0);
		GetLienardWiechertFieldsUniformVelocity(chargePos, fieldPoint, velocity, charge, E_field, B_field);
	}

	// Legacy function - kept for compatibility
	// Given t and r in lab frame of the observer, calculate Lienard-Wiechert potentials 
	// for charge moving along x-axis with given velocity
	Real calcLienardWiechertScalarPotential(Vector3Cartesian r_at_point, Real t, Vector3Cartesian rs_charge_pos,
																					Real q, Real charge_velocity)
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

		Real scalar_potential = q / ((r_at_point - r_ret_charge_pos).NormL2() * (1 - beta * ScalarProduct(ns, charge_v) / c));

		return scalar_potential;
	}
}

#endif // MPL_LIENARDWIECHERT_POTENTIAL_H