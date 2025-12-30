#if !defined MPL_LIENARDWIECHERT_POTENTIAL_H
#define MPL_LIENARDWIECHERT_POTENTIAL_H

#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/VectorTypes.h"

#include "base/BaseUtils.h"

using namespace MML;
using namespace MML::Utils;

namespace MPL
{
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