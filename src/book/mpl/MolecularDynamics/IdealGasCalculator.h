#if !defined MPL_IDEAL_GAS_CALCULATOR_H
#define MPL_IDEAL_GAS_CALCULATOR_H

#include "MMLBase.h"

#include "../Base/PhysicalConstants.h"

using namespace MML;

namespace MPL
{
	class IdealGasCalculator
	{
		// Ideal Gas Law: PV = nRT
	public:
		static Real NumOfMoleculesInMoles(Real Pressure_Pa, Real Volume_m3, Real Temp_K)
		{
			return (Pressure_Pa * Volume_m3) / (PhyConst::MolarGasConstant * Temp_K);
		}
		static Real NumOfMolecules(Real Pressure_Pa, Real Volume_m3, Real Temp_K)
		{
			return NumOfMoleculesInMoles(Pressure_Pa, Volume_m3, Temp_K) * PhyConst::AvogardoConstant;
		}

		static Real Pressure(Real n_mole, Real Volume_m3, Real Temp_K)
		{
			// Ideal Gas Law: PV = nRT => P = nRT / V
			return (n_mole * PhyConst::MolarGasConstant * Temp_K) / Volume_m3;
		}
		static Real Volume(Real n_mole, Real Pressure_Pa, Real Temp_K)
		{
			// Ideal Gas Law: PV = nRT => V = nRT / P
			return (n_mole * PhyConst::MolarGasConstant * Temp_K) / Pressure_Pa;
		}
		static Real Temperature(Real n_mole, Real Pressure_Pa, Real Volume_m3)
		{
			// Ideal Gas Law: PV = nRT => T = PV / (nR)
			return (Pressure_Pa * Volume_m3) / (n_mole * PhyConst::MolarGasConstant);
		}

		static Real AvgSpeed(Real Temp_K, Real molarMass_kg_per_mol)
		{
			// Average speed of gas molecules: v_avg = sqrt((3 * R * T) / M)
			// where M is the molar mass in kg/mol
			return std::sqrt((3 * PhyConst::MolarGasConstant * Temp_K) / molarMass_kg_per_mol);
		}
	};
}

#endif // MPL_IDEAL_GAS_CALCULATOR_H