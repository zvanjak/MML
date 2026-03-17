#ifndef MPL_GASCONSTANTS_H
#define MPL_GASCONSTANTS_H

#include "MMLBase.h"

using namespace MML;

#include "PhysicalConstants.h"

namespace MPL
{
	struct GasData
	{
		std::string name;							// Name of the gas
		std::string formula;					// Chemical formula of the gas 
		
		double molarMass_g_mol;				// Molecular weight in g/mol
		
		double radiusVanDerWaals;			// Van der Waals radius in Angstroms
		double kineticDiameter;				// Kinetic diameter in Angstroms

		double criticalTemperature;		// Critical temperature in K
		double criticalPressure;			// Critical pressure in Pa

		Real MolarMass_kg_mol()	const { return molarMass_g_mol / 1000.0; }
		Real MoleculeMass_kg()	const { return molarMass_g_mol / 1000.0 / PhyConst::AvogardoConstant ; }
		
		Real RadiusVanDerWaals_m()	const { return radiusVanDerWaals * 1e-10; }
		Real RadiusVanDerWaals_nm()	const { return radiusVanDerWaals * 0.1; }
		Real RadiusKinetic_m() 			const { return kineticDiameter / 2 * 1e-10; }
		Real RadiusKinetic_nm() 		const { return kineticDiameter / 2 * 0.1; }
	};

	class GasConstants
	{
	public:
		const static inline GasData _gasses[] = {
			{ "Hydrogen",							"H2",		 2.016,  1.2,		2.9,  33.2,  1.3e6 },
			{ "Helium",								"He",		 4.0026, 1.4,		2.6,   5.2,  2.3e6 },
			{ "Methane",							"CH4",	16.04,	 1.75,	3.8, 190.6,  4.6e6 },
			{ "Ammonia",							"NH3",	17.0305, 1.8,		3.6, 405.5, 11.e6 },
			{ "Neon",									"Ne",		20.1797, 1.54,	2.8,  44.4,	 2.7e6 },
			{ "Nitrogen",							"N2",		28.0134, 1.55,	3.7, 126.2,  3.4e6 },
			{ "Carbon Monoxide",			"CO",		28.0101, 1.76,	3.7, 132.,   3.e6 },
			{ "Oxygen",								"O2",		31.9988, 1.52,	3.5, 154.6,  5.0e6 },
			{ "Hydrogen Chloride",		"HCl",	36.461,  1.75,	3.4, 188.,   8.e6 },
			{ "Argon",								"Ar",		39.948,  1.88,	3.7, 150.8,  4.9e6 },
			{ "Carbon Dioxide",				"CO2",	44.01,	 1.8,		3.3, 304.1,  7.4e6 },
			{ "Nitrous Oxide",				"N2O",	44.0128, 1.8,		3.5, 309.,   7.e6 },
			{ "Krypton",							"Kr",		83.798,  2.,		3.,  209.,   5.e6 },
			{ "Xenon",								"Xe",  131.293,  2.,		4.,  289.,   5.e6},
			{ "Radon",								"Rn",	 222.0,		 2.2,		4.5, 377.0,  6.3e6 },
			{ "Sulfur Hexafluoride",	"SF6", 146.055,  1.9,		5.5, 373.0,  3.7e6 }
		};
		static const GasData& gasByName(const std::string& name)
		{
			for (const auto& gas : _gasses)
			{
				if (gas.name == name)
					return gas;
			}
			throw std::runtime_error("Gas not found: " + name);
		}
		static const GasData& gasByFormula(const std::string& formula)
		{
			for (const auto& gas : _gasses)
			{
				if (gas.formula == formula)
					return gas;
			}
			throw std::runtime_error("Gas not found: " + formula);
		}

		static size_t getGasCount() {
			return sizeof(_gasses) / sizeof(GasData);
		}
	};
}

#endif // MPL_GASCONSTANTS_H