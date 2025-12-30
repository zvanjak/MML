#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/VectorTypes.h"

#include "mpl/Electromagnetism/ElectricChargeDistributionSimulator.h"
#endif

using namespace MML;
using namespace MPL;

// First demo - fill given geometry randomly with electric charges
// and simulate their final distribution

// Second demo - fill given geometry with electric charges one by one
// and simulate their final distribution after each charge addition

// What happens if, after distributing all charges, we attach long wire to the geometry
// and connect it to the ground? Will the charges redistribute?

void Chapter17_electric_charge_distribution()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****          EXAMPLE 17 - electric charge distributions           ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// na neutral solid tube stavljas elektrone jedan po jedan
	// kad ce prvi doci na bocnu stranicu?
	// hoce li ih uopće biti na bocnoj stranici?

	Real radius = 1.0;
	Real length = 10.0;

	// define geometry

	int numCharges = 100;
	int numIntermedSteps = 100;		// number of steps to perform after addition of each charge
	int numSteps = (numCharges + 1) * numIntermedSteps;
	for (int i = 0; i < numSteps; i++)
	{
		// add charge
		Real charge = -1.0;
	Vector3Cartesian pos = Vector3Cartesian(radius, 0.0, 0.0);
	// calculate electric field
	// calculate force on charge
	}
}

