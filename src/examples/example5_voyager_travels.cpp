#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Derivation.h"
#include "core/Function.h"
#endif


using namespace MML;

// SolarSystem gravity simulator

class SolarSystemGravitySimulator
{
	// ima definiranu putanju za sve planete
	ParametricCurve<3> _planetPaths[9];
};

void Example5_Voyager_travels()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 5 - Voyager travels                 ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// to je small body simulator, u zadanom polju tijela, koje se gibaju po zadanim putanjama

}