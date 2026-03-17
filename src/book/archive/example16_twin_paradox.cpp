#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Derivation.h"

#include "mpl/SpecialRelativity/TwinParadoxSimulator.h"
#endif


using namespace MML;

// first demo - twin moves from A to B, distance 10 lightyears, speed 0.8c, then returns to A
// simple case - no acceleration, just constant speed
 

// second demo - twin moves from A to B, distance 10 lightyears, speed 0.8c, then returns to A, 
// but with acceleration at the start and deceleration at the end of the trip


// third demo - twin start at A with constant acceleration (until it reaches speed 0.8c),
// and goes to black hole at B, 10 lightyears away, goes around it, then returns to A
// with deceleration at the end of the trip
// perform claculation within general relativity framework, with Schwarzschild metric

void Example16_Twin_paradox()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 16 - Twin paradox                   ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;
}