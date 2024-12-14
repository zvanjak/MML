#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorTypes.h"

#include "core/Derivation.h"
#endif


using namespace MML;

// Calcualting distibution of charges, for different conductor configurations

struct SimpleCharge
{
	Point3Cartesian _pos;
	Vector3Cartesian _velocity;
	Real _charge;
};

class ConductorGeometry
{
public:
	// entry point - where charges are added
	Vector3Cartesian _entryPoint;

	std::vector<SimpleCharge> _charges;

	ConductorGeometry()
		: _entryPoint(Vector3Cartesian(0,0,0))
	{	}
	ConductorGeometry(Vector3Cartesian entryPoint)
		: _entryPoint(entryPoint)
	{	}

	virtual bool isWithinBoundary(Vector3Cartesian pos) = 0;
	virtual Point3Cartesian projectBackToGeometry(const Point3Cartesian& init_pos, const Point3Cartesian& new_pos) = 0;
};

class ParalelepipedGeometry : public ConductorGeometry
{
public:
	Real _a, _b, _c;

	bool isWithinBoundary(Vector3Cartesian pos) override
	{
		return false;
	}
	Point3Cartesian projectBackToGeometry(const Point3Cartesian& init_pos, const Point3Cartesian& new_pos)
	{
		return Point3Cartesian();
	}

};

class CylinderGeometry : public ConductorGeometry
{
public:
	Real _radius;
	Real _length;

	CylinderGeometry(Real radius, Real length)
		: _radius(radius), _length(length) 
	{	}

	CylinderGeometry(Real radius, Real length, Vector3Cartesian entryPoint)
		: _radius(radius), _length(length), ConductorGeometry(entryPoint)
	{ }

	bool isWithinBoundary(Vector3Cartesian pos) override
	{
		return false;
	}

};


void Example6_electric_charge_distribution()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****          EXAMPLE 6 - electric charge distributions            ****" << std::endl;
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