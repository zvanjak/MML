#if !defined MPL_ELECTRICCHARGEDISTRIBUTIONSIMULATOR_H
#define MPL_ELECTRICCHARGEDISTRIBUTIONSIMULATOR_H

#include "MMLBase.h"

#include "base/VectorTypes.h"

using namespace MML;

namespace MPL
{
	// Calculating distibution of charges, for different conductor configurations

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
			: _entryPoint(Vector3Cartesian(0, 0, 0))
		{
		}
		ConductorGeometry(Vector3Cartesian entryPoint)
			: _entryPoint(entryPoint)
		{
		}

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
		Point3Cartesian projectBackToGeometry(const Point3Cartesian& init_pos, const Point3Cartesian& new_pos) override
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
		{
		}

		CylinderGeometry(Real radius, Real length, Vector3Cartesian entryPoint)
			: _radius(radius), _length(length), ConductorGeometry(entryPoint)
		{
		}

		bool isWithinBoundary(Vector3Cartesian pos) override
		{
			return false;
		}

	};

}

#endif // MPL_ELECTRICCHARGEDISTRIBUTIONSIMULATOR_H