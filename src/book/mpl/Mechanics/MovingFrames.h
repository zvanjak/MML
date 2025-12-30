#if !defined MPL_MOVING_FRAMES_H
#define MPL_MOVING_FRAMES_H

#include "MMLBase.h"

#include "base/VectorTypes.h"
#include "base/Function.h"


using namespace MML;

namespace MPL
{
	class FixedGlobalFrame
	{ };

	// represents a frame that is moving with constant velocity
	class MovingInertialFrame
	{
	public:
		Vec3Cart _velocity;
		Vec3Cart _pos_at_0;
	};

	// represents a frame that is moving with constant acceleration
	class AcceleratingInertialFrame
	{
		public:
		Vec3Cart _velocity;			// velocity at t = 0
		Vec3Cart _acceleration;	// constant acceleration
		Vec3Cart _pos_at_0;			// position at t = 0
	};

	// represents a carousel frame that is rotating with constant angular velocity
	class RotatingCarouselFrame
	{
	public:
		Real _w;						// angular velocity of the frame
		Real _angle_at_t0;	// angle at t = 0
	};
	
	// represents a frame that is orbiting around some point in space
	class CircleOrbitingFrameCartesian 
	{
	public:
		Real _orbitRadius;
		Real _period;
		Real _angle_at_t0;
		
		Pnt3Cart _centerPoint;	// center point of the orbit
		Vec3Cart _axis;					// normal to plane of rotation (axis of rotation)
	};

	// represents a frame that is orbiting around given axis, with polar coordinates
	class RotatingFrameCylindrical
	{
	public:
		Real			_angularVelocity;		// angular velocity of the frame
		Vec3Cart	_axis;							// axis of rotation

		Real _angleAtT0;				// angle at t = 0
	};

	// represents a frame for rotating sphere, with spherical coordinates
	class RotatingFrameSpherical
	{
	public:
		Real			_angularVelocity;		// angular velocity of the frame
		Vec3Cart	_axis;							// axis of rotation (passing through origin)

		Real _phiAtT0;    // azimuthal angle at t = 0
	};

	// represents a local frame on a rotating sphere, with cartesian coordinates
	// x-axis points to the south pole, y-axis points to the east, and z-axis is perpendicular to the surface of the sphere
	class RotatingSphereLocalCartesian
	{
		RotatingFrameSpherical &_parentFrame;

	public:
		// specify position on the sphere where local frame origin is located
		// with latitude and longitude in radians
		// and z-axis pointing outwards from the sphere
		const Real _latitude;
		const Real _longitude;

		RotatingSphereLocalCartesian(RotatingFrameSpherical &parentFrame, Real latitude, Real longitude)
			: _parentFrame(parentFrame), _latitude(latitude), _longitude(longitude) 
		{	}
	};

}

#endif