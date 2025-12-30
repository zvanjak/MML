#if !defined MPL_PROJECTILES_3D_H
#define MPL_PROJECTILES_3D_H

#include "interfaces/IODESystem.h"

#include "base/VectorTypes.h"
#include "base/InterpolatedFunction.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

using namespace MML;

namespace MPL
{
	// Projectile in 3D motion ODE system interface
	class IProjectileODE : public IODESystem
	{
	public:
		virtual std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "x";		// x position
			case 1: return "y";		// y position
			case 2: return "z";		// z position
			case 3: return "vx";	// x velocity
			case 4: return "vy";	// y velocity
			case 5: return "vz";	// z velocity
			default: return IODESystem::getVarName(ind);
			}
		}
		
		Vector<Real> getInitCond(Real initHeight, Vec3Cart velocity) const
		{
			Vector<Real> initCond(6);
			initCond[0] = 0;							// initial x position
			initCond[1] = 0;							// initial y position
			initCond[2] = initHeight;			// initial z position

			initCond[3] = velocity.X();	// initial x velocity (m/s)
			initCond[4] = velocity.Y();	// initial y velocity (m/s)
			initCond[5] = velocity.Z();	// initial z velocity (m/s)

			return initCond;
		}

		Vector<Real> getInitCond(Real x0, Real y0, Real z0, Vec3Cart velocity) const
		{
			Vector<Real> initCond(6);
			initCond[0] = x0;
			initCond[1] = y0;
			initCond[2] = z0;

			initCond[3] = velocity.X();	// initial x velocity (m/s)
			initCond[4] = velocity.Y();	// initial y velocity (m/s)
			initCond[5] = velocity.Z();	// initial z velocity (m/s)

			return initCond;
		}
		Vector<Real> getInitCond(Real initHeight, Real verticalAngle, Real angleToXAxis,Real velocityMagnitude) const
		{
			Vec3Cart velocity{ velocityMagnitude * cos(verticalAngle) * cos(angleToXAxis),
												 velocityMagnitude * cos(verticalAngle) * sin(angleToXAxis),
												 velocityMagnitude * sin(verticalAngle) };

			return getInitCond(initHeight, velocity);
		}
	};

	// ODE system for projectile in 3D in vacuum
	class ProjectileInVacuumODE : public IProjectileODE
	{
	public:
		int  getDim() const override { return 6; }
		virtual void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[3];			// dx/dt = vx
			dxdt[1] = x[4];			// dy/dt = vy
			dxdt[2] = x[5];			// dz/dt = vz
			dxdt[3] = 0;				// dvx/dt = 0 (no air resistance)
			dxdt[4] = 0;				// dvy/dt = 0 (no air resistance)
			dxdt[5] = -9.81;		// dvz/dt = -g (gravity)

			if( x[2] < 0.0 ) // if the projectile is below ground level
			{
				dxdt[0] = 0.0; 
				dxdt[1] = 0.0; 
				dxdt[2] = 0.0; 
				dxdt[3] = 0.0; 
				dxdt[4] = 0.0; 
				dxdt[5] = 0.0; 
			}
		}

		// Calculate time of flight using the formula: t = (v_z + sqrt(v_z^2 + 2 * g * h)) / g
		Real TimeOfFlight(Real initHeight, Vec3Cart velocity) const
		{
			return (velocity.Z() + std::sqrt(velocity.Z() * velocity.Z() + 2 * 9.81 * initHeight)) / 9.81;
		}
		// Calculate maximum height using the formula: h_max = h + (v_z^2) / (2 * g)
		Real MaxHeight(Real initHeight, Vec3Cart velocity) const
		{
			return initHeight + (velocity.Z() * velocity.Z()) / (2 * 9.81);
		}
		// Calculating landing position of the projectile
		Pnt3Cart calculateLandingPos(Real initHeight, Vec3Cart velocity) const
		{
			// Calculate time of flight
			Real timeOfFlight = TimeOfFlight(initHeight, velocity);

			return Pnt3Cart(
				velocity.X() * timeOfFlight,  // x position at landing
				velocity.Y() * timeOfFlight,  // y position at landing
				0.0                           // z position at landing (ground level)
			);
		}
		Pnt3Cart calculateLandingPos(Real x0, Real y0, Real z0, Vec3Cart velocity) const
		{
			// Calculate time of flight
			Real timeOfFlight = TimeOfFlight(z0, velocity);
			return Pnt3Cart(
				x0 + velocity.X() * timeOfFlight,  // x position at landing
				y0 + velocity.Y() * timeOfFlight,  // y position at landing
				0.0                                 // z position at landing (ground level)
			);
		}
	};

	// Class modeling ODE system for projectile launch, on 2D carousel with counterclockwise rotation, embedded in 3D,
	// where Coriolis force is included in equations of motion
	// Launch point is at radius _r from the center of the carousel, and the carousel rotates with angular velocity _w
	// Local x-y-z coordinate system is defined as follows:
	// - x-axis is directed along the radius of the carousel
	// - y-axis is directed tangentially to the carousel in the direction of rotation
	// - z-axis is directed vertically upwards
	// Origin of local coordinate system is at the center of the carousel
	class ProjectileInVacuumLaunchFromCarouselODE : public IProjectileODE
	{
		Real _w;			// angular velocity of the carousel (rad/s)

	public:
		ProjectileInVacuumLaunchFromCarouselODE(Real angularVelocity) : _w(angularVelocity) {}

		int  getDim() const override { return 6; }
		// define equations of motion in local coordinate system
		virtual void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			// x, y, z are positions
			// vx, vy, vz are velocities
			dxdt[0] = x[3];								// dx/dt = vx
			dxdt[1] = x[4];								// dy/dt = vy
			dxdt[2] = x[5];								// dz/dt = vz

			// Coriolis force 
			dxdt[3] = 2 * _w * x[4];      // dvx/dt = 2*w*vy
			dxdt[4] = -2 * _w * x[3];     // dvy/dt = -2*w*vx

			// Centrifugal force
			dxdt[3] += _w * _w * x[0];    // dvx/dt += w^2 * x
			dxdt[4] += _w * _w * x[1];    // dvy/dt += w^2 * y

			dxdt[5] = -9.81;						// dvz/dt = -g (gravity)
			
			if( x[2] < 0.0 ) // if the projectile is below ground level
			{
				dxdt[0] = 0.0; 
				dxdt[1] = 0.0; 
				dxdt[2] = 0.0; 
				dxdt[3] = 0.0; 
				dxdt[4] = 0.0; 
				dxdt[5] = 0.0; 
			}
		}
	};

	// Class modeling ODE system for projectile launch from the surface of the Earth,
	// from given latitude and longitude, with given initial velocity,
	// using the local coordinate system where:
	// - x-axis is directed along the north-south meridian at the launch point
	// - y-axis is directed tangentially to the Earth in the direction of rotation (east-west parallel
	// - z-axis is directed vertically upwards
	// Takes into account the rotation of the Earth (Coriolis and centrifugal forces) and gravity.
	class ProjectileInVacuumLaunchFromEarthLocalTangentPlaneCoordsODE : public IProjectileODE
	{
		// givenparameters
		Real _latitude;		// latitude of the launch point (in radians)
		Real _longitude;	// longitude of the launch point (in radians)
		Vec3Cart _launchVelocity; // initial velocity of the projectile in local coordinate system

		// calculated parameters
		Real _w;			// angular velocity of the Earth (rad/s) at given latitude

	public:
		ProjectileInVacuumLaunchFromEarthLocalTangentPlaneCoordsODE(Real latitudeDeg, Real longitudeDeg, Vec3Cart launchVelocity)
			: _launchVelocity(launchVelocity)
		{
			// convert latitude and longitude from degrees to radians
			_latitude = Utils::DegToRad(latitudeDeg);
			_longitude = Utils::DegToRad(longitudeDeg);

			Real earthW = 2 * Constants::PI / 86400.0; // angular velocity of the Earth (rad/s)
			
			// calculate angular velocity of the Earth at given latitude
			_w = earthW * cos(_latitude);
		}
		
		int getDim() const override { return 6; }

		virtual void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			// x[0]: x (south), x[1]: y (east), x[2]: z (up)
			// x[3]: vx, x[4]: vy, x[5]: vz

			// Earth's angular velocity in local (south, east, up) coordinates
			const Real omegaE = 2 * Constants::PI / 86400.0;
			
			// Calculate components of angular velocity vector in local coordinates
			const Real omega_x = -omegaE * std::cos(_latitude); // south
			const Real omega_y = 0.0;              // east
			const Real omega_z = omegaE * std::sin(_latitude);  // up

			Vec3Cart r(x[0], x[1], x[2]);
			Vec3Cart v(x[3], x[4], x[5]);
			Vec3Cart omega(omega_x, omega_y, omega_z);

			// Coriolis acceleration: -2 * omega x v
			Vec3Cart a_coriolis = -2.0 * VectorProduct(omega, v);

			// Centrifugal acceleration: -omega x (omega x r)
			Vec3Cart a_centrifugal = -VectorProduct(omega, VectorProduct(omega, r));

			// Gravity (downward, -z)
			Vec3Cart a_gravity(0.0, 0.0, -9.81);

			// Total acceleration
			Vec3Cart a_total = a_coriolis + a_centrifugal + a_gravity;

			// Derivatives
			dxdt[0] = x[3]; // dx/dt = vx
			dxdt[1] = x[4]; // dy/dt = vy
			dxdt[2] = x[5]; // dz/dt = vz
			dxdt[3] = a_total.X();
			dxdt[4] = a_total.Y();
			dxdt[5] = a_total.Z();

			// Stop integration if below ground
			if (x[2] < 0.0)
			{
				for (int i = 0; i < 6; ++i)
					dxdt[i] = 0.0;
			}
		}

		Vector<Real> getInitCond() const
		{
			// Initial position: (0, 0, 0) in local (south, east, up) coordinates
			// Initial velocity: as provided in _launchVelocity (already in local frame)
			Vector<Real> initCond(6);
			initCond[0] = 0.0; // x (south)
			initCond[1] = 0.0; // y (east)
			initCond[2] = 0.0; // z (up)
			initCond[3] = _launchVelocity.X();
			initCond[4] = _launchVelocity.Y();
			initCond[5] = _launchVelocity.Z();
			return initCond;
		}
	};

	// Class modeling ODE system for projectile launch from the surface of the Earth,
	// from given latitude and longitude, with given initial velocity,
	// and solving it within global Cartesian system with origin at the center of the Earth,
	//	with x-axis directed towards the Greenwich meridian,
	//	with y-axis directed towards the 90 degrees east longitude,
	//	with the z-axis directed upwards,
	// Earth gravity force is calculated from that global position, relative to the center of the Earth,
	// with the assumption that the Earth is a perfect sphere with radius 6371 km.
	// Influence of the Earth's rotation is taken into account solely by initial conditions,
	// by transforming given (local) initial velocity to the global Cartesian coordinate system,
	// and simulating projectile for the given time in that global coordinate system.
	// But it also provide transformation formulas that take into account rotation of the Earth
	// with transformation from this global Cartesian coordinate system to the (latitude, longitude, altitude) system
	class ProjectileInVacuumLaunchFromEarthGlobalCartesianCoordsODE : public IProjectileODE
	{
		// given parameters
		Real			_latitude;				// latitude of the launch point (in radians)
		Real			_longitude;				// longitude of the launch point (in radians)
		Vec3Cart	_launchVelocity;	// initial velocity of the projectile in local coordinate system
		
		Real _earthRadius;					// radius of the Earth (in meters)
		Real _w;										// angular velocity of the Earth (rad/s) at given latitude
	
	public:
		ProjectileInVacuumLaunchFromEarthGlobalCartesianCoordsODE(Real latitudeDeg, Real longitudeDeg, Vec3Cart launchVelocity, Real earthRadius = 6371000.0)
			: _earthRadius(earthRadius), _launchVelocity(launchVelocity)
		{
			// convert latitude and longitude from degrees to radians
			_latitude = Utils::DegToRad(latitudeDeg);
			_longitude = Utils::DegToRad(longitudeDeg);

			_w = 2 * Constants::PI / 86400.0; // angular velocity of the Earth (rad/s)
		}

		int getDim() const override { return 6; }
		virtual void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			Vec3Cart r(x[0], x[1], x[2]); // position vector in global Cartesian coordinates
			Vec3Cart v(x[3], x[4], x[5]); // velocity vector in global Cartesian coordinates

			// Calculate gravitational acceleration
			Vec3Cart a_gravity = -9.81 * r.Normalized(); // gravity acts towards the center of the Earth

			// Derivatives
			dxdt[0] = v.X(); // dx/dt = vx
			dxdt[1] = v.Y(); // dy/dt = vy
			dxdt[2] = v.Z(); // dz/dt = vz
			dxdt[3] = a_gravity.X();
			dxdt[4] = a_gravity.Y();
			dxdt[5] = a_gravity.Z();

			// Stop integration if below ground
			if (r.NormL2() < _earthRadius)
			{
				for (int i = 0; i < 6; ++i)
					dxdt[i] = 0.0;
			}
		}
		Vector<Real> getInitCond() const
		{
			Vector<Real> initCond(6);

			// Initial position in global Cartesian coordinates
			Real x0 = _earthRadius * std::cos(_latitude) * std::cos(_longitude);
			Real y0 = _earthRadius * std::cos(_latitude) * std::sin(_longitude);
			Real z0 = _earthRadius * std::sin(_latitude);

			initCond[0] = x0;
			initCond[1] = y0;
			initCond[2] = z0;

			// Local (south, east, up) to global Cartesian transformation
			// South: (-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat))
			// East:  (-sin(lon), cos(lon), 0)
			// Up:    (cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
			Real v_south = _launchVelocity.X();
			Real v_east = _launchVelocity.Y();
			Real v_up = _launchVelocity.Z();

			Real vx_launch =
				v_south * (-std::sin(_latitude) * std::cos(_longitude)) +
				v_east * (-std::sin(_longitude)) +
				v_up * (std::cos(_latitude) * std::cos(_longitude));

			Real vy_launch =
				v_south * (-std::sin(_latitude) * std::sin(_longitude)) +
				v_east * (std::cos(_longitude)) +
				v_up * (std::cos(_latitude) * std::sin(_longitude));

			Real vz_launch =
				v_south * (std::cos(_latitude)) +
				v_east * (0.0) +
				v_up * (std::sin(_latitude));

			// Velocity of the launch point due to Earth's rotation: omega x r0
			// omega = (0, 0, _w)
			// r0 = (x0, y0, z0)
			Real vx_earth = -_w * y0;
			Real vy_earth = _w * x0;
			Real vz_earth = 0.0;

			// Total initial velocity in global frame
			initCond[3] = vx_launch + vx_earth;
			initCond[4] = vy_launch + vy_earth;
			initCond[5] = vz_launch + vz_earth;

			return initCond;
		}
	};
}

#endif // MPL_PROJECTILES_3D_H
