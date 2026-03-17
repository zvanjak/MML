#if !defined MPL_RELATIVISTIC_MECHANICS_SIMULATOR_H
#define MPL_RELATIVISTIC_MECHANICS_SIMULATOR_H

#include "MMLBase.h"


using namespace MML;

namespace MPL
{
	// keypoint - create realistic trajectory of a particle under relativistic mechanics
	// making sure that velocity never exceeds speed of light, and its trajectory is timelike

	// first question is in what frame are we observing the particle? and recording its trajectory?
	// and how we are going to represent the trajectory, and parameterize it?
	class RelativisticMechanicsSimulator
	{
		// linear trajectory, along x-axis, constant velocity, from point A to point B, in time T
		ParametricCurveFromStdFunc<4> getLinearTrajectory(Real x0, Real x1, Real velocity)
		{
			Real t1 = 0.0;
			Real t2 = (x1 - x0) / velocity; 

			ParametricCurveFromStdFunc<4> line(t1, t2, [=](Real t) { return VectorN<Real, 4>{ t, x0 + velocity * t, 0.0, 0.0 }; });

			return line;
		}

		// accelerated trajectory with constant acceleration, until it reaches some velocity, then constant velocity
		// from point A to point B, in time T
		ParametricCurveFromStdFunc<4> getAcceleratedTrajectory(Real x0, Real x1, Real acceleration, Real maxVelocity)
		{
			// acceleration and maxVelocity must have the same sign

			Real t1 = 0.0;
			
			// calculate time to reach maxVelocity
			Real t2 = maxVelocity / acceleration;
			Real x2 = x0 + 0.5 * acceleration * t2 * t2;		// position at time t2

			// calculate time to reach point x1
			Real t3 = (x1 - x2) / maxVelocity;

			ParametricCurveFromStdFunc<4> line(t1, t2+t3, [=](Real t) { 
					if (t < t2)
						return VectorN<Real, 4>{ t, x0 + Real(0.5) * acceleration * t * t, Real(0.0), Real(0.0) };
					else
						return VectorN<Real, 4>{ t, x2 + maxVelocity * (t-t2), Real(0.0), Real(0.0) };
				});

			return line;
		}

		// roundabout trip - with semicircular turn at the middle of journey at point B, and then getting back to starting point
	};
}

#endif