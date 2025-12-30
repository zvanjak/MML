#if !defined MPL_MECHANICS_MOVING_FRAMES_TRANSF_H
#define MPL_MECHANICS_MOVING_FRAMES_TRANSF_H

#include "MMLBase.h"

#include "MovingFrames.h"

using namespace MML;

namespace MPL
{
	template<typename VectorFrom, typename VectorTo, int N>
	class IMovingFrameTransf 
	{
	public:
		virtual       VectorTo            transf(Real t, const VectorFrom& posLocal) const = 0;
		virtual       VectorFrom          transfInverse(Real t, const VectorTo& posParent) const = 0;

		// TODO - add function isInertial() that returns true if the frame is inertial, false otherwise)
		// but, calculate it from the transformation functions

		virtual ~IMovingFrameTransf() {}
	};

	// class representing transformation from moving frame to fixed frame
	class MovingFrameToFixedTransf : public IMovingFrameTransf<Vec3Cart, Vec3Cart, 3>
	{
		MovingInertialFrame& _movingFrame;

	public:
		MovingFrameToFixedTransf(MovingInertialFrame& frame)
			: _movingFrame(frame)
		{
		}

		virtual Vec3Cart transf(Real t, const Vec3Cart& posLocal) const override
		{
			Vec3Cart pos = _movingFrame._pos_at_0 + _movingFrame._velocity * t;
			return posLocal + pos; // transform to fixed frame
		}
		virtual Vec3Cart transfInverse(Real t, const Vec3Cart& posFixed) const override
		{
			Vec3Cart pos = _movingFrame._pos_at_0 + _movingFrame._velocity * t;
			return posFixed - pos; // transform back to moving frame
		}
	};

	// class representing transf. from RotatingCarouselFrame to fixed frame
	// using Cartesian coordinates for local frame 
	class RotatingCarouselLocalCartToFixedTransf : public IMovingFrameTransf<Vec3Cart, Vec3Cart, 3>
	{
		RotatingCarouselFrame& _rotatingFrame;
	public:
		RotatingCarouselLocalCartToFixedTransf(RotatingCarouselFrame& frame)
			: _rotatingFrame(frame)
		{
		}
		virtual Vec3Cart transf(Real t, const Vec3Cart& posLocal) const override
		{
			// resulting global coordinates
			Vec3Cart posGlob;

			Real phi_global = _rotatingFrame._w * t;	// global azimuthal angle considering carousel rotation
			
			posGlob.X() = posLocal.X() * cos(phi_global) - posLocal.Y() * sin(phi_global);
			posGlob.Y() = posLocal.X() * sin(phi_global) + posLocal.Y() * cos(phi_global);
			posGlob.Z() = posLocal.Z();

			return posGlob;
		}
		virtual Vec3Cart transfInverse(Real t, const Vec3Cart& posFixed) const override
		{
			// global coordinates
			Real x_global = posFixed.X();
			Real y_global = posFixed.Y();
			Real z_global = posFixed.Z();
			// local coordinates
			Real x_local, y_local, z_local;

			Real phi_global = _rotatingFrame._w * t; // global azimuthal angle considering carousel rotation
			x_local = x_global * cos(-phi_global) - y_global * sin(-phi_global);
			y_local = x_global * sin(-phi_global) + y_global * cos(-phi_global);
			z_local = z_global;

			Vec3Cart pos = Vec3Cart{ x_local, y_local, z_local };
			return pos; // transform back to moving frame
		}
	};

	// class representing transf. from RotatingCarouselFrame to fixed frame
	// using cylindrical coordinates for local frame 
	class RotatingCarouselLocalCylToFixedTransf : public IMovingFrameTransf<Vec3Cyl, Vec3Cart, 3>
	{
		RotatingCarouselFrame& _rotatingFrame;
	public:
		RotatingCarouselLocalCylToFixedTransf(RotatingCarouselFrame& frame)
			: _rotatingFrame(frame)
		{
		}
		virtual Vec3Cart transf(Real t, const Vec3Cyl& posLocal) const override
		{
			// resulting global coordinates
			Vec3Cart posGlob;
			
			Real phi_global = _rotatingFrame._w * t; // global azimuthal angle considering carousel rotation
			
			posGlob.X() = posLocal.R() * cos(posLocal.Phi() + phi_global);
			posGlob.Y() = posLocal.R() * sin(posLocal.Phi() + phi_global);
			posGlob.Z() = posLocal.Z();
			
			return posGlob;
		}
		virtual Vec3Cyl transfInverse(Real t, const Vec3Cart& posFixed) const override
		{
			// global coordinates
			Real x_global = posFixed.X();
			Real y_global = posFixed.Y();
			Real z_global = posFixed.Z();
			
			// local coordinates
			Real r_local, phi_local, z_local;
			Real phi_global = _rotatingFrame._w * t; // global azimuthal angle considering carousel rotation
			
			r_local = std::sqrt(x_global * x_global + y_global * y_global);
			if (r_local == 0.0)
				phi_local = 0.0;
			else
				phi_local = std::atan2(y_global, x_global) - phi_global;
			z_local = z_global;
			
			Vec3Cyl pos = Vec3Cyl{ r_local, phi_local, z_local };
			return pos; // transform back to moving frame
		}

	};
	
	// class representing local transf. for RotatingCarouselFrame 
	// from local Cartesian to local cylindrical coordinates in local frame 
	class RotatingCarouselLocalCartToLocalCylTransf : public IMovingFrameTransf<Vec3Cart, Vec3Cyl, 3>
	{
		RotatingCarouselFrame& _rotatingFrame;
	public:
		RotatingCarouselLocalCartToLocalCylTransf(RotatingCarouselFrame& frame)
			: _rotatingFrame(frame)
		{	}
		
		virtual Vec3Cyl transf(Real t, const Vec3Cart& posLocal) const override
		{
			Vec3Cyl pos;
			pos.R() = std::sqrt(posLocal.X() * posLocal.X() + posLocal.Y() * posLocal.Y());
			if (pos.R() > 0.0)
			{
				pos.Phi() = std::atan2(posLocal.Y(), posLocal.X()); // azimuthal angle in local cylindrical coordinates
			}
			else
			{
				pos.Phi() = 0.0; // if r_local is zero, set phi_local to zero
			}
			pos.Z() = posLocal.Z(); // z-coordinate remains the same in both coordinate systems
			return pos;
		}
		virtual Vec3Cart transfInverse(Real t, const Vec3Cyl& posLocal) const override
		{
			Vec3Cart pos;
			pos.X() = posLocal.R() * cos(posLocal.Phi());
			pos.Y() = posLocal.R() * sin(posLocal.Phi());
			pos.Z() = posLocal.Z();
			return pos; // transform back to moving frame
		}
	};

	// class representing transformation from Circle orbiting frame to fixed frame
	class CircleOrbitingFrameToFixedTransf : public IMovingFrameTransf<Vec3Cart, Vec3Cart, 3>
	{
		CircleOrbitingFrameCartesian& _orbitingFrame;
	public:
		CircleOrbitingFrameToFixedTransf(CircleOrbitingFrameCartesian& frame)
			: _orbitingFrame(frame)
		{
		}

		virtual Vec3Cart transf(Real t, const Vec3Cart& posLocal) const override
		{
			Real angle = _orbitingFrame._angle_at_t0 + 2 * Constants::PI * t / _orbitingFrame._period;

			Vec3Cart pos = Vec3Cart{ _orbitingFrame._orbitRadius * cos(angle), _orbitingFrame._orbitRadius * sin(angle), 0.0 };

			return posLocal + pos;	// transform to fixed frame
		}

		virtual Vec3Cart transfInverse(Real t, const Vec3Cart& posFixed) const override
		{
			Real angle = _orbitingFrame._angle_at_t0 + 2 * Constants::PI * t / _orbitingFrame._period;

			Vec3Cart pos = Vec3Cart{ _orbitingFrame._orbitRadius * cos(angle), _orbitingFrame._orbitRadius * sin(angle), 0.0 };

			return posFixed - pos;				// transform back to moving frame
		}
	};

	// class representing transformation from RotatingCylindricalFrame to fixed frame
	class RotatingCylindricalFrameToFixedTransf : public IMovingFrameTransf<Vec3Cart, Vec3Cart, 3>
	{
		RotatingFrameCylindrical& _rotatingFrame;

	public:
		RotatingCylindricalFrameToFixedTransf(RotatingFrameCylindrical& frame)
			: _rotatingFrame(frame)
		{
		}

		virtual Vec3Cart transf(Real t, const Vec3Cart& posLocal) const override
		{
			// Compute total rotation angle at time t
			Real angle = _rotatingFrame._angleAtT0 + _rotatingFrame._angularVelocity * t;
			Vec3Cart axis = _rotatingFrame._axis.Normalized();

			// Rodrigues' rotation formula
			Vec3Cart v = posLocal;
			Real cosA = std::cos(angle);
			Real sinA = std::sin(angle);
			Vec3Cart rotated =
				v * cosA +
				VectorProduct(axis, v) * sinA +
				axis * (axis * v) * (1 - cosA);

			return rotated;
		}

		virtual Vec3Cart transfInverse(Real t, const Vec3Cart& posFixed) const override
		{
			// Inverse rotation: use -angle
			Real angle = _rotatingFrame._angleAtT0 + _rotatingFrame._angularVelocity * t;
			Vec3Cart axis = _rotatingFrame._axis.Normalized();

			Vec3Cart v = posFixed;
			Real cosA = std::cos(-angle);
			Real sinA = std::sin(-angle);
			Vec3Cart rotated =
				v * cosA +
				VectorProduct(axis, v) * sinA +
				axis * (axis * v) * (1 - cosA);

			return rotated;
		}

	};


	// class representing transformation from RotatingFrameSpherical to fixed frame
	// TODO - FIX!!!
	class RotatingSphericalFrameToFixedTransf : public IMovingFrameTransf<Vec3Sph, Vec3Cart, 3>
	{
		RotatingFrameSpherical& _rotatingFrame;
	public:
		RotatingSphericalFrameToFixedTransf(RotatingFrameSpherical& frame)
			: _rotatingFrame(frame)
		{
		}
		virtual Vec3Cart transf(Real t, const Vec3Sph& posLocal) const override
		{
			// Compute total rotation angle at time t
			Real angle = _rotatingFrame._phiAtT0 + _rotatingFrame._angularVelocity * t;
			Vec3Cart axis = _rotatingFrame._axis.Normalized();

			// Rodrigues' rotation formula
			Vec3Cart v = posLocal;
			Real cosA = std::cos(angle);
			Real sinA = std::sin(angle);
			Vec3Cart rotated =
				v * cosA +
				VectorProduct(axis, v) * sinA +
				axis * (axis * v) * (1 - cosA);
			return rotated;
		}
		virtual Vec3Sph transfInverse(Real t, const Vec3Cart& posFixed) const override
		{
			// Inverse rotation: use -angle
			Real angle = _rotatingFrame._phiAtT0 + _rotatingFrame._angularVelocity * t;
			Vec3Cart axis = _rotatingFrame._axis.Normalized();
			Vec3Cart v = posFixed;
			Real cosA = std::cos(-angle);
			Real sinA = std::sin(-angle);
			Vec3Cart rotated =
				v * cosA +
				VectorProduct(axis, v) * sinA +
				axis * (axis * v) * (1 - cosA);
			return rotated;
		}
	};
}

#endif
