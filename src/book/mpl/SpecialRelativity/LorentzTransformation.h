#if !defined MPL_LORENTZ_TRANSFORMATION_H
#define MPL_LORENTZ_TRANSFORMATION_H

#include "MMLBase.h"

#include "interfaces/ICoordTransf.h"

#include "base/Vector.h"
#include "base/VectorTypes.h"
#include "base/Matrix.h"
#include "base/MatrixNM.h"
#include "base/Tensor.h"
#include "base/Geometry.h"

#include "core/CoordTransf.h"

using namespace MML;

namespace MPL
{
	// class modeling Lorentz transformation along the X axis, for a system moving at velocity v along x
	// transformation from the rest frame to a frame moving with the given velocity along the X axis.
	class CoordTransfLorentzXAxis : public CoordTransfWithInverse<Vector4Minkowski, Vector4Minkowski, 4>
	{
	private:
		Real _velocity; // v in units of c
		Real _gamma;    // Lorentz factor
		Real _x0;				// distance between origins of nertial and moving frame

	public:
		CoordTransfLorentzXAxis(Real velocity) : _velocity(velocity), _x0(0.0)
		{
			if (_velocity <= -1.0 || _velocity >= 1.0)
				throw std::range_error("Invalid velocity for Lorentz transformation: must be in (-1, 1).");

			_gamma = 1.0 / std::sqrt(1.0 - _velocity * _velocity);
		}
		CoordTransfLorentzXAxis(Real x0, Real velocity)
			: _velocity(velocity), _x0(x0)
		{
			if (_velocity <= -1.0 || _velocity >= 1.0)
				throw std::range_error("Invalid velocity for Lorentz transformation: must be in (-1, 1).");

			_gamma = 1.0 / std::sqrt(1.0 - _velocity * _velocity);
		}

		// Forward Lorentz transformation: from rest frame to moving frame
		Vector4Minkowski transf(const Vector4Minkowski& q) const override
		{
			// q = (t, x, y, z)
			Real t = q[0];
			Real x = q[1];
			Real y = q[2];
			Real z = q[3];
			Real t_prime = _gamma * (t - _velocity * (x - _x0));
			Real x_prime = _gamma * ((x - _x0) - _velocity * t);
			return Vector4Minkowski{ t_prime, x_prime, y, z };
		}

		// Inverse Lorentz transformation: from moving frame to rest frame
		Vector4Minkowski transfInverse(const Vector4Minkowski& q) const override
		{
			Real t_prime = q[0];
			Real x_prime = q[1];
			Real y = q[2];
			Real z = q[3];
			Real t = _gamma * (t_prime + _velocity * x_prime);
			Real x = _gamma * (x_prime + _velocity * t_prime) + _x0;
			return Vector4Minkowski{ t, x, y, z };
		}

		// Scalar functions for each coordinate (optional, for interface compatibility)
		Real func1(const VectorN<Real, 4>& q) const { return _gamma * (q[0] - _velocity * q[1]); }
		Real func2(const VectorN<Real, 4>& q) const { return _gamma * (q[1] - _velocity * q[0]); }
		Real func3(const VectorN<Real, 4>& q) const { return q[2]; }
		Real func4(const VectorN<Real, 4>& q) const { return q[3]; }

		Real funcInverse1(const VectorN<Real, 4>& q) const { return _gamma * (q[0] + _velocity * q[1]); }
		Real funcInverse2(const VectorN<Real, 4>& q) const { return _gamma * (q[1] + _velocity * q[0]); }
		Real funcInverse3(const VectorN<Real, 4>& q) const { return q[2]; }
		Real funcInverse4(const VectorN<Real, 4>& q) const { return q[3]; }

		const IScalarFunction<4>& coordTransfFunc(int i) const override
		{
			static ScalarFunctionFromStdFunc<4> f1([this](const VectorN<Real, 4>& q) { return func1(q); });
			static ScalarFunctionFromStdFunc<4> f2([this](const VectorN<Real, 4>& q) { return func2(q); });
			static ScalarFunctionFromStdFunc<4> f3([this](const VectorN<Real, 4>& q) { return func3(q); });
			static ScalarFunctionFromStdFunc<4> f4([this](const VectorN<Real, 4>& q) { return func4(q); });

			switch (i) {
			case 0: return f1;
			case 1: return f2;
			case 2: return f3;
			default: return f4;
			}
		}
		const IScalarFunction<4>& inverseCoordTransfFunc(int i) const override
		{
			static ScalarFunctionFromStdFunc<4> fInv1([this](const VectorN<Real, 4>& q) { return funcInverse1(q); });
			static ScalarFunctionFromStdFunc<4> fInv2([this](const VectorN<Real, 4>& q) { return funcInverse2(q); });
			static ScalarFunctionFromStdFunc<4> fInv3([this](const VectorN<Real, 4>& q) { return funcInverse3(q); });
			static ScalarFunctionFromStdFunc<4> fInv4([this](const VectorN<Real, 4>& q) { return funcInverse4(q); });

			switch (i) {
			case 0: return fInv1;
			case 1: return fInv2;
			case 2: return fInv3;
			default: return fInv4;
			}
		}
	};

	// Class modeling general Lorentz transformation, for a moving system with velocity given by a vector
	// Transformation from the rest frame to a frame moving with the given velocity vector.
	class CoordTransfLorentzGeneral : public CoordTransfWithInverse<Vector4Minkowski, Vector4Minkowski, 4>
	{
	private:
		Vec3Cart _velocity;		// velocity vector in units of c
		Real _gamma;          // Lorentz factor
		Real _betaX;       // velocity component along X axis
		Real _betaY;       // velocity component along Y axis
		Real _betaZ;       // velocity component along Z axis

		MatrixNM<Real, 4, 4> _transfMatrix; // transformation matrix for the Lorentz transformation

	public:
		CoordTransfLorentzGeneral(Real speed, const Vec3Cart& velocityDir) 
		{
			_velocity = speed * velocityDir.GetAsUnitVector();
			if (speed <= 0.0 || speed >= 1.0)
				throw std::range_error("Invalid velocity for Lorentz transformation: must be in (0, 1).");

			Real beta2 = speed * speed;
			_gamma = 1.0 / std::sqrt(1.0 - speed * speed);
			_betaX = _velocity.X();
			_betaY = _velocity.Y();
			_betaZ = _velocity.Z();

			_transfMatrix = MatrixNM<Real, 4, 4>{
				_gamma,                                                  -_gamma * _betaX,                                        -_gamma * _betaY,                                        -_gamma * _betaZ,
				-_gamma * _betaX, Real(1.0) + (_gamma - Real(1.0)) * _betaX * _betaX / beta2,       (_gamma - Real(1.0)) * _betaX * _betaY / beta2,       (_gamma - Real(1.0)) * _betaX * _betaZ / beta2,
				-_gamma * _betaY,       (_gamma - Real(1.0)) * _betaY * _betaX / beta2, Real(1.0) + (_gamma - Real(1.0)) * _betaY * _betaY / beta2,       (_gamma - Real(1.0)) * _betaY * _betaZ / beta2,
				-_gamma * _betaZ,       (_gamma - Real(1.0)) * _betaZ * _betaX / beta2,       (_gamma - Real(1.0)) * _betaZ * _betaY / beta2, Real(1.0) + (_gamma - Real(1.0)) * _betaZ * _betaZ / beta2
			};
		}

		Vector4Minkowski transf(const Vector4Minkowski& q) const override
		{
			VectorN<Real, 4> res = _transfMatrix * q;

			return Vector4Minkowski{ res[0], res[1], res[2], res[3] };
		}
		Vector4Minkowski transfInverse(const Vector4Minkowski& q) const override
		{
			VectorN<Real, 4> res = _transfMatrix.GetTranspose() * q;

			return Vector4Minkowski{ res[0], res[1], res[2], res[3] };
		}

		// Scalar functions for each coordinate (optional, for interface compatibility)
		// FIX THIS!!!
// Forward Lorentz transformation (already correct)
		Real func1(const VectorN<Real, 4>& q) const {
			// t' = gamma * (t - v�x)
			Vec3Cart x_vec(q[1], q[2], q[3]);
			Real vdotx = _velocity * x_vec;
			return _gamma * (q[0] - vdotx);
		}
		Real func2(const VectorN<Real, 4>& q) const {
			// x'_1 = x_1 + ((gamma-1)/v^2) (v�x) v_x - gamma t v_x
			Vec3Cart x_vec(q[1], q[2], q[3]);
			Real v2 = _velocity * _velocity;
			Real vdotx = _velocity * x_vec;
			return x_vec.X() + ((_gamma - 1.0) / v2) * vdotx * _velocity.X() - _gamma * q[0] * _velocity.X();
		}
		Real func3(const VectorN<Real, 4>& q) const {
			Vec3Cart x_vec(q[1], q[2], q[3]);
			Real v2 = _velocity * _velocity;
			Real vdotx = _velocity * x_vec;
			return x_vec.Y() + ((_gamma - 1.0) / v2) * vdotx * _velocity.Y() - _gamma * q[0] * _velocity.Y();
		}
		Real func4(const VectorN<Real, 4>& q) const {
			Vec3Cart x_vec(q[1], q[2], q[3]);
			Real v2 = _velocity * _velocity;
			Real vdotx = _velocity * x_vec;
			return x_vec.Z() + ((_gamma - 1.0) / v2) * vdotx * _velocity.Z() - _gamma * q[0] * _velocity.Z();
		}

		// Inverse Lorentz transformation (fixed)
		Real funcInverse1(const VectorN<Real, 4>& q) const {
			// t = gamma * (t' + v�x')
			Vec3Cart x_vec(q[1], q[2], q[3]);
			Real vdotx = _velocity * x_vec;
			return _gamma * (q[0] + vdotx);
		}
		Real funcInverse2(const VectorN<Real, 4>& q) const {
			// x_1 = x'_1 + ((gamma-1)/v^2) (v�x') v_x + gamma t' v_x
			Vec3Cart x_vec(q[1], q[2], q[3]);
			Real v2 = _velocity * _velocity;
			Real vdotx = _velocity * x_vec;
			return x_vec.X() + ((_gamma - 1.0) / v2) * vdotx * _velocity.X() + _gamma * q[0] * _velocity.X();
		}
		Real funcInverse3(const VectorN<Real, 4>& q) const {
			Vec3Cart x_vec(q[1], q[2], q[3]);
			Real v2 = _velocity * _velocity;
			Real vdotx = _velocity * x_vec;
			return x_vec.Y() + ((_gamma - 1.0) / v2) * vdotx * _velocity.Y() + _gamma * q[0] * _velocity.Y();
		}
		Real funcInverse4(const VectorN<Real, 4>& q) const {
			Vec3Cart x_vec(q[1], q[2], q[3]);
			Real v2 = _velocity * _velocity;
			Real vdotx = _velocity * x_vec;
			return x_vec.Z() + ((_gamma - 1.0) / v2) * vdotx * _velocity.Z() + _gamma * q[0] * _velocity.Z();
		}

		const IScalarFunction<4>& coordTransfFunc(int i) const override
		{
			static ScalarFunctionFromStdFunc<4> f1([this](const VectorN<Real, 4>& q) { return func1(q); });
			static ScalarFunctionFromStdFunc<4> f2([this](const VectorN<Real, 4>& q) { return func2(q); });
			static ScalarFunctionFromStdFunc<4> f3([this](const VectorN<Real, 4>& q) { return func3(q); });
			static ScalarFunctionFromStdFunc<4> f4([this](const VectorN<Real, 4>& q) { return func4(q); });

			switch (i) {
			case 0: return f1;
			case 1: return f2;
			case 2: return f3;
			default: return f4;
			}
		}

		const IScalarFunction<4>& inverseCoordTransfFunc(int i) const override
		{
			static ScalarFunctionFromStdFunc<4> fInv1([this](const VectorN<Real, 4>& q) { return funcInverse1(q); });
			static ScalarFunctionFromStdFunc<4> fInv2([this](const VectorN<Real, 4>& q) { return funcInverse2(q); });
			static ScalarFunctionFromStdFunc<4> fInv3([this](const VectorN<Real, 4>& q) { return funcInverse3(q); });
			static ScalarFunctionFromStdFunc<4> fInv4([this](const VectorN<Real, 4>& q) { return funcInverse4(q); });

			switch (i) {
			case 0: return fInv1;
			case 1: return fInv2;
			case 2: return fInv3;
			default: return fInv4;
			}
		}

	};
}

#endif