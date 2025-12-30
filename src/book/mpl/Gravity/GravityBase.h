#if !defined MPL_GRAVITY_BASE_H
#define MPL_GRAVITY_BASE_H

#include "MMLBase.h"

#include "base/VectorTypes.h"

using namespace MML;

namespace MPL
{
	struct GravityMass
	{
		Real	_mass;
		Real	_radius;
		std::string _color;
	};

	class GravityBodyState
	{
		Real	_mass;
		Real	_radius;
		std::string _color;

		Vec3Cart	_position;
		Vec3Cart	_velocity;

	public:
		GravityBodyState() : _mass(1), _position{ 0.0, 0.0, 0.0 }, _velocity{ 0.0, 0.0, 0.0 }, _color("Black"), _radius(10)
		{	}
		GravityBodyState(const Real& mass, const Vec3Cart& position)
			: _mass(mass), _position(position), _velocity{ 0.0, 0.0, 0.0 }, _color("Black"), _radius(10)
		{	}
		GravityBodyState(const Real& mass, const Vec3Cart& position, const Vec3Cart& velocity)
			: _mass(mass), _position(position), _velocity(velocity), _color("Black"), _radius(10)
		{	}
		GravityBodyState(const Real& mass, const Vec3Cart& position, const Vec3Cart& velocity, const std::string& color, const Real& radius)
			: _mass(mass), _position(position), _velocity(velocity), _color(color), _radius(radius)
		{	}

		Real	Mass()  const { return _mass; }
		Real& Mass()				{ return _mass; }
		Real	Radius() const	{ return _radius; }
		Real& Radius()				{ return _radius; }
		std::string Color() const { return _color; }
		std::string& Color()			{ return _color; }

		Vec3Cart	R() const { return _position; }
		Vec3Cart& R()				{ return _position; }
		Vec3Cart	V() const { return _velocity; }
		Vec3Cart& V()				{ return _velocity; }
	};

	class GravityBodyState2D
	{
	private:
		Real	_mass;
		Real	_radius;
		std::string _color;

		Vec2Cart	_position;
		Vec2Cart	_velocity;

	public:
		GravityBodyState2D() : _mass(1), _position{ 0.0, 0.0 }, _velocity{ 0.0, 0.0 }, _color("Black"), _radius(10)
		{	}
		GravityBodyState2D(const Real& mass, const Vec2Cart& position)
			: _mass(mass), _position(position), _velocity{ 0.0, 0.0 }, _color("Black"), _radius(10)
		{	}
		GravityBodyState2D(const Real& mass, const Vec2Cart& position, const Vec2Cart& velocity)
			: _mass(mass), _position(position), _velocity(velocity), _color("Black"), _radius(10)
		{	}
		GravityBodyState2D(const Real& mass, const Vec2Cart& position, const Vec2Cart& velocity, 
											 const Real& radius, const std::string& color)
			: _mass(mass), _position(position), _velocity(velocity), _color(color), _radius(radius)
		{	}

		Real	Mass()  const { return _mass; }
		Real& Mass()				{ return _mass; }
		Real Radius() const { return _radius; }
		Real& Radius()			{ return _radius; }
		std::string Color() const { return _color; }
		std::string& Color()			{ return _color; }

		Vec2Cart	R() const { return _position; }
		Vec2Cart& R()				{ return _position; }
		Vec2Cart	V() const { return _velocity; }
		Vec2Cart& V()				{ return _velocity; }
	};

	class GravityPotentialField : public IScalarFunction<3>
	{
	private:
		Real _G;
		Real _Mass;
		Vec3Cart _Position;
	public:
		GravityPotentialField(Real G, Real Mass, Vec3Cart Position) : _G(G), _Mass(Mass), _Position(Position) {}

		Real operator()(const VectorN<Real, 3>& x) const
		{
			Real pot = 0.0;
			Vec3Cart r = x - _Position;
			Real rNorm = r.NormL2();

			if (rNorm > 0.0)
				pot = -_G * _Mass / rNorm;
			else
				pot = Constants::PosInf; 

			return pot;
		}
	};

	class GravityForceField : public IVectorFunction<3>
	{
	private:
		Real _G;
		Real _Mass;
		Vec3Cart _Position;
	public:
		GravityForceField(Real G, Real Mass, Vec3Cart Position) : _G(G), _Mass(Mass), _Position(Position) {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			VectorN<Real, 3> force(0.0);

			Vec3Cart r = x - _Position;
			Real rNorm = r.NormL2();

			if (rNorm > 0.0)
				force = -_G * _Mass / (rNorm * rNorm * rNorm) * r;
			else
				force = Vec3Cart(0.0, 0.0, 0.0); // Avoid division by zero

			return force;
		}
	};
} // namespace MPL

#endif