#if !defined MPL_COLLISION_SIMULATOR_3D_BASE_H
#define MPL_COLLISION_SIMULATOR_3D_BASE_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorTypes.h"

#include <string>


using namespace MML;

namespace MPL
{
	struct Ball3D
	{
	private:
		double _mass;
		double _radius;
		std::string _color;
		Pnt3Cart _position;
		Vec3Cart _velocity;

	public:
		Ball3D(double mass, double radius, std::string color, const Pnt3Cart& position,
			const Vec3Cart& velocity)
			: _mass(mass), _radius(radius), _color(color), _position(position), _velocity(velocity)
		{
		}

		double  Mass() const { return _mass; }
		double& Mass() { return _mass; }

		double  Rad() const { return _radius; }
		double& Rad() { return _radius; }

		std::string Color() const { return _color; }
		std::string& Color() { return _color; }

		Pnt3Cart  Pos() const { return _position; }
		Pnt3Cart& Pos() { return _position; }

		Vec3Cart  V() const { return _velocity; }
		Vec3Cart& V() { return _velocity; }
	};

	// represents single collision with wall
	struct WallCollisionEvent3D
	{
		int _ballInd;
		Real _transferedMomentum; // momentum transferred to the wall by the ball

		WallCollisionEvent3D(int ballInd = -1, Real transferedMomentum = 0.0)
			: _ballInd(ballInd), _transferedMomentum(transferedMomentum) {
		}
	};

	class IWallPressureRecorder3D
	{
	public:
		virtual void SetNumberOfWalls(int nWalls) = 0;
		virtual void StartNextStep() = 0;
		virtual void RecordPressureEvent(int wallIndex, const WallCollisionEvent3D& event) = 0;
	};

	class Wall3D
	{
		std::string _name; // name of the wall, e.g. "left", "right", "bottom", "top"
	public:
		Wall3D() : _name("") {}
		Wall3D(const std::string& name) : _name(name) {}

		const std::string& Name() const { return _name; }
	};

	class BoxWallPressureRecorder3D : public IWallPressureRecorder3D
	{
		std::vector<Wall3D> _walls; // walls of the box, left, right, bottom, top

		// collision events for each wall, indexed first by step index, then by wall index
		std::vector<std::vector<std::vector<WallCollisionEvent3D>>> _collisionEvents;

		int _currStep;

	public:
		BoxWallPressureRecorder3D(int expectedSteps = 100)
			: _currStep(0)
		{
			_collisionEvents.resize(expectedSteps);

			_walls.resize(6); // left, right, front, back, bottom, top
			for (int i = 0; i < expectedSteps; i++)
				_collisionEvents[i].resize(6);
		}

		void SetNumberOfWalls(int nWalls)
		{
			_walls.resize(nWalls);
		}
		void StartNextStep()
		{
			_currStep++;
			if (_currStep >= _collisionEvents.size())
				_collisionEvents.resize((int)(_currStep * 1.5));		// increase size by 50%

			// make sure we have space for all walls in the current step
			if (_collisionEvents[_currStep].size() < _walls.size())
				_collisionEvents[_currStep].resize(_walls.size());
		}
		void RecordPressureEvent(int wallIndex, const WallCollisionEvent3D& event)
		{
			// record pressure event for the given ball
			if (wallIndex < 0 || wallIndex >= _walls.size())
				throw std::out_of_range("Wall index out of range");

			// we expect to have a valid step index and wall index
			// with both vectors resized as needed in StartNextStep
			_collisionEvents[_currStep][wallIndex].push_back(event);
		}

		// get all events for wall
		std::vector<std::vector<WallCollisionEvent3D>> getEventsPerStepForWall(int wallIndex)
		{
			if (wallIndex < 0 || wallIndex >= _walls.size())
				throw std::out_of_range("Wall index out of range");

			std::vector<std::vector<WallCollisionEvent3D>> eventsPerStep;

			// resizing to number of recorded steps
			eventsPerStep.resize(_collisionEvents[0].size());

			for (int i = 0; i < _collisionEvents.size(); i++)
			{
				for (int j = 0; j < _collisionEvents[i][wallIndex].size(); j++)
				{
					eventsPerStep[i].push_back(_collisionEvents[i][wallIndex][j]);
				}
			}

			return eventsPerStep;
		}

		Vector<Real> getTranferedMomentumPerStepForWall(int wallIndex)
		{
			Vector<Real> result;
			if (wallIndex < 0 || wallIndex >= _walls.size())
				throw std::out_of_range("Wall index out of range");

			// Only consider steps that have been recorded
			for (int step = 0; step <= _currStep && step < _collisionEvents.size(); ++step)
			{
				Real sumMomentum = 0.0;
				if (wallIndex < _collisionEvents[step].size())
				{
					for (const auto& event : _collisionEvents[step][wallIndex])
					{
						sumMomentum += event._transferedMomentum;
					}
				}
				result.push_back(sumMomentum);
			}
			return result;
		}
	};
}
#endif // MPL_COLLISION_SIMULATOR_2D_BASE_H