#if !defined MPL_COLLISION_SIMULATOR_2D_BASE_H
#define MPL_COLLISION_SIMULATOR_2D_BASE_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorTypes.h"

#include <string>


using namespace MML;

namespace MPL
{
	struct Ball2D
	{
	private:
		double _mass, _radius;
		std::string _color;
		Pnt2Cart _position;
		Vec2Cart _velocity;

	public:
		Ball2D(double mass, double radius, std::string color, const Pnt2Cart& position, const Vec2Cart& velocity)
			: _mass(mass), _radius(radius), _color(color), _position(position), _velocity(velocity) {
		}

		double  Mass() const { return _mass; }
		double& Mass() { return _mass; }
		double  Rad()  const { return _radius; }
		double& Rad() { return _radius; }

		Pnt2Cart  Pos() const { return _position; }
		Pnt2Cart& Pos() { return _position; }
		Vec2Cart  V()		const { return _velocity; }
		Vec2Cart& V() { return _velocity; }

		std::string Color() const { return _color; }
		std::string& Color() { return _color; }
	};

	// represents single collision with wall
	struct WallCollisionEvent2D
	{
		int		_ballInd;
		Real	_transferedMomentum; // momentum transferred to the wall by the ball

		WallCollisionEvent2D(int ballInd = -1, Real transferedMomentum = 0.0)
			: _ballInd(ballInd), _transferedMomentum(transferedMomentum) {
		}
	};

	class Wall
	{
		std::string _name;		// name of the wall, e.g. "left", "right", "bottom", "top"
	public:
		Wall() : _name("") {}
		Wall(const std::string& name) : _name(name) {}

		const std::string& Name() const { return _name; }
	};

	class IWallPressureRecorder2D
	{
	public:
		virtual void SetNumberOfWalls(int nWalls) = 0;
		virtual void StartNextStep() = 0;
		virtual void RecordPressureEvent(int wallIndex, const WallCollisionEvent2D& event) = 0;

		~IWallPressureRecorder2D() {}
	};

	class BoxContainerWallPressureRecorder2D : public IWallPressureRecorder2D
	{
		std::vector<Wall> _walls;

		// collision events for each wall, indexed first by step index, then by wall index
		std::vector<std::vector<std::vector<WallCollisionEvent2D>>> _collisionEvents;

		int _currStep;

	public:
		BoxContainerWallPressureRecorder2D(int expectedSteps = 100)
			: _currStep(0)
		{
			_collisionEvents.resize(expectedSteps);

			_walls.resize(4); // left, right, bottom, top
			for (int i = 0; i < expectedSteps; i++)
				_collisionEvents[i].resize(4);
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
		void RecordPressureEvent(int wallIndex, const WallCollisionEvent2D& event)
		{
			// record pressure event for the given ball
			if (wallIndex < 0 || wallIndex >= _walls.size())
				throw std::out_of_range("Wall index out of range");

			// we expect to have a valid step index and wall index
			// with both vectors resized as needed in StartNextStep
			_collisionEvents[_currStep][wallIndex].push_back(event);
		}

		// get all events for wall
		std::vector<std::vector<WallCollisionEvent2D>> getEventsPerStepForWall(int wallIndex)
		{
			if (wallIndex < 0 || wallIndex >= _walls.size())
				throw std::out_of_range("Wall index out of range");

			std::vector<std::vector<WallCollisionEvent2D>> eventsPerStep;

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

		Vector<Real> getTotalTransferedMomentumPerStep()
		{
			Vector<Real> result;
			// Only consider steps that have been recorded
			for (int step = 0; step <= _currStep && step < _collisionEvents.size(); ++step)
			{
				Real sumMomentum = 0.0;
				for (int wallIndex = 0; wallIndex < _walls.size(); ++wallIndex)
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

		Real getTotalTransferedMomentumInStep(int timeStep)
		{
			if (timeStep < 0 || timeStep >= _collisionEvents.size())
				throw std::out_of_range("Time step index out of range");
			Real sumMomentum = 0.0;
			for (int wallIndex = 0; wallIndex < _walls.size(); ++wallIndex)
			{
				for (const auto& event : _collisionEvents[timeStep][wallIndex])
				{
					sumMomentum += event._transferedMomentum;
				}
			}
			return sumMomentum;
		}
	};

	class BallCollisionEvent
	{
	public:
		int _ballIndex1;
		int _ballIndex2;
		int _timeStep;

		BallCollisionEvent(int ballIndex1 = -1, int ballIndex2 = -1, int timeStep = -1)
			: _ballIndex1(ballIndex1), _ballIndex2(ballIndex2), _timeStep(timeStep) 
		{	}
	};

	class ICollisionRecorder2D
	{
	public:
		virtual void StartNextStep() = 0;
		virtual void RecordCollisionEvent(int timeStep, int ballIndex1, int ballIndex2) = 0;
		virtual ~ICollisionRecorder2D() {}
	};

	class BoxContainerCollisionRecorder2D : public ICollisionRecorder2D
	{
		std::vector<std::vector<BallCollisionEvent>> _collisionEvents;
		int _currStep;
	public:
		BoxContainerCollisionRecorder2D(int expectedSteps = 100)
			: _currStep(0)
		{
			_collisionEvents.resize(expectedSteps);
		}

		void StartNextStep()
		{
			_currStep++;
			if (_currStep >= _collisionEvents.size())
				_collisionEvents.resize((int)(_currStep * 1.5)); // increase size by 50%
			_collisionEvents[_currStep].clear(); // clear events for the new step
		}
		void RecordCollisionEvent(int timeStep, int ballIndex1, int ballIndex2)
		{
			if (timeStep < 0 || timeStep >= _collisionEvents.size())
				throw std::out_of_range("Time step index out of range");

			BallCollisionEvent event(ballIndex1, ballIndex2, timeStep);
			
			_collisionEvents[timeStep].push_back(event);
		}
		const std::vector<BallCollisionEvent>& GetCollisionsForStep(int step) const
		{
			if (step < 0 || step >= _collisionEvents.size())
				throw std::out_of_range("Time step index out of range");
			return _collisionEvents[step];
		}
	};
}
#endif // MPL_COLLISION_SIMULATOR_2D_BASE_H