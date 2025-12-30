#if !defined MPL_CONTAINERS_H
#define MPL_CONTAINERS_H

#include "MMLBase.h"

#include "CollisionSimulator2DBase.h"

using namespace MML;

namespace MPL
{
	struct BoxContainer2D
	{
		double _width;
		double _height;

		std::vector<Ball2D> _balls;

		IWallPressureRecorder2D* _wallPressureRecorder = nullptr; // for recording wall pressure events, if needed

		BoxContainer2D() : _width(1000), _height(1000) {}
		BoxContainer2D(double width, double height) : _width(width), _height(height) {}
		BoxContainer2D(double width, double height, IWallPressureRecorder2D* wallPressureRecorder)
			: _width(width), _height(height), _wallPressureRecorder(wallPressureRecorder)
		{
			if (_wallPressureRecorder)
				_wallPressureRecorder->SetNumberOfWalls(4);
		}

		void setWallPressureRecorder(IWallPressureRecorder2D* wallPressureRecorder)
		{
			_wallPressureRecorder = wallPressureRecorder;
			if (_wallPressureRecorder)
				_wallPressureRecorder->SetNumberOfWalls(4);
		}

		// for inherited class, DynamicBoxContainer, that will have movable wall in middle
		// it will calculate the pressure on the wall, and update its position and velocity
		virtual void UpdateContainer() { }

		void		AddBall(const Ball2D& body) { _balls.push_back(body); }
		Ball2D& Ball(int i) { return _balls[i]; }

		// check if ball is out of bounds and handle it
		void CheckAndHandleOutOfBounds(int ballIndex)
		{
			Ball2D& ball = _balls[ballIndex];

			// left wall collision
			if (ball.Pos().X() < ball.Rad() && ball.V().X() < 0)
			{
				ball.Pos().X() = ball.Rad() + (ball.Rad() - ball.Pos().X()); // Get back to box!
				ball.V().X() *= -1;
				if (_wallPressureRecorder)
				{
					// record pressure event for the left wall
					Real transferedMomentum = ball.Mass() * fabs(ball.V().X());
					_wallPressureRecorder->RecordPressureEvent(0, WallCollisionEvent2D(ballIndex, transferedMomentum));
				}
			}
			// right wall collision
			if (ball.Pos().X() > _width - ball.Rad() && ball.V().X() > 0)
			{
				ball.Pos().X() -= (ball.Pos().X() + ball.Rad()) - _width;
				ball.V().X() *= -1;

				if (_wallPressureRecorder)
				{
					// record pressure event for the right wall
					Real transferedMomentum = ball.Mass() * fabs(ball.V().X());
					_wallPressureRecorder->RecordPressureEvent(1, WallCollisionEvent2D(ballIndex, transferedMomentum));
				}
			}
			// bottom wall collision
			if (ball.Pos().Y() < ball.Rad() && ball.V().Y() < 0)
			{
				ball.Pos().Y() = ball.Rad() + (ball.Rad() - ball.Pos().Y());
				ball.V().Y() *= -1;

				if (_wallPressureRecorder)
				{
					// record pressure event for the bottom wall
					Real transferedMomentum = ball.Mass() * fabs(ball.V().Y());
					_wallPressureRecorder->RecordPressureEvent(2, WallCollisionEvent2D(ballIndex, transferedMomentum));
				}
			}
			// top wall collision
			if (ball.Pos().Y() > _height - ball.Rad() && ball.V().Y() > 0)
			{
				ball.Pos().Y() -= (ball.Pos().Y() + ball.Rad()) - _height;
				ball.V().Y() *= -1;

				if (_wallPressureRecorder)
				{
					// record pressure event for the top wall
					Real transferedMomentum = ball.Mass() * fabs(ball.V().Y());
					_wallPressureRecorder->RecordPressureEvent(3, WallCollisionEvent2D(ballIndex, transferedMomentum));
				}
			}
		}

		// fill given matrix with the number of balls in each subdivided container
		void SetNumBallsInSubdividedContainer(int nRows, int nCols, Matrix<Vector<int>>& M)
		{
			double cellWidth = _width / nCols;
			double cellHeight = _height / nRows;

			for (int i = 0; i < _balls.size(); i++)
			{
				Ball2D& ball = _balls[i];

				int row = static_cast<int>(ball.Pos().Y() / cellHeight);
				int col = static_cast<int>(ball.Pos().X() / cellWidth);

				if (row >= 0 && row < nRows && col >= 0 && col < nCols)
				{
					M(row, col).push_back(i);
				}
			}
		}
	};

	// Container with the small hole in the middle, where balls can go through

	// Container with thi movable wall in the middle, which can be moved left or right

	// Container representing billiard table, with holes in the corners and sides
}
#endif // MPL_CONTAINERS_H
