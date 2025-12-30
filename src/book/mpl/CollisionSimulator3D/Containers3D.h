#if !defined MPL_CONTAINERS_3D_H
#define MPL_CONTAINERS_3D_H

#include "MMLBase.h"

#include "base/Matrix3D.h"

#include "CollisionSimulator3DBase.h"

using namespace MML;

namespace MPL
{
	struct BoxContainer3D
	{
		double _width;
		double _height;
		double _depth;

		std::vector<Ball3D> _balls;

		IWallPressureRecorder3D* _wallPressureRecorder = nullptr; // for recording wall pressure events, if needed

		BoxContainer3D() : _width(1000), _height(1000), _depth(1000) {}
		BoxContainer3D(double width, double height, double depth) : _width(width), _height(height), _depth(depth) {}
		BoxContainer3D(double width, double height, double depth, IWallPressureRecorder3D* wallPressureRecorder)
			: _width(width), _height(height), _depth(depth), _wallPressureRecorder(wallPressureRecorder) { }

		void SetWallPressureRecorder(IWallPressureRecorder3D* wallPressureRecorder)
		{
			_wallPressureRecorder = wallPressureRecorder;
		}

		void		AddBall(const Ball3D& body) { _balls.push_back(body); }
		Ball3D& Ball(int i)									{ return _balls[i]; }

		std::vector<std::string> getBallColors()
		{
			std::vector<std::string> colors;
			for (const auto& ball : _balls)
			{
				colors.push_back(ball.Color());
			}
			return colors;
		}
		std::vector<Real>				 getBallRadii()
		{
			std::vector<Real> radii;
			for (const auto& ball : _balls)
			{
				radii.push_back(ball.Rad());
			}
			return radii;
		}

		void CheckAndHandleOutOfBounds(int ballIndex)
		{
			Ball3D& ball = _balls[ballIndex];

			if (ball.Pos().X() < ball.Rad() && ball.V().X() < 0)		// checking if ball is out of bounds
			{
				ball.Pos().X() = ball.Rad() + (ball.Rad() - ball.Pos().X());		// moving it back to box
				ball.V().X() *= -1;
				if (_wallPressureRecorder)
				{
					// FIX this alculation of transfered momentum!!!
					Real transferedMomentum = ball.Mass() * std::abs(ball.V().X());
					_wallPressureRecorder->RecordPressureEvent(0, WallCollisionEvent3D(ballIndex, transferedMomentum));
				}
			}

			if (ball.Pos().X() > _width - ball.Rad() && ball.V().X() > 0)
			{
				ball.Pos().X() -= (ball.Pos().X() + ball.Rad()) - _width;
				ball.V().X() *= -1;
			}

			if (ball.Pos().Y() < ball.Rad() && ball.V().Y() < 0)
			{
				ball.Pos().Y() = ball.Rad() + (ball.Rad() - ball.Pos().Y());
				ball.V().Y() *= -1;
			}

			if (ball.Pos().Y() > _height - ball.Rad() && ball.V().Y() > 0)
			{
				ball.Pos().Y() -= (ball.Pos().Y() + ball.Rad()) - _height;
				ball.V().Y() *= -1;
			}

			if (ball.Pos().Z() < ball.Rad() && ball.V().Z() < 0)
			{
				ball.Pos().Z() = ball.Rad() + (ball.Rad() - ball.Pos().Z());
				ball.V().Z() *= -1;
			}

			if (ball.Pos().Z() > _depth - ball.Rad() && ball.V().Z() > 0)
			{
				ball.Pos().Z() -= (ball.Pos().Z() + ball.Rad()) - _depth;
				ball.V().Z() *= -1;
			}
		}

		void SetNumBallsInSubdividedContainer(int n1, int n2, int n3, Matrix3D<Vector<int>>& M)
		{
			double cellWidth = _width / n1;
			double cellHeight = _height / n2;
			double cellDepth = _depth / n3;

			for (int i = 0; i < _balls.size(); i++)
			{
				Ball3D& ball = _balls[i];

				int row = static_cast<int>(ball.Pos().Y() / cellHeight);
				int col = static_cast<int>(ball.Pos().X() / cellWidth);
				int depth = static_cast<int>(ball.Pos().Z() / cellDepth);

				if (row >= 0 && row < n1 && col >= 0 && col < n2 && depth >= 0 && depth < n3)
				{
					M(row, col, depth).push_back(i);
				}
			}
		}
	};
}
#endif // MPL_CONTAINERS_H
