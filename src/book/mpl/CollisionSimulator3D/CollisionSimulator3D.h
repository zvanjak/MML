#if !defined MPL_COLLISION_SIMULATOR_3D_H
#define MPL_COLLISION_SIMULATOR_3D_H

#include "MMLBase.h"

#include <thread>

#include "base/Matrix3D.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"

#include "../MPLBase.h"

#include "CollisionSimulator3DBase.h"
#include "Containers3D.h"
#include "SimResultsCollSim3D.h"

using namespace MML;

namespace MPL
{
	class CollisionSimulator3D
	{
		BoxContainer3D _box;

		// 3D matrix for subdividing space into subcontainers, for faster collision detection
		Matrix3D<Vector<int>> M;

		// for recording wall pressure events, if needed
		IWallPressureRecorder3D* _wallPressureRecorder = nullptr;

	public:
		CollisionSimulator3D() : M(2,2,2) {}
		CollisionSimulator3D(const BoxContainer3D& box) : _box(box), M(2, 2, 2) {}
		CollisionSimulator3D(const BoxContainer3D& box, int n1, int n2, int n3) : _box(box), M(n1, n2, n3) {}
		CollisionSimulator3D(const BoxContainer3D& box, IWallPressureRecorder3D* wallPressureRecorder)
			: _box(box), _wallPressureRecorder(wallPressureRecorder) {
		}
		CollisionSimulator3D(const BoxContainer3D& box, int n1, int n2, int n3, IWallPressureRecorder3D* wallPressureRecorder)
			: _box(box), M(n1, n2, n3), _wallPressureRecorder(wallPressureRecorder) {
		}

		void setWallPressureRecorder(IWallPressureRecorder3D* wallPressureRecorder)
		{
			_wallPressureRecorder = wallPressureRecorder;
			if (_wallPressureRecorder)
				_wallPressureRecorder->SetNumberOfWalls(4);
		}

		double DistBalls(int i, int j)
		{
			return _box._balls[i].Pos().Dist(_box._balls[j].Pos());
		}
		bool	 HasBallsCollided(int i, int j)
		{
			// ako je udaljenost izmedju njihovih centara manja od zbroja radijusa
			if (DistBalls(i, j) < _box._balls[i].Rad() + _box._balls[j].Rad())
				return true;
			else
				return false;
		}

		void HandleCollision(int m, int n, double dt)
		{
			//std::cout << "Collision detected between balls " << m << " and " << n << std::endl;

			Ball3D& ball1 = _box._balls[m];
			Ball3D& ball2 = _box._balls[n];

			// calculating point where they were before collision
			Pnt3Cart x10 = ball1.Pos() - ball1.V() * dt;
			Pnt3Cart x20 = ball2.Pos() - ball2.V() * dt;

			Vec3Cart dx0(x10, x20);
			Vec3Cart dv(ball2.V() - ball1.V());

			double A = dv * dv;
			double B = 2 * dx0 * dv;
			double C = dx0 * dx0 - POW2(ball2.Rad() + ball1.Rad());

			double t1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
			double t2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);

			double tCollision = t1 < t2 ? t1 : t2;

			// calculating position of balls at the point of collision (moving them backwards)
			Pnt3Cart x1 = ball1.Pos() + (tCollision - dt) * ball1.V();
			Pnt3Cart x2 = ball2.Pos() + (tCollision - dt) * ball2.V();

			// https://en.wikipedia.org/wiki/Elastic_collision - calculating new velocities after collision
			double   m1 = ball1.Mass(), m2 = ball2.Mass();

			Vec3Cart v1 = ball1.V(), v2 = ball2.V();

			Vec3Cart v1_v2 = v1 - v2;
			Vec3Cart x1_x2(x2, x1);

			Vec3Cart v1_new = v1 - 2 * m2 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec3Cart(x2, x1);
			Vec3Cart v2_new = v2 - 2 * m1 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec3Cart(x1, x2);

			ball1.V() = v1_new;
			ball2.V() = v2_new;

			// adjusting new ball positions
			ball1.Pos() = x1 + ball1.V() * (dt - tCollision);
			ball2.Pos() = x2 + ball2.V() * (dt - tCollision);
		}

		void SimulateOneStepExact(double dt)
		{
			// check collisions at start configuration!!!
			// that is also an end configuration for previous step
			int NumBalls = _box._balls.size();

			// first, update all ball's positions, and handle out of bounds
			for (int i = 0; i < NumBalls; i++)
			{
				_box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;

				_box.CheckAndHandleOutOfBounds(i);
			}

			// check for collisions, and handle it if there is one
			for (int m = 0; m < NumBalls - 1; m++) {
				for (int n = m + 1; n < NumBalls; n++)
				{
					if (HasBallsCollided(m, n))
						HandleCollision(m, n, dt);
				}
			}
		}
		void SimulateOneStepFast(double dt)
		{
			int NumBalls = _box._balls.size();

			// first, update all ball's positions
			for (int i = 0; i < NumBalls; i++)
				_box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;

			_box.SetNumBallsInSubdividedContainer(M.dim1(), M.dim2(), M.dim3(), M);

			// handle out of bounds situations for all balls
			for (int i = 0; i < NumBalls; i++)
				_box.CheckAndHandleOutOfBounds(i);

			// now, we have to check for collisions only in subcontainers where there are balls
			for (int i = 0; i < M.dim1(); i++) {
				for (int j = 0; j < M.dim2(); j++) {
					for (int k = 0; k < M.dim3(); k++)
					{
						if (M(i, j, k).size() > 1) // if there are at least two balls in this subcontainer
						{
							for (int m = 0; m < M(i, j, k).size() - 1; m++) {
								for (int n = m + 1; n < M(i, j, k).size(); n++)
								{
									int ballIndex1 = M(i, j, k)[m];
									int ballIndex2 = M(i, j, k)[n];

									if (HasBallsCollided(ballIndex1, ballIndex2))
										HandleCollision(ballIndex1, ballIndex2, dt);
								}
							}
							M(i, j, k).Clear();	// we are done with this subcontainer, clear it for next iteration
						}
					}
				}
			}
		}
		void SimulateOneStepFastMultithread(double dt)
		{
			_box.SetNumBallsInSubdividedContainer(M.dim1(), M.dim2(), M.dim3(), M);

			// check collisions at start configuration!!!
			// that is also an end configuration for previous step
			int NumBalls = _box._balls.size();

			// first, update all ball's positions, and handle out of bounds
			for (int i = 0; i < NumBalls; i++)
			{
				_box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;

				_box.CheckAndHandleOutOfBounds(i);
			}

			// now, we have to check for collisions only in subcontainers where there are balls
			std::vector<std::thread> threads;
			for (int i = 0; i < M.dim1(); i++) {
				for (int j = 0; j < M.dim2(); j++) {
					for (int k = 0; k < M.dim3(); k++)
					{
						// if there are at least two balls in this subcontainer
						if (M(i, j, k).size() > 1) // if there are at least two balls in this subcontainer
						{
							threads.emplace_back([this, i, j, k, dt]()
								{
									for (int m = 0; m < M(i, j, k).size() - 1; m++) {
										for (int n = m + 1; n < M(i, j, k).size(); n++)
										{
											int ballIndex1 = M(i, j, k)[m];
											int ballIndex2 = M(i, j, k)[n];

											if (HasBallsCollided(ballIndex1, ballIndex2))
												HandleCollision(ballIndex1, ballIndex2, dt);
										}
									}
									M(i, j, k).Clear();
								});
						}
					}
				}
			}
			for (auto& t : threads) t.join();
		}

		SimResultsCollSim3D Simulate(int numSteps, double timeStep, bool useFast = true)
		{
			int numBalls = _box._balls.size();
			double dt = timeStep; 

			SimResultsCollSim3D results(numBalls);

			for (int i = 0; i < numSteps; i++)
			{
				if (_wallPressureRecorder != nullptr)
					_wallPressureRecorder->StartNextStep();

				for (int j = 0; j < numBalls; j++)  // saving current position and velocity of each ball
				{
					results.BallPosList[j].push_back(_box._balls[j].Pos());
					results.BallVelList[j].push_back(_box._balls[j].V()); 
				}

				if(useFast)
					//SimulateOneStepFast(dt);
					SimulateOneStepFastMultithread(dt);
				else
					SimulateOneStepExact(dt);
			}
			return results;
		}

		bool Serialize(std::string fileName, const SimResultsCollSim3D& simResults, double dT, int saveEveryNSteps = 1)
		{
			// preparing data for Serializer class
			int numBalls = _box._balls.size();
			std::vector<std::string> ballColors;
			std::vector<Real> ballRadius;

			for (const auto& ball : _box._balls)
			{
				ballColors.push_back(ball.Color());
				ballRadius.push_back(ball.Rad());
			}

			auto result = Serializer::SaveParticleSimulation3D(fileName, numBalls, _box._width, _box._height, _box._depth, simResults.BallPosList, ballColors, ballRadius, dT, saveEveryNSteps);
			return result.success;
		}
	};
}

#endif