#if !defined MPL_COLLISION_SIMULATOR_2D_H
#define MPL_COLLISION_SIMULATOR_2D_H

#include "MMLBase.h"

#include <thread>

#include "../MPLBase.h"

#include "base/Matrix.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"
#include "tools/ThreadPool.h"

#include "CollisionSimulator2DBase.h"
#include "Containers2D.h"
#include "SimResultsCollSim2D.h"

using namespace MML;

namespace MPL
{
	class CollisionSimulator2D
	{
		BoxContainer2D _box;

		// M is a matrix of vectors, where each vector contains indices of balls in that subcontainer
		Matrix<Vector<int>> M;

		// for recording wall pressure events, if needed
		IWallPressureRecorder2D* _wallPressureRecorder = nullptr;

		// for recording collision events, if needed
		ICollisionRecorder2D* _collisionRecorder = nullptr;

	public:
		CollisionSimulator2D() : M(2, 2) {}
		CollisionSimulator2D(const BoxContainer2D& box) : _box(box), M(2, 2) {}
		CollisionSimulator2D(const BoxContainer2D& box, int nRows, int nCols) : _box(box), M(nRows, nCols) {}
		CollisionSimulator2D(const BoxContainer2D& box, IWallPressureRecorder2D* wallPressureRecorder)
			: _box(box), _wallPressureRecorder(wallPressureRecorder) {
		}
		CollisionSimulator2D(const BoxContainer2D& box, ICollisionRecorder2D* collisionRecorder)
			: _box(box), _collisionRecorder(collisionRecorder) {
		}
		CollisionSimulator2D(const BoxContainer2D& box, int nRows, int nCols, IWallPressureRecorder2D* wallPressureRecorder)
			: _box(box), M(nRows, nCols), _wallPressureRecorder(wallPressureRecorder) {
		}
		CollisionSimulator2D(const BoxContainer2D& box, int nRows, int nCols, ICollisionRecorder2D* collisionRecorder)
			: _box(box), M(nRows, nCols), _collisionRecorder(collisionRecorder) {
		}
		CollisionSimulator2D(const BoxContainer2D& box, int nRows, int nCols, IWallPressureRecorder2D* wallPressureRecorder, ICollisionRecorder2D* collisionRecorder)
			: _box(box), M(nRows, nCols), _wallPressureRecorder(wallPressureRecorder), _collisionRecorder(collisionRecorder) {
		}

		void setWallPressureRecorder(IWallPressureRecorder2D* wallPressureRecorder)
		{
			_wallPressureRecorder = wallPressureRecorder;
			if (_wallPressureRecorder)
				_wallPressureRecorder->SetNumberOfWalls(4);
		}
		void setCollisionRecorder(ICollisionRecorder2D* collisionRecorder)
		{
			_collisionRecorder = collisionRecorder;
		}

		double	DistBalls(int i, int j)
		{
			return _box._balls[i].Pos().Dist(_box._balls[j].Pos());
		}
		bool		HasBallsCollided(int i, int j)
		{
			// ako je udaljenost izmedju njihovih centara manja od zbroja radijusa
			if (DistBalls(i, j) < _box._balls[i].Rad() + _box._balls[j].Rad())
				return true;
			else
				return false;
		}

		void HandleCollision(int stepNum, int m, int n, double dt)
		{
			if (_collisionRecorder != nullptr)
				_collisionRecorder->RecordCollisionEvent(stepNum, m, n);

			Ball2D& ball1 = _box._balls[m];
			Ball2D& ball2 = _box._balls[n];

			// calculating point where they were before collision
			// (if there was collision with the box wall, then calc.pos. will be outside the box
			// but it doesn't matter, since we need only direction, ie. velocity, to calculate exact collision point)
			Pnt2Cart x10 = ball1.Pos() - ball1.V() * dt;
			Pnt2Cart x20 = ball2.Pos() - ball2.V() * dt;

			Vec2Cart dx0(x10, x20);
			Vec2Cart dv(ball2.V() - ball1.V());

			// first, we have to calculate exact moment of collision, and balls positions then
			double A = dv * dv;
			double B = 2 * dx0 * dv;
			double C = dx0 * dx0 - POW2(ball2.Rad() + ball1.Rad());

			double t1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
			double t2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);

			double tCollision = t1 < t2 ? t1 : t2;

			// calculating position of balls at the point of collision (moving them backwards)
			Pnt2Cart x1 = ball1.Pos() + (tCollision - dt) * ball1.V();
			Pnt2Cart x2 = ball2.Pos() + (tCollision - dt) * ball2.V();

			// based on https://en.wikipedia.org/wiki/Elastic_collision - calculating new velocities after collision
			double   m1 = ball1.Mass(), m2 = ball2.Mass();
			Vec2Cart v1 = ball1.V(), v2 = ball2.V();

			Vec2Cart v1_v2 = v1 - v2;
			Vec2Cart x1_x2(x2, x1);

			Vec2Cart v1_new = v1 - 2 * m2 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec2Cart(x2, x1);
			Vec2Cart v2_new = v2 - 2 * m1 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec2Cart(x1, x2);

			ball1.V() = v1_new;
			ball2.V() = v2_new;

			// adjusting new ball positions (taking into account only time after collision)
			ball1.Pos() = x1 + ball1.V() * (dt - tCollision);
			ball2.Pos() = x2 + ball2.V() * (dt - tCollision);
		}

		void SimulateOneStepExact(double dt, int stepNum)
		{
			int NumBalls = _box._balls.size();
			
			for (int i = 0; i < NumBalls; i++)   // first, update all ball's positions, and handle out of bounds
			{
				_box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;

				_box.CheckAndHandleOutOfBounds(i);
			}

			// check for collisions, and handle it if there is one
			for (int m = 0; m < NumBalls - 1; m++) {
				for (int n = m + 1; n < NumBalls; n++)
				{
					if (HasBallsCollided(m, n))
						HandleCollision(stepNum, m, n, dt);
				}
			}
		}
		void SimulateOneStepFast(double dt, int stepNum)
		{
			int NumBalls = _box._balls.size();

			for (int i = 0; i < NumBalls; i++)
				_box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;

			_box.SetNumBallsInSubdividedContainer(M.RowNum(), M.ColNum(), M);
			
			// handle out of bounds situations for all balls
			// TODO - do this ONLY for balls that are in subcontainers that are at the border of the box!!!
			// even though ... this is LINEAR operation, so it is not a big deal
			for (int i = 0; i < NumBalls; i++)
				_box.CheckAndHandleOutOfBounds(i);

			// now, we check for collisions, separately for each subcontainer
			for (int i = 0; i < M.RowNum(); i++) {
				for (int j = 0; j < M.ColNum(); j++)
				{
					if (M(i, j).size() > 1) // if there are at least two balls in this subcontainer
					{
						for (int m = 0; m < M(i, j).size() - 1; m++) {
							for (int n = m + 1; n < M(i, j).size(); n++)
							{
								int ballIndex1 = M(i, j)[m];
								int ballIndex2 = M(i, j)[n];

								if (HasBallsCollided(ballIndex1, ballIndex2))
									HandleCollision(stepNum, ballIndex1, ballIndex2, dt);
							}
						}
						M(i, j).Clear();	// we are done with this subcontainer, clear it for next iteration
					}
				}
			}
		}
		void SimulateOneStepFastMultithread(double dt, int stepNum)
		{
			_box.SetNumBallsInSubdividedContainer(M.RowNum(), M.ColNum(), M);

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
			for (int i = 0; i < M.RowNum(); i++)
			{
				for (int j = 0; j < M.ColNum(); j++)
				{
					if (M(i, j).size() > 1) // if there are at least two balls in this subcontainer
					{
						threads.emplace_back([this, i, j, stepNum, dt]()
							{
								for (int m = 0; m < M(i, j).size() - 1; m++) {
									for (int n = m + 1; n < M(i, j).size(); n++)
									{
										int ballIndex1 = M(i, j)[m];
										int ballIndex2 = M(i, j)[n];

										if (HasBallsCollided(ballIndex1, ballIndex2))
											HandleCollision(stepNum, ballIndex1, ballIndex2, dt);
									}
								}
								M(i, j).Clear();
							});
					}
				}
			}
			for (auto& t : threads) t.join();
		}
		void SimulateOneStepFastMultithreadThreadPool(double dt, int stepNum)
		{
			_box.SetNumBallsInSubdividedContainer(M.RowNum(), M.ColNum(), M);

			int NumBalls = _box._balls.size();

			// Update positions and handle out of bounds
			for (int i = 0; i < NumBalls; i++) {
				_box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;
				_box.CheckAndHandleOutOfBounds(i);
			}

			unsigned int numThreads = std::thread::hardware_concurrency();
			if (numThreads == 0) 
				numThreads = 2; 

			ThreadPool pool(numThreads);

			std::atomic<int> tasksRemaining = 0;

			for (int i = 0; i < M.RowNum(); i++) {
				for (int j = 0; j < M.ColNum(); j++) {
					if (M(i, j).size() > 1) {
						tasksRemaining++;
						pool.enqueue([this, i, j, stepNum, dt, &tasksRemaining]()
							{
								for (int m = 0; m < M(i, j).size() - 1; m++) {
									for (int n = m + 1; n < M(i, j).size(); n++) 
									{
										int ballIndex1 = M(i, j)[m];
										int ballIndex2 = M(i, j)[n];
										if (HasBallsCollided(ballIndex1, ballIndex2))
											HandleCollision(stepNum, ballIndex1, ballIndex2, dt);
									}
								}
								M(i, j).Clear();
								--tasksRemaining;
							});
					}
				}
			}

			// Wait for all tasks to finish
			while (tasksRemaining > 0) {
				std::this_thread::yield();
			}
		}

		SimResultsCollSim2D Simulate(int numSteps, double timeStep, CollisionSimulatorRunType runType = CollisionSimulatorRunType::RunTypeFastMultithread)
		{
			int			numBalls = _box._balls.size();
			double	dt = timeStep;

			SimResultsCollSim2D results(_box._balls);

			for (int i = 0; i < numSteps; i++)
			{
				if (_wallPressureRecorder != nullptr)
					_wallPressureRecorder->StartNextStep();

				for (int j = 0; j < numBalls; j++)
				{
					results.BallPosList[j].push_back(_box._balls[j].Pos());
					results.BallVelList[j].push_back(_box._balls[j].V());
				}

				switch (runType)
				{	
				case MPL::RunTypeExact:
					SimulateOneStepExact(dt, i);
					break;
				case MPL::RunTypeFast:
					SimulateOneStepFast(dt, i);
					break;
				case MPL::RunTypeFastMultithread:
					SimulateOneStepFastMultithread(dt, i);
					break;
				case MPL::RunTypeFastMultithreadTP:
					SimulateOneStepFastMultithreadThreadPool(dt, i);
					break;
				default:
					break;
				}
			}

			return results;
		}

		bool Serialize(std::string fileName, const SimResultsCollSim2D& simResults, Real dT, int saveEveryNSteps = 1)
		{
			// preparing data for Serializer class
			int numBalls = static_cast<int>(_box._balls.size());
			std::vector<std::string> ballColors;
			std::vector<Real> ballRadius;

			for (const auto& ball : _box._balls)
			{
				ballColors.push_back(ball.Color());
				ballRadius.push_back(ball.Rad());
			}

			auto result = Serializer::SaveParticleSimulation2D(fileName, numBalls, _box._width, _box._height, simResults.BallPosList, ballColors, ballRadius, dT, saveEveryNSteps);
			return result.success;
		}
	};
}
#endif