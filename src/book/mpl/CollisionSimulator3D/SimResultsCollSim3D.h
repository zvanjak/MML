#if !defined MPL_SIMRESULTS_COLLSIM_3D_H
#define MPL_SIMRESULTS_COLLSIM_3D_H

#include "MMLBase.h"

#include "CollisionSimulator3DBase.h"
#include "Containers3D.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace MML;

namespace MPL
{
	class SimResultsCollSim3D
	{
		int _numBalls;
	public:
		std::vector<Ball3D>	 Balls;													// properties of balls
		std::vector<std::vector<Pnt3Cart>> BallPosList;			// positions of all balls at each time step
		std::vector<std::vector<Vec3Cart>> BallVelList;			// positions of all balls at each time step

		SimResultsCollSim3D(int numBalls)
			: BallPosList(numBalls), BallVelList(numBalls), _numBalls(numBalls)
		{	}
		SimResultsCollSim3D(const std::vector<Ball3D>& balls)
			: Balls(balls), BallPosList(balls.size()), BallVelList(balls.size()), _numBalls(balls.size())
		{	}
		SimResultsCollSim3D(const std::vector<Ball3D>& balls, const std::vector<std::vector<Pnt3Cart>>& ballPositions)
			: Balls(balls), BallPosList(ballPositions), BallVelList(balls.size()), _numBalls(balls.size())
		{	}
		SimResultsCollSim3D(const std::vector<Ball3D>& balls, const std::vector<std::vector<Pnt3Cart>>& ballPositions,
												const std::vector<std::vector<Vec3Cart>>& ballVelocities)
			: Balls(balls), BallPosList(ballPositions), BallVelList(ballVelocities), _numBalls(balls.size())
		{	}

		int NumBalls() const { return _numBalls; }

		const std::vector<Pnt3Cart>& getBallPositions(int indBall) const
		{
			if (indBall < 0 || indBall >= BallPosList.size())
				throw std::out_of_range("Index out of range for ball positions");
			return BallPosList[indBall];
		}
		const std::vector<Vec3Cart>& getBallVelocities(int indBall) const
		{
			if (indBall < 0 || indBall >= BallVelList.size())
				throw std::out_of_range("Index out of range for ball velocities");
			return BallVelList[indBall];
		}

		std::vector<Pnt3Cart> getPathForBall(int ind) const
		{
			if (ind < 0 || ind >= BallPosList.size())
				throw std::out_of_range("Index out of range for ball positions");

			std::vector<Pnt3Cart> result;
			for (const auto& pos : BallPosList[ind]) {
				result.push_back(pos);
			}
			return result;
		}

		// statistics for balls
		// for balls speed, at time t, calc min, max, avg speed and speed deviation
		void CalcAllBallsStatistic(int timeStep, double& minSpeed, double& maxSpeed, double& avgSpeed, double& speedDev)
		{
			if (timeStep < 0 || timeStep >= BallVelList[0].size())
				throw std::out_of_range("Time step out of range for ball velocities");

			minSpeed = std::numeric_limits<double>::max();
			maxSpeed = std::numeric_limits<double>::lowest();
			double sumSpeed = 0.0, sumSpeedSq = 0.0;
			avgSpeed = 0.0;
			speedDev = 0.0;
			int count = 0;
			for (const auto& ballVelList : BallVelList)
			{
				if (timeStep < ballVelList.size())
				{
					double speed = ballVelList[timeStep].NormL2();
					minSpeed = std::min(minSpeed, speed);
					maxSpeed = std::max(maxSpeed, speed);
					sumSpeed += speed;
					sumSpeedSq += speed * speed;
					count++;
				}
			}
			if (count > 0)
			{
				avgSpeed = sumSpeed / count;
				speedDev = sqrt(sumSpeedSq / count - avgSpeed * avgSpeed);
			}
		}

		Pnt3Cart CentreOfMass(int timeStep)
		{
			if (timeStep < 0 || timeStep >= BallPosList[0].size())
				throw std::out_of_range("Time step out of range for CentreOfMass() call");

			Pnt3Cart com(0.0, 0.0, 0.0);
			double totalMass = 0.0;
			for (int i = 0; i < BallPosList.size(); i++)
			{
				const std::vector<Pnt3Cart>& ballPositions = BallPosList[i];
				const Ball3D& ball = Balls[i];
				if (ballPositions.size() > timeStep)
				{
					com += ball.Mass() * ballPositions[timeStep];
					totalMass += ball.Mass();
				}
			}
			if (totalMass > 0.0)
				com /= totalMass;
			return com;
		}

		Pnt3Cart AvgPos(int timeStep)
		{
			if (timeStep < 0 || timeStep >= BallPosList[0].size())
				throw std::out_of_range("Time step out of range for AvgPos() call");

			Pnt3Cart avgPos(0.0, 0.0, 0.0);
			int count = 0;
			for (const auto& ballPositions : BallPosList)
			{
				if (ballPositions.size() > timeStep)
				{
					avgPos += ballPositions[timeStep];
					count++;
				}
			}
			if (count > 0)
				avgPos /= count;
			return avgPos;
		}
		
		Real AvgSpeed(int timeStep)
		{
			if (timeStep < 0 || timeStep >= BallVelList[0].size())
				throw std::out_of_range("Time step out of range for AvgSpeed() call");
			double totalSpeed = 0.0;
			int count = 0;
			for (const auto& ballVelocities : BallVelList)
			{
				if (ballVelocities.size() > timeStep)
				{
					totalSpeed += ballVelocities[timeStep].NormL2();
					count++;
				}
			}
			return count > 0 ? totalSpeed / count : 0.0;
		}

		Real TotalKineticEnergy(double timeStep)
		{
			Real totalEnergy = 0.0;
			for (int i = 0; i < BallVelList.size(); i++)
			{
				const std::vector<Vec3Cart>& ballVelocities = BallVelList[i];
				const Ball3D& ball = Balls[i];

				if (ballVelocities.size() > timeStep)
				{
					totalEnergy += 0.5 * ball.Mass() * ballVelocities[timeStep].NormL2() * ballVelocities[timeStep].NormL2();
				}
			}
			return totalEnergy;
		}

		// Statistics for balls by color
		std::vector<std::pair<std::string, Pnt3Cart>> CentreOfMassByColor(int timeStep)
		{
			std::map<std::string, std::pair<Pnt3Cart, double>> mapComByColor; // color -> (centre of mass, total mass)
			for (int i = 0; i < BallPosList.size(); i++)
			{
				const std::vector<Pnt3Cart>& ballPositions = BallPosList[i];
				const Ball3D& ball = Balls[i];
				if (ballPositions.size() > timeStep)
				{
					mapComByColor[ball.Color()].first += ball.Mass() * ballPositions[timeStep];
					mapComByColor[ball.Color()].second += ball.Mass();
				}
			}
			std::vector<std::pair<std::string, Pnt3Cart>> result(mapComByColor.size());
			int index = 0;
			for (const auto& pair : mapComByColor)
			{
				if (pair.second.second > 0.0) // avoid division by zero
					result[index++] = { pair.first, pair.second.first / pair.second.second };
				else
					result[index++] = { pair.first, Pnt3Cart(0.0, 0.0, 0.0) };
			}
			return result;
		}

		std::vector<std::pair<std::string, Pnt3Cart>> AvgPosByColor(int timeStep)
		{
			std::map<std::string, std::pair<Pnt3Cart, int>> mapPosByColor; // color -> (sum of positions, count)
			for (int i = 0; i < BallPosList.size(); i++)
			{
				const std::vector<Pnt3Cart>& ballPositions = BallPosList[i];
				const Ball3D& ball = Balls[i];
				if (ballPositions.size() > timeStep)
				{
					mapPosByColor[ball.Color()].first += ballPositions[timeStep];
					mapPosByColor[ball.Color()].second++;
				}
			}
			std::vector<std::pair<std::string, Pnt3Cart>> result(mapPosByColor.size());
			int index = 0;
			for (const auto& pair : mapPosByColor)
			{
				if (pair.second.second > 0) // avoid division by zero
					result[index++] = { pair.first, pair.second.first / pair.second.second };
				else
					result[index++] = { pair.first, Pnt3Cart(0.0, 0.0, 0.0) };
			}
			return result;
		}

		std::vector<std::pair<std::string, Real>> AvgSpeedByColor(double timeStep)
		{
			std::map<std::string, std::pair<Real, int>> mapSpeedByColor; // color -> (sum of speeds, count)
			for (int i = 0; i < BallVelList.size(); i++)
			{
				const std::vector<Vec3Cart>& ballVelocities = BallVelList[i];
				const Ball3D& ball = Balls[i];
				if (ballVelocities.size() > timeStep)
				{
					double speed = ballVelocities[timeStep].NormL2();
					mapSpeedByColor[ball.Color()].first += speed;
					mapSpeedByColor[ball.Color()].second++;
				}
			}
			std::vector<std::pair<std::string, Real>> result(mapSpeedByColor.size());
			int index = 0;
			for (const auto& pair : mapSpeedByColor)
			{
				if (pair.second.second > 0) // avoid division by zero
					result[index++] = { pair.first, pair.second.first / pair.second.second };
				else
					result[index++] = { pair.first, 0.0 };
			}
			return result;
		}

		std::vector<std::pair<std::string, Real>> TotalKineticEnergyByColor(double timeStep)
		{
			std::map<std::string, Real> energyByColor; // color -> total kinetic energy
			for (int i = 0; i < BallVelList.size(); i++)
			{
				const std::vector<Vec3Cart>& ballVelocities = BallVelList[i];
				const Ball3D& ball = Balls[i];
				if (ballVelocities.size() > timeStep)
				{
					double speed = ballVelocities[timeStep].NormL2();
					energyByColor[ball.Color()] += 0.5 * ball.Mass() * speed * speed;
				}
			}
			std::vector<std::pair<std::string, Real>> result(energyByColor.size());
			int index = 0;
			for (const auto& pair : energyByColor)
			{
				result[index++] = pair;
			}
			return result;
		}
	};
}
#endif // MPL_SIMRESULTS_COLLSIM_2D_H

