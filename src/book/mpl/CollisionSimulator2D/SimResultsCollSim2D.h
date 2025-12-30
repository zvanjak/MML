#if !defined MPL_SIMRESULTS_COLLSIM_2D_H
#define MPL_SIMRESULTS_COLLSIM_2D_H

#include "MMLBase.h"

#include "CollisionSimulator2DBase.h"
#include "Containers2D.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace MML;

namespace MPL
{
	class SimResultsCollSim2D
	{
		int _numBalls;
	public:
		std::vector<Ball2D>	 Balls;													// properties of balls
		std::vector<std::vector<Pnt2Cart>> BallPosList;			// positions of all balls at each time step
		std::vector<std::vector<Vec2Cart>> BallVelList;			// positions of all balls at each time step

		SimResultsCollSim2D(int numBalls)
			: BallPosList(numBalls), BallVelList(numBalls), _numBalls(numBalls)
		{
		}
		SimResultsCollSim2D(const std::vector<Ball2D>& balls)
			: Balls(balls), BallPosList(balls.size()), BallVelList(balls.size()), _numBalls(balls.size())
		{
		}
		SimResultsCollSim2D(const std::vector<Ball2D>& balls, const std::vector<std::vector<Pnt2Cart>>& ballPositions)
			: Balls(balls), BallPosList(ballPositions), BallVelList(balls.size()), _numBalls(balls.size())
		{
		}

		SimResultsCollSim2D(const std::vector<Ball2D>& balls, const std::vector<std::vector<Pnt2Cart>>& ballPositions,
												const std::vector<std::vector<Vec2Cart>>& ballVelocities)
			: Balls(balls), BallPosList(ballPositions), BallVelList(ballVelocities), _numBalls(balls.size())
		{
		}

		int NumBalls() const { return _numBalls; }

		const std::vector<Pnt2Cart>& getBallPositions(int indBall) const
		{
			if (indBall < 0 || indBall >= BallPosList.size())
				throw std::out_of_range("Index out of range for ball positions");
			return BallPosList[indBall];
		}
		const std::vector<Vec2Cart>& getBallVelocities(int indBall) const
		{
			if (indBall < 0 || indBall >= BallVelList.size())
				throw std::out_of_range("Index out of range for ball velocities");
			return BallVelList[indBall];
		}

		std::vector<Pnt2Cart> getPathForBall(int ind) const
		{
			if (ind < 0 || ind >= BallPosList.size())
				throw std::out_of_range("Index out of range for ball positions");

			std::vector<Pnt2Cart> result;
			for (const auto& pos : BallPosList[ind]) {
				result.push_back(pos);
			}
			return result;
		}

		// Statistics for all balls
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

		Pnt2Cart CentreOfMass(int timeStep)
		{
			if (timeStep < 0 || timeStep >= BallPosList[0].size())
				throw std::out_of_range("Time step out of range for ball positions");

			Pnt2Cart com(0.0, 0.0);
			double totalMass = 0.0;
			for (int i = 0; i < BallPosList.size(); i++)
			{
				const std::vector<Pnt2Cart>& ballPositions = BallPosList[i];
				const Ball2D& ball = Balls[i];
				if (ballPositions.size() > timeStep)
				{
					com += ball.Mass() * ballPositions[timeStep];
					totalMass += ball.Mass();
				}
			}
			return totalMass > 0.0 ? com / totalMass : Pnt2Cart(0.0, 0.0); // avoid division by zero
		}

		Pnt2Cart AvgPos(int timeStep)
		{
			if (timeStep < 0 || timeStep >= BallPosList[0].size())
				throw std::out_of_range("Time step out of range for AvgPos() call");

			Pnt2Cart avgPos(0.0, 0.0);
			int count = 0;
			for (int i = 0; i < BallPosList.size(); i++)
			{
				const std::vector<Pnt2Cart>& ballPositions = BallPosList[i];
				if (ballPositions.size() > timeStep)
				{
					avgPos += ballPositions[timeStep];
					count++;
				}
			}
			return count > 0 ? avgPos / count : Pnt2Cart(0.0, 0.0);
		}

		Real AvgSpeed(int timeStep)
		{
			if (timeStep < 0 || timeStep >= BallVelList[0].size())
				throw std::out_of_range("Time step out of range for AvgSpeed() call");

			double totalSpeed = 0.0;
			int count = 0;
			for (int i = 0; i < BallVelList.size(); i++)
			{
				const std::vector<Vec2Cart>& ballVelocities = BallVelList[i];
				if (ballVelocities.size() > timeStep)
				{
					totalSpeed += ballVelocities[timeStep].NormL2();
					count++;
				}
			}
			return count > 0 ? totalSpeed / count : 0.0;
		}

		Real TotalKineticEnergy(int timeStep)
		{
			if (timeStep < 0 || timeStep >= BallVelList[0].size())
				throw std::out_of_range("Time step out of range for TotalKineticEnergy() call");

			Real totalEnergy = 0.0;
			for (int i = 0; i < BallVelList.size(); i++)
			{
				const std::vector<Vec2Cart>& ballVelocities = BallVelList[i];
				const Ball2D& ball = Balls[i];

				if (ballVelocities.size() > timeStep)
				{
					totalEnergy += 0.5 * ball.Mass() * ballVelocities[timeStep].NormL2() * ballVelocities[timeStep].NormL2();
				}
			}
			return totalEnergy;
		}

		Real RootMeanSquareSpeed(int timeStep)
		{
			if (timeStep < 0 || timeStep >= BallVelList[0].size())
				throw std::out_of_range("Time step out of range for RootMeanSquareSpeed() call");
			double sumSpeedSq = 0.0;
			int count = 0;
			for (int i = 0; i < BallVelList.size(); i++)
			{
				const std::vector<Vec2Cart>& ballVelocities = BallVelList[i];
				if (ballVelocities.size() > timeStep)
				{
					double speed = ballVelocities[timeStep].NormL2();
					sumSpeedSq += speed * speed;
					count++;
				}
			}
			return count > 0 ? sqrt(sumSpeedSq / count) : 0.0;
		}

		// Statistics by color
		// for given time step, returns vector of pairs (color, center of mass)
		std::vector<std::pair<std::string, Pnt2Cart>> CentreOfMassByColor(int timeStep)
		{
			std::map<std::string, std::pair<Pnt2Cart, int>> mapComByColor; // color -> (sum of positions, count)
			for (int i = 0; i < BallPosList.size(); i++)
			{
				const std::vector<Pnt2Cart>& ballPositions = BallPosList[i];
				const Ball2D& ball = Balls[i];
				if (ballPositions.size() > timeStep)
				{
					const Pnt2Cart& pos = ballPositions[timeStep];
					mapComByColor[ball.Color()].first += ball.Mass() * pos;
					mapComByColor[ball.Color()].second += ball.Mass(); // accumulate mass for center of mass calculation
				}
			}
			std::vector<std::pair<std::string, Pnt2Cart>> result(mapComByColor.size());
			int index = 0;
			for (const auto& pair : mapComByColor)
			{
				if (pair.second.second > 0) // avoid division by zero
					result[index++] = { pair.first, pair.second.first / static_cast<Real>(pair.second.second) };
				else
					result[index++] = { pair.first, Pnt2Cart(0.0, 0.0) }; // or some other default value
			}
			return result;
		}

		std::vector<std::pair<std::string, Pnt2Cart>> AvgPosByColor(int timeStep)
		{
			std::map<std::string, std::pair<Pnt2Cart, int>> mapPosByColor; // color -> (sum of positions, count)
			for (int i = 0; i < BallPosList.size(); i++)
			{
				const std::vector<Pnt2Cart>& ballPositions = BallPosList[i];
				const Ball2D& ball = Balls[i];
				if (ballPositions.size() > timeStep)
				{
					const Pnt2Cart& pos = ballPositions[timeStep];
					mapPosByColor[ball.Color()].first += pos;
					mapPosByColor[ball.Color()].second++;
				}
			}
			std::vector<std::pair<std::string, Pnt2Cart>> result(mapPosByColor.size());
			int index = 0;
			for (const auto& pair : mapPosByColor)
			{
				if (pair.second.second > 0) 
					result[index++] = { pair.first, pair.second.first / static_cast<Real>(pair.second.second) };
				else
					result[index++] = { pair.first, Pnt2Cart(0.0, 0.0) }; // or some other default value
			}
			return result;
		}

		std::vector<std::pair<std::string, Real>> AvgSpeedByColor(int timeStep)
		{
			std::map<std::string, std::pair<Real, int>> mapSpeedByColor; // color -> (sum of speeds, count)
			for (int i = 0; i < BallVelList.size(); i++)
			{
				const std::vector<Vec2Cart>& ballVelocities = BallVelList[i];
				const Ball2D& ball = Balls[i];
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
				const std::vector<Vec2Cart>& ballVelocities = BallVelList[i];
				const Ball2D& ball = Balls[i];
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

