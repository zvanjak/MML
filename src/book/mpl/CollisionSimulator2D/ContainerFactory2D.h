#if !defined MPL_CONTAINER_FACTORY_H
#define MPL_CONTAINER_FACTORY_H

#include "MMLBase.h"
#include "base/Random.h"

#include "CollisionSimulator2DBase.h"
#include "Containers2D.h"

using namespace MML;

namespace MPL
{
	class ContainerFactory2D
	{
	public:
		// container with five different balls, with given masses, radii, colors, positions and velocities
		static BoxContainer2D CreateConfig()
		{
			BoxContainer2D box1(1000, 800);

			box1.AddBall(Ball2D(1.0, 10.0, "red", Pnt2Cart(100, 150), Vec2Cart(3, 1.5)));
			box1.AddBall(Ball2D(2.0, 20.0, "blue", Pnt2Cart(200, 500), Vec2Cart(-2, -3)));
			box1.AddBall(Ball2D(1.0, 10.0, "green", Pnt2Cart(300, 150), Vec2Cart(1, -2)));
			box1.AddBall(Ball2D(3.0, 20.0, "yellow", Pnt2Cart(100, 400), Vec2Cart(-1, 1)));
			box1.AddBall(Ball2D(1.0, 10.0, "purple", Pnt2Cart(200, 600), Vec2Cart(1, -0.6)));

			return box1;
		}

		// container with N same balls, with given mass, radius, color, positions and velocities
		static BoxContainer2D CreateConfigNSameBalls(double width, double height, int N, double mass, double radius, 
																								 double velocityRange, const std::string& color)
		{
			BoxContainer2D box2(width, height);
			for (int i = 0; i < N; i++)
			{
				double x = Random::UniformReal(radius, width - radius);
				double y = Random::UniformReal(radius, height - radius);
				Real vx, vy;
				Random::UniformVecDirection2(vx, vy, velocityRange);
				box2.AddBall(Ball2D(mass, radius, color, Pnt2Cart(x, y), Vec2Cart(vx, vy)));
			}
			return box2;
		}

		// container with two types of balls, with masses, radii, positions and velocities
		// randomized in given intervals
		static BoxContainer2D CreateConfig1(double width, double height, int N,
			const std::string& color1, double minMass1, double maxMass1, double minRad1, double maxRad1, double velocityRange1,
			const std::string& color2, double minMass2, double maxMass2, double minRad2, double maxRad2, double velocityRange2 )
		{
			BoxContainer2D box2(width, height);

			for (int i = 0; i < N / 2; i++)
			{
				double mass = Random::UniformReal(minMass1, maxMass1);
				double radius = Random::UniformReal(minRad1, maxRad1);
				Point2Cartesian position(Random::UniformReal(radius, box2._width - radius), Random::UniformReal(radius, box2._height - radius));
				Vector2Cartesian velocity(Random::UniformReal(-velocityRange1, velocityRange1), Random::UniformReal(-velocityRange1, velocityRange1));

				box2.AddBall(Ball2D(mass, radius, "Red", position, velocity));
			}

			for (int i = 0; i < N / 2; i++)
			{
				double mass = Random::UniformReal(minMass2, maxMass2);
				double radius = Random::UniformReal(minRad2, maxRad2);
				Point2Cartesian position(Random::UniformReal(radius, box2._width - radius), Random::UniformReal(radius, box2._height - radius));
				Vector2Cartesian velocity(Random::UniformReal(-velocityRange2, velocityRange2), Random::UniformReal(-velocityRange2, velocityRange2));

				box2.AddBall(Ball2D(mass, radius, "Blue", position, velocity));
			}

			return box2;
		}

		// container with two types of balls, initially each occupies one half of the container
		// random balls position within each half, with given mass and velocity for each type
		static BoxContainer2D CreateConfig2(double width, double height,
			const std::string& color1, int nBall1, double mass1, double rad1, double vel1,
			const std::string& color2, int nBall2, double mass2, double rad2, double vel2)
		{
			BoxContainer2D box(width, height);

			for (int i = 0; i < nBall1; i++)
			{
				double x = Random::UniformReal(rad1, width / 2 - rad1);
				double y = Random::UniformReal(rad1, height - rad1);

				Real vx, vy;
				Random::UniformVecDirection2(vx, vy, vel1);

				box.AddBall(Ball2D(mass1, rad1, color1, Pnt2Cart(x, y), Vec2Cart(vx, vy)));
			}
			for (int i = 0; i < nBall2; i++)
			{
				double x = Random::UniformReal(width / 2 + rad2, width - rad2);
				double y = Random::UniformReal(rad2, height - rad2);

				Real vx, vy;
				Random::UniformVecDirection2(vx, vy, vel2);

				box.AddBall(Ball2D(mass2, rad2, color2, Pnt2Cart(x, y), Vec2Cart(vx, vy)));
			}

			return box;
		}

		// container with a central big ball and many small balls around it, with random positions and velocities
		static BoxContainer2D CreateBrownMotionConfig(double width, double height, int numSmallBalls)
		{
			BoxContainer2D box2(width, height);

			// adding central big ball
			box2.AddBall(Ball2D(1000, 20, "Red", Pnt2Cart(width / 2, height / 2), Vec2Cart(0, 0)));

			// adding numSmallBalls much smaller balls
			for (int i = 0; i < numSmallBalls; i++)
			{
				// create random position for a small ball
				double radius = 1.0;
				double mass = 1.0;
				double x = Random::UniformReal(radius, width - radius);
				double y = Random::UniformReal(radius, height - radius);

				// check if the position is not too close to the central ball
				while (box2.Ball(0).Pos().Dist(Pnt2Cart(x, y)) < box2.Ball(0).Rad() + radius)
				{
					x = Random::UniformReal(radius, width - radius);
					y = Random::UniformReal(radius, height - radius);
				}

				// create random velocity
				double speed = Random::UniformReal(10.0, 20.0);
				Real vx, vy;
				Random::UniformVecDirection2(vx, vy, speed);
				box2.AddBall(Ball2D(mass, radius, "Blue", Pnt2Cart(x, y), Vec2Cart(vx, vy)));
			}
			return box2;
		}

		// create a container with two types of balls, 'numEnergetic' of them are at the center and very energetic
		// the other 'numBalls' are scattered around the container with low energy
		static BoxContainer2D CreateConfig3(double width, double height, 
																				int numBalls, double mass1, double rad1, double vel1,
																				int numEnergetic, double posRadius, double mass2, double rad2, double vel2 )
		{
			BoxContainer2D box(width, height);

			for (int i = 0; i < numBalls; i++)
			{
				Real vx, vy;
				Random::UniformVecDirection2(vx, vy, vel1);
				
				Pnt2Cart position(Random::UniformReal(rad1, width - rad1), Random::UniformReal(rad1, height - rad1));
				// if ball is too close to the center, generate new position
				while (Pnt2Cart(width / 2, height / 2).Dist(position) < posRadius + rad1)
				{
					position = Pnt2Cart(Random::UniformReal(rad1, width - rad1), Random::UniformReal(rad1, height - rad1));
				}

				box.AddBall(Ball2D(mass1, rad1, "Blue", position, Vec2Cart(vx, vy)));
			}

			for (int i = 0; i < numEnergetic; i++)
			{
				Pnt2Cart position(width / 2 + Random::UniformReal(-posRadius, posRadius), height / 2 + Random::UniformReal(-posRadius, posRadius));

				Real vx, vy;
				Random::UniformVecDirection2(vx, vy, vel2);

				box.AddBall(Ball2D(mass2, rad2, "Red", position, Vec2Cart(vx, vy)));
			}

			return box;
		}

		// container with equal Balls, N in line and evenly ditributed in the container with given width and height
		// total number of Balls is N*N, with given mass, radius and speed (given value, but random direction)
		static BoxContainer2D CreateConfig4(double width, double height, int N, double mass, double radius, double speed)
		{
			BoxContainer2D box(width, height);
			double stepX = width / N;
			double stepY = height / N;
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					double x = i * stepX + stepX / 2;
					double y = j * stepY + stepY / 2;
					Vec2Cart velocity(Random::UniformReal(-speed, speed), Random::UniformReal(-speed, speed));
					box.AddBall(Ball2D(mass, radius, "Green", Pnt2Cart(x, y), velocity));
				}
			}
			return box;
		}
	};
} // namespace MPL
#endif // MPL_CONTAINER_FACTORY_H
