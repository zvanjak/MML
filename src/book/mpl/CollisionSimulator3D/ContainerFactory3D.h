#if !defined MPL_CONTAINER_FACTORY_3D_H
#define MPL_CONTAINER_FACTORY_3D_H

#include "MMLBase.h"
#include "base/Random.h"

#include "CollisionSimulator3DBase.h"
#include "Containers3D.h"

using namespace MML;

namespace MPL
{
	class ContainerFactory3D
	{
	public:
		static BoxContainer3D CreateConfig()
		{
			BoxContainer3D box1(1000, 1000, 1000);

			box1.AddBall(Ball3D(1.0, 10.0, "red", Pnt3Cart(100, 150, 50), Vec3Cart(3, 1.5, 4)));
			box1.AddBall(Ball3D(2.0, 20.0, "blue", Pnt3Cart(200, 500, 200), Vec3Cart(-2, -3, -6)));
			box1.AddBall(Ball3D(1.0, 10.0, "green", Pnt3Cart(300, 150, 300), Vec3Cart(1, -2, 3)));
			box1.AddBall(Ball3D(3.0, 20.0, "yellow", Pnt3Cart(100, 400, 200), Vec3Cart(-1, 1, 3)));
			box1.AddBall(Ball3D(1.0, 10.0, "purple", Pnt3Cart(200, 600, 500), Vec3Cart(1, -0.6, -2)));

			return box1;
		}

		// container with four types of balls, representing gases in air - nitrogen, oxygen, argon and carbon dioxide
		static BoxContainer3D CreateConfigGases(double width, double height, double depth, int totalMolecules,
																						double rad_N2, double rad_O2, double rad_Ar2, double rad_CO2,
																						double vel_N2, double vel_O2, double vel_Ar2, double vel_CO2,
																						std::string color_N2, std::string color_O2, std::string color_Ar2, std::string color_CO2 )
		{
			// Assuming a rough distribution of gases in air, we can use the following proportions:
			// Nitrogen: 78%, Oxygen: 21%, Argon: 0.93%, Carbon Dioxide: 0.04%
			// Total molecules = 1000 (for simplicity, can be adjusted)
			int nitrogenCount = static_cast<int>(totalMolecules * 0.78);
			int oxygenCount = static_cast<int>(totalMolecules * 0.21);
			int argonCount = static_cast<int>(totalMolecules * 0.0093);
			int carbonDioxideCount = static_cast<int>(totalMolecules * 0.0004);

			BoxContainer3D box(width, height, depth);
			// Nitrogen
			for (int i = 0; i < nitrogenCount; i++)
			{
				Pnt3Cart position(Random::UniformReal(rad_N2, width - rad_N2), Random::UniformReal(rad_N2, height - rad_N2), 
													Random::UniformReal(rad_N2, depth - rad_N2));

				Real vx, vy, vz;
				Random::UniformVecDirection3(vx, vy, vz, vel_N2);
				Vec3Cart velocity(vx, vy, vz);

				box.AddBall(Ball3D(0.02802, rad_N2, color_N2, position, velocity));
			}
			// Oxygen
			for (int i = 0; i < oxygenCount; i++)
			{
				Pnt3Cart position(Random::UniformReal(rad_O2, width - rad_O2), Random::UniformReal(rad_O2, height - rad_O2), Random::UniformReal(rad_O2, depth - rad_O2));
				Real vx, vy, vz;
				Random::UniformVecDirection3(vx, vy, vz, vel_O2);
				Vec3Cart velocity(vx, vy, vz);
				box.AddBall(Ball3D(0.03200, rad_O2, color_O2, position, velocity));
			}
			// Argon
			for (int i = 0; i < argonCount; i++)
			{
				Pnt3Cart position(Random::UniformReal(rad_Ar2, width - rad_Ar2), Random::UniformReal(rad_Ar2, height - rad_Ar2), Random::UniformReal(rad_Ar2, depth - rad_Ar2));
				Real vx, vy, vz;
				Random::UniformVecDirection3(vx, vy, vz, vel_Ar2);
				Vec3Cart velocity(vx, vy, vz);
				box.AddBall(Ball3D(0.03995,rad_Ar2, color_Ar2, position, velocity));
			}
			// Carbon Dioxide
			for (int i = 0; i < carbonDioxideCount; i++)
			{
				Pnt3Cart position(Random::UniformReal(rad_CO2, width - rad_CO2), Random::UniformReal(rad_CO2, height - rad_CO2), Random::UniformReal(rad_CO2, depth - rad_CO2));
				Real vx, vy, vz;
				Random::UniformVecDirection3(vx, vy, vz, vel_CO2);
				Vec3Cart velocity(vx, vy, vz);
				box.AddBall(Ball3D(0.04401, rad_CO2, color_CO2, position, velocity));
			}

			return box;
		}

		// container with N same balls, with given mass, radius, color, positions and velocities
		static BoxContainer3D CreateConfigNSameBalls(double width, double height, double depth, int N,
																			 double mass, double radius, double velocityAbs, const std::string& color)
		{
			BoxContainer3D box(width, height, depth);
			for (int i = 0; i < N; i++)
			{
				Pnt3Cart position(Random::UniformReal(radius, width - radius), Random::UniformReal(radius, height - radius), Random::UniformReal(radius, depth - radius));

				Real vx, vy, vz;
				Random::UniformVecDirection3(vx, vy, vz, velocityAbs);
				Vec3Cart velocity(vx, vy, vz);

				box.AddBall(Ball3D(mass, radius, color, position, velocity));
			}
			return box;
		}

		static BoxContainer3D CreateConfig1(double width, double height, double depth, int N,
			const std::string& color1, double minMass1, double maxMass1, double minRad1, double maxRad1, double velocityRange1,
			const std::string& color2, double minMass2, double maxMass2, double minRad2, double maxRad2, double velocityRange2)
		{
			BoxContainer3D box(width, height, depth);

			for (int i = 0; i < N / 2; i++)
			{
				double mass = Random::UniformReal(minMass1, maxMass1);
				double radius = Random::UniformReal(minRad1, maxRad1);
				Pnt3Cart position(Random::UniformReal(radius, box._width - radius), Random::UniformReal(radius, box._height - radius), Random::UniformReal(radius, box._height - radius));
				Vec3Cart velocity(Random::UniformReal(-5.0, 5.0), Random::UniformReal(-5.0, 5.0), Random::UniformReal(-5.0, 5.0));

				box.AddBall(Ball3D(mass, radius, color1, position, velocity));
			}
			for (int i = 0; i < N / 2; i++)
			{
				double mass = Random::UniformReal(minMass2, maxMass2);
				double radius = Random::UniformReal(minRad2, maxRad2);
				Pnt3Cart position(Random::UniformReal(radius, box._width - radius), Random::UniformReal(radius, box._height - radius), Random::UniformReal(radius, box._height - radius));
				Vec3Cart velocity(Random::UniformReal(-5.0, 5.0), Random::UniformReal(-5.0, 5.0), Random::UniformReal(-5.0, 5.0));

				box.AddBall(Ball3D(mass, radius, color2, position, velocity));
			}

			return box;
		}


		// container with two types of balls, initially each occupies one half of the container
		// random balls position within each half, with given mass and velocity for each type
		static BoxContainer3D CreateConfig2(double width, double height, double depth,
														int nBall1, double mass1, double rad1, double vel1,
														int nBall2, double mass2, double rad2, double vel2)
		{
			BoxContainer3D box(width, height, depth);

			for (int i = 0; i < nBall1; i++)
			{
				double x = Random::UniformReal(rad1, width / 2 - rad1);
				double y = Random::UniformReal(rad1, height - rad1);
				double z = Random::UniformReal(rad1, depth - rad1);

				Real vx, vy, vz;
				Random::UniformVecDirection3(vx, vy, vz, vel1);

				box.AddBall(Ball3D(mass1, rad1, "Red", Pnt3Cart(x, y, z), Vec3Cart(vx, vy, vz)));
			}
			for (int i = 0; i < nBall2; i++)
			{
				double x = Random::UniformReal(width / 2 + rad2, width - rad2);
				double y = Random::UniformReal(rad2, height - rad2);
				double z = Random::UniformReal(rad1, depth - rad1);

				Real vx, vy, vz;
				Random::UniformVecDirection3(vx, vy, vz, vel1);

				box.AddBall(Ball3D(mass2, rad2, "Blue", Pnt3Cart(x, y, z), Vec3Cart(vx, vy, vz)));
			}

			return box;
		}
	};
} // namespace MPL
#endif // MPL_CONTAINER_FACTORY_H
