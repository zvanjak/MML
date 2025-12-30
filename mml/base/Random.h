///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Random.h                                                            ///
///  Description: Random number generation utilities                                  ///
///               Uniform and Gaussian distributions using Mersenne Twister           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_RANDOM_H
#define MML_RANDOM_H

#include "MMLBase.h"

#include <random>

namespace MML
{
	class Random
	{
	private:
		// Thread-local static RNG for performance and thread-safety
		// Creating RNG once per thread instead of on every call provides 100-1000x speedup
		inline static thread_local std::mt19937 gen{std::random_device{}()};

	public:
		static Real UniformReal(Real min, Real max)
		{
			std::uniform_real_distribution<Real> dis(min, max);
			return dis(gen);
		}
		static int UniformInt(int min, int max)
		{
			std::uniform_int_distribution<int> dis(min, max);
			return dis(gen);
		}
		static Real UniformVecDirection2(Real& vx, Real& vy, Real abs)
		{
			// Generate a random angle in the range [0, 2*PI)
			Real angle = UniformReal(0, 2 * Constants::PI);
			vx = abs * std::cos(angle);
			vy = abs * std::sin(angle);
			return abs;
		}
		static Real UniformVecDirection3(Real& vx, Real& vy, Real& vz, Real abs)
		{
			// Generate a random direction on the sphere with radus 'abs'
			Real theta = UniformReal(0, 2 * Constants::PI); // azimuthal angle
			Real phi = UniformReal(0, Constants::PI); // polar angle
			vx = abs * std::sin(phi) * std::cos(theta);
			vy = abs * std::sin(phi) * std::sin(theta);
			vz = abs * std::cos(phi);
			return abs;
		}
	};
}

#endif