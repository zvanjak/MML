///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Random.h                                                            ///
///  Description: Random number generation utilities                                  ///
///               Uniform and Gaussian distributions using Mersenne Twister           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
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
		// Seeds the thread-local RNG for reproducible results (affects calling thread only)
		static void SetSeed(unsigned int seed) { gen.seed(seed); }

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
			// Generate a uniformly distributed random direction on the sphere with radius 'abs'
			// Correct method: sample u ~ Uniform(-1, 1) and set cos(phi) = u
			// This ensures uniform distribution over the sphere surface
			Real theta = UniformReal(0, 2 * Constants::PI); // azimuthal angle
			Real u = UniformReal(-1.0, 1.0);               // cos(polar angle)
			Real sin_phi = std::sqrt(1.0 - u * u);         // sin(phi) = sqrt(1 - cos²(phi))
			vx = abs * sin_phi * std::cos(theta);
			vy = abs * sin_phi * std::sin(theta);
			vz = abs * u;  // cos(phi) = u
			return abs;
		}
	};
}

#endif