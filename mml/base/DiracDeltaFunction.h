///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DiracDeltaFunction.h                                                ///
///  Description: Dirac delta function approximations for discrete sampling           ///
///               Gaussian and sinc approximations for numerical integration          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DELTA_FUNCTION_H
#define MML_DELTA_FUNCTION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

namespace MML
{
	/// @brief Base class for Dirac delta function approximations.
	/// @details Provides common interface for various approximations of the Dirac delta
	///          distribution δ(x), which satisfies ∫δ(x)dx = 1 and δ(x) = 0 for x ≠ 0.
	///          As N → ∞, all approximations approach the true delta function.
	class DiracFunction : public IRealFunction
	{
	protected:
		int _N;   ///< Approximation parameter (larger = sharper peak)
	public:
		/// @brief Constructs Dirac delta approximation.
		/// @param N Sharpness parameter (larger values give narrower peaks)
		DiracFunction(int N) : _N(N) {}
	};

	/// @brief Step function (rectangular) approximation of Dirac delta.
	/// @details Returns N for |x| < 1/(2N), otherwise 0. Rectangle has area 1.
	///          δ_N(x) = N for |x| < 1/(2N), 0 otherwise
	class DiracStep : public DiracFunction
	{
	public:
		/// @brief Constructs step function approximation.
		/// @param N Sharpness parameter (rectangle width = 1/N, height = N)
		DiracStep(int N) : DiracFunction(N) {}

		/// @brief Evaluates the step function approximation.
		/// @param x Point to evaluate
		/// @return N if |x| < 1/(2N), otherwise 0
		Real operator()(const Real x) const
		{
			if (x < -1.0 / (2 * _N) || x > 1.0 / (2 * _N))
				return 0.0;
			else
				return _N;
		}
	};

	/// @brief Gaussian approximation of Dirac delta.
	/// @details Uses normalized Gaussian: δ_N(x) = N/√(2π) · exp(-N²x²/2)
	///          This is a smooth approximation with exponential decay.
	///          Integrates to 1 for all N.
	class DiracExp : public DiracFunction
	{
	public:
		/// @brief Constructs Gaussian approximation.
		/// @param N Sharpness parameter (standard deviation = 1/N)
		DiracExp(int N) : DiracFunction(N) {}

		/// @brief Evaluates the Gaussian approximation.
		/// @param x Point to evaluate
		/// @return N/√(2π) · exp(-N²x²/2)
		Real operator()(const Real x) const { return _N / std::sqrt(2 * Constants::PI) * std::exp(-x * x * _N * _N / 2.0); }
	};

	/// @brief Lorentzian (Cauchy) approximation of Dirac delta.
	/// @details Uses Lorentzian/Cauchy distribution: δ_N(x) = N / (π(1 + N²x²))
	///          Has heavier tails than Gaussian but smoother than step.
	class DiracSqr : public DiracFunction
	{
	public:
		/// @brief Constructs Lorentzian approximation.
		/// @param N Sharpness parameter (FWHM ~ 2/N)
		DiracSqr(int N) : DiracFunction(N) {}

		/// @brief Evaluates the Lorentzian approximation.
		/// @param x Point to evaluate
		/// @return N / (π(1 + N²x²))
		Real operator()(const Real x) const { return _N / Constants::PI / (1 + _N * _N * x * x); }
	};

	/// @brief Sinc function approximation of Dirac delta.
	/// @details Uses sinc function: δ_N(x) = sin(Nx) / (πx)
	///          Oscillates and decays slowly. Common in Fourier analysis.
	///          At x = 0, returns the limit value N/π.
	class DiracSin : public DiracFunction
	{
	public:
		/// @brief Constructs sinc approximation.
		/// @param N Frequency parameter
		DiracSin(int N) : DiracFunction(N) {}

		/// @brief Evaluates the sinc approximation.
		/// @param x Point to evaluate
		/// @return sin(Nx) / (πx), or N/π when x = 0
		Real operator()(const Real x) const { 
			if (std::abs(x) < 1e-15)  // Handle x ≈ 0: lim(x→0) sin(Nx)/(πx) = N/π
				return static_cast<Real>(_N) / Constants::PI;
			return std::sin(_N * x) / (Constants::PI * x); 
		}
	};
}

#endif
