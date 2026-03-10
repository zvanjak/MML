///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationBase.h                                                    ///
///  Description: Base definitions for numerical differentiation                      ///
///               Step size control and derivative types                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_BASE_H
#define MML_DERIVATION_BASE_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

namespace MML
{
	namespace Derivation
	{
		/// @brief Optimal step size for 1st order derivative (forward/backward difference)
		/// @note h = 2√ε minimizes roundoff + truncation error for O(h) methods
		static inline const Real NDer1_h = 2 * std::sqrt(Constants::Eps);
		/// @brief Optimal step size for 2nd order derivative (central difference)
		/// @note h = (3ε)^(1/3) for O(h²) methods, balances O(ε/h) roundoff with O(h²) truncation
		static inline const Real NDer2_h = std::pow(3 * Constants::Eps, 1.0 / 3.0);
		/// @brief Optimal step size for 4th order derivative (5-point stencil)
		/// @note h = (11.25ε)^(1/5) ≈ 0.0012 for double precision
		static inline const Real NDer4_h = std::pow(11.25 * Constants::Eps, 1.0 / 5.0);     // 0.0012009323661373839 for double!
		/// @brief Optimal step size for 6th order derivative (7-point stencil)
		/// @note h = (ε/168)^(1/7) for O(h⁶) methods
		static inline const Real NDer6_h = std::pow(Constants::Eps / 168.0, 1.0 / 7.0);
		/// @brief Optimal step size for 8th order derivative (9-point stencil)
		/// @note h = (551.25ε)^(1/9) for O(h⁸) methods
		static inline const Real NDer8_h = std::pow(551.25 * Constants::Eps, 1.0 / 9.0);
	}
}

#endif // MML_DERIVATION_BASE_H
