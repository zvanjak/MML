///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationBase.h                                                    ///
///  Description: Base definitions for numerical differentiation                      ///
///               Step size control and derivative types                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_BASE_H
#define MML_DERIVATION_BASE_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

namespace MML
{
	namespace Derivation
	{
		static inline const Real NDer1_h = 2 * std::sqrt(Constants::Eps);
		static inline const Real NDer2_h = std::pow(3 * Constants::Eps, 1.0 / 3.0);
		static inline const Real NDer4_h = std::pow(11.25 * Constants::Eps, 1.0 / 5.0);     // 0.0012009323661373839 for double!
		static inline const Real NDer6_h = std::pow(Constants::Eps / 168.0, 1.0 / 7.0);
		static inline const Real NDer8_h = std::pow(551.25 * Constants::Eps, 1.0 / 9.0);
	}
}

#endif // MML_DERIVATION_BASE_H
