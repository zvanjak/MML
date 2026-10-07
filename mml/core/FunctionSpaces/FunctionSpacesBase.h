///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/FunctionSpacesBase.h                                 ///
///  Description: Foundational declarations for finite function-space methods         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_SPACES_BASE_H
#define MML_FUNCTION_SPACES_BASE_H

#include <mml/MMLBase.h>

namespace MML::FunctionSpaces
{
	enum class DiscretizationMethod
	{
		SpectralGalerkin,
		Collocation,
		FiniteDifference,
		FiniteElement,
		MatrixFree
	};
}

#endif // MML_FUNCTION_SPACES_BASE_H
