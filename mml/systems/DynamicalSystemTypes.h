///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DynamicalSystemTypes.h                                              ///
///  Description: Type definitions for dynamical systems analysis                     ///
///               Fixed point types, result structures, section definitions           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DYNAMICAL_SYSTEM_TYPES_H
#define MML_DYNAMICAL_SYSTEM_TYPES_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

#include <vector>
#include <string>
#include <complex>

namespace MML::Systems 
{
	//=============================================================================
	// FIXED POINT TYPES
	//=============================================================================

	/// @brief Classification of fixed point stability (2D and higher)
	enum class FixedPointType {
		Unknown,				 ///< Could not classify
		StableNode,			 ///< All eigenvalues negative real
		UnstableNode,		 ///< All eigenvalues positive real
		Saddle,					 ///< Mixed signs (real eigenvalues)
		SaddleFocus,		 ///< Mixed signs with complex eigenvalues (e.g., Lorenz attractor)
		StableFocus,		 ///< Complex with negative real part (spiral in)
		UnstableFocus,	 ///< Complex with positive real part (spiral out)
		Center,					 ///< Purely imaginary eigenvalues
		StableStarNode,	 ///< Degenerate stable node
		UnstableStarNode ///< Degenerate unstable node
	};

	/// @brief Convert FixedPointType to string
	inline std::string ToString(FixedPointType type) {
		switch (type) {
		case FixedPointType::StableNode:
			return "Stable Node";
		case FixedPointType::UnstableNode:
			return "Unstable Node";
		case FixedPointType::Saddle:
			return "Saddle";
		case FixedPointType::SaddleFocus:
			return "Saddle-Focus";
		case FixedPointType::StableFocus:
			return "Stable Focus";
		case FixedPointType::UnstableFocus:
			return "Unstable Focus";
		case FixedPointType::Center:
			return "Center";
		case FixedPointType::StableStarNode:
			return "Stable Star Node";
		case FixedPointType::UnstableStarNode:
			return "Unstable Star Node";
		default:
			return "Unknown";
		}
	}

	//=============================================================================
	// RESULT STRUCTURES
	//=============================================================================

	/// @brief Result of fixed point analysis
	template<typename Type = Real>
	struct FixedPoint {
		Vector<Type> location;											 ///< Position in state space
		std::vector<std::complex<Type>> eigenvalues; ///< Jacobian eigenvalues
		Matrix<Type> jacobian;											 ///< Jacobian at fixed point
		FixedPointType type;												 ///< Stability classification
		bool isStable;															 ///< Overall stability
		Type convergenceResidual;										 ///< ||f(x*)|| at convergence
		int iterations;															 ///< Newton iterations used
	};

	/// @brief Result of Lyapunov exponent computation
	template<typename Type = Real>
	struct LyapunovResult {
		Vector<Type> exponents;			///< Lyapunov exponents (λ₁ ≥ λ₂ ≥ ... ≥ λₙ)
		Type maxExponent;						///< Largest exponent (λ₁)
		Type sum;										///< Sum of all exponents
		Type kaplanYorkeDimension;	///< Kaplan-Yorke dimension
		bool isChaotic;							///< True if λ₁ > 0
		int numOrthonormalizations; ///< Number of QR steps performed
		Type totalTime;							///< Total integration time
	};

	/// @brief Result of bifurcation sweep
	template<typename Type = Real>
	struct BifurcationDiagram {
		std::string parameterName;											///< Name of swept parameter
		std::vector<Type> parameterValues;							///< Parameter values
		std::vector<std::vector<Type>> attractorValues; ///< Attractor points at each param
		int numTransientSteps;													///< Transient removal steps
		int numRecordedPoints;													///< Points recorded per parameter
	};

	/// @brief Poincaré section definition
	template<typename Type = Real>
	struct PoincareSection {
		int variable;	 ///< Which state variable defines section
		Type value;		 ///< Section at x[variable] = value
		int direction; ///< +1 = positive crossing, -1 = negative, 0 = both

		PoincareSection(int var = 0, Type val = 0, int dir = 0)
				: variable(var)
				, value(val)
				, direction(dir) {}
	};

} // namespace MML::Systems
#endif // MML_DYNAMICAL_SYSTEM_TYPES_H
