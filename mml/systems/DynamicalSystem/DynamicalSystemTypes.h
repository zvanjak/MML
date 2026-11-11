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

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

#include <vector>
#include <string>
#include <complex>
#include <limits>

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
		Center,					 ///< Purely imaginary nonzero spectrum; linearized label only
		StableStarNode,	 ///< Degenerate stable node
		UnstableStarNode ///< Degenerate unstable node
	};

	/// @brief Convert FixedPointType to string
	inline std::string to_string(FixedPointType type) {
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
	/// @deprecated Use to_string() (CONTRACT_POLICY section 11: STL spelling).
	[[deprecated("use to_string()")]]
	inline std::string ToString(FixedPointType type) { return to_string(type); }

	//=============================================================================
	// RESULT STRUCTURES
	//=============================================================================

	/// @brief Result of fixed point analysis
	template<typename Type = Real>
	struct FixedPoint {
		Vector<Type> location;											 ///< Position in state space
		std::vector<std::complex<Type>> eigenvalues; ///< Jacobian eigenvalues
		Matrix<Type> jacobian;											 ///< Jacobian at fixed point
		FixedPointType type = FixedPointType::Unknown;												 ///< Stability classification
		bool isStable = false;															 ///< Strictly decaying linear spectrum; false does not prove instability
		Type convergenceResidual = std::numeric_limits<Type>::quiet_NaN(); ///< ||f(location)||, unavailable until evaluated
		int iterations = 0;															 ///< Newton iterations used
		AlgorithmStatus status = AlgorithmStatus::AlgorithmSpecificFailure; ///< Finding/classification outcome
		std::string error_message = "Fixed-point analysis has not run";
		int function_evaluations = 0; ///< Direct Newton derivs calls; excludes calls internal to user Jacobians
	};

	/// @brief Result of Lyapunov exponent computation
	template<typename Type = Real>
	struct LyapunovResult {
		Vector<Type> exponents;			///< Lyapunov exponents (λ₁ ≥ λ₂ ≥ ... ≥ λₙ)
		Type maxExponent = Type(0);						///< Largest exponent (λ₁)
		Type sum = Type(0);										///< Sum of all exponents
		Type kaplanYorkeDimension = Type(0);	///< Kaplan-Yorke dimension
		bool isChaotic = false;							///< Positive-growth heuristic using Precision::LyapunovChaosThreshold
		int numOrthonormalizations = 0; ///< Number of QR steps performed
		Type totalTime = Type(0);							///< Total integration time
		int function_evaluations = 0; ///< Direct trajectory RHS calls, including rejected QR trials; excludes calls inside Jacobians
		int qr_refinements = 0; ///< Internal block halvings needed to resolve the tangent frame
	};

	/// @brief Result of bifurcation sweep
	template<typename Type = Real>
	struct BifurcationDiagram {
		std::string parameterName;											///< Name of swept parameter
		std::vector<Type> parameterValues;							///< Parameter values
		std::vector<std::vector<Type>> attractorValues; ///< Attractor points at each param
		int numTransientSteps = 0;													///< Legacy rounded transientTime/initialStep estimate, not adaptive work
		int numRecordedPoints = 0;													///< Legacy rounded recordTime/initialStep estimate, not actual output count
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
