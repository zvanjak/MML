///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/TrialSpace1D.h                                       ///
///  Description: One-dimensional finite trial-space interface                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_TRIAL_SPACE_1D_H
#define MML_TRIAL_SPACE_1D_H

#include <mml/core/FunctionSpaces/FunctionSpace1D.h>

#include <string>

namespace MML::FunctionSpaces
{
	class TrialSpace1D
	{
	protected:
		void validateBasisIndex(int basisIndex, const char* context) const
		{
			if (basisIndex < 0 || basisIndex >= dimension())
				throw FunctionSpaceInputError(std::string(context) + ": basis index out of range");
		}

		static void validateDerivativeOrder(int order, const char* context)
		{
			if (order < 0)
				throw FunctionSpaceInputError(std::string(context) + ": derivative order must be non-negative");
		}

		void validateNodeIndex(int nodeIndex, const char* context) const
		{
			if (nodeIndex < 0 || nodeIndex >= nodeCount())
				throw FunctionSpaceInputError(std::string(context) + ": node index out of range");
		}

	public:
		virtual ~TrialSpace1D() = default;

		virtual const FunctionSpace1D& functionSpace() const noexcept = 0;
		virtual int dimension() const noexcept = 0;
		virtual Real basisValue(int basisIndex, Real x) const = 0;

		virtual bool hasBasisDerivative(int order = 1) const noexcept
		{
			return order == 0;
		}

		virtual Real basisDerivative(int basisIndex, int order, Real x) const
		{
			validateDerivativeOrder(order, "TrialSpace1D::basisDerivative");
			if (order == 0)
				return basisValue(basisIndex, x);
			throw FunctionSpaceError("TrialSpace1D::basisDerivative: derivatives are not available for this trial space");
		}

		virtual bool hasNodes() const noexcept { return false; }
		virtual int nodeCount() const noexcept { return 0; }

		virtual Real node(int nodeIndex) const
		{
			validateNodeIndex(nodeIndex, "TrialSpace1D::node");
			throw FunctionSpaceError("TrialSpace1D::node: nodes are not available for this trial space");
		}
	};
}

#endif // MML_TRIAL_SPACE_1D_H
