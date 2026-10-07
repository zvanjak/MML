///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/LinearOperator.h                                     ///
///  Description: Matrix-free linear operators for discretized function spaces        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_SPACES_LINEAR_OPERATOR_H
#define MML_FUNCTION_SPACES_LINEAR_OPERATOR_H

#include <mml/core/FunctionSpaces/FunctionSpace1D.h>

#include <mml/base/Vector/Vector.h>

namespace MML::FunctionSpaces
{
	class LinearOperator
	{
	public:
		virtual ~LinearOperator() = default;

		virtual int rows() const noexcept = 0;
		virtual int cols() const noexcept = 0;
		virtual void apply(const Vector<Real>& input, Vector<Real>& output) const = 0;

		Vector<Real> apply(const Vector<Real>& input) const
		{
			Vector<Real> output(rows());
			apply(input, output);
			return output;
		}

	protected:
		void validateInputOutput(const Vector<Real>& input, const Vector<Real>& output, const char* context) const
		{
			if (input.size() != cols())
				throw VectorDimensionError(std::string(context) + ": input dimension must match operator columns", input.size(), cols());
			if (output.size() != rows())
				throw VectorDimensionError(std::string(context) + ": output dimension must match operator rows", output.size(), rows());
		}
	};
}

#endif // MML_FUNCTION_SPACES_LINEAR_OPERATOR_H
