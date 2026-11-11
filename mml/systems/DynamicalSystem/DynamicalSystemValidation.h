///////////////////////////////////////////////////////////////////////////////////////////
///  MinimalMathLibrary - shared dynamical-system request and callback validation
///  Copyright (c) 2024-2026 Zvonimir Vanjak. MIT License (see LICENSE.md).
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DYNAMICAL_SYSTEM_VALIDATION_H
#define MML_DYNAMICAL_SYSTEM_VALIDATION_H

#include <mml/base/Matrix/Matrix.h>
#include <mml/interfaces/IDynamicalSystem.h>
#include <cmath>
#include <exception>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace MML::Systems
{
	class ParameterRestorationError : public std::runtime_error, public MMLException
	{
		std::exception_ptr _operationException;
		std::exception_ptr _restorationException;

		static std::string Describe(const std::exception_ptr& error)
		{
			if (!error) return "none";
			try { std::rethrow_exception(error); }
			catch (const std::exception& exception) { return exception.what(); }
			catch (...) { return "non-standard exception (retained for inspection)"; }
		}

	public:
		ParameterRestorationError(std::exception_ptr operation, std::exception_ptr restoration)
			: std::runtime_error("Parameter restoration failed; model state is unknown. Operation error: " +
				Describe(operation) + ". Restoration error: " + Describe(restoration)),
			  _operationException(operation), _restorationException(restoration) {}

		std::exception_ptr OperationException() const noexcept { return _operationException; }
		std::exception_ptr RestorationException() const noexcept { return _restorationException; }
	};
}

namespace MML::Systems::DynSysDetail
{
	inline void CheckVariationalDimension(int dimension)
	{
		if (dimension <= 0 || dimension > (std::numeric_limits<int>::max() - dimension) / dimension)
			throw ArgumentError("Lyapunov: positive dimension and representable augmented dimension required");
	}

	inline void CheckParameterScan(const IDynamicalSystem& system, int index,
		Real minimum, Real maximum, int count)
	{
		if (index < 0 || index >= system.getNumParam())
			throw ArgumentError("Parameter scan index out of range");
		if (count < 1)
			throw ArgumentError("Parameter scan requires at least one sample");
		if (!std::isfinite(minimum) || !std::isfinite(maximum))
			throw ArgumentError("Parameter scan endpoints must be finite");
	}

	template<typename ComputeFn>
	auto WithRestoredParameter(IDynamicalSystem& system, int index, ComputeFn&& compute)
	{
		const Real original = system.getParam(index);
		if (!std::isfinite(original))
			throw ArgumentError("Original scan parameter must be finite");
		std::optional<std::invoke_result_t<ComputeFn&>> result;
		std::exception_ptr operationError;
		try { result.emplace(compute()); }
		catch (...) { operationError = std::current_exception(); }

		// Restore explicitly: a throwing setter must never escape a noexcept destructor.
		try { system.setParam(index, original); }
		catch (...) { throw ParameterRestorationError(operationError, std::current_exception()); }
		if (operationError) std::rethrow_exception(operationError);
		return std::move(*result);
	}

	template<typename NonfiniteError = NumericalMethodError>
	inline void CheckState(const Vector<Real>& state, int dimension, const char* context, bool checkFinite = true)
	{
		if (state.size() != dimension)
			throw VectorDimensionError(context, state.size(), dimension);
		if (checkFinite)
			for (Real value : state)
				if (!std::isfinite(value))
					throw NonfiniteError(std::string(context) + ": non-finite state");
	}

	inline void CheckJacobian(const Matrix<Real>& matrix, int dimension)
	{
		if (matrix.rows() != dimension || matrix.cols() != dimension)
			throw MatrixDimensionError("Lyapunov: Jacobian must match the system dimension",
				matrix.rows(), matrix.cols(), dimension, dimension);
		for (int row = 0; row < dimension; ++row)
			for (int column = 0; column < dimension; ++column)
				if (!std::isfinite(matrix(row, column)))
					throw NumericalMethodError("Lyapunov: non-finite Jacobian");
	}

	inline void CheckInitialState(const IODESystem& system, const Vector<Real>& state)
	{
		if (system.getDim() <= 0)
			throw ArgumentError("Dynamical analysis requires a positive system dimension");
		CheckState<ArgumentError>(state, system.getDim(), "Dynamical analysis initial state");
	}

	inline void RequirePositive(Real value, const char* name)
	{
		if (!std::isfinite(value) || value <= 0)
			throw ArgumentError(std::string(name) + " must be finite and positive");
	}

	inline void RequireFinite(Real value, const char* name)
	{
		if (!std::isfinite(value))
			throw ArgumentError(std::string(name) + " must be finite");
	}

	inline void RequireNonnegative(Real value, const char* name)
	{
		if (!std::isfinite(value) || value < 0)
			throw ArgumentError(std::string(name) + " must be finite and nonnegative");
	}

	inline int NominalIntervalCount(Real duration, Real initialStep)
	{
		if (initialStep > 0 && initialStep < 1 && duration > std::numeric_limits<Real>::max() * initialStep)
			throw ResourceLimitError("Bifurcation nominal interval estimate exceeds the integer range");
		const Real count = std::round(duration / initialStep);
		if (!std::isfinite(count) || static_cast<long double>(count) > std::numeric_limits<int>::max())
			throw ResourceLimitError("Bifurcation nominal interval estimate exceeds the integer range");
		return static_cast<int>(count);
	}
}
#endif // MML_DYNAMICAL_SYSTEM_VALIDATION_H
