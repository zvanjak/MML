///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Constraints/BoundConstraints.h                         ///
///  Description: Bound constraints utility for box-constrained optimization          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BOUND_CONSTRAINTS_H
#define MML_BOUND_CONSTRAINTS_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace MML::Optimization
{
	class BoundConstraintsError : public std::runtime_error
	{
	public:
		explicit BoundConstraintsError(const std::string& message)
			: std::runtime_error("BoundConstraintsError: " + message) { }
	};

	class BoundConstraints
	{
		Vector<Real> _lower;
		Vector<Real> _upper;

		static bool IsValidBoundValue(Real value)
		{
			return !std::isnan(value);
		}

		void validateBounds() const
		{
			if (_lower.size() != _upper.size())
				throw BoundConstraintsError("lower and upper bounds must have the same dimension");
			if (_lower.size() == 0)
				throw BoundConstraintsError("bounds dimension must be positive");

			for (int index = 0; index < _lower.size(); ++index) {
				if (!IsValidBoundValue(_lower[index]) || !IsValidBoundValue(_upper[index]))
					throw BoundConstraintsError("bounds cannot be NaN");
				if (_lower[index] > _upper[index])
					throw BoundConstraintsError("lower bound cannot exceed upper bound");
			}
		}

		void validateDimension(int dimension, const char* context) const
		{
			if (dimension != dimension_count())
				throw VectorDimensionError(context, dimension_count(), dimension);
		}

	public:
		BoundConstraints() = default;

		BoundConstraints(const Vector<Real>& lower, const Vector<Real>& upper)
			: _lower(lower), _upper(upper)
		{
			validateBounds();
		}

		BoundConstraints(int dimension, Real lower, Real upper)
			: _lower(dimension), _upper(dimension)
		{
			for (int index = 0; index < dimension; ++index) {
				_lower[index] = lower;
				_upper[index] = upper;
			}
			validateBounds();
		}

		template<int N>
		BoundConstraints(const VectorN<Real, N>& lower, const VectorN<Real, N>& upper)
			: _lower(N), _upper(N)
		{
			for (int index = 0; index < N; ++index) {
				_lower[index] = lower[index];
				_upper[index] = upper[index];
			}
			validateBounds();
		}

		static BoundConstraints Unbounded(int dimension)
		{
			return BoundConstraints(dimension, -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity());
		}

		int dimension_count() const noexcept { return _lower.size(); }
		int size() const noexcept { return dimension_count(); }

		const Vector<Real>& lower() const noexcept { return _lower; }
		const Vector<Real>& upper() const noexcept { return _upper; }

		Real lower(int index) const { return _lower[index]; }
		Real upper(int index) const { return _upper[index]; }

		bool HasFiniteLowerBound(int index) const noexcept { return std::isfinite(_lower[index]); }
		bool HasFiniteUpperBound(int index) const noexcept { return std::isfinite(_upper[index]); }
		bool IsFixed(int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			return std::abs(_upper[index] - _lower[index]) <= tolerance;
		}

		Vector<Real> Project(const Vector<Real>& point) const
		{
			validateDimension(point.size(), "BoundConstraints::Project");
			Vector<Real> projected(point.size());
			for (int index = 0; index < point.size(); ++index)
				projected[index] = std::clamp(point[index], _lower[index], _upper[index]);
			return projected;
		}

		template<int N>
		VectorN<Real, N> Project(const VectorN<Real, N>& point) const
		{
			validateDimension(N, "BoundConstraints::Project");
			VectorN<Real, N> projected;
			for (int index = 0; index < N; ++index)
				projected[index] = std::clamp(point[index], _lower[index], _upper[index]);
			return projected;
		}

		void ProjectInPlace(Vector<Real>& point) const
		{
			validateDimension(point.size(), "BoundConstraints::ProjectInPlace");
			for (int index = 0; index < point.size(); ++index)
				point[index] = std::clamp(point[index], _lower[index], _upper[index]);
		}

		template<int N>
		void ProjectInPlace(VectorN<Real, N>& point) const
		{
			validateDimension(N, "BoundConstraints::ProjectInPlace");
			for (int index = 0; index < N; ++index)
				point[index] = std::clamp(point[index], _lower[index], _upper[index]);
		}

		bool IsFeasible(const Vector<Real>& point, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			validateDimension(point.size(), "BoundConstraints::IsFeasible");
			return MaxViolation(point) <= tolerance;
		}

		template<int N>
		bool IsFeasible(const VectorN<Real, N>& point, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			validateDimension(N, "BoundConstraints::IsFeasible");
			return MaxViolation(point) <= tolerance;
		}

		Real MaxViolation(const Vector<Real>& point) const
		{
			validateDimension(point.size(), "BoundConstraints::MaxViolation");
			Real violation = REAL(0.0);
			for (int index = 0; index < point.size(); ++index) {
				if (point[index] < _lower[index])
					violation = std::max(violation, _lower[index] - point[index]);
				if (point[index] > _upper[index])
					violation = std::max(violation, point[index] - _upper[index]);
			}
			return violation;
		}

		template<int N>
		Real MaxViolation(const VectorN<Real, N>& point) const
		{
			validateDimension(N, "BoundConstraints::MaxViolation");
			Real violation = REAL(0.0);
			for (int index = 0; index < N; ++index) {
				if (point[index] < _lower[index])
					violation = std::max(violation, _lower[index] - point[index]);
				if (point[index] > _upper[index])
					violation = std::max(violation, point[index] - _upper[index]);
			}
			return violation;
		}

		bool IsLowerActive(const Vector<Real>& point, int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			validateDimension(point.size(), "BoundConstraints::IsLowerActive");
			return HasFiniteLowerBound(index) && std::abs(point[index] - _lower[index]) <= tolerance;
		}

		bool IsUpperActive(const Vector<Real>& point, int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			validateDimension(point.size(), "BoundConstraints::IsUpperActive");
			return HasFiniteUpperBound(index) && std::abs(point[index] - _upper[index]) <= tolerance;
		}

		bool IsActive(const Vector<Real>& point, int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			return IsLowerActive(point, index, tolerance) || IsUpperActive(point, index, tolerance);
		}

		bool IsFree(const Vector<Real>& point, int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			return !IsActive(point, index, tolerance);
		}

		template<int N>
		bool IsLowerActive(const VectorN<Real, N>& point, int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			validateDimension(N, "BoundConstraints::IsLowerActive");
			return HasFiniteLowerBound(index) && std::abs(point[index] - _lower[index]) <= tolerance;
		}

		template<int N>
		bool IsUpperActive(const VectorN<Real, N>& point, int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			validateDimension(N, "BoundConstraints::IsUpperActive");
			return HasFiniteUpperBound(index) && std::abs(point[index] - _upper[index]) <= tolerance;
		}

		template<int N>
		bool IsActive(const VectorN<Real, N>& point, int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			return IsLowerActive(point, index, tolerance) || IsUpperActive(point, index, tolerance);
		}

		template<int N>
		bool IsFree(const VectorN<Real, N>& point, int index, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			return !IsActive(point, index, tolerance);
		}
	};
}

#endif // MML_BOUND_CONSTRAINTS_H
