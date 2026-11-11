#ifndef MML_GEOMETRY_CHECKS_H
#define MML_GEOMETRY_CHECKS_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Geometry/GeometryBase/GeometryArithmetic.h>
#include <algorithm>
#include <cmath>

namespace MML::GeometryDetail
{
	inline void RequireFinite(Real value, const char* operation)
	{
		if (!std::isfinite(value))
			throw DomainError(operation);
	}

	inline void RequireNonnegative(Real value, const char* operation)
	{
		if (!std::isfinite(value) || value < 0)
			throw DomainError(operation);
	}

	inline void RequireRange(Real minimum, Real maximum, const char* operation)
	{
		if (std::isnan(minimum) || std::isnan(maximum) || minimum > maximum ||
			minimum == Constants::PosInf || maximum == Constants::NegInf)
			throw DomainError(operation);
	}

	template<int N>
	void RequireFinite(const VectorN<Real, N>& value, const char* operation)
	{
		for (int index = 0; index < N; ++index)
			RequireFinite(value[index], operation);
	}

	inline Real Checked(Real value)
	{
		if (!std::isfinite(value))
			throw GeometryError("Curve/surface geometry - non-finite arithmetic result");
		return value;
	}

	template<int N>
	VectorN<Real, N> Checked(const VectorN<Real, N>& value)
	{
		for (int index = 0; index < N; ++index)
			Checked(value[index]);
		return value;
	}

	template<int N>
	Real Norm(const VectorN<Real, N>& value)
	{
		Real norm{};
		for (int index = 0; index < N; ++index)
			norm = std::hypot(norm, Checked(value[index]));
		return Checked(norm);
	}

	template<int N>
	VectorN<Real, N> Unit(const VectorN<Real, N>& value)
	{
		return detail::GeometryUnitDirection(value, "Curve/surface geometry - non-finite arithmetic result",
			"Curve/surface geometry - undefined unit vector at a zero vector");
	}
}

#endif
