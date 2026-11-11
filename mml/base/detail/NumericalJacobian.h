#ifndef MML_BASE_DETAIL_NUMERICAL_JACOBIAN_H
#define MML_BASE_DETAIL_NUMERICAL_JACOBIAN_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace MML::NumericalJacobianDetail
{
	inline void ValidateJacobianStep(Real h)
	{
		if (!std::isfinite(h) || h < 0)
			throw ArgumentError("Jacobian: step must be finite and nonnegative (zero selects automatic)");
	}

	template<typename Point>
	void ValidateJacobianPoint(const Point& point)
	{
		for (int i = 0; i < point.size(); ++i)
			if (!std::isfinite(point[i]))
				throw ArgumentError("Jacobian: non-finite input coordinate");
	}

	inline void ValidateJacobianStencil(Real x, Real h, int reach)
	{
		const Real minus = x - h, plus = x + h;
		if (!std::isfinite(h) || h <= 0 || !std::isfinite(minus) || !std::isfinite(plus)
			|| !(minus < x && x < plus))
			throw NumericalMethodError("Jacobian: unrepresentable stencil");
		if (reach == 2) {
			const Real minus2 = x - 2 * h, plus2 = x + 2 * h;
			if (!std::isfinite(minus2) || !std::isfinite(plus2)
				|| !(minus2 < minus && plus < plus2))
				throw NumericalMethodError("Jacobian: unrepresentable outer stencil");
		}
	}

	template<typename Value>
	void ValidateJacobianSample(const Value& value, int expectedSize)
	{
		if (value.size() != expectedSize)
			throw VectorDimensionError("Jacobian: output dimension mismatch", expectedSize, value.size());
		for (int i = 0; i < expectedSize; ++i)
			if (!std::isfinite(value[i]))
				throw NumericalMethodError("Jacobian: non-finite function value");
	}

	inline Real FinishJacobianQuotient(Real numerator, int exponent, Real h, Real divisor)
	{
		int stepExponent;
		const Real stepMantissa = std::frexp(h, &stepExponent);
		const Real value = std::scalbn(numerator / (divisor * stepMantissa), exponent - stepExponent);
		if (!std::isfinite(value))
			throw NumericalMethodError("Jacobian: non-finite derivative entry");
		return value;
	}

	inline Real JacobianDifference2(Real plus, Real minus, Real h)
	{
		Real difference = plus - minus;
		int exponent;
		if (std::isfinite(difference)) {
			std::frexp(std::abs(difference), &exponent);
			difference = std::scalbn(difference, -exponent);
		}
		else {
			std::frexp(std::max(std::abs(plus), std::abs(minus)), &exponent);
			difference = std::scalbn(plus, -exponent) - std::scalbn(minus, -exponent);
		}
		return FinishJacobianQuotient(difference, exponent, h, Real(2));
	}

	inline Real JacobianCentralHalfSpan(Real x, Real h)
	{
		const Real plus = x + h, minus = x - h;
		const Real span = plus - minus;
		return std::isfinite(span) ? span / 2 : plus / 2 - minus / 2;
	}

	// Kept below Core so interface defaults and explicit ODE differentiation share one kernel.
	template<typename Derivatives>
	void JacobianODE(const Derivatives& derivs, Real t, const Vector<Real>& x, Matrix<Real>& jac, Real h = 0)
	{
		ValidateJacobianStep(h);
		const int n = x.size();
		if (n <= 0 || !std::isfinite(t))
			throw ArgumentError("ODE Jacobian: invalid dimensions or time");
		ValidateJacobianPoint(x);
		Vector<Real> steps(n);
		for (int column = 0; column < n; ++column) {
			steps[column] = (h > 0 ? h : Precision::CbrtEps) * std::max(Real(1), std::abs(x[column]));
			ValidateJacobianStencil(x[column], steps[column], 1);
		}

		Matrix<Real> result(n, n);
		Vector<Real> point(x), plus(n), minus(n);
		for (int column = 0; column < n; ++column) {
			for (int row = 0; row < n; ++row) {
				plus[row] = std::numeric_limits<Real>::quiet_NaN();
				minus[row] = std::numeric_limits<Real>::quiet_NaN();
			}
			point[column] = x[column] + steps[column];
			derivs(t, point, plus);
			ValidateJacobianSample(plus, n);
			point[column] = x[column] - steps[column];
			derivs(t, point, minus);
			ValidateJacobianSample(minus, n);
			point[column] = x[column];
			const Real halfSpan = JacobianCentralHalfSpan(x[column], steps[column]);
			for (int row = 0; row < n; ++row)
				result(row, column) = JacobianDifference2(plus[row], minus[row], halfSpan);
		}
		jac = std::move(result);
	}
}

#endif
