///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/FunctionSpace1D.h                                    ///
///  Description: One-dimensional function-space descriptors                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_SPACE_1D_H
#define MML_FUNCTION_SPACE_1D_H

#include <mml/core/FunctionSpaces/FunctionSpacesBase.h>

#include <mml/core/Integration/Integration1D.h>
#include <mml/interfaces/IFunction.h>

#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

namespace MML::FunctionSpaces
{
	class FunctionSpaceError : public std::runtime_error
	{
	public:
		explicit FunctionSpaceError(const std::string& message)
			: std::runtime_error("FunctionSpaceError: " + message) { }
	};

	class FunctionSpaceInputError : public std::domain_error
	{
	public:
		explicit FunctionSpaceInputError(const std::string& message)
			: std::domain_error("FunctionSpaceInputError: " + message) { }
	};

	class FunctionSpace1D
	{
	protected:
		static void ValidateFiniteDomain(Real min, Real max, const char* context)
		{
			if (!std::isfinite(min) || !std::isfinite(max))
				throw FunctionSpaceInputError(std::string(context) + ": domain bounds must be finite");
			if (min >= max)
				throw FunctionSpaceInputError(std::string(context) + ": domain min must be less than domain max");
		}

	public:
		virtual ~FunctionSpace1D() = default;

		virtual Real domainMin() const noexcept = 0;
		virtual Real domainMax() const noexcept = 0;
		virtual Real weight(Real x) const = 0;

		Real domainLength() const noexcept { return domainMax() - domainMin(); }

		bool contains(Real x, Real tolerance = Defaults::VectorIsEqualTolerance) const noexcept
		{
			return x >= domainMin() - tolerance && x <= domainMax() + tolerance;
		}

		Real innerProduct(const IRealFunction& f, const IRealFunction& g, Real eps = Precision::DefaultToleranceStrict) const
		{
			class WeightedProduct : public IRealFunction
			{
				const FunctionSpace1D& _space;
				const IRealFunction& _f;
				const IRealFunction& _g;

			public:
				WeightedProduct(const FunctionSpace1D& space, const IRealFunction& f, const IRealFunction& g)
					: _space(space), _f(f), _g(g) { }

				Real operator()(Real x) const override
				{
					return _f(x) * _g(x) * _space.weight(x);
				}
			};

			WeightedProduct product(*this, f, g);
			return IntegrateTrap(product, domainMin(), domainMax(), eps).value;
		}
	};

	class L2IntervalSpace : public FunctionSpace1D
	{
		Real _min;
		Real _max;

	public:
		L2IntervalSpace(Real min, Real max)
			: _min(min), _max(max)
		{
			ValidateFiniteDomain(_min, _max, "L2IntervalSpace");
		}

		Real domainMin() const noexcept override { return _min; }
		Real domainMax() const noexcept override { return _max; }
		Real weight(Real) const override { return REAL(1.0); }
	};

	class WeightedL2IntervalSpace : public FunctionSpace1D
	{
		Real _min;
		Real _max;
		std::function<Real(Real)> _weight;

	public:
		WeightedL2IntervalSpace(Real min, Real max, std::function<Real(Real)> weightFunction)
			: _min(min), _max(max), _weight(std::move(weightFunction))
		{
			ValidateFiniteDomain(_min, _max, "WeightedL2IntervalSpace");
			if (!_weight)
				throw FunctionSpaceInputError("WeightedL2IntervalSpace: weight function must be callable");
		}

		Real domainMin() const noexcept override { return _min; }
		Real domainMax() const noexcept override { return _max; }

		Real weight(Real x) const override
		{
			Real value = _weight(x);
			if (!std::isfinite(value))
				throw FunctionSpaceInputError("WeightedL2IntervalSpace: weight function returned non-finite value");
			return value;
		}
	};
}

#endif // MML_FUNCTION_SPACE_1D_H
