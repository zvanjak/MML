#if !defined MML_FIRST_DERIVATIVE_STENCIL_H
#define MML_FIRST_DERIVATIVE_STENCIL_H

#include "DerivationBase.h"

#include <array>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>

namespace MML::Derivation::Detail
{
	enum class FirstDerivativeOrder { One = 1, Two = 2, Four = 4, Six = 6, Eight = 8 };

	template<FirstDerivativeOrder Order>
	struct FirstDerivativeStencil;
	template<FirstDerivativeOrder Order>
	struct FirstDerivativePartialEvaluationOrder;

	template<>
	struct FirstDerivativeStencil<FirstDerivativeOrder::One>
	{
		static constexpr std::array value_offsets{1, 0};
		static constexpr std::array error_offsets{1, 0, -1};
		static constexpr std::array<Real, 2> value_coefficients{REAL(1.0), REAL(-1.0)};
		static constexpr std::array<Real, 3> truncation_error_coefficients{REAL(0.5), REAL(-1.0), REAL(0.5)};
		static constexpr int max_reach = 1;
	};

	template<>
	struct FirstDerivativeStencil<FirstDerivativeOrder::Two>
	{
		static constexpr std::array value_offsets{1, -1};
		static constexpr std::array error_offsets{1, -1, 2, -2};
		static constexpr std::array<Real, 2> value_coefficients{REAL(0.5), REAL(-0.5)};
		static constexpr std::array<Real, 4> truncation_error_coefficients{REAL(-1.0 / 6.0), REAL(1.0 / 6.0), REAL(1.0 / 12.0), REAL(-1.0 / 12.0)};
		static constexpr int max_reach = 2;
	};

	template<>
	struct FirstDerivativeStencil<FirstDerivativeOrder::Four>
	{
		static constexpr std::array value_offsets{1, -1, 2, -2};
		static constexpr std::array error_offsets{1, -1, 2, -2, 3, -3};
		static constexpr std::array<Real, 4> value_coefficients{REAL(2.0 / 3.0), REAL(-2.0 / 3.0), REAL(-1.0 / 12.0), REAL(1.0 / 12.0)};
		static constexpr std::array<Real, 6> truncation_error_coefficients{REAL(1.0 / 12.0), REAL(-1.0 / 12.0), REAL(-1.0 / 15.0), REAL(1.0 / 15.0), REAL(1.0 / 60.0), REAL(-1.0 / 60.0)};
		static constexpr int max_reach = 3;
	};

	template<>
	struct FirstDerivativeStencil<FirstDerivativeOrder::Six>
	{
		static constexpr std::array value_offsets{1, -1, -2, 2, 3, -3};
		static constexpr std::array error_offsets{1, -1, -2, 2, 3, -3, 4, -4};
		static constexpr std::array<Real, 6> value_coefficients{REAL(3.0 / 4.0), REAL(-3.0 / 4.0), REAL(3.0 / 20.0), REAL(-3.0 / 20.0), REAL(1.0 / 60.0), REAL(-1.0 / 60.0)};
		static constexpr std::array<Real, 8> truncation_error_coefficients{REAL(-1.0 / 20.0), REAL(1.0 / 20.0), REAL(-1.0 / 20.0), REAL(1.0 / 20.0), REAL(-3.0 / 140.0), REAL(3.0 / 140.0), REAL(1.0 / 280.0), REAL(-1.0 / 280.0)};
		static constexpr int max_reach = 4;
	};

	template<>
	struct FirstDerivativeStencil<FirstDerivativeOrder::Eight>
	{
		static constexpr std::array value_offsets{1, -1, -2, 2, 3, -3, -4, 4};
		static constexpr std::array error_offsets{1, -1, -2, 2, 3, -3, -4, 4, 5, -5};
		static constexpr std::array<Real, 8> value_coefficients{REAL(4.0 / 5.0), REAL(-4.0 / 5.0), REAL(1.0 / 5.0), REAL(-1.0 / 5.0), REAL(4.0 / 105.0), REAL(-4.0 / 105.0), REAL(1.0 / 280.0), REAL(-1.0 / 280.0)};
		static constexpr std::array<Real, 10> truncation_error_coefficients{REAL(1.0 / 30.0), REAL(-1.0 / 30.0), REAL(4.0 / 105.0), REAL(-4.0 / 105.0), REAL(3.0 / 140.0), REAL(-3.0 / 140.0), REAL(2.0 / 315.0), REAL(-2.0 / 315.0), REAL(1.0 / 1260.0), REAL(-1.0 / 1260.0)};
		static constexpr int max_reach = 5;
	};

	template<>
	struct FirstDerivativePartialEvaluationOrder<FirstDerivativeOrder::One>
	{
		static constexpr std::array value_offsets{0, 1};
		static constexpr std::array error_offsets{0, 1, -1};
	};
	template<>
	struct FirstDerivativePartialEvaluationOrder<FirstDerivativeOrder::Two>
		: FirstDerivativeStencil<FirstDerivativeOrder::Two> {};
	template<>
	struct FirstDerivativePartialEvaluationOrder<FirstDerivativeOrder::Four>
		: FirstDerivativeStencil<FirstDerivativeOrder::Four> {};
	template<>
	struct FirstDerivativePartialEvaluationOrder<FirstDerivativeOrder::Six>
	{
		static constexpr std::array value_offsets{1, -1, 2, -2, 3, -3};
		static constexpr std::array error_offsets{1, -1, 2, -2, 3, -3, 4, -4};
	};
	template<>
	struct FirstDerivativePartialEvaluationOrder<FirstDerivativeOrder::Eight>
	{
		static constexpr std::array value_offsets{1, -1, 2, -2, 3, -3, 4, -4};
		static constexpr std::array error_offsets{1, -1, 2, -2, 3, -3, 4, -4, 5, -5};
	};

	template<typename Value>
	struct FirstDerivativeStencilResult
	{
		Value value;
		Real error = REAL(0.0);
		int function_evaluations = 0;
	};

	template<FirstDerivativeOrder Order, typename Sample, typename Magnitude, typename ValueOffsets, typename ErrorOffsets,
	         typename ErrorEstimate = std::nullptr_t>
	auto EvaluateFirstDerivativeStencilWithOffsets(Sample&& sample, Magnitude&& magnitude, Real h, bool estimate_error,
	                                              const ValueOffsets& value_offsets, const ErrorOffsets& error_offsets,
	                                              ErrorEstimate&& error_estimate = nullptr)
	{
		using Stencil = FirstDerivativeStencil<Order>;
		using Value = std::remove_cvref_t<std::invoke_result_t<Sample&, int>>;
		std::array<std::optional<Value>, 2 * Stencil::max_reach + 1> samples;
		int function_evaluations = 0;

		auto at = [&](int offset) -> const Value& {
			auto& cached = samples[static_cast<std::size_t>(offset + Stencil::max_reach)];
			if (!cached)
			{
				cached.emplace(std::invoke(sample, offset));
				++function_evaluations;
			}
			return *cached;
		};
		auto norm = [&](const Value& value) { return static_cast<Real>(std::invoke(magnitude, value)); };
		for (int offset : value_offsets)
			at(offset);
		if (estimate_error)
			for (int offset : error_offsets)
				at(offset);

		Value value;
		Real error = REAL(0.0);
		if constexpr (Order == FirstDerivativeOrder::One)
		{
			const Value& yh = at(1);
			const Value& y0 = at(0);
			value = (yh - y0) / h;
			if (estimate_error)
			{
				const Value& ym = at(-1);
				error = norm(yh - REAL(2.0) * y0 + ym) / (REAL(2.0) * h) + (norm(yh) + norm(y0)) * Constants::Eps / h;
			}
		}
		else if constexpr (Order == FirstDerivativeOrder::Two)
		{
			const Value& yh = at(1);
			const Value& ymh = at(-1);
			Value diff = yh - ymh;
			value = diff / (REAL(2.0) * h);
			if (estimate_error)
				error = norm((at(2) - at(-2)) / REAL(2.0) - diff) / (REAL(6.0) * h)
				      + Constants::Eps * (norm(yh) + norm(ymh)) / (REAL(2.0) * h);
		}
		else if constexpr (Order == FirstDerivativeOrder::Four)
		{
			const Value& yh = at(1);
			const Value& ymh = at(-1);
			const Value& y2h = at(2);
			const Value& ym2h = at(-2);
			Value y2 = ym2h - y2h;
			Value y1 = yh - ymh;
			value = (y2 + REAL(8.0) * y1) / (REAL(12.0) * h);
			if (estimate_error)
				error = norm((at(3) - at(-3)) / REAL(2.0) + REAL(2.0) * (ym2h - y2h)
				             + REAL(5.0) * (yh - ymh) / REAL(2.0)) / (REAL(30.0) * h)
			      + Constants::Eps * (norm(y2h) + norm(ym2h) + REAL(8.0) * (norm(ymh) + norm(yh))) / (REAL(12.0) * h);
		}
		else if constexpr (Order == FirstDerivativeOrder::Six)
		{
			const Value& yh = at(1);
			const Value& ymh = at(-1);
			Value y1 = yh - ymh;
			Value y2 = at(-2) - at(2);
			Value y3 = at(3) - at(-3);
			value = (y3 + REAL(9.0) * y2 + REAL(45.0) * y1) / (REAL(60.0) * h);
			if (estimate_error)
				error = norm((at(4) - at(-4) - REAL(6.0) * y3 - REAL(14.0) * y1 - REAL(14.0) * y2) / REAL(2.0)) / (REAL(140.0) * h)
			      + REAL(5.0) * (norm(yh) + norm(ymh)) * Constants::Eps / h;
		}
		else
		{
			const Value& yh = at(1);
			const Value& ymh = at(-1);
			Value y1 = yh - ymh;
			Value y2 = at(-2) - at(2);
			Value y3 = at(3) - at(-3);
			Value y4 = at(-4) - at(4);
			Value tmp1 = REAL(3.0) * y4 / REAL(8.0) + REAL(4.0) * y3;
			Value tmp2 = REAL(21.0) * y2 + REAL(84.0) * y1;
			value = (tmp1 + tmp2) / (REAL(105.0) * h);
			if (estimate_error)
				error = norm((at(5) - at(-5)) / REAL(2.0) + REAL(4.0) * y4 + REAL(27.0) * y3 / REAL(2.0)
				             + REAL(24.0) * y2 + REAL(21.0) * y1) / (REAL(630.0) * h)
			      + REAL(7.0) * (norm(yh) + norm(ymh)) * Constants::Eps / h;
		}
		if constexpr (!std::is_same_v<std::remove_cvref_t<ErrorEstimate>, std::nullptr_t>)
			if (estimate_error)
				error = std::invoke(error_estimate, at, norm, h);

		return FirstDerivativeStencilResult<Value>{std::move(value), error, function_evaluations};
	}

	template<FirstDerivativeOrder Order, typename Sample, typename Magnitude>
	auto EvaluateFirstDerivativeStencil(Sample&& sample, Magnitude&& magnitude, Real h, bool estimate_error)
	{
		using Stencil = FirstDerivativeStencil<Order>;
		return EvaluateFirstDerivativeStencilWithOffsets<Order>(
			std::forward<Sample>(sample), std::forward<Magnitude>(magnitude), h, estimate_error,
			Stencil::value_offsets, Stencil::error_offsets);
	}

	template<FirstDerivativeOrder Order, typename Sample>
	auto EvaluateScalarPartialFirstDerivativeStencil(Sample&& sample, Real h, bool estimate_error)
	{
		using EvaluationOrder = FirstDerivativePartialEvaluationOrder<Order>;
		return EvaluateFirstDerivativeStencilWithOffsets<Order>(
			std::forward<Sample>(sample), [](Real value) { return std::abs(value); }, h, estimate_error,
			EvaluationOrder::value_offsets, EvaluationOrder::error_offsets);
	}
}

#endif // MML_FIRST_DERIVATIVE_STENCIL_H