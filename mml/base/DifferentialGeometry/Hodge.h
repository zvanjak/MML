///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/Hodge.h                                        ///
///  Description: Hodge star and exterior-algebra-derived vector operations           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DIFFERENTIAL_GEOMETRY_HODGE_H
#define MML_DIFFERENTIAL_GEOMETRY_HODGE_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/base/Matrix/MatrixNM.h>

#include <array>
#include <cmath>
#include <limits>
#include <type_traits>

namespace MML
{
	template<int N, class FrameTag> class Metric;
	enum class Orientation
	{
		Positive = 1,
		Negative = -1
	};

	namespace Detail
	{
		template<class Number>
		struct HodgeScaled
		{
			Number mantissa{};
			int exponent = 0;

			HodgeScaled() = default;
			explicit HodgeScaled(Number value, int binaryExponent = 0) : mantissa(value), exponent(binaryExponent) {
				if (!std::isfinite(value))
					throw GeometryError("Hodge star requires finite metric and form components");
				Normalize();
			}

			void Normalize() {
				if (mantissa == 0) {
					exponent = 0;
					return;
				}
				int shift;
				mantissa = std::frexp(mantissa, &shift);
				exponent += shift;
			}

			bool AbsLess(const HodgeScaled& other) const {
				if (mantissa == 0) return other.mantissa != 0;
				if (other.mantissa == 0) return false;
				return exponent == other.exponent ? std::abs(mantissa) < std::abs(other.mantissa) : exponent < other.exponent;
			}

			void Multiply(const HodgeScaled& other) {
				mantissa *= other.mantissa;
				exponent += other.exponent;
				Normalize();
			}

			void Divide(const HodgeScaled& other) {
				if (other.mantissa == 0)
					throw GeometryError("Hodge star encountered a singular metric");
				mantissa /= other.mantissa;
				exponent -= other.exponent;
				Normalize();
			}

			void Add(const HodgeScaled& other) {
				if (other.mantissa == 0) return;
				if (mantissa == 0) {
					*this = other;
					return;
				}
				const int commonExponent = std::max(exponent, other.exponent);
				mantissa = std::scalbn(mantissa, exponent - commonExponent)
					+ std::scalbn(other.mantissa, other.exponent - commonExponent);
				exponent = commonExponent;
				Normalize();
			}

			HodgeScaled SqrtAbs() const {
				if (mantissa == 0)
					throw GeometryError("Hodge star requires a nonsingular metric");
				HodgeScaled result(*this);
				result.mantissa = std::abs(result.mantissa);
				if (result.exponent % 2 != 0) {
					result.mantissa *= 2;
					--result.exponent;
				}
				result.mantissa = std::sqrt(result.mantissa);
				result.exponent /= 2;
				result.Normalize();
				return result;
			}

			template<class Target = Number>
			Target Value() const {
				const Number expanded = std::scalbn(mantissa, exponent);
				if (!std::isfinite(expanded) || (expanded == 0 && mantissa != 0)
					|| std::abs(expanded) > Number(std::numeric_limits<Target>::max()))
					throw GeometryError("Hodge star component is not representable");
				const Target value = static_cast<Target>(expanded);
				if (value == 0 && mantissa != 0)
					throw GeometryError("Hodge star component is not representable");
				return value;
			}
		};

		template<int N>
		HodgeScaled<Real> HodgeDeterminant(const MatrixNM<Real, N, N>& matrix)
		{
			// Partial-pivoting LU, retaining binary exponents for entries and pivots.
			std::array<HodgeScaled<Real>, N * N> a;
			for (int row = 0; row < N; ++row)
				for (int col = 0; col < N; ++col)
					a[row * N + col] = HodgeScaled<Real>(matrix(row, col));
			HodgeScaled<Real> det(1);
			for (int col = 0; col < N; col++) {
				int pivot = col;
				for (int row = col + 1; row < N; row++)
					if (a[pivot * N + col].AbsLess(a[row * N + col]))
						pivot = row;
				if (a[pivot * N + col].mantissa == 0)
					throw GeometryError("Hodge star requires a nonsingular metric");
				if (pivot != col) {
					for (int j = 0; j < N; j++)
						std::swap(a[col * N + j], a[pivot * N + j]);
					det.mantissa = -det.mantissa;
				}
				const auto pivotValue = a[col * N + col];
				det.Multiply(pivotValue);
				for (int row = col + 1; row < N; row++) {
					auto factor = a[row * N + col];
					if (factor.mantissa == 0) continue;
					factor.Divide(pivotValue);
					for (int j = col + 1; j < N; j++) {
						auto product = a[col * N + j];
						product.Multiply(factor);
						product.mantissa = -product.mantissa;
						a[row * N + j].Add(product);
					}
				}
			}
			return det;
		}

		template<int N>
		int LeviCivitaSymbol(const std::array<int, N>& indices)
		{
			bool seen[N] = {};
			for (int i = 0; i < N; i++) {
				if (indices[i] < 0 || indices[i] >= N || seen[indices[i]])
					return 0;
				seen[indices[i]] = true;
			}

			int inversions = 0;
			for (int i = 0; i < N; i++)
				for (int j = i + 1; j < N; j++)
					if (indices[i] > indices[j])
						inversions++;
			return inversions % 2 == 0 ? 1 : -1;
		}

		template<int N, int K>
		std::array<int, K == 0 ? 1 : K> UnflattenMultiIndex(int flat)
		{
			std::array<int, K == 0 ? 1 : K> indices{};
			for (int i = K - 1; i >= 0; i--) {
				indices[i] = flat % N;
				flat /= N;
			}
			return indices;
		}

		template<class Scalar, int N, int K, class FrameTag>
		Scalar RaisedComponent(const DifferentialForm<Scalar, N, K, FrameTag>& form,
			const Metric<N, FrameTag>& metric,
			const std::array<int, K == 0 ? 1 : K>& raisedIndices)
		{
			Scalar result{};
			for (int flat = 0; flat < Detail::IntPow(N, K); flat++) {
				auto lowerIndices = UnflattenMultiIndex<N, K>(flat);
				Scalar term = form.ComponentAt(lowerIndices);
				for (int i = 0; i < K; i++)
					term = static_cast<Scalar>(term * metric.contravariantComponents()(raisedIndices[i], lowerIndices[i]));
				result += term;
			}
			return result;
		}
	}

	template<class Scalar, int N, int K, class FrameTag>
	DifferentialForm<Scalar, N, N - K, FrameTag> HodgeStar(
		const DifferentialForm<Scalar, N, K, FrameTag>& form,
		const Metric<N, FrameTag>& metric,
		Orientation orientation = Orientation::Positive)
	{
		constexpr int ComplementDegree = N - K;
		DifferentialForm<Scalar, N, ComplementDegree, FrameTag> result;
		if (orientation != Orientation::Positive && orientation != Orientation::Negative)
			throw GeometryError("Hodge star requires a valid orientation");
		const auto volumeScale = Detail::HodgeDeterminant(metric.covariantComponents()).SqrtAbs();
		for (int row = 0; row < N; ++row)
			for (int col = 0; col < N; ++col)
				if (!std::isfinite(metric.contravariantComponents()(row, col)))
					throw GeometryError("Hodge star requires a finite inverse metric");
		Scalar orientationSign = orientation == Orientation::Positive ? Scalar{ 1 } : Scalar{ -1 };

		for (int outFlat = 0; outFlat < result.ComponentCount; outFlat++) {
			auto outIndices = Detail::UnflattenMultiIndex<N, ComplementDegree>(outFlat);
			if constexpr (MMLScalar<Scalar>) {
				const auto calculatePart = [&](auto selectPart) {
					using Part = decltype(selectPart(Scalar{}));
					using Work = std::common_type_t<Real, Part>;
					Detail::HodgeScaled<Work> component;
					for (int raisedFlat = 0; raisedFlat < Detail::IntPow(N, K); ++raisedFlat) {
						const auto raisedIndices = Detail::UnflattenMultiIndex<N, K>(raisedFlat);
						std::array<int, N> epsilonIndices{};
						for (int i = 0; i < K; ++i) epsilonIndices[i] = raisedIndices[i];
						for (int i = 0; i < ComplementDegree; ++i) epsilonIndices[K + i] = outIndices[i];
						const int epsilon = Detail::LeviCivitaSymbol<N>(epsilonIndices);
						if (epsilon == 0) continue;
						for (int lowerFlat = 0; lowerFlat < form.ComponentCount; ++lowerFlat) {
							const auto lowerIndices = Detail::UnflattenMultiIndex<N, K>(lowerFlat);
							Detail::HodgeScaled<Work> term(selectPart(form.ComponentAt(lowerIndices)));
							if (term.mantissa == 0) continue;
							term.Multiply(Detail::HodgeScaled<Work>(volumeScale.mantissa, volumeScale.exponent));
							for (int i = 0; i < K; ++i)
								term.Multiply(Detail::HodgeScaled<Work>(
									metric.contravariantComponents()(raisedIndices[i], lowerIndices[i])));
							term.mantissa *= epsilon;
							component.Add(term);
						}
					}
					for (int divisor = 2; divisor <= K; ++divisor)
						component.Divide(Detail::HodgeScaled<Work>(Work(divisor)));
					if (orientation == Orientation::Negative) component.mantissa = -component.mantissa;
					return component.template Value<Part>();
				};
				if constexpr (MMLComplex<Scalar>) {
					result.ComponentAt(outIndices) = Scalar(
						calculatePart([](const Scalar& value) { return value.real(); }),
						calculatePart([](const Scalar& value) { return value.imag(); }));
				} else {
					result.ComponentAt(outIndices) = calculatePart([](Scalar value) { return value; });
				}
				continue;
			}

			Scalar component{};

			for (int raisedFlat = 0; raisedFlat < Detail::IntPow(N, K); raisedFlat++) {
				auto raisedIndices = Detail::UnflattenMultiIndex<N, K>(raisedFlat);
				std::array<int, N> epsilonIndices{};
				for (int i = 0; i < K; i++)
					epsilonIndices[i] = raisedIndices[i];
				for (int i = 0; i < ComplementDegree; i++)
					epsilonIndices[K + i] = outIndices[i];

				int epsilon = Detail::LeviCivitaSymbol<N>(epsilonIndices);
				if (epsilon != 0)
					component += ScalarDetail::ScalarFromInteger<Scalar>(epsilon) * Detail::RaisedComponent(form, metric, raisedIndices);
			}

			result.ComponentAt(outIndices) = orientationSign * Scalar(volumeScale.template Value<Real>())
				* component / ScalarDetail::ScalarFromInteger<Scalar>(Detail::Factorial(K));
		}

		return result;
	}

	template<class Scalar, int N, int K, class FrameTag>
	DifferentialForm<Scalar, N, N - K, FrameTag> hodge_star(
		const DifferentialForm<Scalar, N, K, FrameTag>& form,
		const Metric<N, FrameTag>& metric,
		Orientation orientation = Orientation::Positive)
	{
		return HodgeStar(form, metric, orientation);
	}

	template<int N, class FrameTag>
	TangentVector<3, FrameTag> Cross(
		const TangentVector<N, FrameTag>& a,
		const TangentVector<N, FrameTag>& b,
		const Metric<N, FrameTag>& metric,
		Orientation orientation = Orientation::Positive)
		requires (N == 3)
	{
		Form1<3, FrameTag> aFlat = ToForm(Flat(metric, a));
		Form1<3, FrameTag> bFlat = ToForm(Flat(metric, b));
		Form1<3, FrameTag> crossCovector = HodgeStar(Wedge(aFlat, bFlat), metric, orientation);
		return Sharp(metric, ToCovector(crossCovector));
	}

	template<int N, class FrameTag>
	TangentVector<3, FrameTag> cross(
		const TangentVector<N, FrameTag>& a,
		const TangentVector<N, FrameTag>& b,
		const Metric<N, FrameTag>& metric,
		Orientation orientation = Orientation::Positive)
		requires (N == 3)
	{
		return Cross(a, b, metric, orientation);
	}
}

#endif // MML_DIFFERENTIAL_GEOMETRY_HODGE_H
