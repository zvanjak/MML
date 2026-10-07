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
#include <mml/base/DifferentialGeometry/Metric.h>

#include <array>
#include <cmath>
#include <type_traits>

namespace MML
{
	enum class Orientation
	{
		Positive = 1,
		Negative = -1
	};

	namespace Detail
	{
		template<int N>
		Real Determinant(const MatrixNM<Real, N, N>& matrix)
		{
			MatrixNM<Real, N, N> a(matrix);
			Real det = REAL(1.0);

			for (int col = 0; col < N; col++) {
				int pivot = col;
				Real maxAbs = std::abs(a(col, col));
				for (int row = col + 1; row < N; row++) {
					Real value = std::abs(a(row, col));
					if (value > maxAbs) {
						maxAbs = value;
						pivot = row;
					}
				}

				if (maxAbs == REAL(0.0))
					return REAL(0.0);

				if (pivot != col) {
					for (int j = 0; j < N; j++)
						std::swap(a(col, j), a(pivot, j));
					det = -det;
				}

				Real pivotValue = a(col, col);
				det *= pivotValue;
				for (int row = col + 1; row < N; row++) {
					Real factor = a(row, col) / pivotValue;
					for (int j = col + 1; j < N; j++)
						a(row, j) -= factor * a(col, j);
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
					term *= metric.contravariantComponents()(raisedIndices[i], lowerIndices[i]);
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
		Real det = Detail::Determinant(metric.covariantComponents());
		Real volumeScale = std::sqrt(std::abs(det));
		Scalar orientationSign = orientation == Orientation::Positive ? Scalar{ 1 } : Scalar{ -1 };

		for (int outFlat = 0; outFlat < result.ComponentCount; outFlat++) {
			auto outIndices = Detail::UnflattenMultiIndex<N, ComplementDegree>(outFlat);
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
					component += Scalar(epsilon) * Detail::RaisedComponent(form, metric, raisedIndices);
			}

			result.ComponentAt(outIndices) = orientationSign * Scalar(volumeScale) * component / Scalar(Detail::Factorial(K));
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