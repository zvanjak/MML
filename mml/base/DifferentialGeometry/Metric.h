///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/Metric.h                                       ///
///  Description: Typed metric value object and musical rank-1 operations             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DIFFERENTIAL_GEOMETRY_METRIC_H
#define MML_DIFFERENTIAL_GEOMETRY_METRIC_H

#include <mml/MMLBase.h>

#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Tensor/Rank1Tensor.h>
#include <mml/base/Vector/VectorN.h>

namespace MML
{
	template<int N, class FrameTag>
	class Metric
	{
		MatrixNM<Real, N, N> _covariantComponents;
		MatrixNM<Real, N, N> _contravariantComponents;

	public:
		using value_type = Real;
		using frame_type = FrameTag;
		static constexpr int Dimension = N;

		Metric() : Metric(MatrixNM<Real, N, N>(REAL(1.0))) { }
		explicit Metric(const MatrixNM<Real, N, N>& covariantComponents)
			: _covariantComponents(covariantComponents),
			  _contravariantComponents(covariantComponents.GetInverse())
		{
		}

		static Metric Euclidean()
		{
			return Metric(MatrixNM<Real, N, N>(REAL(1.0)));
		}

		static Metric Diagonal(const VectorN<Real, N>& diagonal)
		{
			MatrixNM<Real, N, N> components;
			for (int i = 0; i < N; i++)
				components(i, i) = diagonal[i];
			return Metric(components);
		}

		Real operator()(int i, int j) const { return _covariantComponents(i, j); }

		const MatrixNM<Real, N, N>& covariantComponents() const noexcept { return _covariantComponents; }
		const MatrixNM<Real, N, N>& contravariantComponents() const noexcept { return _contravariantComponents; }
	};

	template<int N, class FrameTag>
	Metric<N, FrameTag> EuclideanMetric()
	{
		return Metric<N, FrameTag>::Euclidean();
	}

	template<int N, class FrameTag>
	Covector<N, FrameTag> Flat(const Metric<N, FrameTag>& metric, const TangentVector<N, FrameTag>& v)
	{
		Covector<N, FrameTag> result;
		for (int i = 0; i < N; i++) {
			result[i] = REAL(0.0);
			for (int j = 0; j < N; j++)
				result[i] += metric.covariantComponents()(i, j) * v[j];
		}
		return result;
	}

	template<int N, class FrameTag>
	Covector<N, FrameTag> flat(const Metric<N, FrameTag>& metric, const TangentVector<N, FrameTag>& v)
	{
		return Flat(metric, v);
	}

	template<int N, class FrameTag>
	TangentVector<N, FrameTag> Sharp(const Metric<N, FrameTag>& metric, const Covector<N, FrameTag>& alpha)
	{
		TangentVector<N, FrameTag> result;
		for (int i = 0; i < N; i++) {
			result[i] = REAL(0.0);
			for (int j = 0; j < N; j++)
				result[i] += metric.contravariantComponents()(i, j) * alpha[j];
		}
		return result;
	}

	template<int N, class FrameTag>
	TangentVector<N, FrameTag> sharp(const Metric<N, FrameTag>& metric, const Covector<N, FrameTag>& alpha)
	{
		return Sharp(metric, alpha);
	}

	template<int N, class FrameTag>
	Real Inner(const Metric<N, FrameTag>& metric, const TangentVector<N, FrameTag>& a,
		const TangentVector<N, FrameTag>& b)
	{
		return Pair(Flat(metric, a), b);
	}

	template<int N, class FrameTag>
	Real inner(const Metric<N, FrameTag>& metric, const TangentVector<N, FrameTag>& a,
		const TangentVector<N, FrameTag>& b)
	{
		return Inner(metric, a, b);
	}
}

#endif // MML_DIFFERENTIAL_GEOMETRY_METRIC_H