///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Rank1Tensor.h                                                       ///
///  Description: Typed rank-1 tensor wrappers over VectorN                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_RANK1_TENSOR_H
#define MML_RANK1_TENSOR_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/FrameTags.h>
#include <mml/base/Vector/VectorN.h>

#include <initializer_list>
#include <type_traits>

namespace MML
{
	enum class Variance
	{
		Contravariant,
		Covariant
	};

	template<class Scalar, int N, Variance V, class FrameTag>
	class Rank1Tensor
	{
		VectorN<Scalar, N> _components;

	public:
		using value_type = Scalar;
		using frame_type = FrameTag;
		static constexpr int Dimension = N;
		static constexpr Variance IndexVariance = V;

		Rank1Tensor() = default;
		Rank1Tensor(std::initializer_list<Scalar> values) : _components(values) { }
		explicit Rank1Tensor(const VectorN<Scalar, N>& components) : _components(components) { }

		int size() const noexcept { return N; }

		Scalar& operator[](int i) noexcept { return _components[i]; }
		const Scalar& operator[](int i) const noexcept { return _components[i]; }

		Scalar& at(int i) { return _components.at(i); }
		const Scalar& at(int i) const { return _components.at(i); }

		VectorN<Scalar, N>& components() noexcept { return _components; }
		const VectorN<Scalar, N>& components() const noexcept { return _components; }

		Scalar operator()(const Rank1Tensor<Scalar, N, Variance::Contravariant, FrameTag>& v) const
			requires (V == Variance::Covariant)
		{
			Scalar result{};
			for (int i = 0; i < N; i++)
				result += _components[i] * v[i];
			return result;
		}
	};

	template<int N, class FrameTag>
	using TangentVector = Rank1Tensor<Real, N, Variance::Contravariant, FrameTag>;

	template<int N, class FrameTag>
	using Covector = Rank1Tensor<Real, N, Variance::Covariant, FrameTag>;

	template<class Scalar, int N, Variance V, class FrameTag>
	Rank1Tensor<Scalar, N, V, FrameTag> operator+(const Rank1Tensor<Scalar, N, V, FrameTag>& a,
		const Rank1Tensor<Scalar, N, V, FrameTag>& b)
	{
		Rank1Tensor<Scalar, N, V, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = a[i] + b[i];
		return result;
	}

	template<class Scalar, int N, Variance V, class FrameTag>
	Rank1Tensor<Scalar, N, V, FrameTag> operator-(const Rank1Tensor<Scalar, N, V, FrameTag>& a,
		const Rank1Tensor<Scalar, N, V, FrameTag>& b)
	{
		Rank1Tensor<Scalar, N, V, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = a[i] - b[i];
		return result;
	}

	template<class Scalar, int N, Variance V, class FrameTag>
	Rank1Tensor<Scalar, N, V, FrameTag> operator-(const Rank1Tensor<Scalar, N, V, FrameTag>& a)
	{
		Rank1Tensor<Scalar, N, V, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = -a[i];
		return result;
	}

	template<class Scalar, int N, Variance V, class FrameTag>
	Rank1Tensor<Scalar, N, V, FrameTag> operator*(const Rank1Tensor<Scalar, N, V, FrameTag>& a, Scalar scalar)
	{
		Rank1Tensor<Scalar, N, V, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = a[i] * scalar;
		return result;
	}

	template<class Scalar, int N, Variance V, class FrameTag>
	Rank1Tensor<Scalar, N, V, FrameTag> operator*(Scalar scalar, const Rank1Tensor<Scalar, N, V, FrameTag>& a)
	{
		return a * scalar;
	}

	template<class Scalar, int N, Variance V, class FrameTag>
	Rank1Tensor<Scalar, N, V, FrameTag> operator/(const Rank1Tensor<Scalar, N, V, FrameTag>& a, Scalar scalar)
	{
		Rank1Tensor<Scalar, N, V, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = a[i] / scalar;
		return result;
	}

	template<class Scalar, int N, class FrameTag>
	Scalar Pair(const Rank1Tensor<Scalar, N, Variance::Covariant, FrameTag>& alpha,
		const Rank1Tensor<Scalar, N, Variance::Contravariant, FrameTag>& v)
	{
		Scalar result{};
		for (int i = 0; i < N; i++)
			result += alpha[i] * v[i];
		return result;
	}

	template<int N, class FrameTag>
	TangentVector<N, FrameTag> MakeTangentVector(const VectorN<Real, N>& components)
	{
		return TangentVector<N, FrameTag>(components);
	}

	template<int N, class FrameTag>
	Covector<N, FrameTag> MakeCovector(const VectorN<Real, N>& components)
	{
		return Covector<N, FrameTag>(components);
	}
}

#endif // MML_RANK1_TENSOR_H