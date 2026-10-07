///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/Point.h                                        ///
///  Description: Typed point wrapper for coordinate-frame aware geometry             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DIFFERENTIAL_GEOMETRY_POINT_H
#define MML_DIFFERENTIAL_GEOMETRY_POINT_H

#include <mml/MMLBase.h>

#include <mml/base/Tensor/Rank1Tensor.h>
#include <mml/base/Vector/VectorN.h>

#include <initializer_list>

namespace MML
{
	template<class Scalar, int N, class FrameTag>
	class Point
	{
		VectorN<Scalar, N> _coordinates;

	public:
		using value_type = Scalar;
		using frame_type = FrameTag;
		static constexpr int Dimension = N;

		Point() = default;
		Point(std::initializer_list<Scalar> values) : _coordinates(values) { }
		explicit Point(const VectorN<Scalar, N>& coordinates) : _coordinates(coordinates) { }

		int size() const noexcept { return N; }

		Scalar& operator[](int i) noexcept { return _coordinates[i]; }
		const Scalar& operator[](int i) const noexcept { return _coordinates[i]; }

		Scalar& at(int i) { return _coordinates.at(i); }
		const Scalar& at(int i) const { return _coordinates.at(i); }

		VectorN<Scalar, N>& coordinates() noexcept { return _coordinates; }
		const VectorN<Scalar, N>& coordinates() const noexcept { return _coordinates; }
	};

	template<class Scalar, int N, class FrameTag>
	Point<Scalar, N, FrameTag> operator+(const Point<Scalar, N, FrameTag>& point,
		const Rank1Tensor<Scalar, N, Variance::Contravariant, FrameTag>& displacement)
	{
		Point<Scalar, N, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = point[i] + displacement[i];
		return result;
	}

	template<class Scalar, int N, class FrameTag>
	Point<Scalar, N, FrameTag> operator+(const Rank1Tensor<Scalar, N, Variance::Contravariant, FrameTag>& displacement,
		const Point<Scalar, N, FrameTag>& point)
	{
		return point + displacement;
	}

	template<class Scalar, int N, class FrameTag>
	Point<Scalar, N, FrameTag> operator-(const Point<Scalar, N, FrameTag>& point,
		const Rank1Tensor<Scalar, N, Variance::Contravariant, FrameTag>& displacement)
	{
		Point<Scalar, N, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = point[i] - displacement[i];
		return result;
	}

	template<class Scalar, int N, class FrameTag>
	Rank1Tensor<Scalar, N, Variance::Contravariant, FrameTag> operator-(const Point<Scalar, N, FrameTag>& a,
		const Point<Scalar, N, FrameTag>& b)
	{
		Rank1Tensor<Scalar, N, Variance::Contravariant, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = a[i] - b[i];
		return result;
	}
}

#endif // MML_DIFFERENTIAL_GEOMETRY_POINT_H