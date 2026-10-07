///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/CoordinateMap.h                                ///
///  Description: Typed coordinate maps and frame-aware push-forward/pull-back        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_COORDINATE_MAP_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_COORDINATE_MAP_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/base/DifferentialGeometry/Point.h>
#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/core/CoordTransf/CoordTransfBase.h>

#include <array>

namespace MML
{
	namespace Detail
	{
		template<int N, int K>
		std::array<int, K == 0 ? 1 : K> CoordinateMapUnflattenIndex(int flat)
		{
			std::array<int, K == 0 ? 1 : K> indices{};
			for (int i = K - 1; i >= 0; i--) {
				indices[i] = flat % N;
				flat /= N;
			}
			return indices;
		}
	}

	template<class FromFrame, class ToFrame, int N>
	class CoordinateMap
	{
	public:
		virtual Point<Real, N, ToFrame> map_point(const Point<Real, N, FromFrame>& point) const = 0;
		virtual TangentVector<N, ToFrame> push_forward(const TangentVector<N, FromFrame>& vector,
			const Point<Real, N, FromFrame>& atPoint) const = 0;
		virtual Covector<N, FromFrame> pull_back(const Covector<N, ToFrame>& covector,
			const Point<Real, N, FromFrame>& atPoint) const = 0;
		virtual MatrixNM<Real, N, N> jacobian(const Point<Real, N, FromFrame>& atPoint) const = 0;

		Form1<N, FromFrame> pull_back(const Form1<N, ToFrame>& form,
			const Point<Real, N, FromFrame>& atPoint) const
		{
			return ToForm(pull_back(ToCovector(form), atPoint));
		}

		virtual ~CoordinateMap() { }
	};

	template<class FrameTag, int N>
	class IdentityCoordinateMap : public CoordinateMap<FrameTag, FrameTag, N>
	{
	public:
		Point<Real, N, FrameTag> map_point(const Point<Real, N, FrameTag>& point) const override
		{
			return point;
		}

		TangentVector<N, FrameTag> push_forward(const TangentVector<N, FrameTag>& vector,
			const Point<Real, N, FrameTag>&) const override
		{
			return vector;
		}

		Covector<N, FrameTag> pull_back(const Covector<N, FrameTag>& covector,
			const Point<Real, N, FrameTag>&) const override
		{
			return covector;
		}

		MatrixNM<Real, N, N> jacobian(const Point<Real, N, FrameTag>&) const override
		{
			return MatrixNM<Real, N, N>(REAL(1.0));
		}
	};

	template<class VectorFrom, class VectorTo, int N, class FromFrame, class ToFrame>
	class CoordTransfCoordinateMap : public CoordinateMap<FromFrame, ToFrame, N>
	{
		const CoordTransfWithInverse<VectorFrom, VectorTo, N>& _transform;

	public:
		using CoordinateMap<FromFrame, ToFrame, N>::pull_back;

		explicit CoordTransfCoordinateMap(const CoordTransfWithInverse<VectorFrom, VectorTo, N>& transform)
			: _transform(transform)
		{
		}

		Point<Real, N, ToFrame> map_point(const Point<Real, N, FromFrame>& point) const override
		{
			return Point<Real, N, ToFrame>(_transform.transf(VectorFrom(point.coordinates())));
		}

		TangentVector<N, ToFrame> push_forward(const TangentVector<N, FromFrame>& vector,
			const Point<Real, N, FromFrame>& atPoint) const override
		{
			return TangentVector<N, ToFrame>(_transform.transfVecContravariant(
				VectorFrom(vector.components()), VectorFrom(atPoint.coordinates())));
		}

		Covector<N, FromFrame> pull_back(const Covector<N, ToFrame>& covector,
			const Point<Real, N, FromFrame>& atPoint) const override
		{
			return Covector<N, FromFrame>(_transform.transfInverseVecCovariant(
				VectorTo(covector.components()), VectorFrom(atPoint.coordinates())));
		}

		/// @brief Jacobian convention: row = target coordinate, column = source coordinate.
		/// @details For x' = map(x), this returns J(row i, col j) = partial x'_i / partial x_j.
		MatrixNM<Real, N, N> jacobian(const Point<Real, N, FromFrame>& atPoint) const override
		{
			return _transform.jacobian(atPoint.coordinates());
		}
	};

	using PolarToCartesian2DMap = CoordTransfCoordinateMap<Vector2Polar, Vector2Cartesian, 2, Polar2, Cartesian2>;
	using CartesianToPolar2DMap = CoordTransfCoordinateMap<Vector2Cartesian, Vector2Polar, 2, Cartesian2, Polar2>;
	using CylindricalToCartesian3DMap = CoordTransfCoordinateMap<Vector3Cylindrical, Vector3Cartesian, 3, Cylindrical3, Cartesian3>;
	using CartesianToCylindrical3DMap = CoordTransfCoordinateMap<Vector3Cartesian, Vector3Cylindrical, 3, Cartesian3, Cylindrical3>;
	using SphericalToCartesian3DMap = CoordTransfCoordinateMap<Vector3Spherical, Vector3Cartesian, 3, Spherical3, Cartesian3>;
	using CartesianToSpherical3DMap = CoordTransfCoordinateMap<Vector3Cartesian, Vector3Spherical, 3, Cartesian3, Spherical3>;

	template<class FrameTag, int N>
	IdentityCoordinateMap<FrameTag, N> MakeIdentityCoordinateMap()
	{
		return IdentityCoordinateMap<FrameTag, N>();
	}

	template<class Transform>
	PolarToCartesian2DMap MakePolarToCartesian2DMap(const Transform& transform)
	{
		return PolarToCartesian2DMap(transform);
	}

	template<class Transform>
	CartesianToPolar2DMap MakeCartesianToPolar2DMap(const Transform& transform)
	{
		return CartesianToPolar2DMap(transform);
	}

	template<class Transform>
	CylindricalToCartesian3DMap MakeCylindricalToCartesian3DMap(const Transform& transform)
	{
		return CylindricalToCartesian3DMap(transform);
	}

	template<class Transform>
	CartesianToCylindrical3DMap MakeCartesianToCylindrical3DMap(const Transform& transform)
	{
		return CartesianToCylindrical3DMap(transform);
	}

	template<class Transform>
	SphericalToCartesian3DMap MakeSphericalToCartesian3DMap(const Transform& transform)
	{
		return SphericalToCartesian3DMap(transform);
	}

	template<class Transform>
	CartesianToSpherical3DMap MakeCartesianToSpherical3DMap(const Transform& transform)
	{
		return CartesianToSpherical3DMap(transform);
	}

	template<class FromFrame, class ToFrame, int N>
	Point<Real, N, ToFrame> map_point(const CoordinateMap<FromFrame, ToFrame, N>& map,
		const Point<Real, N, FromFrame>& point)
	{
		return map.map_point(point);
	}

	template<class FromFrame, class ToFrame, int N>
	TangentVector<N, ToFrame> push_forward(const CoordinateMap<FromFrame, ToFrame, N>& map,
		const TangentVector<N, FromFrame>& vector,
		const Point<Real, N, FromFrame>& atPoint)
	{
		return map.push_forward(vector, atPoint);
	}

	template<class FromFrame, class ToFrame, int N>
	Covector<N, FromFrame> pull_back(const CoordinateMap<FromFrame, ToFrame, N>& map,
		const Covector<N, ToFrame>& covector,
		const Point<Real, N, FromFrame>& atPoint)
	{
		return map.pull_back(covector, atPoint);
	}

	template<class FromFrame, class ToFrame, int N, int K>
	DifferentialForm<Real, N, K, FromFrame> pull_back(const CoordinateMap<FromFrame, ToFrame, N>& map,
		const DifferentialForm<Real, N, K, ToFrame>& form,
		const Point<Real, N, FromFrame>& atPoint)
	{
		DifferentialForm<Real, N, K, FromFrame> result;
		MatrixNM<Real, N, N> jac = map.jacobian(atPoint);

		for (int sourceFlat = 0; sourceFlat < result.ComponentCount; sourceFlat++) {
			auto sourceIndices = Detail::CoordinateMapUnflattenIndex<N, K>(sourceFlat);
			Real component{};

			for (int targetFlat = 0; targetFlat < form.ComponentCount; targetFlat++) {
				auto targetIndices = Detail::CoordinateMapUnflattenIndex<N, K>(targetFlat);
				Real term = form.ComponentAt(targetIndices);
				for (int i = 0; i < K; i++)
					term *= jac(targetIndices[i], sourceIndices[i]);
				component += term;
			}

			result.ComponentAt(sourceIndices) = component;
		}

		return result;
	}
}

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_COORDINATE_MAP_H