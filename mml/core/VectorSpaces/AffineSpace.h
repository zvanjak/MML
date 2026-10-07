///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorSpaces/AffineSpace.h                                          ///
///  Description: Fixed-size affine spaces, frames, points, and affine maps           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_SPACES_AFFINE_SPACE_H
#define MML_VECTOR_SPACES_AFFINE_SPACE_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/Point.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/core/VectorSpaces/Basis.h>
#include <mml/core/VectorSpaces/LinearMap.h>

#include <string>
#include <type_traits>

namespace MML::VectorSpaces
{
	template<class Scalar, int N>
	class PointInSpace
	{
		VectorN<Scalar, N> _coordinates;

	public:
		using scalar_type = Scalar;
		static constexpr int Dimension = N;

		PointInSpace() = default;
		explicit PointInSpace(const VectorN<Scalar, N>& standardCoordinates)
			: _coordinates(standardCoordinates)
		{
		}

		const VectorN<Scalar, N>& coordinates() const noexcept { return _coordinates; }
		VectorN<Scalar, N>& coordinates() noexcept { return _coordinates; }

		Scalar& operator[](int i) noexcept { return _coordinates[i]; }
		const Scalar& operator[](int i) const noexcept { return _coordinates[i]; }

		template<class FrameTag>
		Point<Scalar, N, FrameTag> toPoint() const
		{
			return Point<Scalar, N, FrameTag>(_coordinates);
		}

		template<class FrameTag>
		static PointInSpace fromPoint(const Point<Scalar, N, FrameTag>& point)
		{
			return PointInSpace(point.coordinates());
		}
	};

	template<class Scalar, int N>
	VectorN<Scalar, N> operator-(const PointInSpace<Scalar, N>& a, const PointInSpace<Scalar, N>& b)
	{
		return a.coordinates() - b.coordinates();
	}

	template<class Scalar, int N>
	PointInSpace<Scalar, N> operator+(const PointInSpace<Scalar, N>& point, const VectorN<Scalar, N>& displacement)
	{
		return PointInSpace<Scalar, N>(point.coordinates() + displacement);
	}

	template<class Scalar, int N>
	PointInSpace<Scalar, N> operator+(const VectorN<Scalar, N>& displacement, const PointInSpace<Scalar, N>& point)
	{
		return point + displacement;
	}

	template<class Scalar, int N>
	PointInSpace<Scalar, N> operator-(const PointInSpace<Scalar, N>& point, const VectorN<Scalar, N>& displacement)
	{
		return PointInSpace<Scalar, N>(point.coordinates() - displacement);
	}

	template<class Scalar, int N>
	VectorN<Scalar, N + 1> HomogeneousPoint(const PointInSpace<Scalar, N>& point)
	{
		VectorN<Scalar, N + 1> result;
		for (int i = 0; i < N; i++)
			result[i] = point.coordinates()[i];
		result[N] = Scalar{ 1 };
		return result;
	}

	template<class Scalar, int N>
	VectorN<Scalar, N + 1> HomogeneousVector(const VectorN<Scalar, N>& vector)
	{
		VectorN<Scalar, N + 1> result;
		for (int i = 0; i < N; i++)
			result[i] = vector[i];
		result[N] = Scalar{};
		return result;
	}

	template<class Scalar, int HomogeneousN>
	PointInSpace<Scalar, HomogeneousN - 1> PointFromHomogeneous(const VectorN<Scalar, HomogeneousN>& coordinates,
		Real tolerance = Defaults::VectorIsEqualTolerance) requires (HomogeneousN > 1)
	{
		constexpr int N = HomogeneousN - 1;
		if (std::abs(coordinates[N]) <= Scalar(tolerance))
			throw VectorDimensionError("PointFromHomogeneous - point weight must be non-zero", HomogeneousN, N);

		VectorN<Scalar, N> standard;
		for (int i = 0; i < N; i++)
			standard[i] = coordinates[i] / coordinates[N];
		return PointInSpace<Scalar, N>(standard);
	}

	template<class Scalar, int HomogeneousN>
	VectorN<Scalar, HomogeneousN - 1> VectorFromHomogeneous(const VectorN<Scalar, HomogeneousN>& coordinates,
		Real tolerance = Defaults::VectorIsEqualTolerance) requires (HomogeneousN > 1)
	{
		constexpr int N = HomogeneousN - 1;
		if (std::abs(coordinates[N]) > Scalar(tolerance))
			throw VectorDimensionError("VectorFromHomogeneous - vector weight must be zero", HomogeneousN, N);

		VectorN<Scalar, N> vector;
		for (int i = 0; i < N; i++)
			vector[i] = coordinates[i];
		return vector;
	}

	template<class Scalar, int N>
	class AffineFrame
	{
		PointInSpace<Scalar, N> _origin;
		Basis<Scalar, N> _basis;

	public:
		using scalar_type = Scalar;
		static constexpr int Dimension = N;

		AffineFrame() : _origin(), _basis(Basis<Scalar, N>::Standard()) { }

		AffineFrame(const PointInSpace<Scalar, N>& origin, const Basis<Scalar, N>& basis)
			: _origin(origin), _basis(basis)
		{
		}

		static AffineFrame Standard()
		{
			return AffineFrame();
		}

		const PointInSpace<Scalar, N>& origin() const noexcept { return _origin; }
		const Basis<Scalar, N>& basis() const noexcept { return _basis; }

		PointInSpace<Scalar, N> pointFromCoordinates(const VectorN<Scalar, N>& frameCoordinates) const
		{
			return _origin + _basis.coordinatesInStandard(frameCoordinates);
		}

		VectorN<Scalar, N> coordinatesOf(const PointInSpace<Scalar, N>& point) const
		{
			return _basis.coordinatesFromStandard(point - _origin);
		}

		VectorN<Scalar, N> vectorFromCoordinates(const VectorN<Scalar, N>& frameCoordinates) const
		{
			return _basis.coordinatesInStandard(frameCoordinates);
		}

		VectorN<Scalar, N> coordinatesOfVector(const VectorN<Scalar, N>& vector) const
		{
			return _basis.coordinatesFromStandard(vector);
		}
	};

	template<class Scalar, int N>
	class AffineSpace
	{
		std::string _name;

	public:
		using scalar_type = Scalar;
		static constexpr int Dimension = N;

		explicit AffineSpace(std::string name = {})
			: _name(std::move(name))
		{
		}

		const std::string& name() const noexcept { return _name; }

		PointInSpace<Scalar, N> origin() const
		{
			return PointInSpace<Scalar, N>();
		}

		PointInSpace<Scalar, N> point(const VectorN<Scalar, N>& standardCoordinates) const
		{
			return PointInSpace<Scalar, N>(standardCoordinates);
		}

		VectorN<Scalar, N> vector(const VectorN<Scalar, N>& standardCoordinates) const
		{
			return standardCoordinates;
		}

		AffineFrame<Scalar, N> standardFrame() const
		{
			return AffineFrame<Scalar, N>::Standard();
		}
	};

	template<class Scalar, int DomainN, int CodomainN>
	class AffineMap
	{
		LinearMap<Scalar, DomainN, CodomainN> _linear;
		VectorN<Scalar, CodomainN> _translation;

	public:
		using scalar_type = Scalar;
		static constexpr int DomainDimension = DomainN;
		static constexpr int CodomainDimension = CodomainN;

		AffineMap()
			: _linear(LinearMap<Scalar, DomainN, CodomainN>::Zero()), _translation()
		{
		}

		AffineMap(const LinearMap<Scalar, DomainN, CodomainN>& linearPart,
			const VectorN<Scalar, CodomainN>& translation)
			: _linear(linearPart), _translation(translation)
		{
		}

		static AffineMap Translation(const VectorN<Scalar, CodomainN>& translation) requires (DomainN == CodomainN)
		{
			return AffineMap(LinearMap<Scalar, DomainN, CodomainN>::Identity(), translation);
		}

		static AffineMap Identity() requires (DomainN == CodomainN)
		{
			return AffineMap(LinearMap<Scalar, DomainN, CodomainN>::Identity(), VectorN<Scalar, CodomainN>());
		}

		const LinearMap<Scalar, DomainN, CodomainN>& linearPart() const noexcept { return _linear; }
		const VectorN<Scalar, CodomainN>& translation() const noexcept { return _translation; }

		PointInSpace<Scalar, CodomainN> apply(const PointInSpace<Scalar, DomainN>& point) const
		{
			return PointInSpace<Scalar, CodomainN>(_linear(point.coordinates()) + _translation);
		}

		VectorN<Scalar, CodomainN> applyVector(const VectorN<Scalar, DomainN>& vector) const
		{
			return _linear(vector);
		}

		PointInSpace<Scalar, CodomainN> operator()(const PointInSpace<Scalar, DomainN>& point) const
		{
			return apply(point);
		}

		MatrixNM<Scalar, CodomainN + 1, DomainN + 1> homogeneousMatrix() const
		{
			MatrixNM<Scalar, CodomainN + 1, DomainN + 1> result;
			for (int i = 0; i < CodomainN; i++) {
				for (int j = 0; j < DomainN; j++)
					result(i, j) = _linear.matrix()(i, j);
				result(i, DomainN) = _translation[i];
			}
			result(CodomainN, DomainN) = Scalar{ 1 };
			return result;
		}

		VectorN<Scalar, CodomainN + 1> applyHomogeneous(const VectorN<Scalar, DomainN + 1>& coordinates) const
		{
			return homogeneousMatrix() * coordinates;
		}

		AffineMap inverse() const requires (DomainN == CodomainN)
		{
			LinearMap<Scalar, DomainN, CodomainN> inverseLinear = _linear.inverse();
			return AffineMap(inverseLinear, -(inverseLinear(_translation)));
		}
	};

	template<class Scalar, int DomainN, int MiddleN, int CodomainN>
	AffineMap<Scalar, DomainN, CodomainN> Compose(const AffineMap<Scalar, MiddleN, CodomainN>& g,
		const AffineMap<Scalar, DomainN, MiddleN>& f)
	{
		return AffineMap<Scalar, DomainN, CodomainN>(
			Compose(g.linearPart(), f.linearPart()),
			g.applyVector(f.translation()) + g.translation());
	}

	template<class Scalar, int DomainN, int MiddleN, int CodomainN>
	AffineMap<Scalar, DomainN, CodomainN> compose(const AffineMap<Scalar, MiddleN, CodomainN>& g,
		const AffineMap<Scalar, DomainN, MiddleN>& f)
	{
		return Compose(g, f);
	}
}

#endif // MML_VECTOR_SPACES_AFFINE_SPACE_H