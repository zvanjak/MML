#include <catch2/catch_all.hpp>

#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/VectorSpaces.h>
#endif

#include <type_traits>

namespace
{
	template<class Left, class Right, class = void>
	struct HasPlus : std::false_type { };

	template<class Left, class Right>
	struct HasPlus<Left, Right, std::void_t<decltype(std::declval<Left>() + std::declval<Right>())>> : std::true_type { };

	template<class MapType>
	concept HasIdentity = requires {
		MapType::Identity();
	};

	template<class MapType>
	concept HasInverse = requires(MapType map) {
		map.inverse();
	};

	template<class MapType, class VectorType>
	concept HasTranslation = requires(VectorType vector) {
		MapType::Translation(vector);
	};
}

using namespace MML;
using namespace MML::VectorSpaces;

namespace MML::Tests::Core::VectorSpacesTests
{
	constexpr Real VectorSpaceTolerance = Testing::Tol(REAL(1e-10), REAL(1e-6));

	TEST_CASE("Basis - Standard basis round trips coordinates", "[VectorSpaces][Basis]")
	{
		Basis<Real, 3> standard = Basis<Real, 3>::Standard();
		VectorN<Real, 3> coords{ REAL(1.0), REAL(2.0), REAL(3.0) };

		VectorN<Real, 3> standardCoords = standard.coordinatesInStandard(coords);
		VectorN<Real, 3> roundTrip = standard.coordinatesFromStandard(standardCoords);

		REQUIRE(standard.determinant() == REAL(1.0));
		REQUIRE(standard.orientationSign() == 1);
		REQUIRE(standardCoords.IsEqualTo(coords));
		REQUIRE(roundTrip.IsEqualTo(coords));
	}

	TEST_CASE("Basis - Scaled basis converts to and from standard coordinates", "[VectorSpaces][Basis]")
	{
		Basis<Real, 3> scaled({
			VectorN<Real, 3>{ REAL(2.0), REAL(0.0), REAL(0.0) },
			VectorN<Real, 3>{ REAL(0.0), REAL(3.0), REAL(0.0) },
			VectorN<Real, 3>{ REAL(0.0), REAL(0.0), REAL(4.0) }
		});

		VectorN<Real, 3> coords{ REAL(1.0), REAL(2.0), REAL(3.0) };
		VectorN<Real, 3> standardCoords = scaled.coordinatesInStandard(coords);
		VectorN<Real, 3> roundTrip = scaled.coordinatesFromStandard(standardCoords);

		REQUIRE(scaled.determinant() == Catch::Approx(REAL(24.0)));
		REQUIRE(scaled.orientationSign() == 1);
		REQUIRE(standardCoords[0] == REAL(2.0));
		REQUIRE(standardCoords[1] == REAL(6.0));
		REQUIRE(standardCoords[2] == REAL(12.0));
		REQUIRE(roundTrip.IsEqualTo(coords, VectorSpaceTolerance));
	}

	TEST_CASE("Basis - ChangeBasis converts between rotated bases", "[VectorSpaces][Basis]")
	{
		const Real invSqrt2 = REAL(1.0) / std::sqrt(REAL(2.0));
		Basis<Real, 2> standard = Basis<Real, 2>::Standard();
		Basis<Real, 2> rotated({
			VectorN<Real, 2>{ invSqrt2, invSqrt2 },
			VectorN<Real, 2>{ -invSqrt2, invSqrt2 }
		});

		VectorN<Real, 2> standardCoords{ REAL(1.0), REAL(0.0) };
		VectorN<Real, 2> rotatedCoords = ChangeBasis(standardCoords, standard, rotated);
		VectorN<Real, 2> back = ChangeBasis(rotatedCoords, rotated, standard);

		REQUIRE(rotated.determinant() == Catch::Approx(REAL(1.0)));
		REQUIRE(rotated.orientationSign() == 1);
		REQUIRE(rotatedCoords[0] == Catch::Approx(invSqrt2));
		REQUIRE(rotatedCoords[1] == Catch::Approx(-invSqrt2));
		REQUIRE(back.IsEqualTo(standardCoords, VectorSpaceTolerance));
	}

	TEST_CASE("Basis - Orientation detects reflected bases", "[VectorSpaces][Basis]")
	{
		Basis<Real, 2> reflected({
			VectorN<Real, 2>{ REAL(1.0), REAL(0.0) },
			VectorN<Real, 2>{ REAL(0.0), -REAL(1.0) }
		});

		REQUIRE(reflected.determinant() == Catch::Approx(-REAL(1.0)));
		REQUIRE(reflected.orientationSign() == -1);
	}

	TEST_CASE("Basis - Dependent basis vectors are rejected", "[VectorSpaces][Basis]")
	{
		REQUIRE_THROWS_AS((Basis<Real, 2>({
			VectorN<Real, 2>{ REAL(1.0), REAL(2.0) },
			VectorN<Real, 2>{ REAL(2.0), REAL(4.0) }
		})), SingularMatrixError);
	}

	TEST_CASE("Subspace - NullSpace and ColumnSpace satisfy rank-nullity", "[VectorSpaces][Subspace]")
	{
		Matrix<Real> A(2, 3);
		A(0, 0) = REAL(1.0); A(0, 1) = REAL(2.0); A(0, 2) = REAL(3.0);
		A(1, 0) = REAL(2.0); A(1, 1) = REAL(4.0); A(1, 2) = REAL(6.0);

		auto kernel = Subspace<Real>::NullSpace(A);
		auto image = Subspace<Real>::ColumnSpace(A);

		REQUIRE(kernel.ambientDimension() == 3);
		REQUIRE(image.ambientDimension() == 2);
		REQUIRE(kernel.dimension() == 2);
		REQUIRE(image.dimension() == 1);
		REQUIRE(kernel.dimension() + image.dimension() == A.cols());
	}

	TEST_CASE("Subspace - Membership projection residual and distance", "[VectorSpaces][Subspace]")
	{
		Matrix<Real> basis(3, 1);
		basis(0, 0) = REAL(1.0);
		basis(1, 0) = REAL(0.0);
		basis(2, 0) = REAL(0.0);
		Subspace<Real> xAxis(3, basis);

		Vector<Real> v{ REAL(2.0), REAL(3.0), REAL(4.0) };
		Vector<Real> projection = xAxis.project(v);
		Vector<Real> residual = xAxis.residual(v);

		REQUIRE(projection[0] == Catch::Approx(REAL(2.0)));
		REQUIRE(projection[1] == Catch::Approx(REAL(0.0)));
		REQUIRE(projection[2] == Catch::Approx(REAL(0.0)));
		REQUIRE(residual[0] == Catch::Approx(REAL(0.0)));
		REQUIRE(residual[1] == Catch::Approx(REAL(3.0)));
		REQUIRE(residual[2] == Catch::Approx(REAL(4.0)));
		REQUIRE(xAxis.distance(v) == Catch::Approx(REAL(5.0)));
		REQUIRE(xAxis.contains(projection));
		REQUIRE_FALSE(xAxis.contains(v, REAL(1e-8)));
	}

	TEST_CASE("Subspace - Orthogonal complement is perpendicular", "[VectorSpaces][Subspace]")
	{
		Matrix<Real> basis(3, 1);
		basis(0, 0) = REAL(1.0);
		Subspace<Real> xAxis(3, basis);

		Subspace<Real> complement = xAxis.orthogonalComplement();

		REQUIRE(complement.ambientDimension() == 3);
		REQUIRE(complement.dimension() == 2);

		for (int j = 0; j < complement.dimension(); j++) {
			Vector<Real> q = complement.orthonormalBasis().VectorFromColumn(j);
			REQUIRE(std::abs(q[0]) < REAL(1e-8));
		}
	}

	TEST_CASE("Subspace - Sum combines spans", "[VectorSpaces][Subspace]")
	{
		Matrix<Real> xBasis(3, 1);
		xBasis(0, 0) = REAL(1.0);
		Matrix<Real> yBasis(3, 1);
		yBasis(1, 0) = REAL(1.0);

		Subspace<Real> xAxis(3, xBasis);
		Subspace<Real> yAxis(3, yBasis);
		Subspace<Real> xyPlane = xAxis.sum(yAxis);

		Vector<Real> inPlane{ REAL(2.0), REAL(3.0), REAL(0.0) };
		Vector<Real> outOfPlane{ REAL(2.0), REAL(3.0), REAL(4.0) };

		REQUIRE(xyPlane.dimension() == 2);
		REQUIRE(xyPlane.contains(inPlane, REAL(1e-8)));
		REQUIRE_FALSE(xyPlane.contains(outOfPlane, REAL(1e-8)));
	}

	TEST_CASE("Subspace - Fundamental subspace factories return expected dimensions", "[VectorSpaces][Subspace]")
	{
		Matrix<Real> A(3, 2);
		A(0, 0) = REAL(1.0); A(0, 1) = REAL(0.0);
		A(1, 0) = REAL(0.0); A(1, 1) = REAL(1.0);
		A(2, 0) = REAL(0.0); A(2, 1) = REAL(0.0);

		auto row = Subspace<Real>::RowSpace(A);
		auto leftNull = Subspace<Real>::LeftNullSpace(A);

		REQUIRE(row.dimension() == 2);
		REQUIRE(leftNull.dimension() == 1);
		REQUIRE(row.ambientDimension() == A.cols());
		REQUIRE(leftNull.ambientDimension() == A.rows());
	}

	TEST_CASE("LinearMap - Application identity zero and inverse", "[VectorSpaces][LinearMap]")
	{
		static_assert(HasIdentity<LinearMap<Real, 2, 2>>);
		static_assert(!HasIdentity<LinearMap<Real, 2, 3>>);
		static_assert(HasInverse<LinearMap<Real, 2, 2>>);
		static_assert(!HasInverse<LinearMap<Real, 2, 3>>);

		MatrixNM<Real, 2, 2> scaling{
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(3.0)
		};
		LinearMap<Real, 2, 2> map(scaling);
		VectorN<Real, 2> x{ REAL(4.0), REAL(5.0) };

		VectorN<Real, 2> y = map(x);
		VectorN<Real, 2> back = map.inverse()(y);
		VectorN<Real, 2> identity = LinearMap<Real, 2, 2>::Identity()(x);
		VectorN<Real, 2> zero = LinearMap<Real, 2, 2>::Zero()(x);

		REQUIRE(y[0] == REAL(8.0));
		REQUIRE(y[1] == REAL(15.0));
		REQUIRE(back.IsEqualTo(x, VectorSpaceTolerance));
		REQUIRE(identity.IsEqualTo(x, VectorSpaceTolerance));
		REQUIRE(zero.isZero());
	}

	TEST_CASE("LinearMap - Composition is dimension-safe and associative for compatible maps", "[VectorSpaces][LinearMap]")
	{
		MatrixNM<Real, 2, 3> projection{
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0)
		};
		MatrixNM<Real, 2, 2> scaling{
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(3.0)
		};
		MatrixNM<Real, 1, 2> sum{
			REAL(1.0), REAL(1.0)
		};

		LinearMap<Real, 3, 2> f(projection);
		LinearMap<Real, 2, 2> g(scaling);
		LinearMap<Real, 2, 1> h(sum);
		VectorN<Real, 3> x{ REAL(4.0), REAL(5.0), REAL(6.0) };

		auto hg_f = Compose(Compose(h, g), f);
		auto h_gf = Compose(h, Compose(g, f));

		REQUIRE(hg_f(x)[0] == Catch::Approx(REAL(23.0)));
		REQUIRE(h_gf(x)[0] == Catch::Approx(hg_f(x)[0]));
	}

	TEST_CASE("LinearMap - Kernel image rank and nullity", "[VectorSpaces][LinearMap]")
	{
		MatrixNM<Real, 2, 3> projection{
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0)
		};
		LinearMap<Real, 3, 2> map(projection);

		Subspace<Real> kernel = map.kernel();
		Subspace<Real> image = map.image();
		Vector<Real> zAxis{ REAL(0.0), REAL(0.0), REAL(7.0) };
		Vector<Real> notKernel{ REAL(1.0), REAL(0.0), REAL(7.0) };

		REQUIRE(map.rank() == 2);
		REQUIRE(map.nullity() == 1);
		REQUIRE(kernel.ambientDimension() == 3);
		REQUIRE(kernel.dimension() == 1);
		REQUIRE(image.ambientDimension() == 2);
		REQUIRE(image.dimension() == 2);
		REQUIRE(kernel.contains(zAxis, REAL(1e-8)));
		REQUIRE_FALSE(kernel.contains(notKernel, REAL(1e-8)));
		REQUIRE(map.rank() + map.nullity() == 3);
	}

	TEST_CASE("LinearMap - matrixInBases changes coordinate representation", "[VectorSpaces][LinearMap]")
	{
		MatrixNM<Real, 2, 2> identity = MatrixNM<Real, 2, 2>::Identity();
		LinearMap<Real, 2, 2> map(identity);
		Basis<Real, 2> domain({
			VectorN<Real, 2>{ REAL(2.0), REAL(0.0) },
			VectorN<Real, 2>{ REAL(0.0), REAL(1.0) }
		});
		Basis<Real, 2> codomain = Basis<Real, 2>::Standard();

		MatrixNM<Real, 2, 2> represented = map.matrixInBases(domain, codomain);

		REQUIRE(represented(0, 0) == Catch::Approx(REAL(2.0)));
		REQUIRE(represented(1, 1) == Catch::Approx(REAL(1.0)));
	}

	TEST_CASE("DualSpace - Dual basis evaluates primal basis as Kronecker delta", "[VectorSpaces][DualSpace]")
	{
		Basis<Real, 2> basis({
			VectorN<Real, 2>{ REAL(2.0), REAL(0.0) },
			VectorN<Real, 2>{ REAL(0.0), REAL(3.0) }
		});
		DualSpace<Real, 2> dualSpace(basis);

		VectorInSpace<Real, 2> e0(VectorN<Real, 2>{ REAL(1.0), REAL(0.0) }, basis);
		VectorInSpace<Real, 2> e1(VectorN<Real, 2>{ REAL(0.0), REAL(1.0) }, basis);
		auto eps0 = dualSpace.dualBasisCovector(0);
		auto eps1 = dualSpace.dualBasisCovector(1);

		REQUIRE(eps0(e0) == Catch::Approx(REAL(1.0)));
		REQUIRE(eps0(e1) == Catch::Approx(REAL(0.0)));
		REQUIRE(eps1(e0) == Catch::Approx(REAL(0.0)));
		REQUIRE(eps1(e1) == Catch::Approx(REAL(1.0)));
	}

	TEST_CASE("DualSpace - Pullback satisfies alpha(f(v)) identity", "[VectorSpaces][DualSpace]")
	{
		MatrixNM<Real, 2, 3> projection{
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0)
		};
		LinearMap<Real, 3, 2> map(projection);
		DualSpace<Real, 2> codomainDual;
		auto alpha = codomainDual.covector(VectorN<Real, 2>{ REAL(4.0), REAL(5.0) });
		VectorInSpace<Real, 3> v(VectorN<Real, 3>{ REAL(1.0), REAL(2.0), REAL(3.0) });
		VectorInSpace<Real, 2> fv(map(v.coordinates()));

		auto pulled = codomainDual.pullback(map, alpha);

		REQUIRE(pulled(v) == Catch::Approx(alpha(fv)));
		REQUIRE(pulled.coordinates()[0] == Catch::Approx(REAL(4.0)));
		REQUIRE(pulled.coordinates()[1] == Catch::Approx(REAL(5.0)));
		REQUIRE(pulled.coordinates()[2] == Catch::Approx(REAL(0.0)));
	}

	TEST_CASE("DualSpace - Annihilator covectors vanish on subspace basis", "[VectorSpaces][DualSpace]")
	{
		Matrix<Real> basis(3, 2);
		basis(0, 0) = REAL(1.0);
		basis(1, 1) = REAL(1.0);
		Subspace<Real> xyPlane(3, basis);
		DualSpace<Real, 3> dualSpace;

		Subspace<Real> annihilator = dualSpace.annihilator(xyPlane);

		REQUIRE(annihilator.ambientDimension() == 3);
		REQUIRE(annihilator.dimension() == 1);
		for (int j = 0; j < annihilator.dimension(); j++) {
			Vector<Real> alpha = annihilator.orthonormalBasis().VectorFromColumn(j);
			REQUIRE(std::abs(alpha[0]) < REAL(1e-8));
			REQUIRE(std::abs(alpha[1]) < REAL(1e-8));
			REQUIRE(std::abs(std::abs(alpha[2]) - REAL(1.0)) < REAL(1e-8));
		}
	}

	TEST_CASE("DualSpace - Conversion to and from typed Covector", "[VectorSpaces][DualSpace]")
	{
		Basis<Real, 2> scaled({
			VectorN<Real, 2>{ REAL(2.0), REAL(0.0) },
			VectorN<Real, 2>{ REAL(0.0), REAL(3.0) }
		});
		DualSpace<Real, 2> dualSpace(scaled);
		Covector<2, Cartesian2> standardCovector{ REAL(6.0), REAL(12.0) };

		CovectorInSpace<Real, 2> inDualBasis = dualSpace.fromCovector(standardCovector);
		Covector<2, Cartesian2> roundTrip = inDualBasis.toCovector<Cartesian2>();

		REQUIRE(inDualBasis.coordinates()[0] == Catch::Approx(REAL(12.0)));
		REQUIRE(inDualBasis.coordinates()[1] == Catch::Approx(REAL(36.0)));
		REQUIRE(roundTrip.components().IsEqualTo(standardCovector.components(), VectorSpaceTolerance));
	}

	TEST_CASE("InnerProductSpace - Weighted inner norm and distance match Gram matrix", "[VectorSpaces][InnerProductSpace]")
	{
		MatrixNM<Real, 2, 2> gram{
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(3.0)
		};
		InnerProductSpace<Real, 2> space(gram);
		VectorN<Real, 2> a{ REAL(1.0), REAL(2.0) };
		VectorN<Real, 2> b{ REAL(3.0), REAL(4.0) };

		REQUIRE(space.isSymmetric());
		REQUIRE(space.isPositiveDefinite());
		REQUIRE(space.inner(a, b) == Catch::Approx(REAL(30.0)));
		REQUIRE(space.norm(a) == Catch::Approx(std::sqrt(REAL(14.0))));
		REQUIRE(space.distance(a, b) == Catch::Approx(std::sqrt(REAL(20.0))));
	}

	TEST_CASE("InnerProductSpace - Gram matrix transforms into scaled basis", "[VectorSpaces][InnerProductSpace]")
	{
		MatrixNM<Real, 2, 2> gram{
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(3.0)
		};
		InnerProductSpace<Real, 2> space(gram);
		Basis<Real, 2> scaled({
			VectorN<Real, 2>{ REAL(2.0), REAL(0.0) },
			VectorN<Real, 2>{ REAL(0.0), REAL(1.0) }
		});

		MatrixNM<Real, 2, 2> represented = space.gramMatrixInBasis(scaled);

		REQUIRE(represented(0, 0) == Catch::Approx(REAL(8.0)));
		REQUIRE(represented(0, 1) == Catch::Approx(REAL(0.0)));
		REQUIRE(represented(1, 0) == Catch::Approx(REAL(0.0)));
		REQUIRE(represented(1, 1) == Catch::Approx(REAL(3.0)));
	}

	TEST_CASE("InnerProductSpace - Gram-Schmidt produces orthonormal basis", "[VectorSpaces][InnerProductSpace]")
	{
		MatrixNM<Real, 2, 2> gram{
			REAL(4.0), REAL(0.0),
			REAL(0.0), REAL(9.0)
		};
		InnerProductSpace<Real, 2> space(gram);
		Basis<Real, 2> standard = Basis<Real, 2>::Standard();

		Basis<Real, 2> orthonormal = space.orthonormalize(standard);
		MatrixNM<Real, 2, 2> represented = space.gramMatrixInBasis(orthonormal);

		REQUIRE(orthonormal.matrix()(0, 0) == Catch::Approx(REAL(0.5)));
		REQUIRE(orthonormal.matrix()(1, 1) == Catch::Approx(REAL(1.0) / REAL(3.0)));
		REQUIRE(represented(0, 0) == Catch::Approx(REAL(1.0)));
		REQUIRE(represented(0, 1) == Catch::Approx(REAL(0.0)).margin(VectorSpaceTolerance));
		REQUIRE(represented(1, 0) == Catch::Approx(REAL(0.0)).margin(VectorSpaceTolerance));
		REQUIRE(represented(1, 1) == Catch::Approx(REAL(1.0)));
	}

	TEST_CASE("InnerProductSpace - Projection residual is orthogonal under Gram matrix", "[VectorSpaces][InnerProductSpace]")
	{
		MatrixNM<Real, 2, 2> gram{
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(1.0)
		};
		InnerProductSpace<Real, 2> space(gram);
		const Real invSqrt2 = REAL(1.0) / std::sqrt(REAL(2.0));
		Matrix<Real> basis(2, 1);
		basis(0, 0) = invSqrt2;
		basis(1, 0) = invSqrt2;
		Subspace<Real> diagonal(2, basis);
		Vector<Real> vector{ REAL(1.0), REAL(0.0) };

		Vector<Real> projection = space.orthogonalProjection(diagonal, vector);
		Vector<Real> residual = vector - projection;
		VectorN<Real, 2> basisVector{ invSqrt2, invSqrt2 };
		VectorN<Real, 2> residualFixed{ residual[0], residual[1] };

		REQUIRE(projection[0] == Catch::Approx(REAL(2.0) / REAL(3.0)));
		REQUIRE(projection[1] == Catch::Approx(REAL(2.0) / REAL(3.0)));
		REQUIRE(space.inner(basisVector, residualFixed) == Catch::Approx(REAL(0.0)).margin(VectorSpaceTolerance));
	}

	TEST_CASE("InnerProductSpace - Orthogonal complement uses Gram orthogonality", "[VectorSpaces][InnerProductSpace]")
	{
		MatrixNM<Real, 2, 2> gram{
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(1.0)
		};
		InnerProductSpace<Real, 2> space(gram);
		const Real invSqrt2 = REAL(1.0) / std::sqrt(REAL(2.0));
		Matrix<Real> basis(2, 1);
		basis(0, 0) = invSqrt2;
		basis(1, 0) = invSqrt2;
		Subspace<Real> diagonal(2, basis);

		Subspace<Real> complement = space.orthogonalComplement(diagonal);

		REQUIRE(complement.ambientDimension() == 2);
		REQUIRE(complement.dimension() == 1);
		Vector<Real> complementVector = complement.orthonormalBasis().VectorFromColumn(0);
		VectorN<Real, 2> u{ invSqrt2, invSqrt2 };
		VectorN<Real, 2> w{ complementVector[0], complementVector[1] };
		REQUIRE(space.inner(u, w) == Catch::Approx(REAL(0.0)).margin(VectorSpaceTolerance));
	}

	TEST_CASE("InnerProductSpace - Adjoint satisfies defining identity", "[VectorSpaces][InnerProductSpace]")
	{
		MatrixNM<Real, 2, 2> domainGram{
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(1.0)
		};
		MatrixNM<Real, 2, 2> codomainGram{
			REAL(3.0), REAL(0.0),
			REAL(0.0), REAL(4.0)
		};
		InnerProductSpace<Real, 2> domain(domainGram);
		InnerProductSpace<Real, 2> codomain(codomainGram);
		MatrixNM<Real, 2, 2> matrix{
			REAL(1.0), REAL(2.0),
			REAL(3.0), REAL(4.0)
		};
		LinearMap<Real, 2, 2> map(matrix);
		VectorN<Real, 2> x{ REAL(5.0), REAL(6.0) };
		VectorN<Real, 2> y{ REAL(7.0), REAL(8.0) };

		auto mapAdjoint = adjoint(map, domain, codomain);

		REQUIRE(codomain.inner(map(x), y) == Catch::Approx(domain.inner(x, mapAdjoint(y))));
	}

	TEST_CASE("InnerProductSpace - Indefinite Gram matrices are rejected", "[VectorSpaces][InnerProductSpace]")
	{
		MatrixNM<Real, 2, 2> indefinite{
			REAL(1.0), REAL(0.0),
			REAL(0.0), -REAL(1.0)
		};

		REQUIRE_THROWS_AS((InnerProductSpace<Real, 2>(indefinite)), MatrixNumericalError);
	}

	TEST_CASE("AffineSpace - Point arithmetic distinguishes points and vectors", "[VectorSpaces][AffineSpace]")
	{
		using Point2 = PointInSpace<Real, 2>;
		static_assert(!HasPlus<Point2, Point2>::value, "point + point must not be a valid affine operation");

		Point2 p(VectorN<Real, 2>{ REAL(1.0), REAL(2.0) });
		Point2 q(VectorN<Real, 2>{ REAL(4.0), REAL(6.0) });
		VectorN<Real, 2> displacement = q - p;
		Point2 moved = p + displacement;
		Point2 movedBack = moved - displacement;

		REQUIRE(displacement[0] == Catch::Approx(REAL(3.0)));
		REQUIRE(displacement[1] == Catch::Approx(REAL(4.0)));
		REQUIRE(moved.coordinates().IsEqualTo(q.coordinates(), VectorSpaceTolerance));
		REQUIRE(movedBack.coordinates().IsEqualTo(p.coordinates(), VectorSpaceTolerance));
	}

	TEST_CASE("AffineSpace - Frames convert point and vector coordinates", "[VectorSpaces][AffineSpace]")
	{
		PointInSpace<Real, 2> origin(VectorN<Real, 2>{ REAL(10.0), REAL(20.0) });
		Basis<Real, 2> basis({
			VectorN<Real, 2>{ REAL(2.0), REAL(0.0) },
			VectorN<Real, 2>{ REAL(0.0), REAL(3.0) }
		});
		AffineFrame<Real, 2> frame(origin, basis);
		VectorN<Real, 2> frameCoords{ REAL(4.0), REAL(5.0) };

		PointInSpace<Real, 2> point = frame.pointFromCoordinates(frameCoords);
		VectorN<Real, 2> roundTrip = frame.coordinatesOf(point);
		VectorN<Real, 2> vector = frame.vectorFromCoordinates(frameCoords);
		VectorN<Real, 2> vectorRoundTrip = frame.coordinatesOfVector(vector);

		REQUIRE(point.coordinates()[0] == Catch::Approx(REAL(18.0)));
		REQUIRE(point.coordinates()[1] == Catch::Approx(REAL(35.0)));
		REQUIRE(roundTrip.IsEqualTo(frameCoords, VectorSpaceTolerance));
		REQUIRE(vector[0] == Catch::Approx(REAL(8.0)));
		REQUIRE(vector[1] == Catch::Approx(REAL(15.0)));
		REQUIRE(vectorRoundTrip.IsEqualTo(frameCoords, VectorSpaceTolerance));
	}

	TEST_CASE("AffineSpace - Affine maps compose and invert", "[VectorSpaces][AffineSpace]")
	{
		static_assert(HasIdentity<AffineMap<Real, 2, 2>>);
		static_assert(!HasIdentity<AffineMap<Real, 2, 3>>);
		static_assert(HasTranslation<AffineMap<Real, 2, 2>, VectorN<Real, 2>>);
		static_assert(!HasTranslation<AffineMap<Real, 2, 3>, VectorN<Real, 3>>);
		static_assert(HasInverse<AffineMap<Real, 2, 2>>);
		static_assert(!HasInverse<AffineMap<Real, 2, 3>>);

		const Real angle = REAL(0.5) * Constants::PI;
		MatrixNM<Real, 2, 2> rotation{
			std::cos(angle), -std::sin(angle),
			std::sin(angle), std::cos(angle)
		};
		AffineMap<Real, 2, 2> rigid(LinearMap<Real, 2, 2>(rotation), VectorN<Real, 2>{ REAL(3.0), REAL(4.0) });
		AffineMap<Real, 2, 2> translate = AffineMap<Real, 2, 2>::Translation(VectorN<Real, 2>{ REAL(5.0), REAL(6.0) });
		PointInSpace<Real, 2> point(VectorN<Real, 2>{ REAL(1.0), REAL(2.0) });

		auto composed = Compose(translate, rigid);
		PointInSpace<Real, 2> sequential = translate(rigid(point));
		PointInSpace<Real, 2> direct = composed(point);
		PointInSpace<Real, 2> recovered = rigid.inverse()(rigid(point));

		REQUIRE(direct.coordinates().IsEqualTo(sequential.coordinates(), VectorSpaceTolerance));
		REQUIRE(recovered.coordinates().IsEqualTo(point.coordinates(), VectorSpaceTolerance));
		REQUIRE(rigid.applyVector(VectorN<Real, 2>{ REAL(1.0), REAL(0.0) })[0] == Catch::Approx(REAL(0.0)).margin(VectorSpaceTolerance));
		REQUIRE(rigid.applyVector(VectorN<Real, 2>{ REAL(1.0), REAL(0.0) })[1] == Catch::Approx(REAL(1.0)));
	}

	TEST_CASE("AffineSpace - Homogeneous helpers round trip points and vectors", "[VectorSpaces][AffineSpace]")
	{
		PointInSpace<Real, 2> point(VectorN<Real, 2>{ REAL(7.0), REAL(8.0) });
		VectorN<Real, 2> vector{ REAL(3.0), REAL(4.0) };
		AffineMap<Real, 2, 2> translate = AffineMap<Real, 2, 2>::Translation(VectorN<Real, 2>{ REAL(10.0), REAL(20.0) });

		VectorN<Real, 3> pointH = HomogeneousPoint(point);
		VectorN<Real, 3> vectorH = HomogeneousVector(vector);
		PointInSpace<Real, 2> pointRoundTrip = PointFromHomogeneous(pointH);
		VectorN<Real, 2> vectorRoundTrip = VectorFromHomogeneous(vectorH);
		VectorN<Real, 3> mappedH = translate.applyHomogeneous(pointH);

		REQUIRE(pointH[2] == Catch::Approx(REAL(1.0)));
		REQUIRE(vectorH[2] == Catch::Approx(REAL(0.0)));
		REQUIRE(pointRoundTrip.coordinates().IsEqualTo(point.coordinates(), VectorSpaceTolerance));
		REQUIRE(vectorRoundTrip.IsEqualTo(vector, VectorSpaceTolerance));
		REQUIRE(mappedH[0] == Catch::Approx(REAL(17.0)));
		REQUIRE(mappedH[1] == Catch::Approx(REAL(28.0)));
		REQUIRE(mappedH[2] == Catch::Approx(REAL(1.0)));
	}

	TEST_CASE("AffineSpace - PointInSpace bridges to typed Point", "[VectorSpaces][AffineSpace]")
	{
		PointInSpace<Real, 2> point(VectorN<Real, 2>{ REAL(2.0), REAL(5.0) });
		Point<Real, 2, Cartesian2> typed = point.toPoint<Cartesian2>();
		PointInSpace<Real, 2> roundTrip = PointInSpace<Real, 2>::fromPoint(typed);

		REQUIRE(typed.coordinates().IsEqualTo(point.coordinates(), VectorSpaceTolerance));
		REQUIRE(roundTrip.coordinates().IsEqualTo(point.coordinates(), VectorSpaceTolerance));
	}
}