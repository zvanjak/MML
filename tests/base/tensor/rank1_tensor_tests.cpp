#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/DifferentialGeometry/Point.h>
#include <mml/base/Tensor/Rank1Tensor.h>
#endif

#include <type_traits>
#include <utility>

using namespace MML;

namespace MML::Tests::Base::Rank1TensorTests
{
	template<class, class, class = void>
	struct HasAddition : std::false_type { };

	template<class A, class B>
	struct HasAddition<A, B, std::void_t<decltype(std::declval<A>() + std::declval<B>())>> : std::true_type { };

	template<class, class, class = void>
	struct HasCall : std::false_type { };

	template<class A, class B>
	struct HasCall<A, B, std::void_t<decltype(std::declval<A>()(std::declval<B>()))>> : std::true_type { };

	template<class, class, class = void>
	struct HasPointAddition : std::false_type { };

	template<class A, class B>
	struct HasPointAddition<A, B, std::void_t<decltype(std::declval<A>() + std::declval<B>())>> : std::true_type { };

	TEST_CASE("Rank1Tensor - Construction metadata and component access", "[Rank1Tensor][DifferentialForms]")
	{
		using Vector = TangentVector<3, Cartesian3>;
		using Alpha = Covector<3, Cartesian3>;

		static_assert(Vector::Dimension == 3, "TangentVector dimension metadata is fixed at compile time");
		static_assert(Vector::IndexVariance == Variance::Contravariant, "TangentVector is contravariant");
		static_assert(Alpha::IndexVariance == Variance::Covariant, "Covector is covariant");
		static_assert(std::is_same<typename Vector::frame_type, Cartesian3>::value, "Frame tag is part of the type");

		Vector v{ REAL(1.0), REAL(2.0), REAL(3.0) };
		REQUIRE(v.size() == 3);
		REQUIRE(v[0] == REAL(1.0));
		REQUIRE(v[1] == REAL(2.0));
		REQUIRE(v[2] == REAL(3.0));

		VectorN<Real, 3> raw{ REAL(4.0), REAL(5.0), REAL(6.0) };
		Vector fromRaw(raw);
		REQUIRE(fromRaw.components() == raw);

		fromRaw.components()[0] = REAL(7.0);
		REQUIRE(fromRaw[0] == REAL(7.0));
	}

	TEST_CASE("Rank1Tensor - Same category arithmetic and covector pairing", "[Rank1Tensor][DifferentialForms]")
	{
		TangentVector<3, Cartesian3> a{ REAL(1.0), REAL(2.0), REAL(3.0) };
		TangentVector<3, Cartesian3> b{ REAL(4.0), -REAL(1.0), REAL(0.5) };
		Covector<3, Cartesian3> alpha{ REAL(2.0), -REAL(3.0), REAL(4.0) };

		auto sum = a + b;
		auto diff = a - b;
		auto scaled = a * REAL(2.0);
		auto divided = scaled / REAL(2.0);

		REQUIRE(sum[0] == REAL(5.0));
		REQUIRE(sum[1] == REAL(1.0));
		REQUIRE(sum[2] == REAL(3.5));
		REQUIRE(diff[0] == -REAL(3.0));
		REQUIRE(diff[1] == REAL(3.0));
		REQUIRE(diff[2] == REAL(2.5));
		REQUIRE(divided.components() == a.components());

		Real expectedPairing = REAL(2.0) * REAL(1.0) - REAL(3.0) * REAL(2.0) + REAL(4.0) * REAL(3.0);
		REQUIRE(Pair(alpha, a) == expectedPairing);
		REQUIRE(alpha(a) == expectedPairing);
	}

	TEST_CASE("Rank1Tensor - Invalid category operations are not available", "[Rank1Tensor][DifferentialForms]")
	{
		using CartesianVector = TangentVector<3, Cartesian3>;
		using SphericalVector = TangentVector<3, Spherical3>;
		using CartesianCovector = Covector<3, Cartesian3>;

		static_assert(HasAddition<CartesianVector, CartesianVector>::value, "Same frame tangent vectors can be added");
		static_assert(!HasAddition<CartesianVector, CartesianCovector>::value, "Tangent vectors and covectors cannot be added");
		static_assert(!HasAddition<CartesianVector, SphericalVector>::value, "Different frames cannot be added");
		static_assert(HasCall<CartesianCovector, CartesianVector>::value, "Covectors can evaluate tangent vectors in the same frame");
		static_assert(!HasCall<CartesianCovector, CartesianCovector>::value, "Covectors cannot evaluate covectors");
		static_assert(!HasCall<CartesianVector, CartesianCovector>::value, "Tangent vectors are not linear functionals");

		SUCCEED("Invalid category operations are rejected by the overload set.");
	}

	TEST_CASE("Point - Typed coordinates and displacement operations", "[Point][DifferentialForms]")
	{
		Point<Real, 3, Cartesian3> p{ REAL(1.0), REAL(2.0), REAL(3.0) };
		TangentVector<3, Cartesian3> dx{ REAL(0.5), -REAL(1.0), REAL(2.0) };

		auto q = p + dx;
		auto r = dx + p;
		auto back = q - dx;
		auto displacement = q - p;

		REQUIRE(q[0] == REAL(1.5));
		REQUIRE(q[1] == REAL(1.0));
		REQUIRE(q[2] == REAL(5.0));
		REQUIRE(r.coordinates() == q.coordinates());
		REQUIRE(back.coordinates() == p.coordinates());
		REQUIRE(displacement.components() == dx.components());

		static_assert(HasPointAddition<Point<Real, 3, Cartesian3>, TangentVector<3, Cartesian3>>::value, "Point plus displacement is valid");
		static_assert(!HasPointAddition<Point<Real, 3, Cartesian3>, Point<Real, 3, Cartesian3>>::value, "Point plus point is invalid");
		static_assert(!HasPointAddition<Point<Real, 3, Cartesian3>, TangentVector<3, Spherical3>>::value, "Point plus different-frame vector is invalid");
	}
}