#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#endif

#include <type_traits>
#include <utility>

using namespace MML;

namespace MML::Tests::Base::DifferentialFormTests
{
	template<class, class, class = void>
	struct HasWedge : std::false_type { };

	template<class A, class B>
	struct HasWedge<A, B, std::void_t<decltype(Wedge(std::declval<A>(), std::declval<B>()))>> : std::true_type { };

	template<class, class, class = void>
	struct HasFormCall : std::false_type { };

	template<class Form, class... Args>
	struct HasFormCall<Form, std::tuple<Args...>, std::void_t<decltype(std::declval<Form>()(std::declval<Args>()...))>> : std::true_type { };

	TEST_CASE("DifferentialForm - Metadata component storage and scalar evaluation", "[DifferentialForm][DifferentialForms]")
	{
		using F0 = ScalarForm<3, Cartesian3>;
		using F2 = Form2<3, Cartesian3>;

		static_assert(F0::Degree == 0, "Scalar forms have degree zero");
		static_assert(F2::Degree == 2, "Two-forms have degree two");
		static_assert(F2::Dimension == 3, "Dimension is fixed at compile time");
		static_assert(F2::ComponentCount == 9, "Full component storage uses N^K entries");
		static_assert(std::is_same<typename F2::frame_type, Cartesian3>::value, "Frame tag is part of the form type");

		F0 scalar;
		scalar.Component() = REAL(2.5);
		REQUIRE(scalar() == REAL(2.5));

		F2 form;
		form.Component(0, 1) = REAL(3.0);
		form.Component(1, 0) = -REAL(3.0);
		REQUIRE(form.Component(0, 1) == REAL(3.0));
		REQUIRE(form.Component(1, 0) == -REAL(3.0));
	}

	TEST_CASE("DifferentialForm - Basis one-forms and covector conversion", "[DifferentialForm][DifferentialForms]")
	{
		auto dx = BasisOneForm<0, 3, Cartesian3>();
		auto dy = BasisOneForm<1, 3, Cartesian3>();
		OneForm<3, Cartesian3> dz = BasisOneForm<2, 3, Cartesian3>();
		TangentVector<3, Cartesian3> v{ REAL(4.0), REAL(5.0), REAL(6.0) };

		REQUIRE(dx(v) == REAL(4.0));
		REQUIRE(dy(v) == REAL(5.0));
		REQUIRE(dz(v) == REAL(6.0));
		static_assert(std::is_same<OneForm<3, Cartesian3>, Form1<3, Cartesian3>>::value, "OneForm is a readability alias for Form1");

		Covector<3, Cartesian3> alpha{ REAL(2.0), -REAL(3.0), REAL(7.0) };
		Form1<3, Cartesian3> form = ToForm(alpha);
		Covector<3, Cartesian3> roundTrip = ToCovector(form);

		REQUIRE(form(v) == Pair(alpha, v));
		REQUIRE(roundTrip.components() == alpha.components());
	}

	TEST_CASE("DifferentialForm - Two-form evaluation is antisymmetric", "[DifferentialForm][DifferentialForms]")
	{
		Form2<3, Cartesian3> omega;
		omega.Component(0, 1) = REAL(1.0);
		omega.Component(1, 0) = -REAL(1.0);

		TangentVector<3, Cartesian3> a{ REAL(2.0), REAL(3.0), REAL(0.0) };
		TangentVector<3, Cartesian3> b{ -REAL(1.0), REAL(4.0), REAL(0.0) };

		REQUIRE(omega(a, b) == REAL(11.0));
		REQUIRE(omega(b, a) == -REAL(11.0));
		REQUIRE(omega.IsAlternating());
	}

	TEST_CASE("DifferentialForm - Alternating helpers detect and repair raw component tables", "[DifferentialForm][DifferentialForms]")
	{
		Form2<3, Cartesian3> raw;
		raw.Component(0, 1) = REAL(4.0);
		raw.Component(1, 0) = REAL(1.0);
		raw.Component(2, 2) = REAL(7.0);

		REQUIRE_FALSE(raw.IsAlternating());

		Form2<3, Cartesian3> alternating = raw.AlternatingPart();
		REQUIRE(alternating.IsAlternating());
		REQUIRE(alternating.Component(0, 1) == Catch::Approx(REAL(1.5)));
		REQUIRE(alternating.Component(1, 0) == Catch::Approx(-REAL(1.5)));
		REQUIRE(alternating.Component(2, 2) == REAL(0.0));
		REQUIRE(alternating.AlternatingPart().Component(0, 1) == Catch::Approx(alternating.Component(0, 1)));
	}

	TEST_CASE("DifferentialForm - SetAlternatingComponent fills signed permutations", "[DifferentialForm][DifferentialForms]")
	{
		Form3<3, Cartesian3> volume;
		volume.SetAlternatingComponent(REAL(2.0), 0, 1, 2);

		REQUIRE(volume.IsAlternating());
		REQUIRE(volume.Component(0, 1, 2) == REAL(2.0));
		REQUIRE(volume.Component(1, 0, 2) == -REAL(2.0));
		REQUIRE(volume.Component(2, 1, 0) == -REAL(2.0));
		REQUIRE(volume.Component(0, 0, 1) == REAL(0.0));
	}

	TEST_CASE("DifferentialForm - Wedge product of basis forms is antisymmetric", "[DifferentialForm][DifferentialForms]")
	{
		auto dx = BasisOneForm<0, 3, Cartesian3>();
		auto dy = BasisOneForm<1, 3, Cartesian3>();
		auto dxWedgeDy = Wedge(dx, dy);
		auto dyWedgeDx = wedge(dy, dx);

		TangentVector<3, Cartesian3> a{ REAL(2.0), REAL(3.0), REAL(0.0) };
		TangentVector<3, Cartesian3> b{ -REAL(1.0), REAL(4.0), REAL(0.0) };

		REQUIRE(dxWedgeDy.Component(0, 1) == REAL(1.0));
		REQUIRE(dxWedgeDy.Component(1, 0) == -REAL(1.0));
		REQUIRE(dyWedgeDx.Component(0, 1) == -REAL(1.0));
		REQUIRE(dxWedgeDy.IsAlternating());
		REQUIRE(dyWedgeDx.IsAlternating());
		REQUIRE(dxWedgeDy(a, b) == REAL(11.0));
		REQUIRE(dyWedgeDx(a, b) == -REAL(11.0));
	}

	TEST_CASE("DifferentialForm - Volume form evaluates determinant-style", "[DifferentialForm][DifferentialForms]")
	{
		auto dx = BasisOneForm<0, 3, Cartesian3>();
		auto dy = BasisOneForm<1, 3, Cartesian3>();
		auto dz = BasisOneForm<2, 3, Cartesian3>();
		VolumeForm<3, Cartesian3> volume = Wedge(Wedge(dx, dy), dz);

		TangentVector<3, Cartesian3> e1{ REAL(1.0), REAL(0.0), REAL(0.0) };
		TangentVector<3, Cartesian3> e2{ REAL(0.0), REAL(1.0), REAL(0.0) };
		TangentVector<3, Cartesian3> e3{ REAL(0.0), REAL(0.0), REAL(1.0) };
		TangentVector<3, Cartesian3> sum{ REAL(1.0), REAL(1.0), REAL(0.0) };

		REQUIRE(volume(e1, e2, e3) == REAL(1.0));
		REQUIRE(volume(e2, e1, e3) == -REAL(1.0));
		REQUIRE(volume(e1, sum, e3) == REAL(1.0));
	}

	TEST_CASE("DifferentialForm - Wedge degree overflow is rejected", "[DifferentialForm][DifferentialForms]")
	{
		using OneForm2D = Form1<2, Cartesian2>;
		using TwoForm2D = Form2<2, Cartesian2>;

		static_assert(HasWedge<OneForm2D, OneForm2D>::value, "1+1 degree is valid in 2D");
		static_assert(!HasWedge<TwoForm2D, OneForm2D>::value, "2+1 degree overflows in 2D");
		SUCCEED("Invalid wedge degree overflow is rejected by the overload set.");
	}

	TEST_CASE("DifferentialForm - Form evaluation rejects wrong arity and vector categories", "[DifferentialForm][DifferentialForms]")
	{
		using OneForm3D = Form1<3, Cartesian3>;
		using TwoForm3D = Form2<3, Cartesian3>;
		using CartesianVector = TangentVector<3, Cartesian3>;
		using CartesianCovector = Covector<3, Cartesian3>;
		using SphericalVector = TangentVector<3, Spherical3>;

		static_assert(HasFormCall<OneForm3D, std::tuple<CartesianVector>>::value, "A one-form evaluates on one same-frame tangent vector");
		static_assert(!HasFormCall<OneForm3D, std::tuple<CartesianVector, CartesianVector>>::value, "A one-form rejects two vectors");
		static_assert(!HasFormCall<OneForm3D, std::tuple<CartesianCovector>>::value, "A one-form rejects covectors");
		static_assert(!HasFormCall<OneForm3D, std::tuple<SphericalVector>>::value, "A one-form rejects different-frame tangent vectors");
		static_assert(HasFormCall<TwoForm3D, std::tuple<CartesianVector, CartesianVector>>::value, "A two-form evaluates on two same-frame tangent vectors");
		static_assert(!HasFormCall<TwoForm3D, std::tuple<CartesianVector>>::value, "A two-form rejects one vector");
		SUCCEED("Invalid form evaluations are rejected by the overload set.");
	}
}