// Naming-convention regression tests for interface vocabulary (CONTRACT_POLICY section 11,
// task MinimalMathLibrary-wewi.5.4): everything declared on an mml/interfaces type is lowerCamel
// (numCovar, isContravar, contract, component, transformsAsTensor, getEquidistantCovering,
// getValues, getTrace); the PascalCase 2.0 spellings remain as [[deprecated]] non-virtual forwarders
// on the interface, so callers keep working while implementers override the lowerCamel name.
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/interfaces/ITensor.h>
#include <mml/interfaces/ITensorField.h>
#include <mml/interfaces/IInterval.h>
#include <mml/interfaces/IFunction.h>
#include <mml/base/Tensor/Tensor2.h>
#include <mml/base/Tensor/Tensor3.h>
#include <mml/base/Tensor/Tensor4.h>
#include <mml/base/Intervals.h>
#include <mml/base/Function.h>
#include <mml/core/Fields/TensorFieldAdapters.h>
#endif

#include <type_traits>

using namespace MML;

namespace MML::Tests::Interfaces::InterfaceNamingTests
{
	// An implementer written against the 2.1 vocabulary only.
	class DiagonalField2 : public ITensorField2<3>
	{
	public:
		DiagonalField2() : ITensorField2<3>(0, 2) {}
		Real component(int i, int j, const VectorN<Real, 3>& pos) const override { return i == j ? pos[i] : 0.0; }
		Tensor2<3> operator()(const VectorN<Real, 3>& pos) const override
		{
			Tensor2<3> t(0, 2);
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					t(i, j) = component(i, j, pos);
			return t;
		}
	};

	class ConnectionField3 : public ITensorField3<2>
	{
	public:
		ConnectionField3() : ITensorField3<2>(1, 2) {}
		bool transformsAsTensor() const override { return false; }
		Real component(int, int, int, const VectorN<Real, 2>&) const override { return 0.0; }
		Tensor3<2> operator()(const VectorN<Real, 2>&) const override { return Tensor3<2>(1, 2); }
	};

	TEST_CASE("InterfaceNaming_ITensor_vocabulary_is_lowerCamel", "[ITensor][naming]")
	{
		Tensor2<3> g(2, 0);           // two covariant indices
		g(0, 0) = 1.0; g(1, 1) = 2.0; g(2, 2) = 3.0;
		const ITensor2<3>& ig = g;
		REQUIRE(ig.numCovar() == 2);
		REQUIRE(ig.numContravar() == 0);
		REQUIRE_FALSE(ig.isContravar(0));
		REQUIRE(g.isCovar(1));

		Tensor2<3> mixed(1, 1);       // one covariant, one contravariant: trace is defined
		mixed(0, 0) = 1.0; mixed(1, 1) = 2.0; mixed(2, 2) = 3.0;
		REQUIRE(mixed.contract() == Catch::Approx(6.0));

		Tensor3<2> t3(1, 2);
		REQUIRE(t3.numCovar() + t3.numContravar() == 3);
		Tensor4<2> t4(2, 2);
		REQUIRE(t4.numCovar() == 2);
	}

	TEST_CASE("InterfaceNaming_ITensorField_component_and_transformsAsTensor", "[ITensorField][naming]")
	{
		DiagonalField2 field;
		VectorN<Real, 3> pos{ 1.0, 2.0, 3.0 };
		const ITensorField2<3>& ifield = field;
		REQUIRE(ifield.component(1, 1, pos) == 2.0);
		REQUIRE(ifield.component(0, 1, pos) == 0.0);
		REQUIRE(field(pos)(2, 2) == 3.0);

		ConnectionField3 conn;
		const ITensorField3<2>& iconn = conn;
		REQUIRE_FALSE(iconn.transformsAsTensor());
	}

	TEST_CASE("InterfaceNaming_IInterval_IFunction_helpers", "[IInterval][IFunction][naming]")
	{
		ClosedInterval unit(0.0, 1.0);
		std::vector<Real> pts;
		unit.getEquidistantCovering(5, pts);
		REQUIRE(pts.size() == 5);
		REQUIRE(pts.front() == 0.0);
		REQUIRE(pts.back() == 1.0);
		const IInterval& iunit = unit;
		std::vector<Real> pts2;
		iunit.getEquidistantCovering(3, pts2);
		REQUIRE(pts2[1] == Catch::Approx(0.5));

		RealFunction square(+[](const Real x) { return x * x; });
		Vector<Real> xs, ys;
		square.getValues(0.0, 2.0, 3, xs, ys);
		REQUIRE(xs.size() == 3);
		REQUIRE(ys[2] == Catch::Approx(4.0));

		ParametricCurve<2> circle(+[](Real t) { return VectorN<Real, 2>{ std::cos(t), std::sin(t) }; });
		auto trace = circle.getTrace(0.0, Constants::PI, 3);
		REQUIRE(trace.size() == 3);
		REQUIRE(trace[2][0] == Catch::Approx(-1.0));
	}

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	TEST_CASE("InterfaceNaming_deprecated_aliases_forward", "[ITensor][naming][deprecated]")
	{
		Tensor2<3> mixed(1, 1);
		mixed(0, 0) = 1.0; mixed(1, 1) = 2.0; mixed(2, 2) = 3.0;
		REQUIRE(mixed.NumCovar() == mixed.numCovar());
		REQUIRE(mixed.NumContravar() == mixed.numContravar());
		REQUIRE(mixed.IsContravar(0) == mixed.isContravar(0));
		REQUIRE(mixed.Contract() == mixed.contract());
		Tensor3<2> t3(1, 2);
		REQUIRE(t3.NumCovar() == t3.numCovar());
		REQUIRE(t3.Contract(0, 1).size() == t3.contract(0, 1).size());
		Tensor4<2> t4(2, 2);
		REQUIRE(t4.Contract(0, 1).numCovar() == t4.contract(0, 1).numCovar());

		DiagonalField2 field;
		VectorN<Real, 3> pos{ 1.0, 2.0, 3.0 };
		REQUIRE(field.Component(1, 1, pos) == field.component(1, 1, pos));   // forwarder on the interface
		ConnectionField3 conn;
		REQUIRE(conn.TransformsAsTensor() == conn.transformsAsTensor());

		ClosedInterval unit(0.0, 1.0);
		std::vector<Real> a, b;
		unit.GetEquidistantCovering(4, a);
		unit.getEquidistantCovering(4, b);
		REQUIRE(a == b);

		RealFunction square(+[](const Real x) { return x * x; });
		Vector<Real> xs, ys, xs2, ys2;
		square.GetValues(0.0, 2.0, 3, xs, ys);
		square.getValues(0.0, 2.0, 3, xs2, ys2);
		REQUIRE(ys.IsEqualTo(ys2));
		ParametricCurve<2> circle(+[](Real t) { return VectorN<Real, 2>{ std::cos(t), std::sin(t) }; });
		REQUIRE(circle.GetTrace(0.0, 1.0, 2).size() == circle.getTrace(0.0, 1.0, 2).size());
	}
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
}
