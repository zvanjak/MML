// Naming-convention regression tests for the derivation entry points (CONTRACT_POLICY section 11,
// task MinimalMathLibrary-wewi.5.5): numerical, forward-AD and reverse-AD Jacobians / Hessians /
// gradients share one Calc<Object>[<Method>] vocabulary, argument types select fixed vs dynamic
// size, and the 2.0 spellings remain as [[deprecated]] forwarders.
#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Function.h>
#include <mml/core/Derivation/Jacobians.h>
#include <mml/core/Derivation/Hessians.h>
#include <mml/core/Derivation/ADJacobians.h>
#endif

#include <cmath>
#include <vector>

using namespace MML;

namespace MML::Tests::Core::DerivationNamingTests
{
	namespace
	{
		VectorN<Real, 2> PolarFixed(const VectorN<Real, 2>& q) { return { q[0] * std::cos(q[1]), q[0] * std::sin(q[1]) }; }
		Vector<Real> PolarDyn(const Vector<Real>& q) { return Vector<Real>{ q[0] * std::cos(q[1]), q[0] * std::sin(q[1]) }; }
		auto PolarAD = [](const auto& q) {
			using std::cos; using std::sin;
			return std::vector<std::decay_t<decltype(q[0])>>{ q[0] * cos(q[1]), q[0] * sin(q[1]) };
		};
		auto BowlAD = [](const auto& x) { return x[0] * x[0] + 3.0 * x[0] * x[1]; };
	}

	TEST_CASE("DerivationNaming_Calc_vocabulary_agrees_across_methods", "[Derivation][AD][naming]")
	{
		VectorFunction<2> polar(PolarFixed);
		VectorN<Real, 2> q0{ 2.0, 0.5 };
		const Real r = 2.0, phi = 0.5;
		// exact: d(x)/d(phi) = -r sin(phi), d(y)/d(r) = sin(phi)
		MatrixNM<Real, 2, 2> Jfix = Derivation::CalcJacobian(polar, q0);
		REQUIRE(Jfix(0, 1) == Catch::Approx(-r * std::sin(phi)).epsilon(1e-9));
		REQUIRE(Jfix(1, 0) == Catch::Approx(std::sin(phi)).epsilon(1e-9));

		int evals = 0;
		Matrix<Real> Jdyn = Derivation::CalcJacobian(PolarDyn, Vector<Real>{ 2.0, 0.5 }, 0.0, &evals);
		REQUIRE(evals == 8);
		REQUIRE(Jdyn(0, 1) == Catch::Approx(Jfix(0, 1)).epsilon(1e-9));
		Matrix<Real> Jrect = Derivation::CalcJacobian(PolarDyn, Vector<Real>{ 2.0, 0.5 }, 2);   // explicit outputDim
		REQUIRE(Jrect.rows() == 2);
		Matrix<Real> Jin(2, 2);
		Derivation::CalcJacobianInPlace(PolarDyn, Vector<Real>{ 2.0, 0.5 }, Jin);
		REQUIRE(Jin.IsEqualTo(Jdyn, 1e-12));
		Matrix<Real> J2 = Derivation::CalcJacobianOrder2(PolarDyn, Vector<Real>{ 2.0, 0.5 });
		REQUIRE(J2(0, 1) == Catch::Approx(Jfix(0, 1)).epsilon(1e-6));

		Matrix<Real> Jfwd = AD::CalcJacobianForwardAD(PolarAD, std::vector<Real>{ 2.0, 0.5 });
		Matrix<Real> Jrev = AD::CalcJacobianReverseAD(PolarAD, Vector<Real>{ 2.0, 0.5 });
		REQUIRE(Jfwd(0, 1) == Catch::Approx(-r * std::sin(phi)).epsilon(1e-14));
		REQUIRE(Jrev.IsEqualTo(Jfwd, 1e-14));

		Matrix<Real> Jode;
		Derivation::CalcJacobianODE([](Real, const Vector<Real>& x, Vector<Real>& dxdt) { dxdt = PolarDyn(x); }, 0.0, Vector<Real>{ 2.0, 0.5 }, Jode);
		REQUIRE(Jode(0, 1) == Catch::Approx(Jfix(0, 1)).epsilon(1e-6));
		Matrix<Real> Jmap;
		Derivation::CalcJacobianMap(PolarDyn, Vector<Real>{ 2.0, 0.5 }, Jmap);
		REQUIRE(Jmap.IsEqualTo(Jode, 1e-12));

		ScalarFunction<2> bowl([](const VectorN<Real, 2>& x) { return x[0] * x[0] + 3.0 * x[0] * x[1]; });
		MatrixNM<Real, 2, 2> H = Derivation::CalcHessian(bowl, q0);
		Matrix<Real> Had = AD::CalcHessianForwardAD(BowlAD, std::vector<Real>{ 2.0, 0.5 });
		REQUIRE(H(0, 0) == Catch::Approx(2.0).epsilon(1e-7));
		REQUIRE(H(0, 1) == Catch::Approx(3.0).epsilon(1e-7));
		REQUIRE(Had(0, 1) == Catch::Approx(3.0).epsilon(1e-14));

		std::vector<Real> gf = AD::CalcGradientForwardAD(BowlAD, std::vector<Real>{ 2.0, 0.5 });
		std::vector<Real> gr = AD::CalcGradientReverseAD(BowlAD, std::vector<Real>{ 2.0, 0.5 });
		REQUIRE(gf[0] == Catch::Approx(2 * 2.0 + 3.0 * 0.5));
		REQUIRE(gr[1] == Catch::Approx(3.0 * 2.0));
		auto [v, d] = AD::CalcDerivativeForwardAD([](auto x) { return x * x * x; }, 2.0);
		REQUIRE(v == Catch::Approx(8.0));
		REQUIRE(d == Catch::Approx(12.0));
	}

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	TEST_CASE("DerivationNaming_deprecated_aliases_forward", "[Derivation][AD][naming][deprecated]")
	{
		VectorFunction<2> polar(PolarFixed);
		VectorN<Real, 2> q0{ 2.0, 0.5 };
		Vector<Real> p{ 2.0, 0.5 };
		REQUIRE(Derivation::calcJacobian(polar, q0).IsEqualTo(Derivation::CalcJacobian(polar, q0)));
		REQUIRE(Derivation::calcJacobianDyn(PolarDyn, p).IsEqualTo(Derivation::CalcJacobian(PolarDyn, p)));
		REQUIRE(Derivation::calcJacobianDyn(PolarDyn, p, 2).IsEqualTo(Derivation::CalcJacobian(PolarDyn, p, 2)));
		REQUIRE(Derivation::calcJacobianDyn2(PolarDyn, p).IsEqualTo(Derivation::CalcJacobianOrder2(PolarDyn, p)));
		Matrix<Real> a(2, 2), b(2, 2);
		Derivation::calcJacobianDynInPlace(PolarDyn, p, a);
		Derivation::CalcJacobianInPlace(PolarDyn, p, b);
		REQUIRE(a.IsEqualTo(b));
		Derivation::calcJacobianDynColumns(PolarDyn, p, a, 2, 0.0, nullptr);
		Derivation::CalcJacobianColumns(PolarDyn, p, b, 2, 0.0, nullptr);
		REQUIRE(a.IsEqualTo(b));
		auto derivs = [](Real, const Vector<Real>& x, Vector<Real>& dxdt) { dxdt = PolarDyn(x); };
		REQUIRE(Derivation::calcJacobianODE(derivs, 0.0, p).IsEqualTo(Derivation::CalcJacobianODE(derivs, 0.0, p)));
		Derivation::calcJacobianMap(PolarDyn, p, a);
		Derivation::CalcJacobianMap(PolarDyn, p, b);
		REQUIRE(a.IsEqualTo(b));

		ScalarFunction<2> bowl([](const VectorN<Real, 2>& x) { return x[0] * x[0] + 3.0 * x[0] * x[1]; });
		REQUIRE(Derivation::calcHessian(bowl, q0).IsEqualTo(Derivation::CalcHessian(bowl, q0)));

		REQUIRE(AD::jacobianForwardAD(PolarAD, p).IsEqualTo(AD::CalcJacobianForwardAD(PolarAD, p)));
		REQUIRE(AD::jacobianReverseAD(PolarAD, p).IsEqualTo(AD::CalcJacobianReverseAD(PolarAD, p)));
		REQUIRE(AD::hessianForwardAD(BowlAD, p).IsEqualTo(AD::CalcHessianForwardAD(BowlAD, p)));
		REQUIRE(AD::gradient(BowlAD, std::vector<Real>{ 2.0, 0.5 }) == AD::CalcGradientForwardAD(BowlAD, std::vector<Real>{ 2.0, 0.5 }));
		REQUIRE(AD::gradientReverse(BowlAD, std::vector<Real>{ 2.0, 0.5 }) == AD::CalcGradientReverseAD(BowlAD, std::vector<Real>{ 2.0, 0.5 }));
		REQUIRE(AD::derivative([](auto x) { return x * x; }, 3.0) == AD::CalcDerivativeForwardAD([](auto x) { return x * x; }, 3.0));
	}
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
}
