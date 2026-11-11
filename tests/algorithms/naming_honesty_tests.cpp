// Naming-honesty regression tests (CONTRACT_POLICY section 11, task MinimalMathLibrary-wewi.5.6):
// result-returning root finders are spelled FindRootXDetailed (plain FindRootX returns Real and
// throws), IsAntiSymmetric is plain anti-symmetry while IsSkewHermitian/IsHermitian conjugate,
// to_string() is the one string-conversion spelling, FunctionsAnalyzer predicates are PascalCase,
// and the 2.0 spellings forward with identical meaning.
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Function.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/Matrix3D.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Tensor/Tensor2.h>
#include <mml/base/Geometry/Geometry3DBodies.h>
#include <mml/algorithms/RootFinding.h>
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#endif

#include <type_traits>

using namespace MML;

namespace MML::Tests::Algorithms::NamingHonestyTests
{
	namespace
	{
		Real Shifted(Real x) { return x * x - MML::Real(2.0); }
	}

	TEST_CASE("NamingHonesty_root_finding_plain_vs_Detailed", "[RootFinding][naming]")
	{
		RealFunction f(Shifted);
		RootFinding::RootFindingConfig config;
		config.tolerance = TOL(1e-12, Precision::EqualityTolerance);
		const Real rootTolerance = TOL(1e-11, Precision::EqualityTolerance);
		const Real methodTolerance = TOL(1e-10, Precision::EqualityTolerance);

		static_assert(std::is_same_v<decltype(RootFinding::FindRootBrent(f, MML::Real(1.0), MML::Real(2.0), MML::Real(1e-12))), Real>);
		static_assert(std::is_same_v<decltype(RootFinding::FindRootBrentDetailed(f, MML::Real(1.0), MML::Real(2.0), MML::Real(1e-12))), RootFinding::RootFindingResult>);
		static_assert(std::is_same_v<decltype(RootFinding::FindRootBrentDetailed(f, 1.0, 2.0, config)), RootFinding::RootFindingResult>);

		Real root = RootFinding::FindRootBrent(f, MML::Real(1.0), MML::Real(2.0), config.tolerance);
		REQUIRE(IsWithinAbs(root, std::sqrt(Real(2)), rootTolerance));
		REQUIRE(std::abs(f(root)) <= Real(4) * config.tolerance);

		auto r1 = RootFinding::FindRootBrentDetailed(f, MML::Real(1.0), MML::Real(2.0), config.tolerance);
		auto r2 = RootFinding::FindRootBrentDetailed(f, 1.0, 2.0, config);
		REQUIRE(r1.converged);
		REQUIRE(r2.converged);
		REQUIRE(r1.root == root);
		REQUIRE(r2.root == root);

		// every scalar method has the same three spellings
		REQUIRE(RootFinding::FindRootBisectionDetailed(f, MML::Real(1.0), MML::Real(2.0), methodTolerance).converged);
		REQUIRE(RootFinding::FindRootFalsePositionDetailed(f, MML::Real(1.0), MML::Real(2.0), methodTolerance).converged);
		REQUIRE(RootFinding::FindRootRiddersDetailed(f, MML::Real(1.0), MML::Real(2.0), methodTolerance).converged);
		REQUIRE(RootFinding::FindRootSecantDetailed(f, MML::Real(1.0), MML::Real(2.0), methodTolerance).converged);
		REQUIRE(RootFinding::FindRootNewtonDetailed(f, MML::Real(1.0), MML::Real(2.0), methodTolerance).converged);
		RealFunction df(+[](Real x) { return MML::Real(2.0) * x; });
		REQUIRE(RootFinding::FindRootNewtonDetailed(f, df, MML::Real(1.0), MML::Real(2.0), methodTolerance).converged);
		REQUIRE(IsWithinAbs(RootFinding::FindRootNewtonDetailed(f, df, 1.0, 2.0, config).root, root, rootTolerance));

		// callables
		auto rl = RootFinding::FindRootBrentDetailed([](Real x) { return x * x - MML::Real(2.0); }, MML::Real(1.0), MML::Real(2.0), config.tolerance);
		REQUIRE(rl.converged);
		REQUIRE(rl.root == root);

		// the Detailed form reports failure instead of throwing
		auto bad = RootFinding::FindRootBrentDetailed(f, MML::Real(3.0), MML::Real(4.0), config.tolerance);   // no sign change
		REQUIRE_FALSE(bad.converged);
		REQUIRE_THROWS_AS(RootFinding::FindRootBrent(f, MML::Real(3.0), MML::Real(4.0), config.tolerance), RootFindingError);
	}

	TEST_CASE("NamingHonesty_symmetry_predicates_conjugation", "[Matrix][naming]")
	{
		Matrix<Complex> a(2, 2, { Complex(0, 1), Complex(1, 2), Complex(-1, 2), Complex(0, -3) });
		// a(1,0) = -a(0,1)? (-1,2) vs -(1,2) = (-1,-2): no -> not anti-symmetric
		// a(1,0) = -conj(a(0,1)) = -(1,-2) = (-1,2): yes, diagonal purely imaginary -> skew-Hermitian
		REQUIRE_FALSE(a.IsAntiSymmetric());
		REQUIRE(a.IsSkewHermitian());
		REQUIRE_FALSE(a.IsHermitian());

		Matrix<Complex> h(2, 2, { Complex(1, 0), Complex(2, 3), Complex(2, -3), Complex(4, 0) });
		REQUIRE(h.IsHermitian());
		REQUIRE_FALSE(h.IsSymmetric());
		Matrix<Complex> s(2, 2, { Complex(1, 1), Complex(2, 3), Complex(2, 3), Complex(4, 0) });
		REQUIRE(s.IsSymmetric());
		REQUIRE_FALSE(s.IsHermitian());

		Matrix<Real> r(2, 2, { 0.0, 5.0, -5.0, 0.0 });
		REQUIRE(r.IsAntiSymmetric());
		REQUIRE(r.IsSkewHermitian());          // identical for real types
		REQUIRE(Matrix<Real>(2, 2, { 1.0, 2.0, 2.0, 1.0 }).IsHermitian());
	}

	TEST_CASE("NamingHonesty_to_string_and_accessors", "[naming]")
	{
		REQUIRE(to_string(AlgorithmStatus::Success) == std::string("Success"));
		Matrix3D<Real> m(1, 2, 3);
		REQUIRE(m.to_string().find("1") != std::string::npos);
		REQUIRE(m.to_string(3).size() >= m.to_string().size());
		Cube3D cube(1.0);
		const IBody& body = cube;
		REQUIRE_FALSE(body.to_string().empty());

		Matrix<Real> sym(2, 2, { 2.0, 1.0, 1.0, 2.0 });
		MatrixAnalyzer<Real> analyzer(sym);
		REQUIRE(analyzer.getMatrix()(0, 1) == 1.0);
		Tensor2<2> t(1, 1);
		t(0, 1) = 7.0;
		REQUIRE(t.ToMatrix()(0, 1) == 7.0);

		RealFunction square(+[](Real x) { return x * x; });
		RealFunctionAnalyzer fa(square);
		REQUIRE(fa.IsDefinedAtPoint(1.0));
		REQUIRE(fa.IsMonotonicOnSamples(1.0, 5.0, 50));
		REQUIRE_FALSE(fa.IsMonotonicOnSamples(-1.0, 1.0, 50));
		REQUIRE(fa.AppearsContinuous(-1.0, 1.0, 50));
		REQUIRE(fa.IsLocalOptimum(MML::Real(0.0), MML::Real(1e-6)));
		REQUIRE_FALSE(fa.IsInflectionPoint(MML::Real(0.0), MML::Real(1e-6)));
	}

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	TEST_CASE("NamingHonesty_deprecated_aliases_forward", "[naming][deprecated]")
	{
		RealFunction f(Shifted);
		RootFinding::RootFindingConfig config;
		config.tolerance = TOL(1e-12, Precision::EqualityTolerance);
		REQUIRE(RootFinding::FindRootBrent(f, 1.0, 2.0, config).root == RootFinding::FindRootBrentDetailed(f, 1.0, 2.0, config).root);
		REQUIRE(RootFinding::FindRootBisection(f, 1.0, 2.0, config).root == RootFinding::FindRootBisectionDetailed(f, 1.0, 2.0, config).root);
		const auto newtonFunction = [](Real x) { return x * x - Real(2); };
		const auto legacyNewton = RootFinding::FindRootNewton(newtonFunction, 1.0, 2.0, config);
		const auto detailedNewton = RootFinding::FindRootNewtonDetailed(newtonFunction, 1.0, 2.0, config);
		REQUIRE(legacyNewton.converged);
		REQUIRE(detailedNewton.converged);
		REQUIRE(legacyNewton.status == detailedNewton.status);
		REQUIRE(legacyNewton.root == detailedNewton.root);
		REQUIRE(legacyNewton.function_value == detailedNewton.function_value);
		REQUIRE(legacyNewton.function_evaluations == detailedNewton.function_evaluations);
		REQUIRE(std::abs(legacyNewton.function_value) <= config.tolerance);
		RealFunction df(+[](Real x) { return MML::Real(2.0) * x; });
		REQUIRE(RootFinding::FindRootNewton(f, df, 1.0, 2.0, config).root == RootFinding::FindRootNewtonDetailed(f, df, 1.0, 2.0, config).root);

		// the complex alias keeps its historical conjugated (skew-Hermitian) meaning
		Matrix<Complex> a(2, 2, { Complex(0, 1), Complex(1, 2), Complex(-1, 2), Complex(0, -3) });
		REQUIRE(a.isAntiSymmetric() == a.IsSkewHermitian());
		Matrix<Real> r(2, 2, { 0.0, 5.0, -5.0, 0.0 });
		REQUIRE(r.isAntiSymmetric() == r.IsAntiSymmetric());

		REQUIRE(ToString(AlgorithmStatus::Success) == to_string(AlgorithmStatus::Success));
		Matrix3D<Real> m(1, 2, 3);
		REQUIRE(m.ToString() == m.to_string());
		REQUIRE(m.ToStringFull(3) == m.to_string(3));
		Cube3D cube(1.0);
		REQUIRE(cube.ToString() == cube.to_string());
		Matrix<Real> sym(2, 2, { 2.0, 1.0, 1.0, 2.0 });
		MatrixAnalyzer<Real> analyzer(sym);
		REQUIRE(&analyzer.GetMatrix() == &analyzer.getMatrix());
		Tensor2<2> t(1, 1);
		REQUIRE(t.GetMatrix().IsEqualTo(t.ToMatrix()));
		RealFunction square(+[](Real x) { return x * x; });
		RealFunctionAnalyzer fa(square);
		REQUIRE(fa.isMonotonic(1.0, 5.0, 50) == fa.IsMonotonicOnSamples(1.0, 5.0, 50));
		REQUIRE(fa.isContinuous(-1.0, 1.0, 50) == fa.AppearsContinuous(-1.0, 1.0, 50));
		REQUIRE(fa.IsMonotonic(1.0, 5.0, 50) == fa.IsMonotonicOnSamples(1.0, 5.0, 50));
		REQUIRE(fa.IsContinuous(-1.0, 1.0, 50) == fa.AppearsContinuous(-1.0, 1.0, 50));
		REQUIRE(fa.FindDiscontinuities(-1.0, 1.0, 8).size() == fa.FindDiscontinuityCandidates(-1.0, 1.0, 8).size());
	}
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
}
