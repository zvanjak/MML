#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/ODESolvers/ODESteppers.h>
#endif
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <array>
#include <cmath>
#include <type_traits>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	class ExtrapolationGrowthODE : public IODESystem {
		int _dimension;
	public:
		mutable int evaluations = 0;
		explicit ExtrapolationGrowthODE(int dimension = 1) : _dimension(dimension) {}
		int getDim() const override { return _dimension; }
		void derivs(Real, const Vector<Real>& state, Vector<Real>& rhs) const override {
			++evaluations;
			rhs = state;
		}
	};

	class ExtrapolationPolynomialODE : public IODESystem {
		int _degree;
	public:
		explicit ExtrapolationPolynomialODE(int degree) : _degree(degree) {}
		int getDim() const override { return 1; }
		void derivs(Real time, const Vector<Real>&, Vector<Real>& rhs) const override {
			const Real z = 2 * time - 1;
			Real previous = 1, current = z;
			for (int i = 2; i <= _degree; ++i) {
				const Real next = 2 * z * current - previous;
				previous = current;
				current = next;
			}
			rhs[0] = current;
		}
	};

	struct ColumnReference {
		Real candidate;
		Real epsilon;
		int evaluations; // Includes the supplied initial derivative value.
		Real conditioningMagnitude;
	};

	// Exact-rational modified-midpoint values combined with independent Lagrange weights at zero.
	// Conditioning magnitude is sum |L_i(0) Y_i|, accounting for cancellation in higher columns.
	const std::array<ColumnReference, 8> evenColumns{{
		{Real(39.666666666666666666666666666667L), Real(0.22950819672131148), 8, Real(53.666666666666666666666666666667L)},
		{Real(50.451851851851851851851851851852L), Real(0.046581565729276644), 14, Real(125.11851851851851851851851851852L)},
		{Real(53.838095238095238095238095238095L), Real(0.007718729497124773), 22, Real(277.43809523809523809523809523810L)},
		{Real(54.500895943562610229276895943563L), Real(0.0009553729815696765), 32, Real(603.93017283950617283950617283951L)},
		{Real(54.589043165339461635757932054228L), Real(8.809412066918134e-05), 44, Real(1307.2123813086776049739012701976L)},
		{Real(54.597501676957459270384440452468L), Real(6.209728166729058e-06), 58, Real(2828.7300831695389518518770219450L)},
		{Real(54.598113842261990410138558286706L), Real(3.4407940206159724e-07), 74, Real(6132.8410716745531560346375161190L)},
		{Real(54.598148409147736725132566871183L), Real(1.5351290936523966e-08), 92, Real(13331.865553765511698753582612750L)}
	}};
	const std::array<ColumnReference, 8> bulirschColumns{{
		{Real(39.666666666666666666666666666667L), Real(0.22950819672131148), 8, Real(53.666666666666666666666666666667L)},
		{Real(50.451851851851851851851851851852L), Real(0.046581565729276644), 14, Real(125.11851851851851851851851851852L)},
		{Real(53.838095238095238095238095238095L), Real(0.007718729497124773), 22, Real(277.43809523809523809523809523810L)},
		{Real(54.527829816883314825701656977377L), Real(0.0006900789718743499), 34, Real(302.74878219783569577808260935833L)},
		{Real(54.594285001692409099816507223915L), Real(3.735500016271414e-05), 50, Real(430.35712478781261850044918827988L)},
		{Real(54.598051695389934369855744122555L), Real(9.409536602115212e-07), 74, Real(387.24946524127967073578258623987L)},
		{Real(54.598148599520782624726156021381L), Real(1.3616703816977757e-08), 106, Real(485.68690121063423175789939504184L)},
		{Real(54.598150023725828985701602132541L), Real(8.894462151825779e-11), 154, Real(413.75354995303849310963751651833L)}
	}};

	// H=8 keeps all final corrections resolvable in float; retain H=4 for wider Real types.
	const std::array<ColumnReference, 8> floatEvenColumns{{
		{Real(467.66666666666669), Real(0.38691322901849218), 8, Real(537.66666666666663)},
		{Real(1176.8814814814814), Real(0.13380233466515176), 14, Real(1981.148148148148)},
		{Real(1974.4095238095238), Real(0.050465994058160601), 22, Real(6024.0666666666666)},
		{Real(2550.4790546737213), Real(0.018062293078486252), 32, Real(16344.454363315697)},
		{Real(2838.6025510840327), Real(0.0056370075120453942), 44, Real(41364.407072488553)},
		{Real(2943.9295872985667), Real(0.0014598185033611131), 58, Real(100219.90562508794)},
		{Real(2973.2379799843502), Real(0.00030794014386016026), 74, Real(236151.76729415989)},
		{Real(2979.6445349577862), Real(5.3071253783774595e-05), 92, Real(546490.800361278)}
	}};
	const std::array<ColumnReference, 8> floatBulirschColumns{{
		{Real(467.66666666666669), Real(0.38691322901849218), 8, Real(537.66666666666663)},
		{Real(1176.8814814814814), Real(0.13380233466515176), 14, Real(1981.148148148148)},
		{Real(1974.4095238095238), Real(0.050465994058160601), 22, Real(6024.0666666666666)},
		{Real(2638.516789687983), Real(0.013977879681803496), 34, Real(8856.3377420689358)},
		{Real(2906.9721048224751), Real(0.0028849068338862283), 50, Real(16171.590848978574)},
		{Real(2973.1162615728449), Real(0.00030888800670801701), 74, Real(16842.77348175654)},
		{Real(2980.4772408640815), Real(1.928830779741851e-05), 106, Real(23560.950524873795)},
		{Real(2980.9446091434402), Real(5.4421082154738444e-07), 154, Real(21138.153317878867)}
	}};

	struct ConvergenceReference {
		int columns;
		int steps;
		Real epsilon;
		Real error;
	};
	[[maybe_unused]] const std::array<ConvergenceReference, 9> convergence{{
		{2, 4, Real(0.00017926454770775203), Real(2.8819380513975028e-05)},
		{2, 8, Real(2.14967890056256e-05), Real(1.93371803233617e-06)},
		{2, 16, Real(2.6188677674429928e-06), Real(1.2517741942831266e-07)},
		{3, 2, Real(1.0089723817327902e-05), Real(2.4678715039403465e-06)},
		{3, 4, Real(3.3059627010999028e-07), Real(4.6024589733275756e-08)},
		{3, 8, Real(1.0492660548084885e-08), Real(7.8529855196112324e-10)},
		{4, 1, Real(3.6651232155366207e-06), Real(1.508518569044884e-06)},
		{4, 2, Real(3.5195445827650423e-08), Real(8.6941304236838222e-09)},
		{4, 4, Real(2.9718644780569967e-10), Real(4.1222483957916669e-11)}
	}};

	struct PolynomialConvergenceReference {
		int degree;
		int columns;
		int steps;
		std::array<Real, 16> epsilons;
		Real error;
	};
	// y'=T_degree(2t-1), y(0)=0: exact endpoint is -1/(degree^2-1) for even degree.
	// Independent exact-rational trapezoidal sums and Lagrange weights give these errors.
	// Per-step epsilons separate the desired column from every earlier correction.
	[[maybe_unused]] const std::array<PolynomialConvergenceReference, 9> floatConvergence{{
		{4, 2, 2, {Real(2.5), Real(0.83333333333333337)}, Real(0.0041666666666666666)},
		{4, 2, 4, {
			Real(0.18275862068965518), Real(0.030726256983240222),
			Real(0.04954954954954955), Real(0.18275862068965518)
		}, Real(0.00026041666666666666)},
		{4, 2, 8, {
			Real(0.13216424294268606), Real(0.010067114093959731),
			Real(0.00059958574076092886), Real(0.0066690378801351593),
			Real(0.010904332654841524), Real(0.0011446409989594173),
			Real(0.016844226893175929), Real(0.04829634260706471)
		}, Real(1.6276041666666666e-05)},
		{12, 3, 2, {
			Real(0.74565724959483215), Real(0.51259901694779431)
		}, Real(0.28995029207086603)},
		{12, 3, 8, {
			Real(0.32161076945113154), Real(0.069258075126744159),
			Real(0.01294787599141743), Real(0.0099731911079799663),
			Real(0.022689808916409083), Real(0.019013050107936622),
			Real(0.041451733959265664), Real(0.30703048026832169)
		}, Real(0.0001760046845185957)},
		{12, 3, 16, {
			Real(0.62224736938524239), Real(0.019431411270466083),
			Real(0.011485592413283278), Real(0.0082786318558365091),
			Real(0.007200417655057887), Real(0.0067008218452377314),
			Real(0.0040993521316468762), Real(0.0051290112393034522),
			Real(0.0059915274206130045), Real(0.0055394246062377582),
			Real(0.0067008218452377314), Real(0.0056056746751230501),
			Real(0.0091227492549850776), Real(0.011485592413283278),
			Real(0.015862031763474303), Real(0.32118513186670161)
		}, Real(2.8518817317326145e-06)},
		{14, 4, 2, {
			Real(0.098914578878179227), Real(0.044627592285215956)
		}, Real(0.12364050129852087)},
		{14, 4, 4, {
			Real(0.12525094169743059), Real(0.00653130538046313),
			Real(0.0056275143445669073), Real(0.11912653090946912)
		}, Real(0.00091682692973212791)},
		{14, 4, 8, {
			Real(0.045422422397469846), Real(0.00075521644872808232),
			Real(0.00028673314554907817), Real(0.00014031652939452656),
			Real(0.00019643493912855755), Real(0.00028673314554907817),
			Real(0.00075521644872808232), Real(0.045422422397469846)
		}, Real(4.097597655934174e-06)}
	}};

	template<typename Stepper>
	void CheckFloatConvergence() {
		Real previousError = 0;
		int previousColumns = 0, previousSteps = 0;
		for (const auto& reference : floatConvergence) {
			CAPTURE(reference.degree, reference.columns, reference.steps);
			ExtrapolationPolynomialODE system(reference.degree);
			Stepper stepper(system);
			Vector<Real> state{0}, derivative{1};
			const Real h = Real(1) / reference.steps;
			for (int i = 0; i < reference.steps; ++i) {
				CAPTURE(i, reference.epsilons[i]);
				const auto result = stepper.doStep(i * h, state, derivative, h, reference.epsilons[i]);
				REQUIRE(result.accepted);
				REQUIRE(result.funcEvals == evenColumns[reference.columns - 2].evaluations - 1);
				CHECK(result.hDone == h);
			}
			const Real exact = -Real(1) / (reference.degree * reference.degree - 1);
			const Real error = std::abs(state[0] - exact);
			CHECK_THAT(error, WithinAbs(reference.error, Precision::RoundoffZeroThreshold * reference.steps));
			REQUIRE(error > 16 * Precision::Eps);
			if (reference.columns == previousColumns) {
				// The sixth-order fixture first refines by four to avoid an earlier-column cancellation.
				const Real ratio = Real(reference.steps) / previousSteps;
				const Real observedOrder = std::log(previousError / error) / std::log(ratio);
				CAPTURE(observedOrder);
				CHECK(observedOrder > 2 * reference.columns - 1);
				CHECK(observedOrder < 2 * reference.columns + 1);
			}
			previousError = error;
			previousColumns = reference.columns;
			previousSteps = reference.steps;
		}
	}

	Real TwoColumnGrowth(Real h) {
		return 1 + h * (1 + h * (Real(1) / 2 + h * (Real(1) / 6 + h * (Real(1) / 24 + h / 192))));
	}
}

TEMPLATE_TEST_CASE("Bulirsch Stoer initializes each extrapolated column from its own midpoint estimate",
	"[BulirschStoer][bulirsch-extrapolation]", BulirschStoer_Stepper, BulirschStoerRational_Stepper) {
	ExtrapolationGrowthODE system;
	TestType stepper(system);
	Vector<Real> state{1}, derivative{1};
	const auto result = stepper.doStep(0, state, derivative, Real(0.5), 1);
	const Real expected = Real(10129) / 6144;
	const Real correction = Real(49) / 24576;
	REQUIRE(result.accepted);
	CHECK(result.funcEvals == 7);
	CHECK(result.funcEvals == system.evaluations);
	CHECK(result.hDone == Real(0.5));
	CHECK_THAT(state[0], WithinAbs(expected, Precision::RoundoffZeroThreshold));
	CHECK(derivative[0] == state[0]);
	CHECK_THAT(result.errMax, WithinAbs(correction / (1 + expected + Precision::DivisionSafetyThreshold),
		Precision::RoundoffZeroThreshold));
	CHECK_THAT(stepper.interpolate(Real(0.5))[0], WithinAbs(expected, Precision::RoundoffZeroThreshold));
}

TEMPLATE_TEST_CASE("Bulirsch Stoer accumulates every Neville correction through all nine columns",
	"[BulirschStoer][bulirsch-extrapolation]", BulirschStoer_Stepper, BulirschStoerRational_Stepper) {
	constexpr bool floatPrecision = std::is_same_v<Real, float>;
	const Real h = floatPrecision ? 8 : 4;
	const auto& references = std::is_same_v<TestType, BulirschStoer_Stepper>
		? (floatPrecision ? floatEvenColumns : evenColumns)
		: (floatPrecision ? floatBulirschColumns : bulirschColumns);
	ExtrapolationGrowthODE system;
	for (const auto& reference : references) {
		CAPTURE(reference.evaluations, reference.epsilon);
		TestType stepper(system);
		Vector<Real> state{1}, derivative{1};
		const int before = system.evaluations;
		const auto result = stepper.doStep(0, state, derivative, h, reference.epsilon);
		REQUIRE(result.accepted);
		CHECK(result.funcEvals == reference.evaluations - 1);
		CHECK(result.funcEvals == system.evaluations - before);
		CHECK(result.hDone == h);
		CHECK(result.errMax >= 0);
		CHECK(result.errMax < 1);
		const Real roundoff = Precision::RoundoffZeroThreshold * reference.candidate +
			Precision::Eps * reference.conditioningMagnitude;
		CHECK(IsWithinAbs(state[0], reference.candidate, roundoff));
		CHECK(derivative[0] == state[0]);
	}
	CHECK(TestType(system).order() == 18);
}

TEMPLATE_TEST_CASE("Bulirsch Stoer recurrence is componentwise reversible and reusable",
	"[BulirschStoer][bulirsch-extrapolation]", BulirschStoer_Stepper, BulirschStoerRational_Stepper) {
	ExtrapolationGrowthODE system(3);
	TestType stepper(system);
	Vector<Real> state{1, -3, 0}, derivative = state;
	Real time = 0, expected = 1;
	const Real rejectedStep = std::is_same_v<Real, float> ? 16 : 4;
	const auto rejected = stepper.doStep(time, state, derivative, rejectedStep, Precision::StrictTolerance);
	REQUIRE_FALSE(rejected.accepted);
	CHECK(state[0] == 1);
	CHECK(state[1] == -3);
	CHECK(state[2] == 0);
	CHECK(derivative[0] == 1);
	CHECK(derivative[1] == -3);
	CHECK(derivative[2] == 0);
	for (Real h : {Real(0.5), Real(-0.25), Real(0.125)}) {
		CAPTURE(h);
		const auto result = stepper.doStep(time, state, derivative, h, 1);
		REQUIRE(result.accepted);
		CHECK(result.funcEvals == 7);
		expected *= TwoColumnGrowth(h);
		CHECK_THAT(state[0], WithinAbs(expected, Precision::RoundoffZeroThreshold * expected));
		CHECK_THAT(state[1], WithinAbs(-3 * expected, Precision::RoundoffZeroThreshold * 3 * expected));
		CHECK(state[2] == 0);
		CHECK(derivative[0] == state[0]);
		CHECK(derivative[1] == state[1]);
		CHECK(derivative[2] == 0);
		time += h;
	}
}

TEMPLATE_TEST_CASE("Bulirsch Stoer fixed columns exhibit global fourth sixth and eighth order convergence",
	"[BulirschStoer][bulirsch-extrapolation]", BulirschStoer_Stepper, BulirschStoerRational_Stepper) {
	if constexpr (std::is_same_v<Real, float>) {
		CheckFloatConvergence<TestType>();
	} else {
		ExtrapolationGrowthODE system;
		Real previousError = 0;
		for (const auto& reference : convergence) {
			CAPTURE(reference.columns, reference.steps, reference.epsilon);
			TestType stepper(system);
			Vector<Real> state{1}, derivative{1};
			const Real h = Real(1) / reference.steps;
			for (int i = 0; i < reference.steps; ++i) {
				const auto result = stepper.doStep(i * h, state, derivative, h, reference.epsilon);
				REQUIRE(result.accepted);
				CHECK(result.funcEvals == evenColumns[reference.columns - 2].evaluations - 1);
				CHECK(result.hDone == h);
			}
			const Real error = std::abs(state[0] - std::exp(Real(1)));
			CHECK_THAT(error, WithinAbs(reference.error,
				Precision::RoundoffZeroThreshold * reference.steps * std::exp(Real(1))));
			if (reference.steps != (reference.columns == 2 ? 4 : reference.columns == 3 ? 2 : 1)) {
				const Real observedOrder = std::log2(previousError / error);
				CAPTURE(observedOrder);
				CHECK(observedOrder > 2 * reference.columns - 1);
				CHECK(observedOrder < 2 * reference.columns + 1);
			}
			previousError = error;
		}
	}
}

TEMPLATE_TEST_CASE("Bulirsch Stoer rejection preserves accepted caches and accounts for every midpoint evaluation",
	"[BulirschStoer][bulirsch-extrapolation]", BulirschStoer_Stepper, BulirschStoerRational_Stepper) {
	ExtrapolationGrowthODE system(3);
	TestType stepper(system);
	Vector<Real> state{1, -3, 0}, derivative = state;
	const auto accepted = stepper.doStep(0, state, derivative, Real(0.5), 1);
	REQUIRE(accepted.accepted);
	CHECK(system.evaluations == accepted.funcEvals);
	const Vector<Real> savedState = state, savedDerivative = derivative, savedFinal = stepper.getFinalDeriv();
	const std::array<Real, 3> times{{0, Real(0.25), Real(0.5)}};
	const std::array<Vector<Real>, 3> savedInterpolation{{
		stepper.interpolate(times[0]), stepper.interpolate(times[1]), stepper.interpolate(times[2])
	}};

	system.evaluations = 0;
	const auto rejected = stepper.doStep(Real(0.5), state, derivative, 16, Precision::StrictTolerance);
	REQUIRE_FALSE(rejected.accepted);
	const int expectedEvaluations = std::is_same_v<TestType, BulirschStoer_Stepper> ? 90 : 152;
	CHECK(rejected.funcEvals == expectedEvaluations);
	CHECK(system.evaluations == expectedEvaluations);
	CHECK(rejected.hDone == 16);
	CHECK(rejected.hNext > 0);
	CHECK(rejected.hNext < rejected.hDone);
	for (int component = 0; component < system.getDim(); ++component) {
		CAPTURE(component);
		CHECK(state[component] == savedState[component]);
		CHECK(derivative[component] == savedDerivative[component]);
		CHECK(stepper.getFinalDeriv()[component] == savedFinal[component]);
	}
	for (size_t i = 0; i < times.size(); ++i) {
		const auto interpolation = stepper.interpolate(times[i]);
		for (int component = 0; component < system.getDim(); ++component) {
			CAPTURE(times[i], component);
			CHECK(interpolation[component] == savedInterpolation[i][component]);
		}
	}
	const auto reused = stepper.doStep(Real(0.5), state, derivative, Real(-0.25), 1);
	REQUIRE(reused.accepted);
	CHECK(reused.funcEvals == 7);
	CHECK(system.evaluations == expectedEvaluations + reused.funcEvals);
	for (int component = 0; component < system.getDim(); ++component) {
		const Real expected = savedState[component] * TwoColumnGrowth(Real(-0.25));
		CHECK_THAT(state[component], WithinAbs(expected, Precision::RoundoffZeroThreshold * std::abs(expected)));
		CHECK(derivative[component] == state[component]);
	}
}
