// Keep the infrastructure first so aggregate headers cannot supply its dependencies.
#include <mml/algorithms/ODESolvers/ODEStepperInfrastructure.h>

#include <catch2/catch_test_macros.hpp>

using namespace MML;

namespace {
	class StationaryHeaderODE final : public IODESystem {
	public:
		mutable int calls = 0;
		int getDim() const override { return 2; }
		void derivs(Real, const Vector<Real>&, Vector<Real>& derivative) const override {
			++calls;
			derivative[0] = 0;
			derivative[1] = 0;
		}
	};
}

TEST_CASE("ODE stepper infrastructure supports DOP853 interpolation as the first include",
	"[ode-header][ode-work]") {
	for (Real step : {Real(-1), Real(1)}) {
		StationaryHeaderODE system;
		DormandPrince8Interpolator interpolator;
		const Vector<Real> state{Real(2), Real(-3)};
		std::vector<Vector<Real>> stages(RKCoeff::DormandPrince8::extended_stages, Vector<Real>{0, 0});
		Vector<Real> workspace;

		CHECK_FALSE(interpolator.ready());
		CHECK_THROWS_AS(interpolator.interpolate(0), ODESolverError);
		CHECK(interpolator.setStep(system, 0, step, state, state, stages, workspace) == 3);
		CHECK(interpolator.ready());
		CHECK(system.calls == 3);
		for (Real fraction : {Real(0), Real(0.5), Real(1)}) {
			const auto value = interpolator.interpolate(fraction * step);
			REQUIRE(value.size() == state.size());
			CHECK(value[0] == state[0]);
			CHECK(value[1] == state[1]);
		}
		CHECK(system.calls == 3);
		interpolator.reset();
		CHECK_FALSE(interpolator.ready());
		CHECK_THROWS_AS(interpolator.interpolate(0), ODESolverError);
	}
}
