#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/interfaces/IFunction.h>
#include <mml/interfaces/IDynamicalSystem.h>
#include <mml/interfaces/IODESystemDAEWithEvents.h>
#include <mml/interfaces/ITensorField.h>
#include <mml/core/Fields/TensorFieldAdapters.h>
#include <mml/core/CoordTransf/CoordTransf2D.h>
#include <mml/core/Derivation/Jacobians.h>
#include <mml/algorithms/DAESolvers.h>
#endif
#include <catch2/catch_all.hpp>
#include <climits>
#include <functional>
#include <type_traits>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	class IdentityFunction : public IRealFunction {
	public:
		mutable std::vector<Real> arguments;
		Real operator()(Real x) const override { arguments.push_back(x); return x; }
	};
	class IdentityCurve : public IParametricCurve<2> {
	public:
		mutable std::vector<Real> arguments;
		Real getMinT() const override { return -std::numeric_limits<Real>::max(); }
		Real getMaxT() const override { return std::numeric_limits<Real>::max(); }
		VectorN<Real, 2> operator()(Real t) const override {
			arguments.push_back(t);
			return {t, -t};
		}
	};
	class ParametricLine : public IParametricCurveParametrized<2> {
		Vector<Real> _parameters{2, 3};
	public:
		int getNumParam() const override { return 2; }
		Real getParam(int i) const override { return _parameters.at(i); }
		void setParam(int i, Real x) override { _parameters.at(i) = x; }
		Real getMinT() const override { return -1; }
		Real getMaxT() const override { return 1; }
		VectorN<Real, 2> operator()(Real t) const override { return {t * _parameters[0], _parameters[1]}; }
		VectorN<Real, 2> operator()(Real t, const Vector<Real>& p) const override { return {t * p.at(0), p.at(1)}; }
	};

	template<class View>
	void CheckSingleArgument(const View& view) {
		if constexpr (std::is_invocable_v<const View&, Real>) {
			CHECK(view(Real(0.5))[0] == 1);
			CHECK(view(Real(0.5))[1] == 3);
		} else {
			FAIL("The parametrized abstract view hides its inherited curve call");
		}
	}

	class CombinedDAE : public IODESystemDAEWithJacobian, public IODESystemDAEParametrized {
	protected:
		Real _rate = 2;
	public:
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		int getNumParam() const override { return 1; }
		Real getParam(int) const override { return _rate; }
		void setParam(int, Real rate) override { _rate = rate; }
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& f) const override { f[0] = -y[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - _rate * x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -_rate; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 1; }
	};
	class EventDAE : public IODESystemDAEWithEvents, public IODESystemDAEParametrized {
	public:
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		int getNumParam() const override { return 0; }
		Real getParam(int) const override { throw IndexError("No parameters"); }
		void setParam(int, Real) override { throw IndexError("No parameters"); }
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& f) const override { f[0] = -y[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 1; }
		int getNumEvents() const override { return 1; }
		Real eventFunction(int, Real, const Vector<Real>& x, const Vector<Real>&) const override { return x[0] - Real(0.5); }
	};
	template<class Model>
	void CheckDAEViews(Model& model) {
		if constexpr (std::is_convertible_v<Model*, IODESystemDAE*>) {
			IODESystemDAE& base = model;
			IODESystemDAEWithJacobian& jacobian = model;
			IODESystemDAEParametrized& parameters = model;
			CHECK(&base == static_cast<IODESystemDAE*>(&jacobian));
			CHECK(&base == static_cast<IODESystemDAE*>(&parameters));
			CHECK(base.getDiffDim() == 1);
			CHECK(base.getAlgDim() == 1);
			Vector<Real> f(1), g(1);
			base.diffEqs(0, Vector<Real>{1}, Vector<Real>{2}, f);
			base.algConstraints(0, Vector<Real>{1}, Vector<Real>{2}, g);
			CHECK(f[0] == -2);
			Matrix<Real> fx(1, 1), fy(1, 1), gx(1, 1), gy(1, 1);
			jacobian.allJacobians(0, Vector<Real>{1}, Vector<Real>{2}, fx, fy, gx, gy);
			CHECK(fx(0, 0) == 0);
			CHECK(fy(0, 0) == -1);
			CHECK(gy(0, 0) == 1);
		} else {
			FAIL("DAE capabilities contain an ambiguous common DAE base");
		}
	}

	class NumericalSystem : public IDynamicalSystem {
	public:
		int dimension = 2;
		mutable int evaluations = 0;
		bool stopAtCallback = false;
		int getDim() const override { return dimension; }
		int getNumParam() const override { return 0; }
		Real getParam(int) const override { throw IndexError("No parameters"); }
		void setParam(int, Real) override { throw IndexError("No parameters"); }
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& f) const override {
			++evaluations;
			if (stopAtCallback) throw std::logic_error("Callback reached before input validation");
			f[0] = x[0] * x[0] + t * x[1];
			f[1] = 3 * x[0] - 2 * x[1];
		}
	};
	class CallbackSystem : public NumericalSystem {
	public:
		std::function<void(Real, const Vector<Real>&, Vector<Real>&)> callback;
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& f) const override {
			++evaluations;
			callback(t, x, f);
		}
	};
	class ParameterCounts : public IParametrized {
	public:
		int count = 0;
		mutable int reads = 0;
		int writes = 0;
		int getNumParam() const override { return count; }
		Real getParam(int) const override { ++reads; throw IndexError("No parameters"); }
		void setParam(int, Real) override { ++writes; throw IndexError("No parameters"); }
	};
	class SampleFunction : public IdentityFunction {
	public:
		bool throwAtSecond = false;
		Real operator()(Real x) const override {
			IdentityFunction::operator()(x);
			if (arguments.size() == 2) {
				if (throwAtSecond) throw std::overflow_error("Sample callback witness");
				return std::numeric_limits<Real>::quiet_NaN();
			}
			return x;
		}
	};
	class SampleCurve : public IdentityCurve {
	public:
		bool throwAtSecond = false;
		VectorN<Real, 2> operator()(Real t) const override {
			const auto value = IdentityCurve::operator()(t);
			if (arguments.size() == 2) {
				if (throwAtSecond) throw std::overflow_error("Curve callback witness");
				return {t, std::numeric_limits<Real>::infinity()};
			}
			return value;
		}
	};

	class Field2 : public ITensorField2<2> {
	public:
		bool wrongVariance = false;
		Field2(int contra, int co) : ITensorField2<2>(contra, co) {}
		Tensor2<2> operator()(const VectorN<Real, 2>&) const override {
			return wrongVariance ? Tensor2<2>(getNumContravar(), getNumCovar()) : Tensor2<2>(getNumCovar(), getNumContravar());
		}
		Real component(int, int, const VectorN<Real, 2>&) const override { return 0; }
	};
	class Field3 : public ITensorField3<2> {
	public:
		bool wrongVariance = false;
		Field3() = default;
		Field3(int contra, int co) : ITensorField3<2>(contra, co) {}
		Tensor3<2> operator()(const VectorN<Real, 2>&) const override {
			return wrongVariance ? Tensor3<2>(getNumContravar(), getNumCovar()) : Tensor3<2>(getNumCovar(), getNumContravar());
		}
		Real component(int, int, int, const VectorN<Real, 2>&) const override { return 0; }
	};
	class Field4 : public ITensorField4<2> {
	public:
		bool wrongVariance = false;
		Field4() = default;
		Field4(int contra, int co) : ITensorField4<2>(contra, co) {}
		Tensor4<2> operator()(const VectorN<Real, 2>&) const override {
			return wrongVariance ? Tensor4<2>(getNumContravar(), getNumCovar()) : Tensor4<2>(getNumCovar(), getNumContravar());
		}
		Real component(int, int, int, int, const VectorN<Real, 2>&) const override { return 0; }
	};
	class Field5 : public ITensorField5<2> {
	public:
		bool wrongVariance = false;
		Field5() = default;
		Field5(int contra, int co) : ITensorField5<2>(contra, co) {}
		Tensor5<2> operator()(const VectorN<Real, 2>&) const override {
			return wrongVariance ? Tensor5<2>(getNumContravar(), getNumCovar()) : Tensor5<2>(getNumCovar(), getNumContravar());
		}
		Real component(int, int, int, int, int, const VectorN<Real, 2>&) const override { return 0; }
	};
	class MixedField2 : public ITensorField2<2> {
	public:
		MixedField2() : ITensorField2<2>(1, 1) {}
		Tensor2<2> operator()(const VectorN<Real, 2>&) const override {
			Tensor2<2> value(1, 1);
			value.setContravar(0, true);
			value.setContravar(1, false);
			value(0, 1) = 7;
			return value;
		}
		Real component(int i, int j, const VectorN<Real, 2>& p) const override { return (*this)(p)(i, j); }
	};
	class MixedField3 : public ITensorField3<2> {
	public:
		MixedField3() : ITensorField3<2>(2, 1) {}
		Tensor3<2> operator()(const VectorN<Real, 2>&) const override {
			Tensor3<2> value(1, 2);
			value.setContravar(0, true);
			value.setContravar(1, false);
			value(1, 0, 1) = 7;
			return value;
		}
		Real component(int i, int j, int k, const VectorN<Real, 2>& p) const override { return (*this)(p)(i, j, k); }
	};
	template<class Field, int Rank>
	void CheckCounts() {
		for (const auto& counts : {std::pair{-1, Rank + 1}, std::pair{Rank + 1, -1},
			std::pair{0, Rank - 1}, std::pair{INT_MAX, INT_MAX}}) {
			try {
				Field invalid(counts.first, counts.second);
				FAIL("Expected invalid variance counts to fail");
			} catch (const TensorCovarContravarNumError& error) {
				CHECK(error.num_contravariant() == counts.first);
				CHECK(error.num_covariant() == counts.second);
			}
		}
		for (int contra = 0; contra <= Rank; ++contra) {
			const Field field(contra, Rank - contra);
			CHECK(field.getNumContravar() == contra);
			CHECK(field.getNumCovar() == Rank - contra);
		}
	}
	void CheckDecayConsumer(IODESystemDAEWithJacobian& model, Real rate) {
		DAESolverConfig config;
		config.step_size = Real(1) / 8;
		const auto result = SolveDAEBackwardEuler(model, 0, Vector<Real>{1}, Vector<Real>{rate}, 1, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.total_steps == 8);
		// Independent BE recurrence after eight exact binary time steps.
		const Real expected = std::pow(1 + rate / 8, Real(-8));
		CHECK_THAT(result.solution.getXValuesAtEnd()[0], WithinAbs(expected, Precision::DefaultToleranceStrict));
		CHECK_THAT(result.solution.getYValuesAtEnd()[0], WithinAbs(rate * expected, Precision::DefaultToleranceStrict));
		CHECK(result.max_constraint_violation <= Precision::DefaultToleranceStrict);
	}
}

TEST_CASE("Interface contracts scalar sampling retains finite extreme endpoints and midpoint",
	"[interfaces][interface-contract]") {
	IdentityFunction function;
	const Real max = std::numeric_limits<Real>::max();
	for (const auto& bounds : {std::pair{-max, max}, std::pair{max, -max}, std::pair{max, max}}) {
		Vector<Real> x, y;
		function.getValues(bounds.first, bounds.second, 5, x, y);
		CHECK(x[0] == bounds.first);
		CHECK(x[4] == bounds.second);
		CHECK(x[2] == (bounds.first == bounds.second ? max : 0));
		for (int i = 0; i < 5; ++i) { CHECK(std::isfinite(x[i])); CHECK(y[i] == x[i]); }
	}
}
TEST_CASE("Interface contracts curve traces preserve extreme ranges and legacy empty counts",
	"[interfaces][interface-contract]") {
	IdentityCurve curve;
	const Real max = std::numeric_limits<Real>::max();
	for (const auto& bounds : {std::pair{-max, max}, std::pair{max, -max}}) {
		const auto trace = curve.getTrace(bounds.first, bounds.second, 5);
		REQUIRE(trace.size() == 5);
		CHECK(trace[0][0] == bounds.first);
		CHECK(trace[4][0] == bounds.second);
		CHECK(trace[2][0] == 0);
		for (const auto& point : trace) { CHECK(std::isfinite(point[0])); CHECK(point[1] == -point[0]); }
	}
	for (int count : {INT_MIN, -1, 0}) CHECK(curve.getTrace(0, 1, count).empty());
}
TEST_CASE("Interface contracts parameter dimension metadata states expected before actual",
	"[interfaces][interface-contract]") {
	ParametricLine curve;
	for (const Vector<Real>& parameters : {Vector<Real>{1}, Vector<Real>{1, 2, 3}}) {
		try {
			curve.setParams(parameters);
			FAIL("Expected a parameter dimension error");
		} catch (const VectorDimensionError& error) {
			CHECK(error.expected() == 2);
			CHECK(error.actual() == parameters.size());
		}
		CHECK(curve.getParams() == Vector<Real>{2, 3});
	}
}
TEST_CASE("Interface contracts parametrized curves remain callable through every abstract view",
	"[interfaces][interface-contract]") {
	ParametricLine curve;
	CheckSingleArgument(static_cast<const IParametricCurveParametrized<2>&>(curve));
	CheckSingleArgument(static_cast<const IParametricCurve<2>&>(curve));
	CheckSingleArgument(static_cast<const IRealToVectorFunction<2>&>(curve));
	CheckSingleArgument(static_cast<const IFunction<VectorN<Real, 2>, Real>&>(curve));
	const IParametricCurveParametrized<2>& view = curve;
	CHECK(view(Real(0.5), Vector<Real>{4, 5}) == VectorN<Real, 2>{2, 5});
}
TEST_CASE("Interface contracts DAE Jacobian parameter and event capabilities share one base",
	"[interfaces][DAE][interface-contract]") {
	CombinedDAE model;
	EventDAE events;
	CheckDAEViews(model);
	CheckDAEViews(events);
}
TEST_CASE("Interface contracts numerical Jacobian rejects short state before any callback",
	"[interfaces][jacobian][interface-contract]") {
	NumericalSystem system;
	system.stopAtCallback = true;
	Matrix<Real> jac(1, 1, {99});
	// The old code reaches the throwing first callback, never its later unchecked short-state access.
	CHECK_THROWS_AS(system.jacobian(0, Vector<Real>{1}, jac), VectorDimensionError);
	CHECK(system.evaluations == 0);
	CHECK(jac.rows() == 1);
	CHECK(jac(0, 0) == 99);
}
TEST_CASE("Interface contracts tensor fields reject impossible variance counts for every rank",
	"[interfaces][TensorField][interface-contract]") {
	CheckCounts<Field2, 2>();
	CheckCounts<Field3, 3>();
	CheckCounts<Field4, 4>();
	CheckCounts<Field5, 5>();
	CHECK(Field3().getNumCovar() == 3);
	CHECK(Field4().getNumCovar() == 4);
	CHECK(Field5().getNumCovar() == 5);
}

TEST_CASE("Interface contracts sampling validates used inputs and preserves outputs on callback failure",
	"[interfaces][interface-contract]") {
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const Real inf = std::numeric_limits<Real>::infinity();
	IdentityFunction function;
	Vector<Real> x{7, 8}, y{9};
	for (int count : {INT_MIN, -1, 0, 1})
		CHECK_THROWS_AS(function.getValues(0, 1, count, x, y), ArgumentError);
	for (Real invalid : {nan, inf, -inf}) {
		CHECK_THROWS_AS(function.getValues(invalid, 1, 2, x, y), NumericInputError);
		CHECK_THROWS_AS(function.getValues(0, invalid, 2, x, y), NumericInputError);
	}
	CHECK(function.arguments.empty());
	for (bool throws : {false, true}) {
		SampleFunction failing;
		failing.throwAtSecond = throws;
		if (throws) CHECK_THROWS_AS(failing.getValues(0, 1, 5, x, y), std::overflow_error);
		else CHECK_THROWS_AS(failing.getValues(0, 1, 5, x, y), NumericInputError);
		CHECK(failing.arguments.size() == 2);
		CHECK(x == Vector<Real>{7, 8});
		CHECK(y == Vector<Real>{9});
		SampleCurve curve;
		curve.throwAtSecond = throws;
		if (throws) CHECK_THROWS_AS(curve.getTrace(0, 1, 5), std::overflow_error);
		else CHECK_THROWS_AS(curve.getTrace(0, 1, 5), NumericInputError);
		CHECK(curve.arguments.size() == 2);
	}
	IdentityCurve curve;
	for (int count : {INT_MIN, -1, 0}) CHECK(curve.getTrace(nan, inf, count).empty());
	CHECK(curve.arguments.empty());
	CHECK(curve.getTrace(2, nan, 1)[0] == VectorN<Real, 2>{2, -2});
	for (Real invalid : {nan, inf, -inf}) {
		CHECK_THROWS_AS(curve.getTrace(invalid, 0, 1), NumericInputError);
		CHECK_THROWS_AS(curve.getTrace(0, invalid, 2), NumericInputError);
	}
	CHECK(curve.arguments == std::vector<Real>{2});
	const Real adjacent = std::nextafter(Real(1), Real(2));
	for (const auto& bounds : {std::pair{Real(1), adjacent}, std::pair{adjacent, Real(1)}, std::pair{Real(1), Real(1)}}) {
		const auto trace = curve.getTrace(bounds.first, bounds.second, 5);
		function.getValues(bounds.first, bounds.second, 5, x, y);
		CHECK(trace.front()[0] == bounds.first);
		CHECK(trace.back()[0] == bounds.second);
		for (int i = 0; i < 5; ++i) {
			CHECK(trace[i][0] >= std::min(bounds.first, bounds.second));
			CHECK(trace[i][0] <= std::max(bounds.first, bounds.second));
			CHECK(trace[i][0] == x[i]);
			CHECK(y[i] == x[i]);
		}
	}
}

TEST_CASE("Interface contracts bulk parameters reject negative counts and permit empty parameter sets",
	"[interfaces][interface-contract]") {
	ParameterCounts parameters;
	CHECK(parameters.getParams().size() == 0);
	CHECK_NOTHROW(parameters.setParams(Vector<Real>{}));
	CHECK_THROWS_AS(parameters.setParams(Vector<Real>{1}), VectorDimensionError);
	parameters.count = -1;
	CHECK_THROWS_AS(parameters.getParams(), ArgumentError);
	CHECK_THROWS_AS(parameters.setParams(Vector<Real>{}), ArgumentError);
	CHECK(parameters.reads == 0);
	CHECK(parameters.writes == 0);
}

TEST_CASE("Interface contracts numerical Jacobians preflight all input and stencil failures",
	"[interfaces][jacobian][interface-contract]") {
	NumericalSystem system;
	Matrix<Real> jac(1, 1, {99});
	const Real inf = std::numeric_limits<Real>::infinity(), nan = std::numeric_limits<Real>::quiet_NaN();
	for (const auto& state : {Vector<Real>{}, Vector<Real>{1}, Vector<Real>{1, 2, 3}}) {
		try {
			system.jacobian(0, state, jac);
			FAIL("Expected state dimension failure");
		} catch (const VectorDimensionError& error) {
			CHECK(error.expected() == 2);
			CHECK(error.actual() == state.size());
		}
	}
	for (int dim : {0, -1}) {
		system.dimension = dim;
		CHECK_THROWS_AS(system.jacobian(0, Vector<Real>{}, jac), ArgumentError);
	}
	system.dimension = 2;
	for (Real invalid : {inf, -inf, nan}) {
		CHECK_THROWS_AS(system.jacobian(invalid, Vector<Real>{1, 2}, jac), ArgumentError);
		CHECK_THROWS_AS(system.jacobian(0, Vector<Real>{1, invalid}, jac), ArgumentError);
	}
	CHECK_THROWS_AS(system.jacobian(0, Vector<Real>{1, std::numeric_limits<Real>::max()}, jac), NumericalMethodError);
	CHECK(system.evaluations == 0);
	CHECK(jac.rows() == 1);
	CHECK(jac.cols() == 1);
	CHECK(jac(0, 0) == 99);
}

TEST_CASE("Interface contracts numerical Jacobians reject malformed callbacks without publishing partial results",
	"[interfaces][jacobian][interface-contract]") {
	const Vector<Real> state{2, 3};
	for (int mode = 0; mode < 5; ++mode) {
		CallbackSystem system;
		Matrix<Real> jac(1, 1, {99});
		system.callback = [mode, &system](Real, const Vector<Real>&, Vector<Real>& f) {
			f[0] = 1;
			f[1] = 2;
			if (system.evaluations != 3) return; // Fail only after a complete first column.
			if (mode == 0) f.resize(1);
			if (mode == 1) f.resize(3);
			if (mode == 2) f[1] = std::numeric_limits<Real>::infinity();
			if (mode == 3) f[1] = std::numeric_limits<Real>::quiet_NaN();
			if (mode == 4) throw std::overflow_error("Derivative callback witness");
		};
		if (mode < 2) CHECK_THROWS_AS(system.jacobian(0, state, jac), VectorDimensionError);
		else if (mode < 4) CHECK_THROWS_AS(system.jacobian(0, state, jac), NumericalMethodError);
		else CHECK_THROWS_AS(system.jacobian(0, state, jac), std::overflow_error);
		CHECK(system.evaluations == 3);
		CHECK(jac.rows() == 1);
		CHECK(jac.cols() == 1);
		CHECK(jac(0, 0) == 99);
		CHECK(state == Vector<Real>{2, 3});
	}
	CallbackSystem omitted;
	omitted.callback = [](Real, const Vector<Real>&, Vector<Real>& f) { f[0] = 1; };
	Matrix<Real> jac(1, 1, {99});
	CHECK_THROWS_AS(omitted.jacobian(0, state, jac), NumericalMethodError);
	CHECK(omitted.evaluations == 1);
	CHECK(jac(0, 0) == 99);
}

TEST_CASE("Interface contracts numerical Jacobians match independent mathematics and canonical ODE differentiation",
	"[interfaces][jacobian][interface-contract]") {
	NumericalSystem system;
	const Vector<Real> state{2, 3};
	Matrix<Real> jac;
	system.jacobian(Real(0.5), state, jac);
	CHECK(system.evaluations == 4);
	const Matrix<Real> expected(2, 2, {4, Real(0.5), 3, -2});
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			CHECK_THAT(jac(i, j), WithinAbs(expected(i, j), Precision::CbrtEps));
	const std::function<void(Real, const Vector<Real>&, Vector<Real>&)> callback =
		[&system](Real t, const Vector<Real>& x, Vector<Real>& f) { system.derivs(t, x, f); };
	const auto canonical = Derivation::CalcJacobianODE(callback, Real(0.5), state);
	CHECK(system.evaluations == 8);
	CHECK(canonical == jac);
	CHECK(state == Vector<Real>{2, 3});

	CallbackSystem huge;
	huge.dimension = 1;
	const Real max = std::numeric_limits<Real>::max();
	huge.callback = [max](Real, const Vector<Real>& x, Vector<Real>& f) { f[0] = std::copysign(max, x[0]); };
	// Raw subtraction max-(-max) overflows, but division by a large central span is finite.
	const auto scaled = Derivation::CalcJacobianODE(
		[&huge](Real t, const Vector<Real>& x, Vector<Real>& f) { huge.derivs(t, x, f); },
		0, Vector<Real>{0}, max / 2);
	CHECK_THAT(scaled(0, 0), WithinAbs(2, Precision::RoundoffZeroThreshold));
	CHECK(huge.evaluations == 2);
	// Nonfinite final quotient is not published either.
	Matrix<Real> preserved(1, 1, {99});
	CHECK_THROWS_AS(huge.jacobian(0, Vector<Real>{0}, preserved), NumericalMethodError);
	CHECK(preserved(0, 0) == 99);
}

TEST_CASE("Interface contracts composed index one DAE capabilities work in a real solver",
	"[interfaces][DAE][interface-contract]") {
	CombinedDAE model;
	CheckDecayConsumer(model, 2);
	IODESystemDAEParametrized& parameters = model;
	parameters.setParams(Vector<Real>{3});
	CHECK(parameters.getParams() == Vector<Real>{3});
	Matrix<Real> fx(1, 1), fy(1, 1), gx(1, 1), gy(1, 1);
	model.allJacobians(0, Vector<Real>{1}, Vector<Real>{3}, fx, fy, gx, gy);
	CHECK(gx(0, 0) == -3);
	CHECK(gy(0, 0) == 1);
	Vector<Real> constraints(1), derivatives(1);
	model.algConstraints(0, Vector<Real>{1}, Vector<Real>{3}, constraints);
	model.diffEqs(0, Vector<Real>{1}, Vector<Real>{3}, derivatives);
	CHECK(constraints[0] == 0);
	CHECK(derivatives[0] == -3);
	CheckDecayConsumer(model, 3);
	EventDAE events;
	const IODESystemDAEWithEvents& eventView = events;
	CHECK(eventView.eventFunction(0, 0, Vector<Real>{1}, Vector<Real>{1}) == Real(0.5));
	CHECK(static_cast<const IODESystemDAE*>(&eventView) == static_cast<const IODESystemDAE*>(&events));
	CheckDecayConsumer(events, 1);
}

TEST_CASE("Interface contracts tensor adapters reject returned variance inconsistent with field declaration",
	"[interfaces][TensorField][interface-contract]") {
	CoordTransfCart2DRotation identity(0);
	Field2 field2(2, 0);
	Field3 field3(3, 0);
	Field4 field4(4, 0);
	Field5 field5(5, 0);
	field2.wrongVariance = field3.wrongVariance = field4.wrongVariance = field5.wrongVariance = true;
	TransformedTensorField2<Vector2Cartesian, Vector2Cartesian, 2> adapted2(field2, identity);
	TransformedTensorField3<Vector2Cartesian, Vector2Cartesian, 2> adapted3(field3, identity);
	TransformedTensorField4<Vector2Cartesian, Vector2Cartesian, 2> adapted4(field4, identity);
	TransformedTensorField5<Vector2Cartesian, Vector2Cartesian, 2> adapted5(field5, identity);
	const VectorN<Real, 2> point{1, 2};
	CHECK_THROWS_AS(adapted2(point), TensorCovarContravarArithmeticError);
	CHECK_THROWS_AS(adapted3.component(0, 0, 0, point), TensorCovarContravarArithmeticError);
	CHECK_THROWS_AS(adapted4(point), TensorCovarContravarArithmeticError);
	CHECK_THROWS_AS(adapted5.component(0, 0, 0, 0, 0, point), TensorCovarContravarArithmeticError);
}

TEST_CASE("Interface contracts tensor field totals do not impose a contiguous variance layout",
	"[interfaces][TensorField][interface-contract]") {
	MixedField2 field2;
	MixedField3 field3;
	CoordTransfCart2DRotation identity(0);
	TransformedTensorField2<Vector2Cartesian, Vector2Cartesian, 2> adapted2(field2, identity);
	TransformedTensorField3<Vector2Cartesian, Vector2Cartesian, 2> adapted3(field3, identity);
	const VectorN<Real, 2> point{1, 2};
	const auto result2 = adapted2(point);
	const auto result3 = adapted3(point);
	for (int slot = 0; slot < 2; ++slot) CHECK(result2.isContravar(slot) == (slot == 0));
	for (int slot = 0; slot < 3; ++slot) CHECK(result3.isContravar(slot) == (slot != 1));
	CHECK(result2.numContravar() == 1);
	CHECK(result3.numContravar() == 2);
	CHECK_THAT(adapted2.component(0, 1, point), WithinAbs(7, Precision::DefaultToleranceStrict));
	CHECK_THAT(adapted3.component(1, 0, 1, point), WithinAbs(7, Precision::DefaultToleranceStrict));
}
