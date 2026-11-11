#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/DynamicalSystem.h>
#endif
#include <catch2/catch_all.hpp>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

using namespace MML;
using namespace MML::Systems;

namespace {
	const std::array<Real, 3> nonfiniteParameters = {
		std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()
	};

	template<typename Model, std::size_t... I>
	Model ConstructModel(const Vector<Real>& values, std::index_sequence<I...>) {
		return Model(values[static_cast<int>(I)]...);
	}

	template<typename Model, int P>
	void CheckFiniteConstruction() {
		const auto original = Model().getParams();
		for (int index = 0; index < P; ++index) {
			for (Real invalid : nonfiniteParameters) {
				CAPTURE(index, invalid);
				auto candidate = original;
				candidate[index] = invalid;
				CHECK_THROWS_AS((ConstructModel<Model>(candidate, std::make_index_sequence<P>{})), ArgumentError);
			}
		}
	}

	template<typename Model>
	void CheckFiniteUpdates() {
		Model model;
		IParametrized& parameters = model;
		const auto original = parameters.getParams();
		for (int index = 0; index < parameters.getNumParam(); ++index) {
			for (Real invalid : nonfiniteParameters) {
				CAPTURE(index, invalid);
				CHECK_THROWS_AS(parameters.setParam(index, invalid), ArgumentError);
				CHECK(parameters.getParams() == original);
				auto candidate = original;
				candidate[0] += 1;
				candidate[index] = invalid;
				CHECK_THROWS_AS(parameters.setParams(candidate), ArgumentError);
				CHECK(parameters.getParams() == original);
			}
		}
		CHECK_THROWS_AS(parameters.setParams(Vector<Real>(original.size() + 1)), VectorDimensionError);
		CHECK(parameters.getParams() == original);
		CHECK_THROWS_AS(parameters.setParams(Vector<Real>(original.size() - 1)), VectorDimensionError);
		CHECK(parameters.getParams() == original);
	}

	class EmptyParameterFlow : public DynamicalSystemBase<1, 0> {
	public:
		void derivs(Real, const Vector<Real>&, Vector<Real>& result) const override { result[0] = 0; }
	};

	class InterfaceMetadataFlow : public IDynamicalSystem {
	public:
		int getDim() const override { return 1; }
		int getNumParam() const override { return 1; }
		Real getParam(int) const override { return 0; }
		void setParam(int, Real) override {}
		void derivs(Real, const Vector<Real>&, Vector<Real>& result) const override { result[0] = 0; }
	};

	class InterfaceMetadataMap : public MML::IDiscreteMap<1> {
	public:
		int getNumParam() const override { return 0; }
		Real getParam(int) const override { throw IndexError("fixture has no parameters"); }
		void setParam(int, Real) override { throw IndexError("fixture has no parameters"); }
		void map(const Vector<Real>& state, Vector<Real>& result) const override { result = state; }
		void jacobian(const Vector<Real>&, Matrix<Real>& result) const override { result = Matrix<Real>::Identity(1); }
	};

	template<typename Model>
	void CheckContinuousMetadata() {
		Model model;
		const IDynamicalSystem& interface = model;
		for (int index : {-1, interface.getDim(), std::numeric_limits<int>::max()})
			CHECK_THROWS_AS(interface.getStateName(index), IndexError);
		for (int index : {-1, interface.getNumParam(), std::numeric_limits<int>::max()}) {
			CHECK_THROWS_AS(interface.getParamName(index), IndexError);
			CHECK_THROWS_AS(interface.getParamRange(index), IndexError);
			CHECK_THROWS_AS(interface.getParam(index), IndexError);
			CHECK_THROWS_AS(model.setParam(index, 0), IndexError);
		}
		for (int index = 0; index < interface.getDim(); ++index)
			CHECK_FALSE(interface.getStateName(index).empty());
		for (int index = 0; index < interface.getNumParam(); ++index) {
			CHECK_FALSE(interface.getParamName(index).empty());
			CHECK(interface.getParamRange(index).first <= interface.getParamRange(index).second);
		}
	}

	template<typename Model>
	void CheckMapMetadata() {
		Model model;
		for (int index : {-1, model.getDim(), std::numeric_limits<int>::max()})
			CHECK_THROWS_AS(model.getStateName(index), IndexError);
		for (int index : {-1, model.getNumParam(), std::numeric_limits<int>::max()}) {
			CHECK_THROWS_AS(model.getParamName(index), IndexError);
			CHECK_THROWS_AS(model.getParam(index), IndexError);
			CHECK_THROWS_AS(model.setParam(index, 0), IndexError);
		}
		for (int index = 0; index < model.getDim(); ++index)
			CHECK_FALSE(model.getStateName(index).empty());
		for (int index = 0; index < model.getNumParam(); ++index)
			CHECK_FALSE(model.getParamName(index).empty());
	}
}

TEST_CASE("Continuous model constructors reject every nonfinite parameter", "[ModelDomainContracts][ModelDomainBaseline]") {
	CheckFiniteConstruction<LorenzSystem, 3>();
	CheckFiniteConstruction<RosslerSystem, 3>();
	CheckFiniteConstruction<VanDerPolSystem, 1>();
	CheckFiniteConstruction<DuffingSystem, 5>();
	CheckFiniteConstruction<ChuaCircuit, 4>();
	CheckFiniteConstruction<HenonHeilesSystem, 1>();
	CheckFiniteConstruction<DoublePendulumSystem, 5>();
}

TEST_CASE("Continuous model parameter updates reject nonfinite values without mutation", "[ModelDomainContracts][ModelDomainBaseline]") {
	CheckFiniteUpdates<LorenzSystem>();
	CheckFiniteUpdates<RosslerSystem>();
	CheckFiniteUpdates<VanDerPolSystem>();
	CheckFiniteUpdates<DuffingSystem>();
	CheckFiniteUpdates<ChuaCircuit>();
	CheckFiniteUpdates<HenonHeilesSystem>();
	CheckFiniteUpdates<DoublePendulumSystem>();
}

TEST_CASE("Double pendulum construction rejects nonpositive masses and lengths", "[ModelDomainContracts][ModelDomainBaseline]") {
	for (int index = 0; index < 4; ++index) {
		for (Real invalid : {Real(0), Real(-1)}) {
			auto values = DoublePendulumSystem().getParams();
			values[index] = invalid;
			CHECK_THROWS_AS((ConstructModel<DoublePendulumSystem>(values, std::make_index_sequence<5>{})), ArgumentError);
		}
	}
}

TEST_CASE("Discrete models reject nonfinite construction and updates atomically", "[ModelDomainContracts][ModelDomainBaseline]") {
	CheckFiniteConstruction<LogisticMap, 1>();
	CheckFiniteConstruction<HenonMap, 2>();
	CheckFiniteConstruction<StandardMap, 1>();
	CheckFiniteConstruction<TentMap, 1>();
	CheckFiniteUpdates<LogisticMap>();
	CheckFiniteUpdates<HenonMap>();
	CheckFiniteUpdates<StandardMap>();
	CheckFiniteUpdates<TentMap>();
}

TEST_CASE("Continuous models guard all state and parameter metadata indices", "[ModelDomainContracts]") {
	CheckContinuousMetadata<LorenzSystem>();
	CheckContinuousMetadata<RosslerSystem>();
	CheckContinuousMetadata<VanDerPolSystem>();
	CheckContinuousMetadata<DuffingSystem>();
	CheckContinuousMetadata<ChuaCircuit>();
	CheckContinuousMetadata<HenonHeilesSystem>();
	CheckContinuousMetadata<DoublePendulumSystem>();
	CheckContinuousMetadata<EmptyParameterFlow>();
}

TEST_CASE("Default dynamical interface metadata validates indices", "[ModelDomainContracts]") {
	InterfaceMetadataFlow model;
	const IDynamicalSystem& interface = model;
	CHECK(interface.getStateName(0) == "x0");
	CHECK(interface.getParamName(0) == "p0");
	CHECK(interface.getParamRange(0).first == -Precision::UnboundedIntervalSamplingLimit);
	CHECK(interface.getParamRange(0).second == Precision::UnboundedIntervalSamplingLimit);
	for (int index : {-1, 1, std::numeric_limits<int>::max()}) {
		CHECK_THROWS_AS(interface.getStateName(index), IndexError);
		CHECK_THROWS_AS(interface.getParamName(index), IndexError);
		CHECK_THROWS_AS(interface.getParamRange(index), IndexError);
	}
}

TEST_CASE("Discrete models and interface defaults guard metadata indices", "[ModelDomainContracts]") {
	CheckMapMetadata<LogisticMap>();
	CheckMapMetadata<HenonMap>();
	CheckMapMetadata<StandardMap>();
	CheckMapMetadata<TentMap>();
	CheckMapMetadata<InterfaceMetadataMap>();
	CHECK(InterfaceMetadataMap().getStateName(0) == "x0");
}

TEST_CASE("Double pendulum scalar and bulk domain failures preserve all parameters", "[ModelDomainContracts]") {
	DoublePendulumSystem model;
	IParametrized& interface = model;
	const auto original = interface.getParams();
	for (int index = 0; index < 4; ++index) {
		for (Real invalid : {Real(0), Real(-1)}) {
			CHECK_THROWS_AS(interface.setParam(index, invalid), ArgumentError);
			CHECK(interface.getParams() == original);
			auto candidate = original;
			candidate[0] = 2;
			candidate[index] = invalid;
			CHECK_THROWS_AS(interface.setParams(candidate), ArgumentError);
			CHECK(interface.getParams() == original);
		}
	}
}

TEST_CASE("Double pendulum accepts signed and zero gravity and positive values outside suggested ranges", "[ModelDomainContracts]") {
	DoublePendulumSystem model(20, 30, 10, 20, -4);
	const Vector<Real> state{Constants::PI / 2, Constants::PI / 2, 0, 0};
	Vector<Real> derivative(4);
	for (Real gravity : {Real(-4), Real(0), Real(40)}) {
		model.setParam(4, gravity);
		model.derivs(0, state, derivative);
		CHECK(IsWithinAbs(derivative[2], -gravity / 10, Precision::LooseTolerance));
		CHECK(IsWithinAbs(derivative[3], Real(0), Precision::LooseTolerance));
		CHECK(std::isfinite(model.getEnergy(state)));
	}
	model.setParams(Vector<Real>{2, 3, 4, 5, -6});
	CHECK(model.getParams() == Vector<Real>{2, 3, 4, 5, -6});
}

TEST_CASE("Double pendulum aligned rods retain a nonsingular denominator at large mass ratio", "[ModelDomainContracts]") {
	const Real largeMass = std::ldexp(Real(1), std::numeric_limits<Real>::digits + 1);
	DoublePendulumSystem model(1, largeMass, 1, 1, 0);
	Vector<Real> derivative(4);
	model.derivs(0, Vector<Real>{0, 0, 0, 0}, derivative);
	CHECK(derivative == Vector<Real>{0, 0, 0, 0});
}

TEST_CASE("Double pendulum positive domain has no arbitrary small-value cutoff", "[ModelDomainContracts]") {
	const Real small = std::numeric_limits<Real>::min();
	DoublePendulumSystem model(small, small, small, small, 0);
	CHECK(model.getParams() == Vector<Real>{small, small, small, small, 0});
	model.setParam(0, std::numeric_limits<Real>::denorm_min());
	CHECK(model.getParam(0) == std::numeric_limits<Real>::denorm_min());
}

TEST_CASE("Finite coefficients are not bounded by an undocumented magnitude limit", "[ModelDomainContracts]") {
	const Real large = std::numeric_limits<Real>::max();
	LorenzSystem lorenz(large, -large, large);
	HenonMap map(large, -large);
	CHECK(lorenz.getParams() == Vector<Real>{large, -large, large});
	CHECK(map.getParams() == Vector<Real>{large, -large});
	lorenz.setParams(Vector<Real>{-large, large, -large});
	map.setParams(Vector<Real>{-large, large});
	CHECK(lorenz.getParam(0) == -large);
	CHECK(map.getParam(0) == -large);
}

TEST_CASE("Suggested continuous model parameter ranges remain advisory", "[ModelDomainContracts]") {
	LorenzSystem lorenz(-2, -3, -4);
	RosslerSystem rossler(-2, -3, -4);
	VanDerPolSystem vanDerPol(-2);
	DuffingSystem duffing(-2, 4, -3, -4, -5);
	ChuaCircuit chua(-2, -3, 4, 5);
	HenonHeilesSystem henonHeiles(-2);
	for (IDynamicalSystem* model : std::array<IDynamicalSystem*, 6>{&lorenz, &rossler, &vanDerPol, &duffing, &chua, &henonHeiles}) {
		auto candidate = model->getParams();
		for (int index = 0; index < candidate.size(); ++index)
			candidate[index] = model->getParamRange(index).second + 1;
		model->setParams(candidate);
		CHECK(model->getParams() == candidate);
		for (int index = 0; index < candidate.size(); ++index)
			model->setParam(index, model->getParamRange(index).first - 1);
	}
}

TEST_CASE("Suggested discrete map regimes do not restrict finite mathematical parameters", "[ModelDomainContracts]") {
	LogisticMap logistic(-2);
	HenonMap henon(-2, 0);
	StandardMap standard(-2);
	TentMap tent(-2);
	for (IParametrized* model : std::array<IParametrized*, 4>{&logistic, &henon, &standard, &tent}) {
		for (int index = 0; index < model->getNumParam(); ++index) {
			model->setParam(index, -3);
			CHECK(model->getParam(index) == -3);
		}
	}
	CHECK(henon.jacobianDeterminant() == 3);
}

TEST_CASE("Zero-parameter continuous models retain valid empty bulk access", "[ModelDomainContracts]") {
	EmptyParameterFlow model;
	CHECK(model.getDim() == 1);
	CHECK(model.getNumParam() == 0);
	CHECK(model.getParams().empty());
	CHECK_NOTHROW(model.setParams(Vector<Real>{}));
	CHECK_THROWS_AS(model.setParams(Vector<Real>{1}), VectorDimensionError);
}

TEST_CASE("Energy invariant access validates its index before evaluation", "[ModelDomainContracts]") {
	HenonHeilesSystem henon;
	DoublePendulumSystem pendulum;
	for (const IDynamicalSystem* model : std::array<const IDynamicalSystem*, 2>{&henon, &pendulum}) {
		CHECK(model->getInvariantName(0) == "Energy");
		CHECK(std::isfinite(model->computeInvariant(0, model->getDefaultInitialCondition())));
		for (int index : {-1, 1, std::numeric_limits<int>::max()}) {
			CHECK_THROWS_AS(model->getInvariantName(index), IndexError);
			CHECK_THROWS_AS(model->computeInvariant(index, Vector<Real>{}), IndexError);
		}
	}
	InterfaceMetadataFlow withoutInvariants;
	CHECK_THROWS_AS(withoutInvariants.getInvariantName(0), IndexError);
	CHECK_THROWS_AS(withoutInvariants.computeInvariant(0, Vector<Real>{}), IndexError);
}
