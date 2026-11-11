#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/MetricTensor.h>
#include <mml/core/DifferentialGeometry/MetricAdapters.h>
#include <mml/core/DifferentialGeometry/Frames.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/CoordTransf/CoordTransf3D.h>
#endif
#include <limits>
#include <type_traits>

using namespace MML;
namespace MetricFieldContracts {
	Real PolynomialFirstDifferenceBudget(int order) {
		// Only the first two stencil orders have truncation error for these polynomials.
		if (order == 1) return 2 * Derivation::NDer1_h;
		if (order == 2) return 2 * Derivation::NDer2_h * Derivation::NDer2_h;
		return Precision::RoundoffZeroThreshold;
	}
	Real FrameDifferenceBudget(Real angularStep) {
		// Twice the five-point second-derivative stencil's coefficient L1 norm is 32/3.
		return std::max(Precision::LooseTolerance, (Real(32) / 3) * Precision::Eps / (angularStep * angularStep)
			+ POW4(angularStep));
	}
	struct StoredMetric : MetricTensorField<2> {
		bool inverse;
		bool asymmetric = false;
		explicit StoredMetric(bool inv = false) : MetricTensorField<2>(inv ? 2 : 0, inv ? 0 : 2), inverse(inv) {}
		Real component(int i, int j, const VectorN<Real, 2>&) const override {
			if (i != j) return asymmetric && i == 0 ? Real(1) : Real(0);
			return inverse ? Real(1) / (i == 0 ? 2 : 3) : Real(i == 0 ? 2 : 3);
		}
	};
	struct Embed : CoordTransf<VectorN<Real, 2>, VectorN<Real, 3>> {
		ScalarFunction<2> u{+[](const VectorN<Real, 2>& p) { return p[0]; }};
		ScalarFunction<2> v{+[](const VectorN<Real, 2>& p) { return p[1]; }};
		VectorN<Real, 3> transf(const VectorN<Real, 2>& p) const override { return {p[0], p[1], p[1]}; }
		const IScalarFunction<2>& coordTransfFunc(int i) const override { return i == 0 ? u : v; }
		JacobianType jacobian(const VectorN<Real, 2>&) const override { return {1, 0, 0, 1, 0, 1}; }
	};
	struct Reduce : CoordTransf<VectorN<Real, 3>, VectorN<Real, 2>> {
		ScalarFunction<3> u{+[](const VectorN<Real,3>& p) { return p[0]; }};
		ScalarFunction<3> v{+[](const VectorN<Real,3>& p) { return p[1]; }};
		VectorN<Real,2> transf(const VectorN<Real,3>& p) const override { return {p[0],p[1]}; }
		const IScalarFunction<3>& coordTransfFunc(int i) const override { return i == 0 ? u : v; }
		JacobianType jacobian(const VectorN<Real,3>&) const override { return {1,0,0,0,1,0}; }
	};
	struct Quadratic : IScalarFunction<3> {
		mutable int calls = 0;
		int fail_after = 0;
		bool poison = false;
		Real operator()(const VectorN<Real, 3>& p) const override {
			++calls;
			if (fail_after && calls >= fail_after) throw DomainError("bounded callback failure");
			return poison ? std::numeric_limits<Real>::quiet_NaN() : p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
		}
	};
	struct ProvidedQuadratic : Quadratic, IScalarDerivativeProvider<3> {
		mutable int derivatives = 0;
		Real EvaluatePartial(const VectorN<Real, 3>& p, int a) const override { ++derivatives; return 2*p[a]; }
		Real EvaluateSecondPartial(const VectorN<Real, 3>&, int a, int b) const override { ++derivatives; return a == b ? Real(2) : Real(0); }
	};
	struct LinearField : IVectorFunction<3> {
		mutable int calls = 0;
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& p) const override { ++calls; return {p[0], 2*p[1], 3*p[2]}; }
	};
	struct ProvidedLinear : LinearField, IVectorDerivativeProvider<3> {
		mutable int derivatives = 0;
		Real EvaluateComponentPartial(const VectorN<Real, 3>&, int c, int a) const override { ++derivatives; return c == a ? Real(c + 1) : 0; }
	};
	struct ParameterizedCircle : IParametricCurve<2>, ICurveDerivativeProvider<2> {
		Real frequency = 1, radius = 1;
		Real getMinT() const override { return -Constants::PosInf; }
		Real getMaxT() const override { return Constants::PosInf; }
		VectorN<Real, 2> operator()(Real t) const override { return {radius*std::cos(frequency*t), radius*std::sin(frequency*t)}; }
		VectorN<Real, 2> EvaluateDerivative(Real t, int order) const override {
			const Real a = frequency*t;
			if (order == 1) return {-radius*frequency*std::sin(a), radius*frequency*std::cos(a)};
			return {-radius*frequency*frequency*std::cos(a), -radius*frequency*frequency*std::sin(a)};
		}
	};
	template<class F> constexpr bool has_normal_flag = requires(const F& f) { f.normal_defined; };
	template<class F> bool NormalDefined(const F& f) {
		if constexpr (has_normal_flag<F>) return f.normal_defined;
		else return f.normal.NormL2() != 0;
	}
	template<class Result>
	void CheckCount(const Result& r, int actual) {
		CHECK(r.IsSuccess());
		CHECK(r.function_evaluations == actual);
		CHECK_FALSE(r.algorithm_name.empty());
	}
}
using namespace MetricFieldContracts;

TEST_CASE("Metric musical operations preserve exact index slots and inverse pairing", "[metric-field-contract]") {
	StoredMetric metric;
	VectorN<Real, 2> pos{1, 1};
	Tensor2<2> t(2, 0);
	t(0,0)=2; t(0,1)=4; t(1,0)=6; t(1,1)=12;
	for (int slot : {0, 1}) {
		const auto raised = metric.RaiseIndex(t, slot, pos);
		CHECK(raised.isContravar(slot));
		CHECK(raised.isCovar(1-slot));
		const auto both = metric.RaiseIndex(raised, 1-slot, pos);
		CHECK(both.numContravar() == 2);
		const auto lower = metric.LowerIndex(raised, slot, pos);
		CHECK(lower.numCovar() == 2);
		for (int i=0;i<2;++i) for (int j=0;j<2;++j)
			CHECK_THAT(lower(i,j), Catch::Matchers::WithinAbs(t(i,j), Real(1e-12)));
		CHECK_THROWS_AS(metric.RaiseIndex(raised, slot, pos), TensorCovarContravarNumError);
		CHECK_THROWS_AS(metric.LowerIndex(t, slot, pos), TensorCovarContravarNumError);
	}
	CHECK_THROWS_AS(metric.RaiseIndex(t, -1, pos), ArgumentError);
	CHECK_THROWS_AS(metric.LowerIndex(t, 2, pos), ArgumentError);
}

TEST_CASE("Covariant access and typed adapters respect inverse metric storage", "[metric-field-contract]") {
	MetricTensorSphericalContravar inverse;
	const VectorN<Real, 3> pos{2, Constants::PI/2, 0};
	const auto typed = MetricAt<3, Spherical3>(inverse, pos);
	CHECK_THAT(typed.covariantComponents()(1,1), Catch::Matchers::WithinAbs(Real(4), Real(1e-12)));
	const auto cov = inverse.GetCovariantMetric(pos);
	const auto inv = inverse.GetContravariantMetric(pos);
	CHECK_THAT(cov(1,1)*inv(1,1), Catch::Matchers::WithinAbs(Real(1), Real(1e-12)));
	const TangentVector<3, Spherical3> v{1, 2, 3};
	const auto back = Sharp(typed, Flat(typed, v));
	for (int i=0;i<3;++i) CHECK_THAT(back[i], Catch::Matchers::WithinAbs(v[i], Real(1e-12)));
	StoredMetric asymmetric;
	asymmetric.asymmetric = true;
	CHECK_THROWS_AS(asymmetric.GetCovariantMetric({0,0}), DomainError);
}

TEST_CASE("Legacy induced metrics retain every ambient Jacobian row and borrow safely", "[metric-field-contract]") {
	Embed embedding;
	using Adapter = MetricTensorFromCoordTransf<VectorN<Real,2>, VectorN<Real,3>, 2>;
	const Adapter metric(embedding);
	const auto g = metric.GetCovariantMetric({0,0});
	CHECK(g(0,0) == 1);
	CHECK(g(1,1) == 2);
	CHECK(metric.component(1,1,{0,0}) == 2);
	CHECK_FALSE((std::is_constructible_v<Adapter, Embed&&>));
}

TEST_CASE("Lorentzian classification rejects nonfinite data and invalid tolerances", "[metric-field-contract]") {
	MetricTensorMinkowski metric;
	const VectorN<Real,4> origin{0,0,0,0};
	CHECK(metric.IsTimelike({2,1,0,0}, origin));
	CHECK(metric.IsSpacelike({1,2,0,0}, origin));
	CHECK(metric.IsNull({1,1,0,0}, origin));
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	CHECK_THROWS_AS(metric.ClassifyInterval({nan,0,0,0}, origin), DomainError);
	CHECK_THROWS_AS(metric.ClassifyInterval({1,0,0,0}, origin, -1), NumericInputError);
	CHECK_THROWS_AS(metric.ClassifyInterval({1,0,0,0}, origin, nan), NumericInputError);
	struct LocalFrame {};
	CHECK_NOTHROW((MetricAt<4, LocalFrame>(metric, origin)));
}

TEST_CASE("Frenet normal is invariant under regular slow and small reparameterizations", "[metric-field-contract]") {
	const Real angularStep = TOL(1e-3, Precision::SixthRootEps);
	for (Real frequency : {Real(1), Real(1e-8)}) {
		ParameterizedCircle circle; circle.frequency = frequency;
		const auto frame = DifferentialGeometry::ComputeFrenetFrame(circle, 0, angularStep/frequency);
		CHECK_THAT(frame.normal[0], Catch::Matchers::WithinAbs(Real(-1), Real(1e-6)));
		CHECK_THAT(frame.curvature, Catch::Matchers::WithinAbs(Real(1), FrameDifferenceBudget(angularStep)));
		CHECK(NormalDefined(frame));
	}
	ParameterizedCircle small; small.radius = TOL(1e-200, 1e-25);
	const auto frame = DifferentialGeometry::ComputeFrenetFrame(small, 0);
	CHECK_THAT(frame.normal[0], Catch::Matchers::WithinAbs(Real(-1), Real(1e-8)));
	CHECK_THAT(frame.speed, Catch::Matchers::WithinRel(small.radius, Precision::EqualityTolerance));
	CHECK_THAT(frame.curvature, Catch::Matchers::WithinRel(1 / small.radius, Precision::EqualityTolerance));
	CHECK((has_normal_flag<decltype(frame)>));
	CHECK_THROWS_AS(DifferentialGeometry::ComputeFrenetFrame(small, 0, Real(-1)), NumericInputError);
}

TEST_CASE("Detailed finite policy and failed callback counts are truthful", "[metric-field-contract]") {
	Quadratic f; f.poison = true;
	const VectorN<Real,3> p{2,1,Real(.5)};
	FieldOperationConfig config; config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto bad = ScalarFieldOperations::GradientCartDetailed(f,p,config);
	CHECK_FALSE(bad.IsSuccess());
	CHECK(bad.status == AlgorithmStatus::NumericalInstability);
	CHECK(bad.function_evaluations == f.calls);
	CHECK_FALSE(bad.error_message.empty());
	f.poison = false; f.calls = 0; f.fail_after = 3;
	const auto failed = ScalarFieldOperations::GradientCartDetailed(f,p,config);
	CHECK(failed.status == AlgorithmStatus::InvalidInput);
	CHECK(failed.function_evaluations == 3);
	CHECK(std::isnan(failed.value[0]));
	CHECK(std::isnan(failed.error[0]));
	CHECK(failed.elapsed_time_ms >= 0);
	f.calls = 0; f.fail_after = 0; config.derivative_order = 3;
	const auto invalid = ScalarFieldOperations::GradientCartDetailed(f,p,config);
	CHECK(invalid.status == AlgorithmStatus::InvalidInput);
	CHECK(f.calls == 0);
}

TEST_CASE("Every detailed field operator dispatches orders and counts actual callbacks", "[metric-field-contract]") {
	const VectorN<Real,3> p{2,1,Real(.5)};
	Quadratic f;
	LinearField field;
	for (int order : {0,1,2,4,6,8}) for (bool error : {false,true}) {
		FieldOperationConfig config; config.derivative_order=order; config.estimate_error=error;
		auto scalarCheck = [&](auto operation) { f.calls=0; auto r=operation(); CheckCount(r,f.calls); return r; };
		auto vectorCheck = [&](auto operation) { field.calls=0; auto r=operation(); CheckCount(r,field.calls); return r; };
		const auto gradient=scalarCheck([&] { return ScalarFieldOperations::GradientCartDetailed(f,p,config); });
		CHECK_THAT(gradient.value[0], Catch::Matchers::WithinAbs(Real(4),
			std::max(Real(1e-5), 5 * PolynomialFirstDifferenceBudget(order))));
		scalarCheck([&] { return ScalarFieldOperations::GradientSpherDetailed(f,Vec3Sph(p),config); });
		scalarCheck([&] { return ScalarFieldOperations::GradientCylDetailed(f,Vec3Cyl(p),config); });
		const auto div=vectorCheck([&] { return VectorFieldOperations::DivCartDetailed(field,p,config); });
		CHECK_THAT(div.value,Catch::Matchers::WithinAbs(Real(6),
			std::max(Real(1e-6), 7 * PolynomialFirstDifferenceBudget(order))));
		vectorCheck([&] { return VectorFieldOperations::DivSpherDetailed(field,p,config); });
		vectorCheck([&] { return VectorFieldOperations::DivCylDetailed(field,p,config); });
		vectorCheck([&] { return VectorFieldOperations::CurlCartDetailed(field,p,config); });
		vectorCheck([&] { return VectorFieldOperations::CurlSpherDetailed(field,p,config); });
		vectorCheck([&] { return VectorFieldOperations::CurlCylDetailed(field,p,config); });
	}
	for(int order : {0,2,4}) {
		FieldOperationConfig config; config.derivative_order=order; config.estimate_error=true;
		f.calls=0;
		const auto lap=ScalarFieldOperations::LaplacianCartDetailed(f,p,config);
		CheckCount(lap,f.calls);
		CHECK_THAT(lap.value,Catch::Matchers::WithinAbs(Real(6),Real(1e-3)));
		f.calls=0; const auto spher=ScalarFieldOperations::LaplacianSpherDetailed(f,Vec3Sph(p),config); CheckCount(spher,f.calls);
		f.calls=0; const auto cyl=ScalarFieldOperations::LaplacianCylDetailed(f,Vec3Cyl(p),config); CheckCount(cyl,f.calls);
	}
}

TEST_CASE("First-order field derivatives retain forward bias while provider dispatch stays exact", "[metric-field-contract]") {
	ProvidedQuadratic f;
	const VectorN<Real, 3> origin{0, 0, 0};
	const auto exact = ScalarFieldOperations::GradientCartDetailed(f, origin);
	REQUIRE(exact.IsSuccess());
	CHECK(exact.function_evaluations == 0);
	CHECK(f.derivatives == 3);
	for (int axis = 0; axis < 3; ++axis) CHECK(exact.value[axis] == 0);

	FieldOperationConfig config;
	config.derivative_order = 1;
	config.estimate_error = false;
	const auto forward = ScalarFieldOperations::GradientCartDetailed(f, origin, config);
	REQUIRE(forward.IsSuccess());
	CHECK(forward.function_evaluations == 6);
	CHECK(f.calls == 6);
	CHECK(f.derivatives == 3);
	for (int axis = 0; axis < 3; ++axis)
		CHECK_THAT(forward.value[axis], Catch::Matchers::WithinRel(Derivation::NDer1_h, Precision::RoundoffZeroThreshold));
}

TEST_CASE("Unsupported Laplacian orders fail instead of being silently ignored", "[metric-field-contract]") {
	Quadratic f;
	FieldOperationConfig config; config.exception_policy=EvaluationExceptionPolicy::ConvertToStatus;
	for (int order : {1,3,6,8}) {
		config.derivative_order=order;
		const auto r=ScalarFieldOperations::LaplacianCartDetailed(f,{2,1,0},config);
		CHECK(r.status == AlgorithmStatus::InvalidInput);
	}
	CHECK(f.calls==0);
}

TEST_CASE("Raised metric slots reach the actual anisotropic tensor transformation", "[metric-field-contract]") {
	struct AnisotropicMetric : MetricTensorField<3> {
		Real component(int i,int j,const VectorN<Real,3>&) const override { return i == j ? Real(i+2) : 0; }
	} metric;
	const CoordTransf3DCartGeneral transform(MatrixNM<Real,3,3>{2,0,0,0,3,0,0,0,4});
	Tensor2<3> t(2,0); t(0,1)=2;
	const auto raised=metric.RaiseIndex(t,0,{1,1,1});
	const auto transformed=transform.transfTensor2(raised,Vec3Cart{1,1,1});
	CHECK_THAT(transformed(0,1),Catch::Matchers::WithinAbs(Real(2)/3,TOL(1e-10, Precision::RoundoffZeroThreshold)));
	CHECK(transformed.isContravar(0));
	CHECK(transformed.isCovar(1));
	Reduce reduction;
	const MetricTensorFromCoordTransf<VectorN<Real,3>,VectorN<Real,2>,3> reduced(reduction);
	const auto gram=reduced.GetCovariantMetric({0,0,0});
	CHECK(gram(0,0)==1);
	CHECK(gram(1,1)==1);
	CHECK(gram(2,2)==0);
	CHECK(reduced.component(2,2,{0,0,0})==0);
	CHECK_THROWS_AS(reduced.GetContravariantMetric({0,0,0}),SingularMatrixError);
}

TEST_CASE("Zero curvature is explicit without rejecting a regular straight frame", "[metric-field-contract]") {
	struct Line : IParametricCurve<2>, ICurveDerivativeProvider<2> {
		Real getMinT() const override { return -1; }
		Real getMaxT() const override { return 1; }
		VectorN<Real,2> operator()(Real t) const override { return {t,0}; }
		VectorN<Real,2> EvaluateDerivative(Real, int order) const override { return order==1 ? VectorN<Real,2>{1,0} : VectorN<Real,2>{0,0}; }
	} line;
	const auto frame=DifferentialGeometry::ComputeFrenetFrame(line,0);
	CHECK_FALSE(NormalDefined(frame));
	CHECK(frame.curvature==0);
	CHECK(frame.normal[0]==0);
	CHECK(frame.normal[1]==0);
	CHECK(frame.tangent[0]==1);
	CHECK_THROWS_AS(DifferentialGeometry::ComputeFrenetFrame(line,2),DomainError);
}

TEST_CASE("Detailed numerical orders have their advertised stencil costs", "[metric-field-contract]") {
	const VectorN<Real,3> p{2,1,Real(.5)};
	for (int order : {0,1,2,4,6,8}) for (bool error : {false,true}) {
		Quadratic f; LinearField field;
		FieldOperationConfig config; config.derivative_order=order; config.estimate_error=error;
		const int effective=order==0 ? 4 : order;
		const int firstCost=effective==1 ? (error ? 3 : 2) : effective+(error ? 2 : 0);
		CHECK(ScalarFieldOperations::GradientCartDetailed(f,p,config).function_evaluations==3*firstCost);
		CHECK(VectorFieldOperations::DivCartDetailed(field,p,config).function_evaluations==3*firstCost);
		CHECK(VectorFieldOperations::DivSpherDetailed(field,p,config).function_evaluations==1+3*firstCost);
		CHECK(VectorFieldOperations::CurlCylDetailed(field,p,config).function_evaluations==1+6*firstCost);
		MetricTensorCartesian3D metric;
		const int before=field.calls;
		const auto det=VectorFieldOperations::DivergenceViaDetDetailed(field,p,metric,config);
		CHECK(det.function_evaluations==field.calls-before);
		CHECK(det.function_evaluations==3*firstCost);
	}
	for (int order : {0,2,4}) for(bool error : {false,true}) {
		Quadratic f; FieldOperationConfig config; config.derivative_order=order; config.estimate_error=error;
		const int first=order==2 ? 2 : 4;
		const int second=order==2 ? 3 : 5;
		const int fc=first+(error ? 2 : 0), sc=second+(error ? 2 : 0);
		CHECK(ScalarFieldOperations::LaplacianCartDetailed(f,p,config).function_evaluations==3*sc);
		CHECK(ScalarFieldOperations::LaplacianSpherDetailed(f,Vec3Sph(p),config).function_evaluations==3*sc+2*fc);
		CHECK(ScalarFieldOperations::LaplacianCylDetailed(f,Vec3Cyl(p),config).function_evaluations==3*sc+fc);
	}
}

TEST_CASE("Detailed controls apply to every scalar and vector result surface", "[metric-field-contract]") {
	const VectorN<Real,3> p{2,1,Real(.5)};
	Quadratic f; f.poison=true; LinearField field;
	FieldOperationConfig config; config.exception_policy=EvaluationExceptionPolicy::ConvertToStatus;
	auto checkBad=[&](const auto& r) {
		CHECK_FALSE(r.IsSuccess());
		CHECK(r.status==AlgorithmStatus::NumericalInstability);
	};
	checkBad(ScalarFieldOperations::GradientSpherDetailed(f,Vec3Sph(p),config));
	checkBad(ScalarFieldOperations::GradientCylDetailed(f,Vec3Cyl(p),config));
	checkBad(ScalarFieldOperations::LaplacianCartDetailed(f,p,config));
	checkBad(ScalarFieldOperations::LaplacianSpherDetailed(f,Vec3Sph(p),config));
	checkBad(ScalarFieldOperations::LaplacianCylDetailed(f,Vec3Cyl(p),config));
	config.check_finite=false; config.estimate_error=false;
	const auto unchecked=ScalarFieldOperations::GradientCartDetailed(f,p,config);
	CHECK(unchecked.IsSuccess());
	CHECK(std::isnan(unchecked.value[0]));
	CHECK(std::isnan(unchecked.error[0]));
	config.derivative_order=3;
	CHECK(VectorFieldOperations::DivSpherDetailed(field,p,config).status==AlgorithmStatus::InvalidInput);
	CHECK(VectorFieldOperations::DivCylDetailed(field,p,config).status==AlgorithmStatus::InvalidInput);
	CHECK(VectorFieldOperations::CurlCartDetailed(field,p,config).status==AlgorithmStatus::InvalidInput);
	CHECK(VectorFieldOperations::CurlSpherDetailed(field,p,config).status==AlgorithmStatus::InvalidInput);
	CHECK(VectorFieldOperations::CurlCylDetailed(field,p,config).status==AlgorithmStatus::InvalidInput);
	CHECK(field.calls==0);
	config.derivative_order=6;
	CHECK(ScalarFieldOperations::LaplacianSpherDetailed(f,Vec3Sph(p),config).status==AlgorithmStatus::InvalidInput);
	CHECK(ScalarFieldOperations::LaplacianCylDetailed(f,Vec3Cyl(p),config).status==AlgorithmStatus::InvalidInput);
	config.derivative_order=4; f.poison=false; f.calls=0;
	const Real nan=std::numeric_limits<Real>::quiet_NaN();
	CHECK(ScalarFieldOperations::GradientCartDetailed(f,{nan,1,0},config).status==AlgorithmStatus::InvalidInput);
	CHECK(f.calls==0);
}

TEST_CASE("Typed and legacy metrics reject malformed storage without assuming positivity", "[metric-field-contract]") {
	struct Frame {};
	for (Real scale : {TOL(1e-200, 1e-20),Real(1),TOL(1e200, 1e20)}) {
		CHECK_NOTHROW((Metric<2,Frame>(MatrixNM<Real,2,2>{-2*scale,0,0,3*scale})));
		CHECK_THROWS_AS((Metric<2,Frame>(MatrixNM<Real,2,2>{2*scale,scale,0,3*scale})),DomainError);
	}
	struct Mixed : MetricTensorField<2> {
		Mixed() : MetricTensorField<2>(1,1) {}
		Real component(int i,int j,const VectorN<Real,2>&) const override { return i==j ? Real(1) : Real(0); }
	} mixed;
	CHECK_THROWS_AS(mixed.GetCovariantMetric({0,0}),DomainError);
	CHECK_THROWS_AS(mixed.GetContravariantMetric({0,0}),DomainError);
	const Real nan=std::numeric_limits<Real>::quiet_NaN();
	CHECK_THROWS_AS((Metric<2,Frame>(MatrixNM<Real,2,2>{1,0,0,nan})),DomainError);
}

TEST_CASE("Detailed vector and provider routes enforce finite policy with exact failed work", "[metric-field-contract]") {
	struct BadVector : IVectorFunction<3> {
		mutable int calls=0;
		VectorN<Real,3> operator()(const VectorN<Real,3>& p) const override {
			++calls; return {p[0],p[1],std::numeric_limits<Real>::quiet_NaN()};
		}
	} field;
	const VectorN<Real,3> p{2,1,Real(.5)};
	FieldOperationConfig config; config.exception_policy=EvaluationExceptionPolicy::ConvertToStatus;
	auto check=[&](auto operation) {
		field.calls=0;
		const auto result=operation();
		CHECK(result.status==AlgorithmStatus::NumericalInstability);
		CHECK(result.function_evaluations==field.calls);
		CHECK(result.function_evaluations==1);
		CHECK_FALSE(result.error_message.empty());
	};
	check([&] {return VectorFieldOperations::DivCartDetailed(field,p,config);});
	check([&] {return VectorFieldOperations::DivSpherDetailed(field,p,config);});
	check([&] {return VectorFieldOperations::DivCylDetailed(field,p,config);});
	check([&] {return VectorFieldOperations::CurlCartDetailed(field,p,config);});
	check([&] {return VectorFieldOperations::CurlSpherDetailed(field,p,config);});
	check([&] {return VectorFieldOperations::CurlCylDetailed(field,p,config);});
	MetricTensorCartesian3D flat;
	check([&] {return VectorFieldOperations::DivergenceViaDetDetailed(field,p,flat,config);});
	struct BadProvider : ProvidedLinear {
		Real EvaluateComponentPartial(const VectorN<Real,3>&,int,int) const override {return Constants::PosInf;}
	} provider;
	const auto providerFailure=VectorFieldOperations::CurlCartDetailed(provider,p,config);
	CHECK(providerFailure.status==AlgorithmStatus::NumericalInstability);
	CHECK(providerFailure.function_evaluations==0);
	CHECK(provider.calls==0);
	struct BadScalarProvider : ProvidedQuadratic {
		Real EvaluatePartial(const VectorN<Real,3>&,int) const override {return Constants::PosInf;}
		Real EvaluateSecondPartial(const VectorN<Real,3>&,int,int) const override {return Constants::PosInf;}
	} scalar;
	CHECK(ScalarFieldOperations::GradientCartDetailed(scalar,p,config).status==AlgorithmStatus::NumericalInstability);
	CHECK(ScalarFieldOperations::LaplacianSpherDetailed(scalar,Vec3Sph(p),config).status==AlgorithmStatus::NumericalInstability);
	CHECK(scalar.calls==0);
}

TEST_CASE("Darboux curvature survives slow reversed surface parameterizations", "[metric-field-contract]") {
	struct Equator : IParametricCurve<2> {
		Real frequency=Real(1e-8);
		Real getMinT() const override {return Constants::NegInf;}
		Real getMaxT() const override {return Constants::PosInf;}
		VectorN<Real,2> operator()(Real t) const override {return {Constants::PI/2,1+frequency*t};}
	} curve;
	Surfaces::Sphere sphere(1);
	const Real angularStep = TOL(1e-3, Precision::SixthRootEps);
	for(Real frequency : {Real(1e-8),Real(-1e-8),TOL(1e-100, 1e-15),-TOL(1e-100, 1e-15)}) {
		curve.frequency=frequency;
		const auto frame=DifferentialGeometry::ComputeDarbouxFrame(sphere,curve,0,angularStep/std::abs(frequency));
		CHECK_THAT(frame.curvature,Catch::Matchers::WithinAbs(Real(1),FrameDifferenceBudget(angularStep)));
		CHECK_THAT(frame.normalCurvature,Catch::Matchers::WithinAbs(Real(-1),FrameDifferenceBudget(angularStep)));
		// Sphere's inherited normal also contains first-order surface-difference error.
		const Real normalStep = Derivation::ScaleStep(Derivation::NDer1_h, Constants::PI / 2);
		CHECK_THAT(frame.geodesicCurvature,Catch::Matchers::WithinAbs(Real(0),
			TOL(1e-6, FrameDifferenceBudget(angularStep) + normalStep)));
		CHECK_THAT(frame.tangentNormal.NormL2(),Catch::Matchers::WithinAbs(Real(1),TOL(1e-12, Precision::RoundoffZeroThreshold)));
	}
}

TEST_CASE("Physical coordinate formulas match independent polynomial identities", "[metric-field-contract]") {
	struct Quartic : IScalarFunction<3>, IScalarDerivativeProvider<3> {
		Real operator()(const VectorN<Real,3>& p) const override {
			Real f=0; for(int i=0;i<3;++i) f += p[i]*p[i]*p[i]*p[i]; return f;
		}
		Real EvaluatePartial(const VectorN<Real,3>& p, int a) const override { return 4*p[a]*p[a]*p[a]; }
		Real EvaluateSecondPartial(const VectorN<Real,3>& p, int a, int b) const override {
			return a == b ? 12*p[a]*p[a] : Real(0);
		}
	} f;
	struct CubicVector : IVectorFunction<3>, IVectorDerivativeProvider<3> {
		VectorN<Real,3> operator()(const VectorN<Real,3>& p) const override {
			return {p[0]*p[0]*p[0],2*p[1]*p[1]*p[1],3*p[2]*p[2]*p[2]};
		}
		Real EvaluateComponentPartial(const VectorN<Real,3>& p, int c, int a) const override {
			return c == a ? Real(3*(c+1))*p[a]*p[a] : Real(0);
		}
	} field;
	for (Real r : {Real(.75),Real(1.25),Real(2)}) for(Real theta : {Real(.4),Real(1.2),Real(2.4)}) {
		const Real phi=Real(.7), sine=std::sin(theta), cot=std::cos(theta)/sine;
		const VectorN<Real,3> p{r,theta,phi};
		bool analytic = false;
		auto check=[&](Real value,Real expected,Real estimate,Real tolerance=Real(1e-5)) {
			Real budget = tolerance*(1+std::abs(expected));
			if (!analytic) {
				REQUIRE(std::isfinite(estimate));
				REQUIRE(estimate >= 0);
				// Coordinate-weighted estimates include truncation and sample cancellation.
				budget = std::max(budget, 2*estimate);
			}
			CHECK_THAT(value,Catch::Matchers::WithinAbs(expected,budget));
		};
		for (int order : {0,1,2,4,6,8}) {
			CAPTURE(order, r, theta);
			analytic = order == 0;
			FieldOperationConfig config; config.derivative_order=order; config.estimate_error=true;
			const auto s=ScalarFieldOperations::GradientSpherDetailed(f,Vec3Sph(p),config);
			check(s.value[0],4*r*r*r,s.error[0]); check(s.value[1],4*theta*theta*theta/r,s.error[1]); check(s.value[2],4*phi*phi*phi/(r*sine),s.error[2]);
			const auto c=ScalarFieldOperations::GradientCylDetailed(f,Vec3Cyl(p),config);
			check(c.value[0],4*r*r*r,c.error[0]); check(c.value[1],4*theta*theta*theta/r,c.error[1]); check(c.value[2],4*phi*phi*phi,c.error[2]);
			const auto dc=VectorFieldOperations::DivCartDetailed(field,p,config);
			check(dc.value,3*r*r+6*theta*theta+9*phi*phi,dc.error);
			const auto ds=VectorFieldOperations::DivSpherDetailed(field,p,config);
			check(ds.value,5*r*r+(6*theta*theta+2*theta*theta*theta*cot)/r+9*phi*phi/(r*sine),ds.error);
			const auto dy=VectorFieldOperations::DivCylDetailed(field,p,config);
			check(dy.value,4*r*r+6*theta*theta/r+9*phi*phi,dy.error);
			const auto cs=VectorFieldOperations::CurlSpherDetailed(field,p,config);
			check(cs.value[0],3*phi*phi*phi*cot/r,cs.error[0]); check(cs.value[1],-3*phi*phi*phi/r,cs.error[1]); check(cs.value[2],2*theta*theta*theta/r,cs.error[2]);
			const auto cc=VectorFieldOperations::CurlCylDetailed(field,p,config);
			check(cc.value[0],0,cc.error[0]); check(cc.value[1],0,cc.error[1]); check(cc.value[2],2*theta*theta*theta/r,cc.error[2]);
		}
		for(int order : {0,2,4}) {
			CAPTURE(order, r, theta);
			analytic = order == 0;
			FieldOperationConfig config; config.derivative_order=order; config.estimate_error=true;
			const auto lc=ScalarFieldOperations::LaplacianCartDetailed(f,p,config);
			check(lc.value,12*(r*r+theta*theta+phi*phi),lc.error,Real(1e-4));
			const auto ls=ScalarFieldOperations::LaplacianSpherDetailed(f,Vec3Sph(p),config);
			check(ls.value,20*r*r+(4*theta*theta*theta*cot+12*theta*theta+12*phi*phi/(sine*sine))/(r*r),ls.error,Real(1e-4));
			const auto ly=ScalarFieldOperations::LaplacianCylDetailed(f,Vec3Cyl(p),config);
			check(ly.value,16*r*r+12*theta*theta/(r*r)+12*phi*phi,ly.error,Real(1e-4));
		}
	}
}

TEST_CASE("Detailed estimates carry the actual coordinate derivative weights", "[metric-field-contract]") {
	ScalarFunction<3> f(+[](const VectorN<Real,3>& p) {return std::exp(p[0])+std::sin(p[1])+std::pow(p[2],Real(7));});
	VectorFunction<3> field(+[](const VectorN<Real,3>& p) {
		return VectorN<Real,3>{std::exp(p[0])*std::sin(p[1]),std::sin(p[0])*std::exp(p[2]),std::pow(p[1],Real(7))*p[2]*p[2]*p[2]};
	});
	const VectorN<Real,3> p{Real(3.5),Real(.7),Real(.3)};
	const Real r=p[0], sine=std::sin(p[1]), cot=std::cos(p[1])/sine;
	Real first[3]{},second[3]{};
	MatrixNM<Real,3,3> e;
	for(int i=0;i<3;++i) {
		Derivation::NDer4Partial(f,i,p,&first[i]);
		Derivation::NSecDer4Partial(f,i,i,p,&second[i]);
		for(int j=0;j<3;++j) Derivation::NDer4Partial(field,i,j,p,&e(i,j));
	}
	FieldOperationConfig config; config.derivative_order=4; config.estimate_error=true;
	auto check=[](Real actual,Real expected) {
		CHECK_THAT(actual,Catch::Matchers::WithinRel(expected,TOL(1e-12, Precision::EqualityTolerance)));
	};
	const auto gs=ScalarFieldOperations::GradientSpherDetailed(f,Vec3Sph(p),config);
	check(gs.error[0],first[0]); check(gs.error[1],first[1]/r); check(gs.error[2],first[2]/(r*sine));
	const auto gc=ScalarFieldOperations::GradientCylDetailed(f,Vec3Cyl(p),config);
	check(gc.error[1],first[1]/r);
	check(VectorFieldOperations::DivSpherDetailed(field,p,config).error,e(0,0)+e(1,1)/r+e(2,2)/(r*sine));
	check(VectorFieldOperations::DivCylDetailed(field,p,config).error,e(0,0)+e(1,1)/r+e(2,2));
	const auto cs=VectorFieldOperations::CurlSpherDetailed(field,p,config);
	check(cs.error[0],e(2,1)/r+e(1,2)/(r*sine));
	check(cs.error[1],e(0,2)/(r*sine)+e(2,0));
	check(cs.error[2],e(1,0)+e(0,1)/r);
	const auto cc=VectorFieldOperations::CurlCylDetailed(field,p,config);
	check(cc.error[0],e(2,1)/r+e(1,2));
	check(cc.error[1],e(0,2)+e(2,0));
	check(cc.error[2],e(1,0)+e(0,1)/r);
	check(ScalarFieldOperations::LaplacianSpherDetailed(f,Vec3Sph(p),config).error,
		second[0]+2*first[0]/r+std::abs(cot)*first[1]/(r*r)+second[1]/(r*r)+second[2]/(r*r*sine*sine));
	check(ScalarFieldOperations::LaplacianCylDetailed(f,Vec3Cyl(p),config).error,
		second[0]+first[0]/r+second[1]/(r*r)+second[2]);
}

TEST_CASE("Provider derivative work is not mislabeled as function evaluation", "[metric-field-contract]") {
	ProvidedQuadratic f;
	ProvidedLinear field;
	const auto grad=ScalarFieldOperations::GradientCartDetailed(f,{2,1,0});
	CHECK(grad.function_evaluations==0);
	CHECK(f.calls==0);
	CHECK(f.derivatives==3);
	const auto div=VectorFieldOperations::DivCartDetailed(field,{2,1,0});
	CHECK(div.function_evaluations==0);
	CHECK(field.calls==0);
	CHECK(field.derivatives==3);
}
