#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/Function.h>
#include <mml/core/FunctionHelpers.h>
#include <mml/core/Derivation.h>
#include <mml/core/Derivation/ADJacobians.h>
#endif

#include <iostream>

using namespace MML;

void Docs_Demo_Derivation_JacobianVocabulary();

void Docs_Demo_Derivation()
{
	RealFunction       f1{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };

	// numerical derivation of real function (available orders - 1, 2, 4, 6, 8)
	double der_f1 = Derivation::NDer1(f1, 0.5);
	double der_f2 = Derivation::NDer2(f1, 0.5, 1e-6);   // setting explicit step size
	Real err;
	double der_f6 = Derivation::NDer6(f1, 0.5, &err);   // if we need error estimate    
	double der_f8 = Derivation::NDer8(f1, 0.5);
	// we can use default Derive routine (set to NDer4)
	double der_f4 = Derivation::Derive(f1, 0.5);
	// available also with error estimate
	double der_f41 = Derivation::DeriveErr(f1, 0.5, &err);

	// second and third derivatives
	double sec_der_f1   = Derivation::NSecDer2(f1, 0.5);
	double third_der_f1 = Derivation::NThirdDer2(f1, 0.5);

	// creating new function that is derivation of existing function
	RealFuncDerived4    f1_der4(f1);        // 4th order derivation

	// scalar and vector functions
	ScalarFunction<3>   f2Scal([](const VectorN<Real, 3>& x) { return (Real)(1.0 / pow(x.NormL2(), 2)); });
	VectorFunction<3>   f3Vec([](const VectorN<Real, 3>& x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });

	VectorN<Real, 3>    der_point{ 1.0, 1.0, 1.0 };
	double der_f2Scal = Derivation::NDer1Partial(f2Scal, 1, der_point);
	VectorN<Real, 3> der_f2Scal_all = Derivation::NDer1PartialByAll(f2Scal, der_point);

	double der_f3Vec = Derivation::NDer1Partial(f3Vec, 1, 1, der_point);
	VectorN<Real, 3>     der_f3Vec_by1    = Derivation::NDer2PartialByAll(f3Vec, 1, der_point);
	MatrixNM<Real, 3, 3> der_f3Vec_by_all = Derivation::NDer4PartialAllByAll(f3Vec, der_point);

	ParametricCurve<3>  f4Curve([](Real x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });
	VectorN<Real, 3>    der_f4Curve   = Derivation::NDer1(f4Curve, 1.0);
	// Second and third derivatives use optimized direct formulas
	VectorN<Real, 3>    sec_f4Curve   = Derivation::NSecDer2(f4Curve, 1.0);   // O(h²), 3 evals
	VectorN<Real, 3>    third_f4Curve = Derivation::NThirdDer2(f4Curve, 1.0); // O(h²), 4 evals

	Docs_Demo_Derivation_JacobianVocabulary();
}

// Derivation.md "Jacobians, Hessians and Gradients — one vocabulary": the numerical and the two
// AD entry points share the Calc<Object>[<Method>] scheme; argument types select fixed vs dynamic size.
void Docs_Demo_Derivation_JacobianVocabulary()
{
	// Numerical: fixed size through the interface, dynamic size through std::function
	VectorFunction<2> polar([](const VectorN<Real, 2>& q) {
		return VectorN<Real, 2>{ q[0] * std::cos(q[1]), q[0] * std::sin(q[1]) };
	});
	VectorN<Real, 2> q0{ 2.0, 0.5 };
	MatrixNM<Real, 2, 2> Jfix = Derivation::CalcJacobian(polar, q0);

	int evals = 0;
	Matrix<Real> Jdyn = Derivation::CalcJacobian(
		[](const Vector<Real>& q) { return Vector<Real>{ q[0] * std::cos(q[1]), q[0] * std::sin(q[1]) }; },
		Vector<Real>{ 2.0, 0.5 }, 0.0, &evals);              // evals == 8

	// Forward and reverse AD on the same map
	auto polarAD = [](const auto& q) {
		using std::cos; using std::sin;
		return std::vector<std::decay_t<decltype(q[0])>>{ q[0] * cos(q[1]), q[0] * sin(q[1]) };
	};
	Matrix<Real> Jfwd = AD::CalcJacobianForwardAD(polarAD, std::vector<Real>{ 2.0, 0.5 });
	Matrix<Real> Jrev = AD::CalcJacobianReverseAD(polarAD, std::vector<Real>{ 2.0, 0.5 });
	// Jfix, Jdyn, Jfwd, Jrev agree to ~1e-10 (finite differences) / rounding (AD)

	std::cout << "CalcJacobian fixed/dyn/fwdAD/revAD (0,1): " << Jfix(0, 1) << " " << Jdyn(0, 1)
	          << " " << Jfwd(0, 1) << " " << Jrev(0, 1) << "  (evals = " << evals << ")\n";

	ScalarFunction<2> bowl([](const VectorN<Real, 2>& x) { return x[0] * x[0] + 3.0 * x[0] * x[1]; });
	MatrixNM<Real, 2, 2> H = Derivation::CalcHessian(bowl, q0);
	auto bowlAD = [](const auto& x) { return x[0] * x[0] + 3.0 * x[0] * x[1]; };
	Matrix<Real> Had = AD::CalcHessianForwardAD(bowlAD, std::vector<Real>{ 2.0, 0.5 });
	std::vector<Real> g = AD::CalcGradientReverseAD(bowlAD, std::vector<Real>{ 2.0, 0.5 });
	std::cout << "CalcHessian(0,1) = " << H(0, 1) << ", forward-AD " << Had(0, 1)
	          << ", reverse-AD gradient = (" << g[0] << ", " << g[1] << ")\n";
}