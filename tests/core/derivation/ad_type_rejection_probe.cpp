// Compile with MML_AD_REJECT_CASE=1..16 to check intentional compile failures.
#ifdef MML_AD_PROBE_SINGLE_HEADER
#include <mml/single_header/MML.h>
#else
#include <mml/core/Derivation/ADJacobians.h>
#endif

using namespace MML;
using namespace MML::AD;

void ADTypeRejectionProbe() {
	const auto scalar = [](auto x) { using std::sin; return sin(x); };
	const auto field = [](const auto& x) { return x[0] * x[0]; };
	const auto vector = [](const auto& x) { return x; };
#if MML_AD_REJECT_CASE == 1
	Dual<int> x(2, 1);
#elif MML_AD_REJECT_CASE == 2
	Dual<bool> x(true, true);
#elif MML_AD_REJECT_CASE == 3
	Dual<Dual<int>> x;
#elif MML_AD_REJECT_CASE == 4
	Dual x(2, 1);
#elif MML_AD_REJECT_CASE == 5
	CalcDerivativeForwardAD(scalar, 2);
#elif MML_AD_REJECT_CASE == 6
	CalcGradientForwardAD(field, std::vector<int>{2});
#elif MML_AD_REJECT_CASE == 7
	CalcJacobianForwardAD(vector, std::vector<int>{2});
#elif MML_AD_REJECT_CASE == 8
	CalcJacobianForwardAD(vector, Vector<int>{2});
#elif MML_AD_REJECT_CASE == 9
	CalcHessianForwardAD(field, std::vector<int>{2});
#elif MML_AD_REJECT_CASE == 10
	CalcHessianForwardAD(field, Vector<int>{2});
#elif MML_AD_REJECT_CASE == 11
	derivative(scalar, 2);
#elif MML_AD_REJECT_CASE == 12
	gradient(field, std::vector<int>{2});
#elif MML_AD_REJECT_CASE == 13
	jacobianForwardAD(vector, std::vector<int>{2});
#elif MML_AD_REJECT_CASE == 14
	jacobianForwardAD(vector, Vector<int>{2});
#elif MML_AD_REJECT_CASE == 15
	hessianForwardAD(field, std::vector<int>{2});
#elif MML_AD_REJECT_CASE == 16
	hessianForwardAD(field, Vector<int>{2});
#else
	Dual<Real> x(2, 1);
	Dual<Dual<Real>> nested;
	CalcDerivativeForwardAD(scalar, Real(2));
	CalcGradientForwardAD(field, std::vector<Real>{2});
	CalcJacobianForwardAD(vector, std::vector<Real>{2});
	CalcJacobianForwardAD(vector, Vector<Real>{2});
	CalcHessianForwardAD(field, std::vector<Real>{2});
	CalcHessianForwardAD(field, Vector<Real>{2});
#endif
}
