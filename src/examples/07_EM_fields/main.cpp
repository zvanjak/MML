/******************************************************************************
 * MML Example: Electromagnetism in Special Relativity
 * ============================================================================
 *
 * Demonstrates the tensor-based SR/EM helpers:
 *
 *   1. Boosted Coulomb field from a uniformly moving charge
 *   2. Retarded potentials for an oscillating electric dipole proxy
 *   3. Faraday tensor field extraction and electromagnetic invariants
 *
 * Units: natural units c = 1 and 1/(4*pi*epsilon0) = 1.
 * Coordinates: spacetime points are ordered as (ct, x, y, z).
 *
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/Function.h>
#endif

#include "Electromagnetism/Electromagnetism.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace MML;
using namespace MPL;

VectorN<Real, 3> Cross(const VectorN<Real, 3>& a, const VectorN<Real, 3>& b)
{
	return VectorN<Real, 3>{
		a[1] * b[2] - a[2] * b[1],
		a[2] * b[0] - a[0] * b[2],
		a[0] * b[1] - a[1] * b[0]
	};
}

Real Dot(const VectorN<Real, 3>& a, const VectorN<Real, 3>& b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Real Norm(const VectorN<Real, 3>& v)
{
	return std::sqrt(Dot(v, v));
}

void PrintVector(const std::string& label, const VectorN<Real, 3>& v)
{
	std::cout << "  " << std::left << std::setw(26) << label << " = ("
		<< std::right << std::setw(12) << v[0] << ", "
		<< std::setw(12) << v[1] << ", "
		<< std::setw(12) << v[2] << ")\n";
}

void DemoBoostedCoulombField()
{
	std::cout << "\n" << std::string(76, '=') << "\n";
	std::cout << "  SCENARIO 1: Boosted Coulomb field from Lienard-Wiechert potentials\n";
	std::cout << std::string(76, '=') << "\n\n";

	const Real charge = static_cast<MML::Real>(1.5);
	const VectorN<Real, 3> beta{ static_cast<MML::Real>(0.45), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0) };
	const VectorN<Real, 3> observation{ static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.25), static_cast<MML::Real>(0.4) };
	const Real observationTime = static_cast<MML::Real>(4.0);

	LienardWiechertOptions options;
	options.speedOfLight = static_cast<MML::Real>(1.0);
	options.tolerance = static_cast<MML::Real>(1e-11);
	options.velocityDerivativeStep = static_cast<MML::Real>(1e-4);
	options.maxIterations = 120;

	ParametricCurveFromStdFunc<3> source([beta](Real t) {
		return VectorN<Real, 3>{ beta[0] * t, beta[1] * t, beta[2] * t };
	});

	VectorFunctionFromStdFunc<4> covariantPotential([&](const VectorN<Real, 4>& x) {
		VectorN<Real, 3> spatial{ x[1], x[2], x[3] };
		return LienardWiechertPotential(charge, source, spatial, x[0], -static_cast<MML::Real>(20.0), options).CovariantComponents();
	});

	VectorN<Real, 4> event{ observationTime, observation[0], observation[1], observation[2] };
	Tensor2<4> faraday = FaradayTensor(covariantPotential, event, static_cast<MML::Real>(2e-4));
	VectorN<Real, 3> electric = ElectricFieldFromFaradayTensor(faraday);
	VectorN<Real, 3> magnetic = MagneticFieldFromFaradayTensor(faraday);

	auto retarded = LienardWiechertState(source, observation, observationTime, -static_cast<MML::Real>(20.0), options);
	VectorN<Real, 3> R = observation - source(observationTime);
	Real betaSquared = Dot(beta, beta);
	Real denominator = std::pow(Dot(R, R) - Dot(Cross(beta, R), Cross(beta, R)), static_cast<MML::Real>(1.5));
	VectorN<Real, 3> analyticElectric = R * (charge * (static_cast<MML::Real>(1.0) - betaSquared) / denominator);
	VectorN<Real, 3> analyticMagnetic = Cross(beta, analyticElectric);

	std::cout << std::fixed << std::setprecision(7);
	std::cout << "Uniform source velocity beta = (" << beta[0] << ", " << beta[1] << ", " << beta[2] << ")\n";
	std::cout << "Observation event (ct,x,y,z) = (" << event[0] << ", " << event[1] << ", " << event[2] << ", " << event[3] << ")\n";
	std::cout << "Retarded time = " << retarded.retardedTime << ", retarded denominator = " << retarded.denominator << "\n\n";

	PrintVector("E from Faraday tensor", electric);
	PrintVector("E analytic boosted", analyticElectric);
	PrintVector("B from Faraday tensor", magnetic);
	PrintVector("B = beta x E", analyticMagnetic);

	std::cout << "\n  |E_numeric - E_analytic| = " << Norm(electric - analyticElectric) << "\n";
	std::cout << "  |B_numeric - B_analytic| = " << Norm(magnetic - analyticMagnetic) << "\n";
	std::cout << "  F_munu F^munu          = " << FaradayTensorInvariant(faraday) << "\n";
	std::cout << "  epsilon F F            = " << FaradayTensorLeviCivitaInvariant(faraday) << "\n";
}

void DemoRadiatingDipoleProxy()
{
	std::cout << "\n" << std::string(76, '=') << "\n";
	std::cout << "  SCENARIO 2: Oscillating electric dipole from retarded point charges\n";
	std::cout << std::string(76, '=') << "\n\n";

	const Real charge = static_cast<MML::Real>(1.0);
	const Real amplitude = static_cast<MML::Real>(0.08);
	const Real omega = static_cast<MML::Real>(1.6);
	const VectorN<Real, 3> observation{ static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.6), static_cast<MML::Real>(1.4) };
	const Real observationTime = static_cast<MML::Real>(5.0);

	LienardWiechertOptions options;
	options.speedOfLight = static_cast<MML::Real>(1.0);
	options.tolerance = static_cast<MML::Real>(1e-11);
	options.velocityDerivativeStep = static_cast<MML::Real>(1e-4);
	options.maxIterations = 140;

	ParametricCurveFromStdFunc<3> positiveCharge([=](Real t) {
		return VectorN<Real, 3>{ static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), amplitude * std::sin(omega * t) };
	});
	ParametricCurveFromStdFunc<3> negativeCharge([=](Real t) {
		return VectorN<Real, 3>{ static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), -amplitude * std::sin(omega * t) };
	});

	VectorFunctionFromStdFunc<4> dipolePotential([&](const VectorN<Real, 4>& x) {
		VectorN<Real, 3> spatial{ x[1], x[2], x[3] };
		auto plus = LienardWiechertPotential(charge, positiveCharge, spatial, x[0], -static_cast<MML::Real>(20.0), options).CovariantComponents();
		auto minus = LienardWiechertPotential(-charge, negativeCharge, spatial, x[0], -static_cast<MML::Real>(20.0), options).CovariantComponents();

		VectorN<Real, 4> total;
		for (int i = 0; i < 4; i++)
			total[i] = plus[i] + minus[i];
		return total;
	});

	VectorN<Real, 4> event{ observationTime, observation[0], observation[1], observation[2] };
	Tensor2<4> faraday = FaradayTensor(dipolePotential, event, static_cast<MML::Real>(3e-4));
	VectorN<Real, 3> electric = ElectricFieldFromFaradayTensor(faraday);
	VectorN<Real, 3> magnetic = MagneticFieldFromFaradayTensor(faraday);

	auto plusState = LienardWiechertState(positiveCharge, observation, observationTime, -static_cast<MML::Real>(20.0), options);
	auto minusState = LienardWiechertState(negativeCharge, observation, observationTime, -static_cast<MML::Real>(20.0), options);

	std::cout << std::fixed << std::setprecision(7);
	std::cout << "Dipole half-separation z(t) = " << amplitude << " sin(" << omega << " t)\n";
	std::cout << "Observation event (ct,x,y,z) = (" << event[0] << ", " << event[1] << ", " << event[2] << ", " << event[3] << ")\n";
	std::cout << "Retarded times: positive charge " << plusState.retardedTime
		<< ", negative charge " << minusState.retardedTime << "\n\n";

	PrintVector("E dipole", electric);
	PrintVector("B dipole", magnetic);

	std::cout << "\n  |E|                    = " << Norm(electric) << "\n";
	std::cout << "  |B|                    = " << Norm(magnetic) << "\n";
	std::cout << "  E dot B                = " << Dot(electric, magnetic) << "\n";
	std::cout << "  F_munu F^munu          = " << FaradayTensorInvariant(faraday) << "\n";
	std::cout << "  epsilon F F            = " << FaradayTensorLeviCivitaInvariant(faraday) << "\n";
}

int main()
{
	std::cout << "\n";
	std::cout << "MinimalMathLibrary Example 07: Electromagnetism in Special Relativity\n";
	std::cout << "Natural units: c = 1, spacetime order (ct,x,y,z), signature (-,+,+,+)\n";

	DemoBoostedCoulombField();
	DemoRadiatingDipoleProxy();

	std::cout << "\nSR/EM examples complete.\n";
	return 0;
}
