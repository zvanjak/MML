/******************************************************************************
 * MML Example: General Relativity - Curvature and Geodesics
 * ============================================================================
 *
 * Demonstrates MML's differential geometry and GR helpers:
 *
 *   1. Intrinsic curvature of the unit 2-sphere
 *   2. Schwarzschild vacuum curvature diagnostics
 *   3. Fixed-step integration of an equatorial Schwarzschild geodesic
 *
 * Coordinates:
 *   - Unit sphere: (theta, phi)
 *   - Schwarzschild: (t, r, theta, phi), signature (-,+,+,+)
 *
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/algorithms/Geodesic.h>
#include <mml/core/MetricTensor.h>
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace MML;

class UnitSphereMetric : public MetricTensorField<2>
{
public:
	UnitSphereMetric() : MetricTensorField<2>(0, 2) { }

	Real Component(int i, int j, const VectorN<Real, 2>& pos) const override
	{
		if (i == 0 && j == 0)
			return REAL(1.0);
		if (i == 1 && j == 1)
			return std::sin(pos[0]) * std::sin(pos[0]);
		return REAL(0.0);
	}
};

class SchwarzschildMetric : public LorentzianMetric<4>
{
	Real _schwarzschildRadius;

public:
	SchwarzschildMetric(Real schwarzschildRadius)
		: LorentzianMetric<4>(0, 2), _schwarzschildRadius(schwarzschildRadius) { }

	Real Component(int i, int j, const VectorN<Real, 4>& pos) const override
	{
		Real r = pos[1];
		Real theta = pos[2];
		Real f = REAL(1.0) - _schwarzschildRadius / r;

		if (i == 0 && j == 0)
			return -f;
		if (i == 1 && j == 1)
			return REAL(1.0) / f;
		if (i == 2 && j == 2)
			return r * r;
		if (i == 3 && j == 3)
			return r * r * std::sin(theta) * std::sin(theta);
		return REAL(0.0);
	}
};

void PrintSection(const std::string& title)
{
	std::cout << "\n" << std::string(76, '=') << "\n";
	std::cout << "  " << title << "\n";
	std::cout << std::string(76, '=') << "\n\n";
}

void PrintFourVector(const std::string& label, const VectorN<Real, 4>& v)
{
	std::cout << "  " << std::left << std::setw(24) << label << " = ("
		<< std::right << std::setw(12) << v[0] << ", "
		<< std::setw(12) << v[1] << ", "
		<< std::setw(12) << v[2] << ", "
		<< std::setw(12) << v[3] << ")\n";
}

Real MetricNormSquared(const MetricTensorField<4>& metric, const VectorN<Real, 4>& pos, const VectorN<Real, 4>& v)
{
	Real norm = REAL(0.0);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			norm += metric.Component(i, j, pos) * v[i] * v[j];
	return norm;
}

void DemoSphereCurvature()
{
	PrintSection("SCENARIO 1: Curvature of the unit 2-sphere");

	UnitSphereMetric metric;
	VectorN<Real, 2> pos{ Constants::PI / REAL(4.0), Constants::PI / REAL(3.0) };
	Real sinTheta = std::sin(pos[0]);

	std::cout << std::fixed << std::setprecision(8);
	std::cout << "Point (theta, phi) = (" << pos[0] << ", " << pos[1] << ")\n";
	std::cout << "Metric components: g_theta_theta = " << metric.Component(0, 0, pos)
		<< ", g_phi_phi = " << metric.Component(1, 1, pos) << "\n\n";

	std::cout << "  R^theta_{phi theta phi} = " << metric.GetRiemannCurvatureTensor(0, 1, 0, 1, pos) << "\n";
	std::cout << "  expected sin^2(theta)    = " << sinTheta * sinTheta << "\n";
	std::cout << "  Ricci(theta,theta)       = " << metric.GetRicciTensor(0, 0, pos) << "\n";
	std::cout << "  Ricci(phi,phi)           = " << metric.GetRicciTensor(1, 1, pos) << "\n";
	std::cout << "  Ricci scalar             = " << metric.GetRicciScalar(pos) << "\n";
	std::cout << "  Einstein tensor G_00     = " << metric.GetEinsteinTensor(0, 0, pos) << "\n";

	std::cout << "\nFor a unit 2-sphere, scalar curvature R = 2 and the 2D Einstein tensor vanishes.\n";
}

void DemoSchwarzschildGeodesic()
{
	PrintSection("SCENARIO 2: Schwarzschild vacuum and equatorial geodesic");

	const Real schwarzschildRadius = REAL(2.0);
	const Real mass = schwarzschildRadius / REAL(2.0);
	const Real orbitRadius = REAL(10.0);
	SchwarzschildMetric metric(schwarzschildRadius);

	VectorN<Real, 4> pos{ REAL(0.0), orbitRadius, Constants::PI / REAL(2.0), REAL(0.0) };
	Real dt_dTau = REAL(1.0) / std::sqrt(REAL(1.0) - REAL(3.0) * mass / orbitRadius);
	Real dphi_dTau = std::sqrt(mass / (orbitRadius * orbitRadius * orbitRadius))
		/ std::sqrt(REAL(1.0) - REAL(3.0) * mass / orbitRadius);
	VectorN<Real, 4> velocity{ dt_dTau, REAL(0.0), REAL(0.0), dphi_dTau };

	std::cout << std::fixed << std::setprecision(8);
	std::cout << "Schwarzschild radius r_s = " << schwarzschildRadius << ", circular orbit r = " << orbitRadius << "\n";
	std::cout << "Metric is Lorentzian: " << (metric.IsLorentzian() ? "yes" : "no") << "\n";
	std::cout << "Vacuum diagnostics at (t,r,theta,phi) = (0,10,pi/2,0):\n";
	std::cout << "  Ricci scalar             = " << metric.GetRicciScalar(pos) << "\n";
	std::cout << "  Einstein tensor G_tt     = " << metric.GetEinsteinTensor(0, 0, pos) << "\n";
	std::cout << "  R^r_{theta r theta}      = " << metric.GetRiemannCurvatureTensor(1, 2, 1, 2, pos) << "\n\n";

	PrintFourVector("initial position", pos);
	PrintFourVector("initial velocity", velocity);
	std::cout << "  g(u,u) initial           = " << MetricNormSquared(metric, pos, velocity) << "\n\n";

	Real lambdaEnd = REAL(120.0);
	auto solution = IntegrateGeodesicFixedStep<4>(metric, pos, velocity, REAL(0.0), lambdaEnd, 1200);
	Vector<Real> finalState = solution.getXValuesAtEnd();
	VectorN<Real, 4> finalPosition = GeodesicPositionFromState<4>(finalState);
	VectorN<Real, 4> finalVelocity = GeodesicVelocityFromState<4>(finalState);

	PrintFourVector("final position", finalPosition);
	PrintFourVector("final velocity", finalVelocity);
	std::cout << "  g(u,u) final             = " << MetricNormSquared(metric, finalPosition, finalVelocity) << "\n";
	std::cout << "  radial drift             = " << finalPosition[1] - orbitRadius << "\n";
	std::cout << "  angular advance          = " << finalPosition[3] << " rad\n";
}

int main()
{
	std::cout << "\n";
	std::cout << "MinimalMathLibrary Example 08: General Relativity\n";
	std::cout << "Curvature tensors, Ricci/Einstein diagnostics, and geodesic integration\n";

	DemoSphereCurvature();
	DemoSchwarzschildGeodesic();

	std::cout << "\nGR examples complete.\n";
	return 0;
}