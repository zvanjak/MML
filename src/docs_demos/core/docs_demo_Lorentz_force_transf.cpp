#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace MML;

namespace
{
	// CARTESIAN: E, v, B and F have components along (e_x, e_y, e_z).
	Vector3Cartesian LorentzForceCartesian(Real q, const Vector3Cartesian& E,
		const Vector3Cartesian& v, const Vector3Cartesian& B)
	{
		return q * (E + VectorProduct(v, B)); // F = q(E + v x B)
	}

	// CYLINDRICAL: E, v, B and F are PHYSICAL components along
	// (e_rho, e_phi, e_z) at the SAME point, not coordinate components.
	Vector3Cylindrical LorentzForceCylindrical(Real q, const Vector3Cylindrical& E,
		const Vector3Cylindrical& v, const Vector3Cylindrical& B)
	{
		const VectorN<Real, 3> vCrossB{
			v[1] * B[2] - v[2] * B[1],
			v[2] * B[0] - v[0] * B[2],
			v[0] * B[1] - v[1] * B[0]
		};
		// Use componentwise arithmetic, not Vector3Cylindrical's
		// position-style addition and scaling.
		const VectorN<Real, 3> electricComponents{E[0], E[1], E[2]};
		const auto F = q * (electricComponents + vCrossB); // F = q(E + v x B)
		return Vector3Cylindrical{F};
	}

	// The Jacobian returns coordinate components (v^rho, v^phi, v^z),
	// not physical components (v_rho, v_phi, v_z). Here h_phi = rho.
	Vector3Cylindrical PhysicalCylindricalComponents(const Vector3Cylindrical& coordinateComponents, Real rho)
	{
		return {coordinateComponents[0], rho * coordinateComponents[1], coordinateComponents[2]};
	}

	void RequireLorentzAgreement(const VectorN<Real, 3>& actual, const VectorN<Real, 3>& expected,
		Real scale, const char* caseName, const char* quantity)
	{
		const Real normalization = scale == 0 ? Real{1} : scale;
		for (int i = 0; i < 3; ++i)
			if (!IsNearlyEqual(actual[i] / normalization, expected[i] / normalization))
				throw std::runtime_error(std::string(caseName) + ": " + quantity + " comparison failed");
	}
}

void Docs_Demo_LorentzForceTransf()
{
	std::cout << "\n***********************************************************************\n"
		<< "****       LORENTZ FORCE: CARTESIAN VS CYLINDRICAL                  ****\n"
		<< "***********************************************************************\n";

	// Ideal infinite solenoid: Ampere's law gives B = mu0*n*I*e_z inside,
	// B = 0 outside. The ideal current sheet at rho = radius is excluded.
	const Real radius = Real{0.5};              // m
	const Real turnsPerMetre = Real{1000};      // 1/m
	const Real current = Real{2};               // A, producing B along +z
	const Real mu0 = Real{4e-7} * Constants::PI; // H/m, textbook vacuum approximation
	const Real insideB = mu0 * turnsPerMetre * current;
	const Real charge = Real{2e-6};             // C, illustrative positive test charge

	struct Case
	{
		const char* name;
		Vector3Cartesian position; // m
		Vector3Cartesian velocity; // m/s
		Vector3Cartesian electric; // V/m
		Real q;                    // C
	};
	const std::array<Case, 11> cases{{
		{"Cartesian x velocity", {0.12, 0.16, 0.3}, {10000, 0, 0}, {0, 0, 0}, charge},
		{"Cartesian y velocity", {0.12, 0.16, 0.3}, {0, 10000, 0}, {0, 0, 0}, charge},
		{"Radial velocity", {0.12, 0.16, 0.3}, {6000, 8000, 0}, {0, 0, 0}, charge},
		{"Azimuthal velocity", {0.12, 0.16, 0.3}, {-8000, 6000, 0}, {0, 0, 0}, charge},
		{"Parallel to B", {0.12, 0.16, 0.3}, {0, 0, 10000}, {0, 0, 0}, charge},
		{"Mixed velocity, second quadrant", {-0.15, 0.2, -0.4}, {3000, -4000, 5000}, {0, 0, 0}, charge},
		{"Mixed velocity, fourth quadrant", {0.24, -0.18, 0.1}, {-2000, 7000, -3000}, {0, 0, 0}, charge},
		{"Outside solenoid", {0.6, 0.8, 0.3}, {3000, -4000, 5000}, {0, 0, 0}, charge},
		{"Applied E, negative charge", {0.12, 0.16, 0.3}, {3000, -4000, 5000}, {100, -50, 25}, -charge},
		{"At rest in applied E", {-0.15, 0.2, -0.4}, {0, 0, 0}, {100, -50, 25}, charge},
		{"Neutral particle", {0.12, 0.16, 0.3}, {3000, -4000, 5000}, {100, -50, 25}, 0}
	}};

	std::cout << "Infinite solenoid along z: radius = " << radius << " m, n = "
		<< turnsPerMetre << " turns/m, I = " << current << " A\n"
		<< "B(rho < radius) = (0, 0, " << insideB << ") T; B(rho > radius) = 0\n"
		<< "E = 0 unless an additional applied electric field is specified.\n"
		<< "Cylindrical chart excludes rho = 0; no samples lie on the current sheet.\n"
		<< "Physical components use (e_rho, e_phi, e_z), a right-handed orthonormal basis.\n"
		<< "v_phi = rho*v^phi; force coordinate component F^phi = F_phi/rho.\n";

	for (const auto& sample : cases)
	{
		const auto positionCyl = CoordTransfCartToCyl.transf(sample.position);
		const Real rho = positionCyl[0];
		const auto velocityCoordinate = CoordTransfCartToCyl.transfVecContravariant(sample.velocity, sample.position);
		const auto velocityPhysical = PhysicalCylindricalComponents(velocityCoordinate, rho);
		const auto electricCoordinate = CoordTransfCartToCyl.transfVecContravariant(sample.electric, sample.position);
		const auto electricPhysical = PhysicalCylindricalComponents(electricCoordinate, rho);
		const Vector3Cylindrical magneticPhysical{0, 0, rho < radius ? insideB : Real{0}};
		const Vector3Cylindrical magneticCoordinate{magneticPhysical[0], magneticPhysical[1] / rho, magneticPhysical[2]};
		const auto magneticCart = CoordTransfCylToCart.transfVecContravariant(magneticCoordinate, positionCyl);

		const auto magneticForceCart = VectorProduct(sample.velocity, magneticCart) * sample.q;
		const auto forceCart = LorentzForceCartesian(sample.q, sample.electric, sample.velocity, magneticCart);
		const auto forcePhysical = LorentzForceCylindrical(sample.q, electricPhysical, velocityPhysical, magneticPhysical);
		// Vector3Cylindrical arithmetic treats (rho, phi, z) as a position
		// representation. Construct component triples explicitly instead.
		const Vector3Cylindrical forceCoordinate{forcePhysical[0], forcePhysical[1] / rho, forcePhysical[2]};
		const auto forceBack = CoordTransfCylToCart.transfVecContravariant(forceCoordinate, positionCyl);
		const auto cartForceInCyl = PhysicalCylindricalComponents(
			CoordTransfCartToCyl.transfVecContravariant(forceCart, sample.position), rho);

		const Real cosPhi = std::cos(positionCyl[1]), sinPhi = std::sin(positionCyl[1]);
		const VectorN<Real, 3> expectedVelocity{
			sample.velocity[0] * cosPhi + sample.velocity[1] * sinPhi,
			-sample.velocity[0] * sinPhi + sample.velocity[1] * cosPhi,
			sample.velocity[2]
		};
		const Real velocityScale = std::max(Real{1}, sample.velocity.NormL2());
		RequireLorentzAgreement(velocityPhysical, expectedVelocity, velocityScale, sample.name, "Velocity projection");
		RequireLorentzAgreement(CoordTransfCylToCart.transfVecContravariant(velocityCoordinate, positionCyl),
			sample.velocity, velocityScale, sample.name, "Velocity round trip");

		const Real expectedB = std::hypot(sample.position[0], sample.position[1]) < radius ? insideB : Real{0};
		const Vector3Cartesian expectedForce{
			sample.q * (sample.electric[0] + sample.velocity[1] * expectedB),
			sample.q * (sample.electric[1] - sample.velocity[0] * expectedB),
			sample.q * sample.electric[2]
		};
		const Real forceScale = std::abs(sample.q) *
			(sample.electric.NormL2() + sample.velocity.NormL2() * std::abs(expectedB));
		RequireLorentzAgreement(forceCart, expectedForce, forceScale, sample.name, "Analytic Cartesian force");
		RequireLorentzAgreement(forceBack, forceCart, forceScale, sample.name, "Cartesian force");
		RequireLorentzAgreement(forcePhysical, cartForceInCyl, forceScale, sample.name, "Cylindrical force");
		if (magneticForceCart.NormL2() > 0 && sample.velocity.NormL2() > 0 &&
			!IsNearlyZero(ScalarProduct(magneticForceCart.Normalized(), sample.velocity.Normalized())))
			throw std::runtime_error(std::string(sample.name) + ": Magnetic force is not perpendicular to velocity");

		std::cout << "\n" << sample.name << " (q = " << sample.q << " C)\n"
			<< "Position Cartesian [m]           : " << sample.position << '\n'
			<< "Position cylindrical [m,rad,m]   : " << positionCyl << '\n'
			<< "Velocity Cartesian [m/s]         : " << sample.velocity << '\n'
			<< "Velocity coordinate [m/s,rad/s,m/s]: " << velocityCoordinate << '\n'
			<< "Velocity physical [m/s]          : " << velocityPhysical << '\n'
			<< "E physical cylindrical [V/m]     : " << electricPhysical << '\n'
			<< "B physical cylindrical [T]       : " << magneticPhysical << '\n'
			<< "Force Cartesian [N]              : " << forceCart << '\n'
			<< "Force physical cylindrical [N]   : " << forcePhysical << '\n'
			<< "Force transformed back [N]       : " << forceBack << '\n'
			<< "Cartesian difference norm [N]    : " << (forceBack - forceCart).NormL2() << '\n';
	}
	std::cout << "\nAll 11 Lorentz-force cases agree in both coordinate systems.\n";
}
