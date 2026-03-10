///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ContinuousSystems.h                                                 ///
///  Description: Classic continuous dynamical systems                                ///
///               Lorenz, Rössler, Van der Pol, Duffing, Chua, Hénon-Heiles, etc.     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CONTINUOUS_SYSTEMS_H
#define MML_CONTINUOUS_SYSTEMS_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "systems/DynamicalSystemBase.h"

#include <cmath>

namespace MML::Systems 
{
	//=============================================================================
	// LORENZ SYSTEM
	//=============================================================================

	/// @brief The classic Lorenz attractor
	///
	/// dx/dt = σ(y - x)
	/// dy/dt = x(ρ - z) - y
	/// dz/dt = xy - βz
	///
	/// Classic parameters: σ = 10, ρ = 28, β = 8/3
	class LorenzSystem : public DynamicalSystemBase<3, 3> {
	public:
		/// @brief Construct with default chaotic parameters
		LorenzSystem(Real sigma = 10.0, Real rho = 28.0, Real beta = 8.0 / 3.0) {
			_params = {sigma, rho, beta};
			_stateNames = {"x", "y", "z"};
			_paramNames = {"sigma", "rho", "beta"};
			_paramRanges = {{0.0, 50.0}, {0.0, 100.0}, {0.0, 10.0}};
		}

		void derivs(Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override {
			Real sigma = _params[0];
			Real rho = _params[1];
			Real beta = _params[2];

			dydt[0] = sigma * (y[1] - y[0]);
			dydt[1] = y[0] * (rho - y[2]) - y[1];
			dydt[2] = y[0] * y[1] - beta * y[2];
		}

		void jacobian(Real /*t*/, const Vector<Real>& y, Matrix<Real>& J) const override {
			Real sigma = _params[0];
			Real rho = _params[1];
			Real beta = _params[2];

			J.Resize(3, 3);
			J(0, 0) = -sigma;
			J(0, 1) = sigma;
			J(0, 2) = 0;
			J(1, 0) = rho - y[2];
			J(1, 1) = -1;
			J(1, 2) = -y[0];
			J(2, 0) = y[1];
			J(2, 1) = y[0];
			J(2, 2) = -beta;
		}

		bool hasAnalyticalJacobian() const override { return true; }
		bool isDissipative() const override { return true; }
		Vector<Real> getDefaultInitialCondition() const override { return Vector<Real>({1.0, 1.0, 1.0}); }
		
		/// @brief Divergence of flow (always negative for dissipative)
		Real getDivergence() const { return -(_params[0] + 1 + _params[2]); }
	};

	//=============================================================================
	// RÖSSLER SYSTEM
	//=============================================================================

	/// @brief Rössler attractor
	///
	/// dx/dt = -y - z
	/// dy/dt = x + ay
	/// dz/dt = b + z(x - c)
	class RosslerSystem : public DynamicalSystemBase<3, 3> {
	public:
		/// @brief Construct with default parameters
		RosslerSystem(Real a = 0.2, Real b = 0.2, Real c = 5.7) {
			_params = {a, b, c};
			_stateNames = {"x", "y", "z"};
			_paramNames = {"a", "b", "c"};
			_paramRanges = {{0.0, 1.0}, {0.0, 1.0}, {0.0, 20.0}};
		}

		void derivs(Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override {
			Real a = _params[0];
			Real b = _params[1];
			Real c = _params[2];

			dydt[0] = -y[1] - y[2];
			dydt[1] = y[0] + a * y[1];
			dydt[2] = b + y[2] * (y[0] - c);
		}

		void jacobian(Real /*t*/, const Vector<Real>& y, Matrix<Real>& J) const override {
			Real a = _params[0];
			Real c = _params[2];

			J.Resize(3, 3);
			J(0, 0) = 0;
			J(0, 1) = -1;
			J(0, 2) = -1;
			J(1, 0) = 1;
			J(1, 1) = a;
			J(1, 2) = 0;
			J(2, 0) = y[2];
			J(2, 1) = 0;
			J(2, 2) = y[0] - c;
		}
	};

	//=============================================================================
	// VAN DER POL OSCILLATOR
	//=============================================================================

	/// @brief Van der Pol oscillator (2D)
	///
	/// x'' - μ(1 - x²)x' + x = 0
	/// Written as:
	/// dx/dt = y
	/// dy/dt = μ(1 - x²)y - x
	class VanDerPolSystem : public DynamicalSystemBase<2, 1> {
	public:
		/// @brief Construct with nonlinearity parameter
		VanDerPolSystem(Real mu = 1.0) {
			_params = {mu};
			_stateNames = {"x", "v"};
			_paramNames = {"mu"};
			_paramRanges = {{0.0, 10.0}};
		}

		void derivs(Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override {
			Real mu = _params[0];
			dydt[0] = y[1];
			dydt[1] = mu * (1 - y[0] * y[0]) * y[1] - y[0];
		}

		void jacobian(Real /*t*/, const Vector<Real>& y, Matrix<Real>& J) const override {
			Real mu = _params[0];
			J.Resize(2, 2);
			J(0, 0) = 0;
			J(0, 1) = 1;
			J(1, 0) = -2 * mu * y[0] * y[1] - 1;
			J(1, 1) = mu * (1 - y[0] * y[0]);
		}
	};

	//=============================================================================
	// DUFFING OSCILLATOR
	//=============================================================================

	/// @brief Duffing oscillator with forcing
	///
	/// x'' + δx' + αx + βx³ = γcos(ωt)
	/// State: [x, v, θ] where θ = ωt mod 2π
	class DuffingSystem : public DynamicalSystemBase<3, 5> {
	public:
		/// @brief Construct with parameters
		DuffingSystem(Real delta = 0.3, Real alpha = -1.0, Real beta = 1.0, Real gamma = 0.5, Real omega = 1.2) {
			_params = {delta, alpha, beta, gamma, omega};
			_stateNames = {"x", "v", "theta"};
			_paramNames = {"delta", "alpha", "beta", "gamma", "omega"};
			_paramRanges = {{0.0, 1.0}, {-2.0, 2.0}, {0.0, 2.0}, {0.0, 1.0}, {0.0, 3.0}};
		}

		void derivs(Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override {
			Real delta = _params[0];
			Real alpha = _params[1];
			Real beta = _params[2];
			Real gamma = _params[3];
			Real omega = _params[4];

			dydt[0] = y[1];
			dydt[1] = -delta * y[1] - alpha * y[0] - beta * y[0] * y[0] * y[0] + gamma * std::cos(y[2]);
			dydt[2] = omega;
		}

		void jacobian(Real /*t*/, const Vector<Real>& y, Matrix<Real>& J) const override {
			Real delta = _params[0];
			Real alpha = _params[1];
			Real beta = _params[2];
			Real gamma = _params[3];

			J.Resize(3, 3);
			J(0, 0) = 0;
			J(0, 1) = 1;
			J(0, 2) = 0;
			J(1, 0) = -alpha - 3 * beta * y[0] * y[0];
			J(1, 1) = -delta;
			J(1, 2) = -gamma * std::sin(y[2]);
			J(2, 0) = 0;
			J(2, 1) = 0;
			J(2, 2) = 0;
		}

		bool hasAnalyticalJacobian() const override { return true; }
		bool isAutonomous() const override { return false; }  // Has periodic forcing
	};

	//=============================================================================
	// CHUA'S CIRCUIT
	//=============================================================================

	/// @brief Chua's circuit - electronic chaos
	///
	/// dx/dt = α(y - x - f(x))
	/// dy/dt = x - y + z
	/// dz/dt = -βy
	///
	/// f(x) = bx + 0.5(a-b)(|x+1| - |x-1|)  (piecewise-linear)
	class ChuaCircuit : public DynamicalSystemBase<3, 4> {
	public:
		ChuaCircuit(Real alpha = 15.6, Real beta = 28.0, Real a = -1.143, Real b = -0.714) {
			_params = {alpha, beta, a, b};
			_stateNames = {"x", "y", "z"};
			_paramNames = {"alpha", "beta", "a", "b"};
			_paramRanges = {{0.0, 30.0}, {0.0, 50.0}, {-2.0, 0.0}, {-1.0, 0.0}};
		}

		void derivs(Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override {
			Real alpha = _params[0];
			Real beta = _params[1];
			Real a = _params[2];
			Real b = _params[3];

			// Chua's nonlinearity
			Real x = y[0];
			Real fx = b * x + 0.5 * (a - b) * (std::abs(x + 1) - std::abs(x - 1));

			dydt[0] = alpha * (y[1] - x - fx);
			dydt[1] = y[0] - y[1] + y[2];
			dydt[2] = -beta * y[1];
		}

		void jacobian(Real /*t*/, const Vector<Real>& y, Matrix<Real>& J) const override {
			Real alpha = _params[0];
			Real beta = _params[1];
			Real a = _params[2];
			Real b = _params[3];

			// Derivative of Chua nonlinearity
			Real x = y[0];
			Real dfx;
			if (x < -1)
				dfx = b;
			else if (x > 1)
				dfx = b;
			else
				dfx = a;

			J.Resize(3, 3);
			J(0, 0) = alpha * (-1 - dfx);
			J(0, 1) = alpha;
			J(0, 2) = 0;
			J(1, 0) = 1;
			J(1, 1) = -1;
			J(1, 2) = 1;
			J(2, 0) = 0;
			J(2, 1) = -beta;
			J(2, 2) = 0;
		}
	};

	//=============================================================================
	// HÉNON-HEILES SYSTEM
	//=============================================================================

	/// @brief Hénon-Heiles Hamiltonian system
	///
	/// H = 0.5(px² + py² + x² + y²) + x²y - y³/3
	/// 4D phase space: [x, y, px, py]
	class HenonHeilesSystem : public DynamicalSystemBase<4, 1> {
	public:
		/// @brief Construct with energy level
		HenonHeilesSystem(Real energy = 0.1) {
			_params = {energy};
			_stateNames = {"x", "y", "px", "py"};
			_paramNames = {"E"};
			_paramRanges = {{0.0, 0.166}};
		}

		void derivs(Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override {
			// Hamilton's equations
			dydt[0] = y[2];                          // dx/dt = ∂H/∂px = px
			dydt[1] = y[3];                          // dy/dt = ∂H/∂py = py
			dydt[2] = -y[0] - 2 * y[0] * y[1];       // dpx/dt = -∂H/∂x
			dydt[3] = -y[1] - y[0] * y[0] + y[1] * y[1]; // dpy/dt = -∂H/∂y
		}

		void jacobian(Real /*t*/, const Vector<Real>& y, Matrix<Real>& J) const override {
			J.Resize(4, 4);
			// Zero out
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					J(i, j) = 0;

			J(0, 2) = 1; // dx/dpx
			J(1, 3) = 1; // dy/dpy
			J(2, 0) = -1 - 2 * y[1]; // dpx/dx
			J(2, 1) = -2 * y[0];     // dpx/dy
			J(3, 0) = -2 * y[0];     // dpy/dx
			J(3, 1) = -1 + 2 * y[1]; // dpy/dy
		}

		bool hasAnalyticalJacobian() const override { return true; }
		bool isHamiltonian() const override { return true; }
		bool isDissipative() const override { return false; }
		
		Vector<Real> getDefaultInitialCondition() const override { 
			return Vector<Real>({0.3, 0.0, 0.2, 0.0}); // x=0.3, y=0, px=0.2, py=0
		}

		int getNumInvariants() const override { return 1; }
		std::string getInvariantName(int /*i*/) const override { return "Energy"; }
		Real computeInvariant(int /*i*/, const Vector<Real>& x) const override { return getEnergy(x); }

		/// @brief Get Hamiltonian value
		Real getEnergy(const Vector<Real>& y) const {
			Real x = y[0], yc = y[1], px = y[2], py = y[3];
			return 0.5 * (px * px + py * py + x * x + yc * yc) + x * x * yc - yc * yc * yc / 3.0;
		}
	};

	//=============================================================================
	// DOUBLE PENDULUM
	//=============================================================================

	/// @brief Double pendulum - classic example of chaos
	///
	/// State: [θ1, θ2, ω1, ω2] (angles and angular velocities)
	/// Parameters: [m1, m2, L1, L2, g]
	class DoublePendulumSystem : public DynamicalSystemBase<4, 5> {
	public:
		DoublePendulumSystem(Real m1 = 1.0, Real m2 = 1.0, Real L1 = 1.0, Real L2 = 1.0, Real gravity = 9.80665) {
			_params = {m1, m2, L1, L2, gravity};
			_stateNames = {"theta1", "theta2", "omega1", "omega2"};
			_paramNames = {"m1", "m2", "L1", "L2", "g"};
			_paramRanges = {{0.1, 10.0}, {0.1, 10.0}, {0.1, 5.0}, {0.1, 5.0}, {0.0, 20.0}};
		}

		void derivs(Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override {
			Real m1 = _params[0], m2 = _params[1];
			Real L1 = _params[2], L2 = _params[3];
			Real g = _params[4];

			Real th1 = y[0], th2 = y[1];
			Real w1 = y[2], w2 = y[3];
			Real dth = th1 - th2;

			Real c = std::cos(dth);
			Real s = std::sin(dth);

			// Denominators
			Real M = m1 + m2;
			Real d1 = L1 * (M - m2 * c * c);
			Real d2 = L2 * (M - m2 * c * c);

			// Numerators
			Real n1 = m2 * L1 * w1 * w1 * s * c + m2 * g * std::sin(th2) * c + m2 * L2 * w2 * w2 * s - M * g * std::sin(th1);

			Real n2 = -m2 * L2 * w2 * w2 * s * c + M * (g * std::sin(th1) * c - L1 * w1 * w1 * s - g * std::sin(th2));

			dydt[0] = w1;
			dydt[1] = w2;
			dydt[2] = n1 / d1;
			dydt[3] = n2 / d2;
		}

		void jacobian(Real /*t*/, const Vector<Real>& y, Matrix<Real>& J) const override {
			J.Resize(4, 4);
			// Numerical approximation for simplicity (analytical is very complex)
			Real eps = 1e-6;
			Vector<Real> f0(4), f1(4);
			Vector<Real> yp = y;

			derivs(0, y, f0);

			for (int j = 0; j < 4; ++j) {
				yp = y;
				yp[j] += eps;
				derivs(0, yp, f1);
				for (int i = 0; i < 4; ++i)
					J(i, j) = (f1[i] - f0[i]) / eps;
			}
		}

		bool isDissipative() const override { return false; }  // Conservative
		
		Vector<Real> getDefaultInitialCondition() const override { 
			return Vector<Real>({0.5, 0.5, 0.0, 0.0}); // Small angles, at rest
		}

		int getNumInvariants() const override { return 1; }
		std::string getInvariantName(int /*i*/) const override { return "Energy"; }
		Real computeInvariant(int /*i*/, const Vector<Real>& x) const override { return getEnergy(x); }

		/// @brief Total energy (should be conserved)
		Real getEnergy(const Vector<Real>& y) const {
			Real m1 = _params[0], m2 = _params[1];
			Real L1 = _params[2], L2 = _params[3];
			Real g = _params[4];

			Real th1 = y[0], th2 = y[1];
			Real w1 = y[2], w2 = y[3];

			// Kinetic energy
			Real T =
				0.5 * (m1 + m2) * L1 * L1 * w1 * w1 + 0.5 * m2 * L2 * L2 * w2 * w2 + m2 * L1 * L2 * w1 * w2 * std::cos(th1 - th2);

			// Potential energy (reference at y = 0)
			Real V = -(m1 + m2) * g * L1 * std::cos(th1) - m2 * g * L2 * std::cos(th2);

			return T + V;
		}
	};

} // namespace MML::Systems
#endif // MML_CONTINUOUS_SYSTEMS_H
