///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IDynamicalSystem.h                                                  ///
///  Description: Interface for dynamical systems analysis                            ///
///               Extends IODESystemParametrized with metadata and Jacobian           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IDynamicalSystem.h
 * @brief Interface for dynamical systems with advanced analysis capabilities.
 * 
 * Defines the contract for dynamical systems that support:
 * - Parameter and state variable metadata (names, ranges)
 * - Analytical or numerical Jacobian computation
 * - System property flags (autonomous, Hamiltonian, dissipative)
 * - Conserved quantities (invariants)
 * 
 * This interface extends IODESystemParametrized to enable advanced analysis
 * features like fixed point finding, Lyapunov exponent computation, and
 * bifurcation analysis.
 * 
 * @see IODESystemParametrized, DynamicalSystemBase, LorenzSystem
 */

#if !defined MML_IDYNAMICALSYSTEM_H
#define MML_IDYNAMICALSYSTEM_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "interfaces/IODESystem.h"

#include <string>
#include <utility>
#include <limits>
#include <cmath>

namespace MML
{
	//=============================================================================
	// IDYNAMICALSYSTEM INTERFACE
	//=============================================================================

	/// @brief Extended interface for dynamical systems analysis
	///
	/// Extends IODESystemParametrized with:
	/// - Parameter metadata (names, ranges)
	/// - State variable metadata
	/// - Optional analytical Jacobian
	/// - System property flags
	class IDynamicalSystem : public IODESystemParametrized {
	public:
		virtual ~IDynamicalSystem() = default;

		//=========================================================================
		// METADATA
		//=========================================================================

		/** @brief Get name of state variable */
		virtual std::string getStateName(int i) const { return "x" + std::to_string(i); }

		/** @brief Get name of parameter */
		virtual std::string getParamName(int i) const { return "p" + std::to_string(i); }

		/** @brief Get valid range for parameter */
		virtual std::pair<Real, Real> getParamRange(int i) const {
			return {-1e10, 1e10}; // Default: no restriction
		}

		/** @brief Get default initial condition */
		virtual Vector<Real> getDefaultInitialCondition() const { return Vector<Real>(getDim(), 0.0); }

		//=========================================================================
		// JACOBIAN
		//=========================================================================

		/** @brief Does system provide analytical Jacobian? */
		virtual bool hasAnalyticalJacobian() const { return false; }

		/** @brief Compute Jacobian matrix at point (x, t) */
		virtual void jacobian(Real t, const Vector<Real>& x, Matrix<Real>& J) const {
			// Default: numerical Jacobian via finite differences
			int n = getDim();
			J.Resize(n, n);
			Vector<Real> f1(n), f2(n), xp(x);

			Real eps = std::sqrt(std::numeric_limits<Real>::epsilon());

			for (int j = 0; j < n; ++j) {
				Real h = eps * std::max<Real>(std::abs(x[j]), Real(1.0));
				xp[j] = x[j] + h;
				derivs(t, xp, f1);
				xp[j] = x[j] - h;
				derivs(t, xp, f2);
				xp[j] = x[j];

				for (int i = 0; i < n; ++i)
					J(i, j) = (f1[i] - f2[i]) / (2.0 * h);
			}
		}

		//=========================================================================
		// SYSTEM PROPERTIES
		//=========================================================================

		/** @brief Is the system autonomous (no explicit time dependence)? */
		virtual bool isAutonomous() const { return true; }

		/** @brief Is the system Hamiltonian (energy conserving)? */
		virtual bool isHamiltonian() const { return false; }

		/** @brief Is the system dissipative (contracting phase space)? */
		virtual bool isDissipative() const { return false; }

		//=========================================================================
		// INVARIANTS
		//=========================================================================

		/** @brief Number of conserved quantities */
		virtual int getNumInvariants() const { return 0; }

		/** @brief Name of invariant */
		virtual std::string getInvariantName(int i) const { return ""; }

		/** @brief Compute invariant value at state x */
		virtual Real computeInvariant(int i, const Vector<Real>& x) const { return 0.0; }
	};

} // namespace MML

#endif // MML_IDYNAMICALSYSTEM_H
