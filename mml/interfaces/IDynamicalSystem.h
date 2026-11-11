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

#include <mml/MMLBase.h>
#include <algorithm>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/interfaces/IODESystem.h>
#include <mml/base/detail/NumericalJacobian.h>

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
	class IDynamicalSystem : public IODESystemParametrized, public IODESystemWithJacobian {
	public:
		virtual ~IDynamicalSystem() = default;

		//=========================================================================
		// METADATA
		//=========================================================================

		/** @brief Get name of state variable */
		virtual std::string getStateName(int i) const {
			if (i < 0 || i >= getDim())
				throw IndexError("IDynamicalSystem::getStateName - state index out of range");
			return "x" + std::to_string(i);
		}

		/** @brief Get name of parameter */
		virtual std::string getParamName(int i) const {
			if (i < 0 || i >= getNumParam())
				throw IndexError("IDynamicalSystem::getParamName - parameter index out of range");
			return "p" + std::to_string(i);
		}

		/** @brief Get a suggested exploration range, not a model admissibility constraint. */
		virtual std::pair<Real, Real> getParamRange(int i) const {
			if (i < 0 || i >= getNumParam())
				throw IndexError("IDynamicalSystem::getParamRange - parameter index out of range");
			return {-Precision::UnboundedIntervalSamplingLimit, Precision::UnboundedIntervalSamplingLimit};
		}

		/** @brief Get default initial condition */
		virtual Vector<Real> getDefaultInitialCondition() const { return Vector<Real>(getDim(), 0.0); }

		//=========================================================================
		// JACOBIAN
		//=========================================================================

		/** @brief Does system provide analytical Jacobian? Default: no - jacobian() below is finite differences. */
		bool hasAnalyticalJacobian() const override { return false; }

		/** @brief Compute Jacobian matrix at point (x, t) - the IODESystemWithJacobian contract.
		 *  Default: central finite differences of derivs(); override with the analytic form and
		 *  return true from hasAnalyticalJacobian().
		 *  Validates the dimension, state, time and all stencils before callbacks. Each derivs()
		 *  call must write exactly getDim() finite entries; omitted entries are rejected.
		 *  Unrepresentable stencils/results throw NumericalMethodError. J is unchanged on failure. */
		void jacobian(Real t, const Vector<Real>& x, Matrix<Real>& J) const override {
			const int n = getDim();
			if (n <= 0)
				throw ArgumentError("IDynamicalSystem::jacobian - dimension must be positive");
			if (x.size() != n)
				throw VectorDimensionError("IDynamicalSystem::jacobian - state dimension mismatch", n, x.size());
			NumericalJacobianDetail::JacobianODE(
				[this](Real time, const Vector<Real>& state, Vector<Real>& output) { derivs(time, state, output); },
				t, x, J);
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
		virtual std::string getInvariantName(int i) const {
			if (i < 0 || i >= getNumInvariants())
				throw IndexError("IDynamicalSystem::getInvariantName - invariant index out of range");
			return "";
		}

		/** @brief Compute invariant value at state x */
		virtual Real computeInvariant(int i, const Vector<Real>& x) const {
			if (i < 0 || i >= getNumInvariants())
				throw IndexError("IDynamicalSystem::computeInvariant - invariant index out of range");
			return 0.0;
		}
	};

} // namespace MML

#endif // MML_IDYNAMICALSYSTEM_H
