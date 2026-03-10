///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IODESystemDAE.h                                                     ///
///  Description: Interface for differential-algebraic equation (DAE) systems        ///
///               Defines evaluation protocol for DAE solvers                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IODESystemDAE.h
 * @brief Interfaces for differential-algebraic equation (DAE) systems.
 * 
 * Defines the contract for DAE systems that can be solved by MML's DAE solvers.
 * DAEs combine differential equations with algebraic constraints:
 * 
 * **Semi-Explicit Index-1 Form (supported):**
 * @f[
 *   \frac{d\mathbf{x}}{dt} = \mathbf{f}(t, \mathbf{x}, \mathbf{y})
 * @f]
 * @f[
 *   \mathbf{0} = \mathbf{g}(t, \mathbf{x}, \mathbf{y})
 * @f]
 * 
 * Where:
 * - **x** (n-dimensional): Differential variables (have time derivatives)
 * - **y** (m-dimensional): Algebraic variables (satisfy constraints)
 * - **f**: Differential equations computing dx/dt
 * - **g**: Algebraic constraint equations (must equal zero)
 * 
 * **Index-1 Requirement:**
 * The Jacobian ∂g/∂y must be nonsingular for the DAE to be index-1.
 * This ensures the algebraic variables y can be uniquely determined
 * from the constraints given x and t.
 * 
 * **Common Use Cases:**
 * - Constrained mechanical systems (pendulum, rigid bodies)
 * - Electrical circuits (Kirchhoff's laws as constraints)
 * - Chemical reaction networks (mass balance constraints)
 * - Power systems (load flow equations)
 * 
 * @see DAESolver, DAESolution, IODESystem
 */

#if !defined MML_IODESYSTEM_DAE_H
#define MML_IODESYSTEM_DAE_H

#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "interfaces/IParametrized.h"

namespace MML
{
	/**
	 * @brief Interface for semi-explicit index-1 DAE systems.
	 * 
	 * Represents a differential-algebraic system in semi-explicit form:
	 * @f[
	 *   \frac{d\mathbf{x}}{dt} = \mathbf{f}(t, \mathbf{x}, \mathbf{y})
	 * @f]
	 * @f[
	 *   \mathbf{0} = \mathbf{g}(t, \mathbf{x}, \mathbf{y})
	 * @f]
	 * 
	 * **Implementation Requirements:**
	 * - getDiffDim() returns the dimension of differential variables x
	 * - getAlgDim() returns the dimension of algebraic variables y
	 * - diffEqs() computes the differential equations f(t, x, y)
	 * - algConstraints() computes the constraint residuals g(t, x, y)
	 * 
	 * **Consistent Initial Conditions:**
	 * Initial values (x₀, y₀) must satisfy g(t₀, x₀, y₀) = 0.
	 * Use ComputeConsistentIC() utility to find consistent y₀ given x₀.
	 * 
	 * @note This interface supports index-1 DAEs only. Higher-index DAEs
	 *       require index reduction techniques before using these solvers.
	 * 
	 * @example
	 * @code
	 * // Simple pendulum in Cartesian coordinates
	 * // Differential: dx/dt = vx, dy/dt = vy, dvx/dt = λx, dvy/dt = λy - g
	 * // Constraint: x² + y² = L²
	 * class PendulumDAE : public IODESystemDAE {
	 *     Real L = 1.0, g = 9.81;
	 * public:
	 *     int getDiffDim() const override { return 4; }  // x, y, vx, vy
	 *     int getAlgDim() const override { return 1; }   // λ (tension/L)
	 *     
	 *     void diffEqs(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	 *                  Vector<Real>& dxdt) const override {
	 *         Real x = diff[0], y = diff[1], vx = diff[2], vy = diff[3];
	 *         Real lambda = alg[0];
	 *         dxdt[0] = vx;
	 *         dxdt[1] = vy;
	 *         dxdt[2] = lambda * x;
	 *         dxdt[3] = lambda * y - g;
	 *     }
	 *     
	 *     void algConstraints(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	 *                         Vector<Real>& constraints) const override {
	 *         Real x = diff[0], y = diff[1];
	 *         constraints[0] = x*x + y*y - L*L;
	 *     }
	 * };
	 * @endcode
	 */
	class IODESystemDAE
	{
	public:
		virtual ~IODESystemDAE() = default;
		
		//=========================================================================
		//                           Dimensions
		//=========================================================================
		
		/**
		 * @brief Get the dimension of differential variables.
		 * @return Number of differential variables (x) in the system
		 */
		virtual int getDiffDim() const = 0;
		
		/**
		 * @brief Get the dimension of algebraic variables.
		 * @return Number of algebraic variables (y) in the system
		 */
		virtual int getAlgDim() const = 0;
		
		/**
		 * @brief Get the total dimension (differential + algebraic).
		 * @return Total number of variables in the system
		 */
		int getTotalDim() const { return getDiffDim() + getAlgDim(); }

		//=========================================================================
		//                           Equations
		//=========================================================================
		
		/**
		 * @brief Compute the differential equations dx/dt = f(t, x, y).
		 * 
		 * Evaluates the right-hand side of the differential equations.
		 * 
		 * @param t Current time/independent variable
		 * @param x Current differential state vector (size = getDiffDim())
		 * @param y Current algebraic state vector (size = getAlgDim())
		 * @param[out] dxdt Output derivative vector dx/dt (size = getDiffDim())
		 */
		virtual void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                     Vector<Real>& dxdt) const = 0;

		/**
		 * @brief Compute the algebraic constraint residuals g(t, x, y).
		 * 
		 * Evaluates the constraint equations. For a valid solution,
		 * the output should be zero (or within numerical tolerance).
		 * 
		 * @param t Current time/independent variable
		 * @param x Current differential state vector (size = getDiffDim())
		 * @param y Current algebraic state vector (size = getAlgDim())
		 * @param[out] g Output constraint residual vector (size = getAlgDim())
		 */
		virtual void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                            Vector<Real>& g) const = 0;

		//=========================================================================
		//                           Variable Names
		//=========================================================================
		
		/**
		 * @brief Get a human-readable name for a differential variable.
		 * 
		 * Override to provide meaningful names for output and debugging
		 * (e.g., "position_x", "velocity", "angle").
		 * 
		 * @param i Index of the differential variable (0 to getDiffDim()-1)
		 * @return Name string for the variable
		 */
		virtual std::string getDiffVarName(int i) const
		{
			if (i < 0 || i >= getDiffDim())
				return "x?";
			return "x" + std::to_string(i + 1);
		}
		
		/**
		 * @brief Get a human-readable name for an algebraic variable.
		 * 
		 * Override to provide meaningful names for output and debugging
		 * (e.g., "tension", "current", "lagrange_mult").
		 * 
		 * @param i Index of the algebraic variable (0 to getAlgDim()-1)
		 * @return Name string for the variable
		 */
		virtual std::string getAlgVarName(int i) const
		{
			if (i < 0 || i >= getAlgDim())
				return "y?";
			return "y" + std::to_string(i + 1);
		}
	};

	/**
	 * @brief DAE system with analytic Jacobian matrices.
	 * 
	 * Extends IODESystemDAE to provide the four Jacobian matrices needed
	 * for efficient implicit integration:
	 * 
	 * - **∂f/∂x**: How differential equations depend on differential variables
	 * - **∂f/∂y**: How differential equations depend on algebraic variables
	 * - **∂g/∂x**: How constraints depend on differential variables
	 * - **∂g/∂y**: How constraints depend on algebraic variables (MUST be nonsingular)
	 * 
	 * The combined system Jacobian for Newton iteration is:
	 * @f[
	 *   \begin{bmatrix}
	 *     I - h\frac{\partial f}{\partial x} & -h\frac{\partial f}{\partial y} \\
	 *     \frac{\partial g}{\partial x} & \frac{\partial g}{\partial y}
	 *   \end{bmatrix}
	 * @f]
	 * 
	 * **Index-1 Requirement:**
	 * For the DAE to be index-1, the Jacobian ∂g/∂y (returned by jacobian_gy)
	 * must be nonsingular. If this matrix is singular, the DAE has index > 1
	 * and requires index reduction before solving.
	 * 
	 * @note Providing analytic Jacobians significantly improves both accuracy
	 *       and performance compared to numerical differentiation.
	 */
	class IODESystemDAEWithJacobian : public IODESystemDAE
	{
	public:
		/**
		 * @brief Compute Jacobian ∂f/∂x (differential eqs w.r.t. differential vars).
		 * 
		 * @param t Current time
		 * @param x Differential state vector
		 * @param y Algebraic state vector
		 * @param[out] df_dx Output Jacobian matrix (getDiffDim() × getDiffDim())
		 */
		virtual void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                         Matrix<Real>& df_dx) const = 0;
		
		/**
		 * @brief Compute Jacobian ∂f/∂y (differential eqs w.r.t. algebraic vars).
		 * 
		 * @param t Current time
		 * @param x Differential state vector
		 * @param y Algebraic state vector
		 * @param[out] df_dy Output Jacobian matrix (getDiffDim() × getAlgDim())
		 */
		virtual void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                         Matrix<Real>& df_dy) const = 0;
		
		/**
		 * @brief Compute Jacobian ∂g/∂x (constraints w.r.t. differential vars).
		 * 
		 * @param t Current time
		 * @param x Differential state vector
		 * @param y Algebraic state vector
		 * @param[out] dg_dx Output Jacobian matrix (getAlgDim() × getDiffDim())
		 */
		virtual void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                         Matrix<Real>& dg_dx) const = 0;
		
		/**
		 * @brief Compute Jacobian ∂g/∂y (constraints w.r.t. algebraic vars).
		 * 
		 * **CRITICAL:** This matrix MUST be nonsingular for index-1 DAEs.
		 * If singular, the DAE has index > 1 and cannot be solved directly.
		 * 
		 * @param t Current time
		 * @param x Differential state vector
		 * @param y Algebraic state vector
		 * @param[out] dg_dy Output Jacobian matrix (getAlgDim() × getAlgDim())
		 */
		virtual void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                         Matrix<Real>& dg_dy) const = 0;
		
		/**
		 * @brief Compute all Jacobians simultaneously (optional optimization).
		 * 
		 * Override this if computing multiple Jacobians shares significant
		 * computation. Default implementation calls individual jacobian_* methods.
		 * 
		 * @param t Current time
		 * @param x Differential state vector
		 * @param y Algebraic state vector
		 * @param[out] df_dx Jacobian ∂f/∂x
		 * @param[out] df_dy Jacobian ∂f/∂y
		 * @param[out] dg_dx Jacobian ∂g/∂x
		 * @param[out] dg_dy Jacobian ∂g/∂y
		 */
		virtual void allJacobians(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                          Matrix<Real>& df_dx, Matrix<Real>& df_dy,
		                          Matrix<Real>& dg_dx, Matrix<Real>& dg_dy) const
		{
			jacobian_fx(t, x, y, df_dx);
			jacobian_fy(t, x, y, df_dy);
			jacobian_gx(t, x, y, dg_dx);
			jacobian_gy(t, x, y, dg_dy);
		}
	};

	/**
	 * @brief DAE system with adjustable parameters.
	 * 
	 * Extends IODESystemDAE to support systems where coefficients or parameters
	 * can be modified at runtime. Useful for parameter sweeps, sensitivity
	 * analysis, and fitting DAE models to data.
	 * 
	 * @example
	 * @code
	 * // Pendulum with adjustable length and gravity
	 * class ParametrizedPendulum : public IODESystemDAEParametrized {
	 *     Real L = 1.0, g = 9.81;  // Default values
	 * public:
	 *     int getNumParam() const override { return 2; }
	 *     Real getParam(int i) const override { return i == 0 ? L : g; }
	 *     void setParam(int i, Real val) override { 
	 *         if (i == 0) L = val; else g = val; 
	 *     }
	 *     // ... rest of DAE implementation
	 * };
	 * @endcode
	 */
	class IODESystemDAEParametrized : public IODESystemDAE, public IParametrized
	{
	public:
		virtual ~IODESystemDAEParametrized() = default;
	};

} // namespace MML

#endif // MML_IODESYSTEM_DAE_H
