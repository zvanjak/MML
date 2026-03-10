///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESystem.h                                                         ///
///  Description: Ordinary differential equation system representations               ///
///               Initial value problem definitions for ODE solvers                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_ODE_SYSTEM_H
#define MML_ODE_SYSTEM_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Matrix/Matrix.h"


namespace MML
{
	/// @brief Basic ODE system representation using function pointers
	/// @details Wraps a C-style function pointer for use with ODE solvers.
	///          The function computes derivatives dx/dt = f(t, x).
	class ODESystem : public IODESystem
	{
	protected:
		int _dim;
		void (*_func)(Real, const Vector<Real>&, Vector<Real>&);

	public:
		/// @brief Default constructor (creates empty system)
		ODESystem() 
				: _dim(0), _func(nullptr) { }
		
		/// @brief Construct ODE system with function pointer
		/// @param n System dimension (number of equations)
		/// @param inFunc Function pointer computing derivatives: f(t, x, dxdt)
		ODESystem(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&)) 
				: _dim(n), _func(inFunc) { }
		
		/// @brief Virtual destructor for proper cleanup in derived classes
		virtual ~ODESystem() = default;

		/// @brief Get system dimension (number of equations)
		int		getDim() const { return _dim; }
		
		/// @brief Compute derivatives at given point
		/// @param t Current time
		/// @param x Current state vector
		/// @param dxdt Output: computed derivatives dx/dt
		void	derivs(const Real t, const Vector<Real> &x, Vector<Real> &dxdt) const
		{
			_func(t, x, dxdt);
		}

		/// @brief Function call operator for computing derivatives
		/// @param t Current time
		/// @param x Current state vector  
		/// @param dxdt Output: computed derivatives dx/dt
		void  operator()(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
		{
			derivs(t, x, dxdt);
		}
	};

	/// @brief ODE system with analytical Jacobian for stiff solvers
	/// @details Extends ODESystem with Jacobian computation capability.
	///          Required for implicit methods like BDF or Rosenbrock.
	///          The Jacobian matrix J[i][j] = ∂f_i/∂x_j.
	class ODESystemWithJacobian : public ODESystem
	{
	private:
		void (*_funcJac)(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&);

	public:
		/// @brief Default constructor (creates empty system)
		ODESystemWithJacobian() : _funcJac(nullptr) { }
		
		/// @brief Construct ODE system with derivatives and Jacobian
		/// @param n System dimension (number of equations)
		/// @param inFunc Function pointer computing derivatives: f(t, x, dxdt)
		/// @param inFuncJac Function pointer computing Jacobian: J(t, x, dxdt, dydx)
		ODESystemWithJacobian(int n,
			void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&),
			void (*inFuncJac)(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		) : ODESystem(n, inFunc), _funcJac(inFuncJac) { }

		/// @brief Compute derivatives and Jacobian at given point
		/// @param t Current time
		/// @param x Current state vector
		/// @param dxdt Output: computed derivatives dx/dt
		/// @param dydx Output: Jacobian matrix J[i][j] = ∂f_i/∂x_j
		void jacobian(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx) const
		{
			if (_funcJac != nullptr)
				_funcJac(t, x, dxdt, dydx);
		}
	};

}
#endif