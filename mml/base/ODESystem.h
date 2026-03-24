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

#include <functional>


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
			if (_func == nullptr)
				throw std::runtime_error("ODESystem::derivs() - system function is null (default-constructed ODESystem)");
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
			if (_funcJac == nullptr)
				throw NotImplementedError("ODESystemWithJacobian::jacobian() - no Jacobian function provided");
			_funcJac(t, x, dxdt, dydx);
		}
	};

	/// @brief ODE system using std::function for lambdas with captured state
	/// @details Wraps a std::function callable, enabling lambdas, functors, and
	///          std::bind expressions as ODE right-hand sides.
	class ODESystemFromStdFunc : public IODESystem
	{
	private:
		int _dim;
		std::function<void(Real, const Vector<Real>&, Vector<Real>&)> _func;

	public:
		ODESystemFromStdFunc(int n, std::function<void(Real, const Vector<Real>&, Vector<Real>&)> inFunc)
			: _dim(n), _func(std::move(inFunc)) { }

		int  getDim() const { return _dim; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const { _func(t, x, dxdt); }
		void operator()(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const { _func(t, x, dxdt); }
	};

	/// @brief ODE system with Jacobian using std::function callables
	class ODESystemWithJacobianFromStdFunc : public ODESystemFromStdFunc, public IODESystemWithJacobian
	{
	private:
		std::function<void(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&)> _funcJac;

	public:
		ODESystemWithJacobianFromStdFunc(int n,
			std::function<void(Real, const Vector<Real>&, Vector<Real>&)> inFunc,
			std::function<void(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&)> inFuncJac)
			: ODESystemFromStdFunc(n, std::move(inFunc)), _funcJac(std::move(inFuncJac)) { }

		// Disambiguate IODESystem methods inherited from both bases
		int  getDim() const override { return ODESystemFromStdFunc::getDim(); }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override { ODESystemFromStdFunc::derivs(t, x, dxdt); }

		void jacobian(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx) const override
		{
			_funcJac(t, x, dxdt, dydx);
		}
	};

}
#endif