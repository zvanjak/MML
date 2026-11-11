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

#include <mml/MMLBase.h>

#include <mml/interfaces/IODESystem.h>

#include <mml/base/Matrix/Matrix.h>

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
				throw NotImplementedError("ODESystem::derivs() - system function is null (default-constructed ODESystem)");
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
	/// @details Function-pointer counterpart of ODESystem that also provides the Jacobian
	///          J[i][j] = ∂f_i/∂x_j. Implements IODESystemWithJacobian (and therefore IODESystem)
	///          through a single base chain, so instances can be passed both to the explicit
	///          solvers (IODESystem&) and to the implicit ones (IODESystemWithJacobian&:
	///          Backward Euler, BDF2, Rosenbrock, Radau).
	class ODESystemWithJacobian : public IODESystemWithJacobian
	{
	public:
		using DerivsFunc        = void (*)(Real, const Vector<Real>&, Vector<Real>&);
		using JacobianFunc      = void (*)(Real, const Vector<Real>&, Matrix<Real>&);                 ///< J(t, x, J) - the canonical form
		using DerivsJacobianFunc = void (*)(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&); ///< jacobn(t, x, dxdt, J) - NR style, both at once

	private:
		int _dim;
		DerivsFunc _func;
		JacobianFunc _funcJacOnly;
		DerivsJacobianFunc _funcJac;

	public:
		/// @brief Default constructor (creates empty system)
		ODESystemWithJacobian() : _dim(0), _func(nullptr), _funcJacOnly(nullptr), _funcJac(nullptr) { }
		
		/// @brief Construct ODE system with derivatives and a combined derivatives+Jacobian callback
		/// @param n System dimension (number of equations)
		/// @param inFunc Function pointer computing derivatives: f(t, x, dxdt)
		/// @param inFuncJac Function pointer computing derivatives and Jacobian together: jacobn(t, x, dxdt, J)
		ODESystemWithJacobian(int n, DerivsFunc inFunc, DerivsJacobianFunc inFuncJac)
			: _dim(n), _func(inFunc), _funcJacOnly(nullptr), _funcJac(inFuncJac) { }

		/// @brief Construct ODE system with derivatives and a Jacobian-only callback (canonical form)
		/// @param inFuncJac Function pointer computing the Jacobian: J(t, x, J)
		ODESystemWithJacobian(int n, DerivsFunc inFunc, JacobianFunc inFuncJac)
			: _dim(n), _func(inFunc), _funcJacOnly(inFuncJac), _funcJac(nullptr) { }

		/// @brief Get system dimension (number of equations)
		int		getDim() const override { return _dim; }

		/// @brief Compute derivatives at given point
		void	derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			if (_func == nullptr)
				throw NotImplementedError("ODESystemWithJacobian::derivs() - system function is null (default-constructed system)");
			_func(t, x, dxdt);
		}

		/// @brief Function call operator for computing derivatives
		void  operator()(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
		{
			derivs(t, x, dxdt);
		}

		/// @brief Compute the Jacobian J[i][j] = ∂f_i/∂x_j at (t, x)
		void jacobian(Real t, const Vector<Real>& x, Matrix<Real>& J) const override
		{
			if (_funcJacOnly != nullptr) {
				_funcJacOnly(t, x, J);
				return;
			}
			if (_funcJac == nullptr)
				throw NotImplementedError("ODESystemWithJacobian::jacobian() - no Jacobian function provided");
			Vector<Real> dxdt(_dim);   // the combined callback also produces the derivatives; discard them here
			_funcJac(t, x, dxdt, J);
		}

		/// @brief Derivatives and Jacobian together; one callback invocation when a combined callback was given
		void derivsAndJacobian(Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& J) const override
		{
			if (_funcJac != nullptr) {
				_funcJac(t, x, dxdt, J);
				return;
			}
			IODESystemWithJacobian::derivsAndJacobian(t, x, dxdt, J);   // derivs() + jacobian()
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
	/// @details Single base chain (IODESystemWithJacobian -> IODESystem): the previous dual
	///          inheritance from ODESystemFromStdFunc created two IODESystem subobjects and made
	///          conversion to IODESystem& ambiguous, so the system could not be handed to explicit solvers.
	class ODESystemWithJacobianFromStdFunc : public IODESystemWithJacobian
	{
	public:
		using DerivsFunc         = std::function<void(Real, const Vector<Real>&, Vector<Real>&)>;
		using JacobianFunc       = std::function<void(Real, const Vector<Real>&, Matrix<Real>&)>;                 ///< J(t, x, J) - canonical
		using DerivsJacobianFunc = std::function<void(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&)>; ///< jacobn(t, x, dxdt, J)

	private:
		int _dim;
		DerivsFunc _func;
		JacobianFunc _funcJacOnly;
		DerivsJacobianFunc _funcJac;

	public:
		/// @brief Derivatives plus a combined derivatives+Jacobian callable
		ODESystemWithJacobianFromStdFunc(int n, DerivsFunc inFunc, DerivsJacobianFunc inFuncJac)
			: _dim(n), _func(std::move(inFunc)), _funcJacOnly(), _funcJac(std::move(inFuncJac)) { }

		/// @brief Derivatives plus a Jacobian-only callable (canonical form)
		ODESystemWithJacobianFromStdFunc(int n, DerivsFunc inFunc, JacobianFunc inFuncJac)
			: _dim(n), _func(std::move(inFunc)), _funcJacOnly(std::move(inFuncJac)), _funcJac() { }

		int  getDim() const override { return _dim; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override { _func(t, x, dxdt); }
		void operator()(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const { _func(t, x, dxdt); }

		void jacobian(Real t, const Vector<Real>& x, Matrix<Real>& J) const override
		{
			if (_funcJacOnly) {
				_funcJacOnly(t, x, J);
				return;
			}
			if (!_funcJac)
				throw NotImplementedError("ODESystemWithJacobianFromStdFunc::jacobian() - no Jacobian function provided");
			Vector<Real> dxdt(_dim);
			_funcJac(t, x, dxdt, J);
		}

		void derivsAndJacobian(Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& J) const override
		{
			if (_funcJac) {
				_funcJac(t, x, dxdt, J);
				return;
			}
			IODESystemWithJacobian::derivsAndJacobian(t, x, dxdt, J);
		}
	};

}
#endif