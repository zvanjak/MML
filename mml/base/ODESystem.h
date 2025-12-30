///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESystem.h                                                         ///
///  Description: Ordinary differential equation system representations               ///
///               Initial value problem definitions for ODE solvers                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_ODE_SYSTEM_H
#define MML_ODE_SYSTEM_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Matrix.h"


namespace MML
{
	class ODESystem : public IODESystem
	{
	protected:
		int _dim;
		void (*_func)(Real, const Vector<Real>&, Vector<Real>&);

	public:
		ODESystem() 
				: _dim(0), _func(nullptr) { }
		ODESystem(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&)) 
				: _dim(n), _func(inFunc) { }
		virtual ~ODESystem() = default;

		int		getDim() const { return _dim; }
		void	derivs(const Real t, const Vector<Real> &x, Vector<Real> &dxdt) const
		{
			_func(t, x, dxdt);
		}

		void  operator()(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
		{
			derivs(t, x, dxdt);
		}
	};

	class ODESystemWithJacobian : public ODESystem
	{
	private:
		void (*_funcJac)(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&);

	public:
		ODESystemWithJacobian() : _funcJac(nullptr) { }
		ODESystemWithJacobian(int n,
			void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&),
			void (*inFuncJac)(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		) : ODESystem(n, inFunc), _funcJac(inFuncJac) { }

		void jacobian(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx) const
		{
			if (_funcJac != nullptr)
				_funcJac(t, x, dxdt, dydx);
		}
	};

}
#endif