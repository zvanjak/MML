#if !defined MML_ODE_SYSTEM_DEFS_H
#define MML_ODE_SYSTEM_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "core/ODESystem.h"
#endif

namespace MML::TestBeds
{
    // TODO - add Chebyshev diff eq system
    // TODO - add couple of (complex) linear diff systems
    // TODO - add at least 5 stiff systems
    class LegandreODE : public IODESystem
    {
        int _n; // order
    public:
        LegandreODE(int n) : _n(n) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) const override
        {
            dydx[0] = y[1];
            dydx[1] = 2 * x / (1 - x * x) * y[1] - _n * (_n + 1) * y[0];
        }
    };

    class LaguerreODE : public IODESystem
    {
        int _n; // order
    public:
        LaguerreODE(int n) : _n(n) {}
        
        int  getDim() const override { return 2; }
        void derivs(const Real x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) const override
        {
            dydx[0] = y[1];
            dydx[1] = (x - 1) / x * y[1] - _n / x * y[0];
        }
    };

    class HermiteODE : public IODESystem
    {
        int _n; // order
    public:
        HermiteODE(int n) : _n(n) {}
        
        int  getDim() const override { return 2; }
        void derivs(const Real x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) const override
        {
            dydx[0] = y[1];
            dydx[1] = 2 * x * y[1] - 2 * _n * y[0];
        }
    };

    static void VanDerPol(Real mju, const Real x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) {
        dydx[0] = y[1];
        dydx[1] = mju * (1.0-y[0]*y[0])*y[1]-y[0];
    }
    static void VanDerPolMju0_1(const Real x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) { return VanDerPol(0.1, x, y, dydx); }

    class VanDerPolODE : public IODESystem
    {
        Real _eps;
    public:
        VanDerPolODE(Real eps) : _eps(eps) {}
        
        int  getDim() const override { return 2; }
        void derivs(const Real x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) const override
        {
            dydx[0] = y[1];
            dydx[1] = ( (1.0 - y[0]*y[0]) * y[1] - y[0] ) / _eps;
        }
    };

    class LorenzSystemODE : public IODESystem
    {
        Real _sigma, _rho, _beta;
    public:
        LorenzSystemODE(Real sigma, Real rho, Real beta) : _sigma(sigma), _rho(rho), _beta(beta) {}
        
        int  getDim() const override { return 3; }
        void derivs(const Real x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) const override
        {
            dydx[0] = _sigma * (y[1] - y[0]);
            dydx[1] = y[0] * (_rho - y[2]) - y[1];
            dydx[2] = y[0] * y[1] - _beta * y[2];
        }
    };

    static void  fnc12(Real t, const Vector<Real> &x, Vector<Real> &ret)
    {
        ret[0] =  x[0] +   x[1] - x[2];
        ret[1] = -x[0] + 3*x[1] - x[2];
        ret[2] = -x[0] +   x[1] + x[2];
    }
    static Vector<Real>  fnc12_sol(Real t)
    {
        Vector<Real> ret(3);
        ret[0] = exp(t);
        ret[1] = exp(t) + exp(2*t);
        ret[2] = exp(t) + exp(2*t);
        return ret;
    }

    static void  stiff_sys1_derivs(Real t, const Vector<Real> &x, Vector<Real> &dydx)
    {
        dydx[0] = -0.013*x[0]-1000.0*x[0]*x[2];
        dydx[1] = -2500.0*x[1]*x[2];
        dydx[2] = -0.013*x[0]-1000.0*x[0]*x[2]-2500.0*x[1]*x[2];
    }

    static void  stiff_sys1_jac(const Real t, const Vector<Real> &x, Vector<Real> &dxdt, Matrix<Real> &dydx)
    {
        int n= (int) x.size();

        for (int i=0;i<n;i++) 
            dxdt[i]=0.0;

        dydx[0][0] = -0.013-1000.0*x[2];
        dydx[0][1] = 0.0;
        dydx[0][2] = -1000.0*x[0];
        dydx[1][0] = 0.0;
        dydx[1][1] = -2500.0*x[2];
        dydx[1][2] = -2500.0*x[1];
        dydx[2][0] = -0.013-1000.0*x[2];
        dydx[2][1] = -2500.0*x[2];
        dydx[2][2] = -1000.0*x[0]-2500.0*x[1];
    }   
}
#endif