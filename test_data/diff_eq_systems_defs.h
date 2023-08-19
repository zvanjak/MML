#if !defined MML_ODE_SYSTEM_DEFS_H
#define MML_ODE_SYSTEM_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Vector.h"
#include "core/Matrix.h"
#include "basic_types/ODESystem.h"
#endif

namespace MML::TestBeds
{
    // TODO - add Chebyshev diff eq system
    // TODO - finish implementation of Legandre, Hermite, Laguerre
    // TODO - add couple of (complex) linear diff systems
    // TODO - add at least 5 stiff systems
    class LegandreODE : public IODESystem
    {
        int _n; // order
    public:
        LegandreODE(int n) : _n(n) {}
        
        int  getDim() { return 2; }
        void derivs(const double x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) override
        {
            // dydx[0]= y[1];
            // dydx[1]=((1.0-y[0]*y[0])*y[1]-y[0])/_eps;
        }
    };
    class HermiteODE : public IODESystem
    {
        int _n; // order
    public:
        HermiteODE(int n) : _n(n) {}
        
        int  getDim() { return 2; }
        void derivs(const double x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) override
        {
            // dydx[0]= y[1];
            // dydx[1]=((1.0-y[0]*y[0])*y[1]-y[0])/_eps;
        }
    };

    class VanDerPolODE : IODESystem
    {
        double _eps;
    public:
        VanDerPolODE(double eps) : _eps(eps) {}
        
        int  getDim() { return 2; }
        void derivs(const double x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) override
        {
            dydx[0] = y[1];
            dydx[1] = ( (1.0 - y[0]*y[0]) * y[1] - y[0] ) / _eps;
        }
    };

    void  fnc12(double t, const Vector<Real> &x, Vector<Real> &ret)
    {
        ret[0] =  x[0] +   x[1] - x[2];
        ret[1] = -x[0] + 3*x[1] - x[2];
        ret[2] = -x[0] +   x[1] + x[2];
    }
    Vector<Real>  fnc12_sol(double t)
    {
        Vector<Real> ret(3);
        ret[0] = exp(t);
        ret[1] = exp(t) + exp(2*t);
        ret[2] = exp(t) + exp(2*t);
        return ret;
    }
    void VanDerPol(double eps, const double x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) {
        dydx[0]= y[1];
        dydx[1]=((1.0-y[0]*y[0])*y[1]-y[0])/eps;
    }
    void VanDerPolEps0_1(const double x, const MML::Vector<Real> &y, MML::Vector<Real> &dydx) { return VanDerPol(0.1, x, y, dydx); }

    void  stiff_sys1_derivs(double t, const Vector<Real> &x, Vector<Real> &dydx)
    {
        dydx[0] = -0.013*x[0]-1000.0*x[0]*x[2];
        dydx[1] = -2500.0*x[1]*x[2];
        dydx[2] = -0.013*x[0]-1000.0*x[0]*x[2]-2500.0*x[1]*x[2];
    }

    void  stiff_sys1_jac(const double t, const Vector<Real> &x, Vector<Real> &dxdt, Matrix<Real> &dydx)
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