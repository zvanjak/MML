#if !defined MPL_OSCILLATORS_H
#define MPL_OSCILLATORS_H


namespace MPL
{
	class DuffingOscillator : public IODESystem
	{
	private:
		Real _alpha;  // linear stiffness
		Real _beta;   // nonlinear stiffness
		Real _gamma;  // damping coefficient
		
		Real _F;			// amplitude of the periodic driving force
		Real _omega;  // frequency of the periodic driving force

	public:
		DuffingOscillator(Real alpha, Real beta, Real gamma, Real F, Real omega)
			: _alpha(alpha), _beta(beta), _gamma(gamma), _F(F), _omega(omega)
		{	}
		
		int  getDim() const override { return 2; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];
			dxdt[1] = -2.0*_gamma * x[1] - _alpha * x[0] - _beta * POW3(x[0]) + _F * cos(_omega * t);
		}
	};

	class VanDerPolOscillator : public IODESystem
	{
		private:
			Real _mu;  // nonlinearity parameter

	public:
		VanDerPolOscillator(Real mu) : _mu(mu) { }
		
		int  getDim() const override { return 2; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];
			dxdt[1] = _mu * (1.0 - x[0] * x[0]) * x[1] - x[0];
		}
	};
}

#endif 