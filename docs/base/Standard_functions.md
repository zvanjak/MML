# Standard functions

Defined set of standard functions

~~~ c++
    namespace StdFunctions
    {
        static inline Real Sin(Real x) { return sin(x); }
        static inline Real Cos(Real x) { return cos(x); }
        static inline Real Tan(Real x) { return tan(x); }
        static inline Real Exp(Real x) { return exp(x); }
        static inline Real Log(Real x) { return log(x); }
        static inline Real Sqrt(Real x) { return sqrt(x); }
        static inline Real Pow(Real x, Real y) { return pow(x, y); }
        static inline Real Sinh(Real x) { return sinh(x); }
        static inline Real Cosh(Real x) { return cosh(x); }
        static inline Real Tanh(Real x) { return tanh(x); }
        static inline Real Asin(Real x) { return asin(x); }
        static inline Real Acos(Real x) { return acos(x); }
        static inline Real Atan(Real x) { return atan(x); }
        static inline Real Asinh(Real x) { return asinh(x); }
        static inline Real Acosh(Real x) { return acosh(x); }
        static inline Real Atanh(Real x) { return atanh(x); }

        static inline Real Erf(Real x)  { return std::erf(x); }
        static inline Real Erfc(Real x) { return std::erfc(x); }

        static inline Real TGamma(Real x) { return std::tgamma(x); }
        static inline Real LGamma(Real x) { return std::lgamma(x); }
        static inline Real RiemannZeta(Real x) { return std::riemann_zeta(x); }
        static inline Real CompEllint_1(Real x) { return std::comp_ellint_1(x); }
        static inline Real CompEllint_2(Real x) { return std::comp_ellint_2(x); }

        static inline Real Hermite(unsigned int n, Real x) { return std::hermite(n, x); }
        static inline Real Legendre(unsigned int n, Real x) { return std::legendre(n, x); }
        static inline Real Laguerre(unsigned int n, Real x) { return std::laguerre(n, x); }
        static inline Real SphBessel(unsigned int n, Real x) { return std::sph_bessel(n, x); }
        static inline Real SphLegendre(int n1, int n2, Real x) { return std::sph_legendre(n1, n2, x); }
    }
~~~



