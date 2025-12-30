#if !defined __MML_SCALAR_FUNCTIONS_TEST_BED_H
#define __MML_SCALAR_FUNCTIONS_TEST_BED_H

#include <string>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Function.h"
#endif

namespace MML::TestBeds
{
    template<int N>
    struct TestFunctionScalar
    {
        std::string _funcName;

        ScalarFunction<N> _func;
        Real (*_funcDerived)(const VectorN<Real, N> &, int ind);

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        // gradijent

        TestFunctionScalar(std::string funcName,
                            Real (*f1)(const VectorN<Real, N> &), std::string funcExpr, 
                            Real (*f2)(const VectorN<Real, N> &, int ind), std::string funcDerivedExpr
                            ) : _funcName(funcName),
                                _func(f1), _funcDerived(f2), 
                                _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}
    };    

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                        ORIGINAL TEST FUNCTIONS (2)                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    // Function 1: Separable trigonometric + exponential
    // f(x,y,z) = cos(x) + sin(y) + exp(z)
    static Real TestScalarFunc1(const VectorN<Real, 3> &x) { return cos(x[0]) + sin(x[1]) + exp(x[2]); }
    static Real TestScalarFunc1_derived(const VectorN<Real, 3> &x, int ind) 
    { 
        if( ind == 0 ) return -sin(x[0]);
        else if( ind == 1 ) return cos(x[1]);
        else return exp(x[2]);
    }

    // Function 2: Complex mixed function
    // f(x,y,z) = sin(xy)*exp(z/(y²+1))/(1+x²)
    static Real TestScalarFunc2(const VectorN<Real, 3> &x) { return sin(x[0] * x[1]) * exp(x[2] / (x[1] * x[1] +1 )) / (1 + x[0] * x[0]); }
    static Real TestScalarFunc2_derived(const VectorN<Real, 3> &xVal, int ind) 
    { 
        Real x = xVal[0];
        Real y = xVal[1];
        Real z = xVal[2];
        if( ind == 0 ) return (exp(z / (y*y + 1)) * ((x*x + 1) * y * cos(x * y) - 2 * x * sin(x * y)))/pow((x*x + 1),2);
        else if( ind == 1 ) return (exp(z/(y*y + 1)) * (x * pow((y*y + 1),2) * cos(x * y) - 2 * y * z * sin(x * y)))/((x*x + 1) * pow((y*y + 1),2));
        else return (exp(z/(y*y + 1)) * sin(x * y))/((x*x + 1) * (y*y + 1));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    QUADRATIC FORMS AND BASIC FUNCTIONS (3-5)                          //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    // Function 3: Sphere function (simple quadratic)
    // f(x,y,z) = x² + y² + z²
    // Classic optimization test - convex, smooth, single minimum at origin
    static Real TestScalarFunc_Sphere(const VectorN<Real, 3> &v) 
    { 
        return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; 
    }
    static Real TestScalarFunc_Sphere_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        return 2 * v[ind];
    }

    // Function 4: Saddle function
    // f(x,y,z) = x² - y² + z²
    // Has a saddle point at origin - tests optimization in non-convex settings
    static Real TestScalarFunc_Saddle(const VectorN<Real, 3> &v) 
    { 
        return v[0]*v[0] - v[1]*v[1] + v[2]*v[2]; 
    }
    static Real TestScalarFunc_Saddle_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        if (ind == 0) return 2 * v[0];
        else if (ind == 1) return -2 * v[1];
        else return 2 * v[2];
    }

    // Function 5: Weighted quadratic form
    // f(x,y,z) = x² + 4y² + 9z²  (ellipsoidal bowl)
    // Tests optimization with different scaling in different directions
    static Real TestScalarFunc_WeightedQuadratic(const VectorN<Real, 3> &v) 
    { 
        return v[0]*v[0] + 4*v[1]*v[1] + 9*v[2]*v[2]; 
    }
    static Real TestScalarFunc_WeightedQuadratic_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        if (ind == 0) return 2 * v[0];
        else if (ind == 1) return 8 * v[1];
        else return 18 * v[2];
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CLASSIC OPTIMIZATION TEST FUNCTIONS (6-8)                          //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Function 6: 3D Rosenbrock function (extended)
    // f(x,y,z) = (1-x)² + 100(y-x²)² + (1-y)² + 100(z-y²)²
    // Famous "banana" function - narrow curved valley, very hard to optimize
    static Real TestScalarFunc_Rosenbrock(const VectorN<Real, 3> &v) 
    { 
        Real x = v[0], y = v[1], z = v[2];
        return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x) + (1-y)*(1-y) + 100*(z-y*y)*(z-y*y);
    }
    static Real TestScalarFunc_Rosenbrock_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        Real x = v[0], y = v[1], z = v[2];
        if (ind == 0) return -2*(1-x) - 400*x*(y-x*x);
        else if (ind == 1) return 200*(y-x*x) - 2*(1-y) - 400*y*(z-y*y);
        else return 200*(z-y*y);
    }

    // Function 7: 3D Rastrigin function
    // f(x,y,z) = 30 + (x²-10cos(2πx)) + (y²-10cos(2πy)) + (z²-10cos(2πz))
    // Highly multimodal - many local minima, tests global optimization
    static Real TestScalarFunc_Rastrigin(const VectorN<Real, 3> &v) 
    { 
        const Real A = 10.0;
        const Real pi2 = 2 * 3.14159265358979323846;
        return 3*A + (v[0]*v[0] - A*cos(pi2*v[0])) 
                   + (v[1]*v[1] - A*cos(pi2*v[1])) 
                   + (v[2]*v[2] - A*cos(pi2*v[2]));
    }
    static Real TestScalarFunc_Rastrigin_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        const Real A = 10.0;
        const Real pi2 = 2 * 3.14159265358979323846;
        return 2*v[ind] + A*pi2*sin(pi2*v[ind]);
    }

    // Function 8: 3D Ackley function
    // Complex multimodal with large nearly flat outer region
    // Tests optimizer ability to escape nearly flat regions
    static Real TestScalarFunc_Ackley(const VectorN<Real, 3> &v) 
    { 
        const Real a = 20.0;
        const Real b = 0.2;
        const Real c = 2 * 3.14159265358979323846;
        Real sum1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        Real sum2 = cos(c*v[0]) + cos(c*v[1]) + cos(c*v[2]);
        return -a * exp(-b * sqrt(sum1/3.0)) - exp(sum2/3.0) + a + exp(1.0);
    }
    static Real TestScalarFunc_Ackley_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        const Real a = 20.0;
        const Real b = 0.2;
        const Real c = 2 * 3.14159265358979323846;
        Real sum1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        Real sum2 = cos(c*v[0]) + cos(c*v[1]) + cos(c*v[2]);
        Real r = sqrt(sum1/3.0);
        // df/dxi = (a*b*xi)/(3*r) * exp(-b*r) + (c/3)*sin(c*xi)*exp(sum2/3)
        if (r < 1e-12) return (c/3.0)*sin(c*v[ind])*exp(sum2/3.0);
        return (a*b*v[ind])/(3.0*r) * exp(-b*r) + (c/3.0)*sin(c*v[ind])*exp(sum2/3.0);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         PHYSICAL POTENTIALS (9-12)                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Function 9: Gravitational/Coulomb potential (1/r form)
    // f(x,y,z) = -1/sqrt(x²+y²+z²)
    // Central force potential - singular at origin
    static Real TestScalarFunc_CoulombPotential(const VectorN<Real, 3> &v) 
    { 
        Real r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        if (r2 < 1e-20) return -1e10;  // Regularize singularity
        return -1.0 / sqrt(r2);
    }
    static Real TestScalarFunc_CoulombPotential_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        Real r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        if (r2 < 1e-20) return 0;  // Regularize singularity
        Real r3 = pow(r2, 1.5);
        return v[ind] / r3;
    }

    // Function 10: Harmonic oscillator potential
    // f(x,y,z) = 0.5*(x²+y²+z²)  (quantum harmonic oscillator)
    // Simple quadratic well - used in physics as basis for perturbation theory
    static Real TestScalarFunc_HarmonicPotential(const VectorN<Real, 3> &v) 
    { 
        return 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    static Real TestScalarFunc_HarmonicPotential_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        return v[ind];
    }

    // Function 11: Lennard-Jones potential (12-6 form)
    // f(r) = 4*[(1/r)^12 - (1/r)^6] where r = |x|
    // Models van der Waals interaction between neutral atoms
    static Real TestScalarFunc_LennardJones(const VectorN<Real, 3> &v) 
    { 
        Real r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        if (r2 < 0.01) return 1e10;  // Regularize singularity
        Real r6_inv = 1.0 / (r2*r2*r2);
        Real r12_inv = r6_inv * r6_inv;
        return 4.0 * (r12_inv - r6_inv);
    }
    static Real TestScalarFunc_LennardJones_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        Real r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        if (r2 < 0.01) return 0;  // Regularize singularity
        Real r = sqrt(r2);
        Real r7_inv = 1.0 / pow(r, 7);
        Real r13_inv = 1.0 / pow(r, 13);
        // df/dxi = 4*(-12/r^13 + 6/r^7) * (xi/r) = 4*xi*(-12/r^14 + 6/r^8)
        return 4.0 * v[ind] * (-12.0*r13_inv/r + 6.0*r7_inv/r);
    }

    // Function 12: Morse potential
    // f(r) = D*(1 - exp(-a*(r-r0)))²  where r = |x|, D=1, a=1, r0=1
    // Models diatomic molecular bond stretching
    static Real TestScalarFunc_Morse(const VectorN<Real, 3> &v) 
    { 
        const Real D = 1.0;   // Well depth
        const Real a = 1.0;   // Width parameter
        const Real r0 = 1.0;  // Equilibrium distance
        Real r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        Real term = 1.0 - exp(-a*(r - r0));
        return D * term * term;
    }
    static Real TestScalarFunc_Morse_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        const Real D = 1.0;
        const Real a = 1.0;
        const Real r0 = 1.0;
        Real r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        if (r2 < 1e-20) return 0;
        Real r = sqrt(r2);
        Real exp_term = exp(-a*(r - r0));
        // df/dxi = 2*D*a*(1-exp(-a(r-r0)))*exp(-a(r-r0)) * (xi/r)
        return 2.0*D*a*(1.0 - exp_term)*exp_term * (v[ind]/r);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                 ADDITIONAL MATHEMATICAL FUNCTIONS (13-15)                             //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Function 13: Gaussian function (bell curve in 3D)
    // f(x,y,z) = exp(-(x²+y²+z²))
    // Smooth, bounded, single maximum at origin
    static Real TestScalarFunc_Gaussian(const VectorN<Real, 3> &v) 
    { 
        return exp(-(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
    }
    static Real TestScalarFunc_Gaussian_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        return -2.0 * v[ind] * exp(-(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
    }

    // Function 14: Sinusoidal product
    // f(x,y,z) = sin(x)*sin(y)*sin(z)
    // Separable, periodic, many critical points
    static Real TestScalarFunc_SinProduct(const VectorN<Real, 3> &v) 
    { 
        return sin(v[0]) * sin(v[1]) * sin(v[2]);
    }
    static Real TestScalarFunc_SinProduct_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        if (ind == 0) return cos(v[0]) * sin(v[1]) * sin(v[2]);
        else if (ind == 1) return sin(v[0]) * cos(v[1]) * sin(v[2]);
        else return sin(v[0]) * sin(v[1]) * cos(v[2]);
    }

    // Function 15: Polynomial with cross terms
    // f(x,y,z) = x³ + y³ + z³ - 3xyz
    // Has interesting factorization: (x+y+z)(x²+y²+z²-xy-yz-zx)
    static Real TestScalarFunc_Polynomial(const VectorN<Real, 3> &v) 
    { 
        Real x = v[0], y = v[1], z = v[2];
        return x*x*x + y*y*y + z*z*z - 3*x*y*z;
    }
    static Real TestScalarFunc_Polynomial_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        Real x = v[0], y = v[1], z = v[2];
        if (ind == 0) return 3*x*x - 3*y*z;
        else if (ind == 1) return 3*y*y - 3*x*z;
        else return 3*z*z - 3*x*y;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                             TEST BED CLASS                                            //
    ///////////////////////////////////////////////////////////////////////////////////////////

    class ScalarFunctionsTestBed
    {
    public:
        static int getNumTestFunctionScalar3() { return 15; }

        const static TestFunctionScalar<3>& getTestFunctionScalar3(int i)  { return _listFuncScalar3[i]; }

        const static TestFunctionScalar<3>& getTestFunctionScalar3(const std::string &funcName)
        {
            for (int i = 0; i < getNumTestFunctionScalar3(); i++)
            {
                if (_listFuncScalar3[i]._funcName == funcName)
                    return _listFuncScalar3[i];
            }
            throw std::runtime_error("TestFunctionScalar " + funcName + " not found!");
        }
    
    private:
        const static inline TestFunctionScalar<3> _listFuncScalar3[] = { 
            // Original functions (1-2)
            { "Scalar func 1", 
              TestScalarFunc1, "cos(x) + sin(y) + exp(z)", 
              TestScalarFunc1_derived, "[-sin(x), cos(y), exp(z)]" },
            { "Scalar func 2", 
              TestScalarFunc2, "sin(x*y)*exp(z/(y^2+1))/(1+x^2)", 
              TestScalarFunc2_derived, "complex partial derivatives" },

            // Quadratic forms (3-5)
            { "Sphere", 
              TestScalarFunc_Sphere, "x^2 + y^2 + z^2", 
              TestScalarFunc_Sphere_derived, "[2x, 2y, 2z]" },
            { "Saddle", 
              TestScalarFunc_Saddle, "x^2 - y^2 + z^2", 
              TestScalarFunc_Saddle_derived, "[2x, -2y, 2z]" },
            { "Weighted Quadratic", 
              TestScalarFunc_WeightedQuadratic, "x^2 + 4y^2 + 9z^2", 
              TestScalarFunc_WeightedQuadratic_derived, "[2x, 8y, 18z]" },

            // Classic optimization test functions (6-8)
            { "Rosenbrock 3D", 
              TestScalarFunc_Rosenbrock, "(1-x)^2 + 100(y-x^2)^2 + (1-y)^2 + 100(z-y^2)^2", 
              TestScalarFunc_Rosenbrock_derived, "gradient of Rosenbrock" },
            { "Rastrigin 3D", 
              TestScalarFunc_Rastrigin, "30 + sum_i(x_i^2 - 10*cos(2*pi*x_i))", 
              TestScalarFunc_Rastrigin_derived, "[2x + 20*pi*sin(2*pi*x), ...]" },
            { "Ackley 3D", 
              TestScalarFunc_Ackley, "-20*exp(-0.2*sqrt(sum/3)) - exp(cos_sum/3) + 20 + e", 
              TestScalarFunc_Ackley_derived, "gradient of Ackley" },

            // Physical potentials (9-12)
            { "Coulomb Potential", 
              TestScalarFunc_CoulombPotential, "-1/sqrt(x^2+y^2+z^2)", 
              TestScalarFunc_CoulombPotential_derived, "[x/r^3, y/r^3, z/r^3]" },
            { "Harmonic Potential", 
              TestScalarFunc_HarmonicPotential, "0.5*(x^2+y^2+z^2)", 
              TestScalarFunc_HarmonicPotential_derived, "[x, y, z]" },
            { "Lennard-Jones", 
              TestScalarFunc_LennardJones, "4*[(1/r)^12 - (1/r)^6]", 
              TestScalarFunc_LennardJones_derived, "gradient of LJ potential" },
            { "Morse Potential", 
              TestScalarFunc_Morse, "D*(1-exp(-a*(r-r0)))^2", 
              TestScalarFunc_Morse_derived, "gradient of Morse potential" },

            // Additional mathematical functions (13-15)
            { "Gaussian 3D", 
              TestScalarFunc_Gaussian, "exp(-(x^2+y^2+z^2))", 
              TestScalarFunc_Gaussian_derived, "[-2x*f, -2y*f, -2z*f]" },
            { "Sin Product", 
              TestScalarFunc_SinProduct, "sin(x)*sin(y)*sin(z)", 
              TestScalarFunc_SinProduct_derived, "[cos(x)*sin(y)*sin(z), ...]" },
            { "Polynomial Cross Terms", 
              TestScalarFunc_Polynomial, "x^3 + y^3 + z^3 - 3xyz", 
              TestScalarFunc_Polynomial_derived, "[3x^2-3yz, 3y^2-3xz, 3z^2-3xy]" }
        };
    };
}

#endif