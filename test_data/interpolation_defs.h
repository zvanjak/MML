#if !defined __MML_INTERPOLATION_DEFS_H
#define __MML_INTERPOLATION_DEFS_H

#include <string>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                              TRUE FUNCTIONS                                            //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // 1. Runge function: Classic example of polynomial interpolation failure
    // f(x) = 1/(1+25x²) on [-1, 1]
    // Exhibits severe Runge phenomenon with high-degree polynomial interpolation at uniform nodes
    static Real Runge_TrueFunc(Real x) 
    { 
        return static_cast<Real>(1.0) / (static_cast<Real>(1.0) + static_cast<Real>(25.0) * x * x); 
    }

    // 2. Smooth polynomial: f(x) = x³ - 2x² + x - 1
    // Should be exactly reproducible by polynomial interpolation
    static Real SmoothPoly_TrueFunc(Real x) 
    { 
        return x*x*x - 2*x*x + x - 1; 
    }

    // 3. Smooth transcendental: f(x) = exp(-x²) (Gaussian)
    // Very smooth, ideal for testing convergence rates
    static Real Gaussian_TrueFunc(Real x) 
    { 
        return std::exp(-x * x); 
    }

    // 4. Smooth sinusoidal: f(x) = sin(πx)
    // Tests periodic-like functions (one period)
    static Real SinePi_TrueFunc(Real x) 
    { 
        return std::sin(static_cast<Real>(Constants::PI) * x); 
    }

    // 5. Oscillatory: f(x) = sin(5πx)
    // Rapidly oscillating - needs many points or low order
    static Real OscillatoryHighFreq_TrueFunc(Real x) 
    { 
        return std::sin(static_cast<Real>(5.0 * Constants::PI) * x); 
    }

    // 6. Absolute value: f(x) = |x|
    // Discontinuous first derivative at x=0
    static Real AbsValue_TrueFunc(Real x) 
    { 
        return std::abs(x); 
    }

    // 7. Square root: f(x) = sqrt(x+1)
    // Discontinuous second derivative as x→-1
    static Real SqrtShift_TrueFunc(Real x) 
    { 
        return std::sqrt(x + static_cast<Real>(1.0)); 
    }

    // 8. Exponential: f(x) = exp(x)
    // Classic smooth function, should converge quickly
    static Real Exponential_TrueFunc(Real x) 
    { 
        return std::exp(x); 
    }

    // 9. Logarithm: f(x) = ln(x+2)
    // Smooth but curved, slower convergence near singularity
    static Real LogShift_TrueFunc(Real x) 
    { 
        return std::log(x + static_cast<Real>(2.0)); 
    }

    // 10. Chebyshev polynomial: T_5(x)
    // Tests if interpolation reproduces polynomials exactly
    static Real Chebyshev5_TrueFunc(Real x) 
    { 
        // T_5(x) = 16x^5 - 20x^3 + 5x
        return 16*x*x*x*x*x - 20*x*x*x + 5*x; 
    }

    // 11. Witch of Agnesi: f(x) = 1/(1+x²)
    // Similar to Runge but milder - still exhibits phenomenon
    static Real WitchOfAgnesi_TrueFunc(Real x) 
    { 
        return static_cast<Real>(1.0) / (static_cast<Real>(1.0) + x * x); 
    }

    // 12. Step-like function: f(x) = tanh(10x)
    // Very steep but continuous - challenges polynomial interpolation
    static Real StepLike_TrueFunc(Real x) 
    { 
        return std::tanh(static_cast<Real>(10.0) * x); 
    }

    // 13. Combined oscillatory: f(x) = sin(x) + 0.5*sin(3x) + 0.25*sin(5x)
    // Multi-frequency content
    static Real MultiFreq_TrueFunc(Real x) 
    { 
        return std::sin(x) + static_cast<Real>(0.5)*std::sin(static_cast<Real>(3.0)*x) 
             + static_cast<Real>(0.25)*std::sin(static_cast<Real>(5.0)*x); 
    }

    // 14. Rational function: f(x) = x/(1+x²)
    // Smooth with asymptotic behavior
    static Real RationalSmooth_TrueFunc(Real x) 
    { 
        return x / (static_cast<Real>(1.0) + x * x); 
    }

    // 15. Bessel-like: f(x) = sin(x)/x (sinc function, defined as 1 at x=0)
    static inline Real Sinc_TrueFunc(Real x)
    { 
        if (std::abs(x) < static_cast<Real>(1e-10)) return static_cast<Real>(1.0);
        return std::sin(x) / x; 
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    PRE-COMPUTED DATA SETS (STATIC)                                     //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Runge function on uniform nodes [-1, 1], n=11
    const static inline Vector<Real> runge_uniform_11_x{-MML::Real(1.0), -MML::Real(0.8), -MML::Real(0.6), -MML::Real(0.4), -MML::Real(0.2), MML::Real(0.0), MML::Real(0.2), MML::Real(0.4), MML::Real(0.6), MML::Real(0.8), MML::Real(1.0)};
    const static inline Vector<Real> runge_uniform_11_y{
        MML::Real(0.038461538461538464), MML::Real(0.058823529411764705), MML::Real(0.1), MML::Real(0.2), MML::Real(0.5), MML::Real(1.0),
        MML::Real(0.5), MML::Real(0.2), MML::Real(0.1), MML::Real(0.058823529411764705), MML::Real(0.038461538461538464)
    };

    // Smooth polynomial x³-2x²+x-1 on [-2,2], n=5
    const static inline Vector<Real> smoothpoly_5_x{-2.0, -1.0, 0.0, 1.0, 2.0};
    const static inline Vector<Real> smoothpoly_5_y{-19.0, -5.0, -1.0, -1.0, 1.0};

    // Gaussian exp(-x²) on [-3,3], n=15
    const static inline Vector<Real> gaussian_15_x{
        -MML::Real(3.0), -MML::Real(2.571428571), -MML::Real(2.142857143), -MML::Real(1.714285714), -MML::Real(1.285714286),
        -MML::Real(0.857142857), -MML::Real(0.428571429), MML::Real(0.0), MML::Real(0.428571429), MML::Real(0.857142857),
        MML::Real(1.285714286), MML::Real(1.714285714), MML::Real(2.142857143), MML::Real(2.571428571), MML::Real(3.0)
    };
    const static inline Vector<Real> gaussian_15_y{
        MML::Real(0.000123409804), MML::Real(0.001350455), MML::Real(0.010064502), MML::Real(0.051065857), MML::Real(0.192107629),
        MML::Real(0.480968586), MML::Real(0.832115086), MML::Real(1.0), MML::Real(0.832115086), MML::Real(0.480968586),
        MML::Real(0.192107629), MML::Real(0.051065857), MML::Real(0.010064502), MML::Real(0.001350455), MML::Real(0.000123409804)
    };

    // Sine(πx) on [-1,1], n=9
    const static inline Vector<Real> sinepi_9_x{-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0};
    const static inline Vector<Real> sinepi_9_y{
        MML::Real(0.0), -MML::Real(0.707106781), -MML::Real(1.0), -MML::Real(0.707106781), MML::Real(0.0), MML::Real(0.707106781), MML::Real(1.0), MML::Real(0.707106781), MML::Real(0.0)
    };

    // Exponential exp(x) on [-2,2], n=9
    const static inline Vector<Real> exp_9_x{-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
    const static inline Vector<Real> exp_9_y{
        MML::Real(0.135335283), MML::Real(0.223130160), MML::Real(0.367879441), MML::Real(0.606530660), MML::Real(1.0),
        MML::Real(1.648721271), MML::Real(2.718281828), MML::Real(4.481689070), MML::Real(7.389056099)
    };

    // Absolute value |x| on [-2,2], n=11
    const static inline Vector<Real> absval_11_x{-MML::Real(2.0), -MML::Real(1.6), -MML::Real(1.2), -MML::Real(0.8), -MML::Real(0.4), MML::Real(0.0), MML::Real(0.4), MML::Real(0.8), MML::Real(1.2), MML::Real(1.6), MML::Real(2.0)};
    const static inline Vector<Real> absval_11_y{MML::Real(2.0), MML::Real(1.6), MML::Real(1.2), MML::Real(0.8), MML::Real(0.4), MML::Real(0.0), MML::Real(0.4), MML::Real(0.8), MML::Real(1.2), MML::Real(1.6), MML::Real(2.0)};

    // Step-like tanh(10x) on [-1,1], n=21
    const static inline Vector<Real> steplike_21_x{
        -MML::Real(1.0), -MML::Real(0.9), -MML::Real(0.8), -MML::Real(0.7), -MML::Real(0.6), -MML::Real(0.5), -MML::Real(0.4), -MML::Real(0.3), -MML::Real(0.2), -MML::Real(0.1), MML::Real(0.0),
        MML::Real(0.1), MML::Real(0.2), MML::Real(0.3), MML::Real(0.4), MML::Real(0.5), MML::Real(0.6), MML::Real(0.7), MML::Real(0.8), MML::Real(0.9), MML::Real(1.0)
    };
    const static inline Vector<Real> steplike_21_y{
        -MML::Real(0.9999999959), -MML::Real(0.9999999588), -MML::Real(0.9999996829), -MML::Real(0.9999983370), -MML::Real(0.9999877117),
        -MML::Real(0.9999092043), -MML::Real(0.9993292997), -MML::Real(0.9950547537), -MML::Real(0.9640275801), -MML::Real(0.7615941560), MML::Real(0.0),
        MML::Real(0.7615941560), MML::Real(0.9640275801), MML::Real(0.9950547537), MML::Real(0.9993292997), MML::Real(0.9999092043),
        MML::Real(0.9999877117), MML::Real(0.9999983370), MML::Real(0.9999996829), MML::Real(0.9999999588), MML::Real(0.9999999959)
    };

    // Chebyshev T_5(x) on Chebyshev nodes [-1,1], n=7
    const static inline Vector<Real> cheb5_chebNodes_7_x{
        MML::Real(0.9749279122), MML::Real(0.7818314825), MML::Real(0.4338837391), MML::Real(0.0), -MML::Real(0.4338837391), -MML::Real(0.7818314825), -MML::Real(0.9749279122)
    };
    const static inline Vector<Real> cheb5_chebNodes_7_y{
        MML::Real(0.7818314825), -MML::Real(0.9749279122), -MML::Real(0.4338837391), MML::Real(0.0), MML::Real(0.4338837391), MML::Real(0.9749279122), -MML::Real(0.7818314825)
    };

    // High-frequency sin(5πx) on [-1,1], n=21
    const static inline Vector<Real> highfreq_21_x{
        -MML::Real(1.0), -MML::Real(0.9), -MML::Real(0.8), -MML::Real(0.7), -MML::Real(0.6), -MML::Real(0.5), -MML::Real(0.4), -MML::Real(0.3), -MML::Real(0.2), -MML::Real(0.1), MML::Real(0.0),
        MML::Real(0.1), MML::Real(0.2), MML::Real(0.3), MML::Real(0.4), MML::Real(0.5), MML::Real(0.6), MML::Real(0.7), MML::Real(0.8), MML::Real(0.9), MML::Real(1.0)
    };
    const static inline Vector<Real> highfreq_21_y{
        MML::Real(0.0), -MML::Real(0.951056516), MML::Real(0.587785252), MML::Real(0.587785252), -MML::Real(0.951056516), MML::Real(0.0),
        MML::Real(0.951056516), -MML::Real(0.587785252), -MML::Real(0.587785252), MML::Real(0.951056516), MML::Real(0.0),
        -MML::Real(0.951056516), MML::Real(0.587785252), MML::Real(0.587785252), -MML::Real(0.951056516), MML::Real(0.0),
        MML::Real(0.951056516), -MML::Real(0.587785252), -MML::Real(0.587785252), MML::Real(0.951056516), MML::Real(0.0)
    };

    // Witch of Agnesi on [-5,5], n=15
    const static inline Vector<Real> witch_15_x{
        -MML::Real(5.0), -MML::Real(4.285714286), -MML::Real(3.571428571), -MML::Real(2.857142857), -MML::Real(2.142857143),
        -MML::Real(1.428571429), -MML::Real(0.714285714), MML::Real(0.0), MML::Real(0.714285714), MML::Real(1.428571429),
        MML::Real(2.142857143), MML::Real(2.857142857), MML::Real(3.571428571), MML::Real(4.285714286), MML::Real(5.0)
    };
    const static inline Vector<Real> witch_15_y{
        MML::Real(0.038461538), MML::Real(0.051652893), MML::Real(0.072727273), MML::Real(0.109090909), MML::Real(0.178947368),
        MML::Real(0.329411765), MML::Real(0.662162162), MML::Real(1.0), MML::Real(0.662162162), MML::Real(0.329411765),
        MML::Real(0.178947368), MML::Real(0.109090909), MML::Real(0.072727273), MML::Real(0.051652893), MML::Real(0.038461538)
    };

    // Logarithm ln(x+2) on [-1,3], n=15
    const static inline Vector<Real> log_15_x{
        -MML::Real(1.0), -MML::Real(0.714285714), -MML::Real(0.428571429), -MML::Real(0.142857143), MML::Real(0.142857143),
        MML::Real(0.428571429), MML::Real(0.714285714), MML::Real(1.0), MML::Real(1.285714286), MML::Real(1.571428571),
        MML::Real(1.857142857), MML::Real(2.142857143), MML::Real(2.428571429), MML::Real(2.714285714), MML::Real(3.0)
    };
    const static inline Vector<Real> log_15_y{
        MML::Real(0.0), MML::Real(0.251314428), MML::Real(0.451076387), MML::Real(0.617000854), MML::Real(0.758317636),
        MML::Real(0.880534220), MML::Real(0.987565958), MML::Real(1.098612289), MML::Real(1.190619774), MML::Real(1.274397839),
        MML::Real(1.351355906), MML::Real(1.422467814), MML::Real(1.488489739), MML::Real(1.550015831), MML::Real(1.609437912)
    };

    // Square root sqrt(x+1) on [0,4], n=11
    const static inline Vector<Real> sqrt_11_x{MML::Real(0.0), MML::Real(0.4), MML::Real(0.8), MML::Real(1.2), MML::Real(1.6), MML::Real(2.0), MML::Real(2.4), MML::Real(2.8), MML::Real(3.2), MML::Real(3.6), MML::Real(4.0)};
    const static inline Vector<Real> sqrt_11_y{
        MML::Real(1.0), MML::Real(1.183215957), MML::Real(1.341640786), MML::Real(1.483239697), MML::Real(1.612451550),
        MML::Real(1.732050808), MML::Real(1.843908891), MML::Real(1.949358869), MML::Real(2.049390153), MML::Real(2.144761059), MML::Real(2.236067977)
    };

    // Multi-frequency on [-π, π], n=25
    const static inline Vector<Real> multifreq_25_x{
        -MML::Real(3.141592654), -MML::Real(2.879793266), -MML::Real(2.617993878), -MML::Real(2.356194490), -MML::Real(2.094395102),
        -MML::Real(1.832595714), -MML::Real(1.570796327), -MML::Real(1.308996939), -MML::Real(1.047197551), -MML::Real(0.785398163),
        -MML::Real(0.523598776), -MML::Real(0.261799388), MML::Real(0.0), MML::Real(0.261799388), MML::Real(0.523598776),
        MML::Real(0.785398163), MML::Real(1.047197551), MML::Real(1.308996939), MML::Real(1.570796327), MML::Real(1.832595714),
        MML::Real(2.094395102), MML::Real(2.356194490), MML::Real(2.617993878), MML::Real(2.879793266), MML::Real(3.141592654)
    };
    const static inline Vector<Real> multifreq_25_y{
        MML::Real(0.0), -MML::Real(0.465011316), -MML::Real(0.771329079), -MML::Real(0.866025404), -MML::Real(0.767175251),
        -MML::Real(0.537109750), -MML::Real(0.250000000), MML::Real(0.022096262), MML::Real(0.224014303), MML::Real(0.340536025),
        MML::Real(0.387595562), MML::Real(0.395240461), MML::Real(0.0), -MML::Real(0.395240461), -MML::Real(0.387595562),
        -MML::Real(0.340536025), -MML::Real(0.224014303), -MML::Real(0.022096262), MML::Real(0.250000000), MML::Real(0.537109750),
        MML::Real(0.767175251), MML::Real(0.866025404), MML::Real(0.771329079), MML::Real(0.465011316), MML::Real(0.0)
    };

} // namespace MML::TestBeds

#endif // __MML_INTERPOLATION_DEFS_H
