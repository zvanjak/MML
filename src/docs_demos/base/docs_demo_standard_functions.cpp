///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_standard_functions.cpp                                    ///
///  Description: Brief demonstration of StandardFunctions.h - special math functions ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/StandardFunctions.h>
#include <mml/base/SpecialFunctions.h>
#include <mml/base/DiracDeltaFunction.h>
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Docs_Demo_ForwardShapes()
{
    const Real gammaP = SpecialFunctions::RegularizedGammaP(Real(1e6), Real(1e6));
    const Real betaP = SpecialFunctions::RegularizedBetaI(Real(0.5), Real(1e12), Real(1e12));
    const Real gammaMedian = SpecialFunctions::RegularizedGammaPInv(Real(1e6), Real(0.5));
    const Real finiteBeta = SpecialFunctions::IncompleteBeta(
        Real(0.5), 2, std::numeric_limits<Real>::denorm_min());
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout << gammaP << " " << betaP << " " << gammaMedian << " " << finiteBeta << '\n';
    if (!IsWithinRel(gammaP, Real(0.50013298076087259124432281750320671L), Precision::EqualityTolerance) ||
        !IsWithinRel(betaP, Real(0.5), Precision::EqualityTolerance) ||
        !(gammaMedian < Real(1e6) && gammaMedian > Real(1e6) - 1) ||
        !IsWithinRel(finiteBeta, std::log(Real(2)) - Real(0.5), Precision::EqualityTolerance))
        throw NumericalMethodError("Forward shape documentation example failed");
}

void Docs_Demo_BetaInverseTails()
{
    const Real lower = SpecialFunctions::RegularizedBetaIInv(Real(1e-100), 2, 1);
    const Real uniform = SpecialFunctions::RegularizedBetaIInv(Real(1e-300), 1, 1);
    const Real tiny = SpecialFunctions::RegularizedBetaIInv(Real(1e-200), Real(0.5), Real(0.5));
    const Real rounded = SpecialFunctions::RegularizedBetaIInv(Real(0.99), 2, Real(0.1));
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout << lower << '\n' << uniform << '\n' << tiny << '\n' << rounded << '\n';

    bool rejected = false;
    try {
        SpecialFunctions::RegularizedBetaIInv(Real(-0.1), 2, 3);
    }
    catch (const DomainError&) {
        rejected = true;
        std::cout << "Out-of-range probability rejected\n";
    }
    if (!IsWithinRel(lower, Real(1e-50), Precision::EqualityTolerance) ||
        uniform != Real(1e-300) || tiny != 0 || std::signbit(tiny) || rounded != 1 || !rejected)
        throw NumericalMethodError("Beta inverse documentation example failed");
}

void Docs_Demo_GammaInverseTails()
{
    const Real lower = SpecialFunctions::RegularizedGammaPInv(2, Real(1e-100));
    const Real tiny = SpecialFunctions::RegularizedGammaPInv(Real(0.5), Real(1e-200));
    const Real p = std::nextafter(Real(1), Real(0));
    const Real upper = SpecialFunctions::RegularizedGammaPInv(1, p);
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout << lower << '\n' << tiny << '\n' << upper << '\n';

    bool rejected = false;
    try {
        SpecialFunctions::RegularizedGammaPInv(2, Real(-0.1));
    }
    catch (const DomainError&) {
        rejected = true;
        std::cout << "Out-of-range probability rejected\n";
    }
    if (!IsNearlyEqual(lower / Real(1.414213562373095048801688724209698079e-50L),
                       Real(1), Precision::EqualityTolerance) ||
        tiny != 0 || std::signbit(tiny) || upper != -std::log1p(-p) || !rejected)
        throw NumericalMethodError("Gamma inverse documentation example failed");
}

void Docs_Demo_SineCosineEndpoints()
{
    const Real infinity = std::numeric_limits<Real>::infinity();
    Real si, ci;
    SpecialFunctions::SinCosIntegral(-Real(0), si, ci);
    std::cout << std::signbit(si) << ' ' << ci << '\n';  // 1 -inf
    const Real sineLimit = SpecialFunctions::SineIntegral(-infinity);
    const Real cosineLimit = SpecialFunctions::CosineIntegral(infinity);
    std::cout << sineLimit << ' ' << cosineLimit << '\n';  // -pi/2 0

    if (si != 0 || !std::signbit(si) || ci != -infinity ||
        sineLimit != -Constants::PI / 2 || cosineLimit != 0 || std::signbit(cosineLimit))
        throw NumericalMethodError("Sine cosine documentation example failed its endpoint limits");

    bool rejected = false;
    try {
        SpecialFunctions::SinCosIntegral(std::numeric_limits<Real>::quiet_NaN(), si, ci);
    }
    catch (const DomainError&) {
        rejected = true;
        std::cout << "NaN rejected before modifying either output\n";
    }
    if (!rejected || si != 0 || !std::signbit(si) || ci != -infinity)
        throw NumericalMethodError("Sine cosine documentation example failed its NaN contract");
}

void Docs_Demo_ComplexReciprocals()
{
    const Complex z(1, 1);
    const Complex values[] = {
        Functions::Sec(z), Functions::Csc(z), Functions::Ctg(z),
        Functions::Sech(z), Functions::Csch(z), Functions::Ctgh(z)
    };
    for (const Complex value : values)
        std::cout << value << '\n';

    bool rejected = false;
    try {
        Functions::Csc(Complex(0));
    }
    catch (const DomainError&) {
        rejected = true;
        std::cout << "Complex cosecant rejected its zero denominator\n";
    }

    const Complex denominators[] = {
        std::cos(z), std::sin(z), std::tan(z),
        std::cosh(z), std::sinh(z), std::tanh(z)
    };
    for (int i = 0; i < 6; ++i) {
        if (!IsNearlyZero(std::abs(values[i] * denominators[i] - Complex(1)), Precision::RoundoffZeroThreshold))
            throw NumericalMethodError("Complex reciprocal documentation example failed its identity");
    }
    if (!rejected)
        throw NumericalMethodError("Complex cosecant documentation example did not reject its pole");
}

void Docs_Demo_IntegerBoundaries()
{
    const long long factorial = Functions::FactorialInt(20);
    std::cout << factorial << '\n';  // 2432902008176640000
    bool rejected = false;
    try {
        Functions::FactorialInt(21);
    }
    catch (const DomainError&) {
        rejected = true;
        std::cout << "21! exceeds the exact signed 64-bit range\n";
    }
    const int n = std::numeric_limits<int>::max();
    const DiracStep step(n);
    const DiracSqr cauchy(n);
    const Real x = Real(1) / static_cast<Real>(n);
    std::cout << step(0) << '\n';
    std::cout << cauchy(x) << '\n';  // n / (2*pi)

    if (factorial != 2432902008176640000LL || !rejected ||
        step(0) != static_cast<Real>(n) || step(x) != 0 ||
        !IsNearlyEqual(cauchy(x) / (static_cast<Real>(n) / (2 * Constants::PI)),
                       Real(1), Precision::RoundoffZeroThreshold))
        throw NumericalMethodError("Integer boundary documentation example failed");

    const Real ordinary = Functions::Factorial(20);
    const Real huge = Functions::Factorial(std::numeric_limits<int>::max());
    std::cout << ordinary << '\n';
    std::cout << huge << '\n';  // positive infinity
    if (!IsNearlyEqual(ordinary / static_cast<Real>(factorial),
                       Real(1), Precision::RoundoffZeroThreshold) ||
        huge != Constants::PosInf)
        throw NumericalMethodError("Floating factorial documentation example failed");
}

void Docs_Demo_IncompleteGamma()
{
    using namespace MML;
    const Real lower = SpecialFunctions::IncompleteGammaLower(200, 1);
    const Real upper = SpecialFunctions::IncompleteGammaUpper(200, 2000);
    std::cout << lower << '\n';  // approximately 0.001848593963122709
    std::cout << upper << '\n';  // approximately 2.298767256466363e-212

    if (!std::isfinite(lower) || !std::isfinite(upper) || lower <= 0 || upper <= 0 ||
        !IsNearlyEqual(lower, Real(0.0018485939631227090831), Precision::EqualityTolerance) ||
        !IsNearlyEqual(upper / Real(2.2987672564663627783e-212), Real(1), Precision::EqualityTolerance))
        throw NumericalMethodError("Incomplete gamma documentation example did not retain its finite tails");
}

void Docs_Demo_LambertW()
{
    using namespace MML;
    const Real endpoint = Real(-0.367879441171442321595523770161460867L);
    const Real x = std::nextafter(endpoint, Real(0));
    const Real principal = SpecialFunctions::LambertW0(x);
    const Real secondary = SpecialFunctions::LambertWm1(x);
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout << principal << '\n' << secondary << '\n';

    if (!std::isfinite(principal) || !std::isfinite(secondary) || principal <= -1 || secondary >= -1 ||
        SpecialFunctions::LambertW0(endpoint) != -1 || SpecialFunctions::LambertWm1(endpoint) != -1)
        throw NumericalMethodError("Lambert W documentation example did not retain its distinct real branches");
    if constexpr (std::numeric_limits<Real>::digits == 53) {
        if (!IsWithinAbs(principal, Real(-0.9999999846957458715), Precision::RoundoffZeroThreshold) ||
            !IsWithinAbs(secondary, Real(-1.0000000153042542846), Precision::RoundoffZeroThreshold))
            throw NumericalMethodError("Lambert W documentation example disagrees with independent reference values");
    }
}

void Docs_Demo_SpecialFunctionRanges()
{
    using namespace MML;
    const Real tail = SpecialFunctions::Dawson(std::numeric_limits<Real>::max());
    const Real ei = SpecialFunctions::ExpIntegralEi(710);
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout << tail << '\n' << ei << '\n';

    if (!std::isfinite(tail) || tail <= 0 ||
        !IsNearlyEqual((tail * std::numeric_limits<Real>::max()) * 2, Real(1), Precision::EqualityTolerance))
        throw NumericalMethodError("Dawson documentation example lost its finite reciprocal tail");
    const Real expectedEi = Real(3.150915688206201214938641192e305L);
    if (std::isfinite(expectedEi)) {
        if (!std::isfinite(ei) || !IsNearlyEqual(ei / expectedEi, Real(1), Precision::EqualityTolerance))
            throw NumericalMethodError("Ei documentation example lost a finite representable result");
    }
    else if (ei != std::numeric_limits<Real>::infinity())
        throw NumericalMethodError("Ei documentation example did not report genuine range overflow");
}

void Docs_Demo_StandardFunctions()
{
    std::cout << "=== StandardFunctions Demo ===" << std::endl;
    
    Real x = 0.5;
    
    // Basic trig functions
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Sin(" << x << ") = " << Functions::Sin(x) << std::endl;
    std::cout << "Cos(" << x << ") = " << Functions::Cos(x) << std::endl;
    std::cout << "Tan(" << x << ") = " << Functions::Tan(x) << std::endl;
    
    // Hyperbolic functions
    std::cout << "Sinh(" << x << ") = " << Functions::Sinh(x) << std::endl;
    std::cout << "Cosh(" << x << ") = " << Functions::Cosh(x) << std::endl;
    
    // Exponential and log
    std::cout << "Exp(" << x << ") = " << Functions::Exp(x) << std::endl;
    std::cout << "Log(" << x << ") = " << Functions::Log(x) << std::endl;
    
    // Special functions
    std::cout << "Erf(" << x << ") = " << Functions::Erf(x) << std::endl;
    std::cout << "TGamma(5) = " << Functions::TGamma(5.0) << " (= 4! = 24)" << std::endl;
    
#if MML_HAS_STD_SPECIAL_FUNCTIONS
    std::cout << "Bessel J0(1) = " << Functions::CylBesselJ(0, 1.0) << std::endl;
    std::cout << "Legendre P2(0.5) = " << Functions::Legendre(2, 0.5) << std::endl;
#else
    std::cout << "(Special functions not available on this platform)" << std::endl;
#endif
    
    Docs_Demo_IncompleteGamma();
    Docs_Demo_ForwardShapes();
    Docs_Demo_LambertW();
    Docs_Demo_SpecialFunctionRanges();
    std::cout << "=== Demo Complete ===" << std::endl;
}
