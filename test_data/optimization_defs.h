#if !defined __MML_OPTIMIZATION_DEFS_H
#define __MML_OPTIMIZATION_DEFS_H

#include <string>
#include <cmath>
#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                      1D OPTIMIZATION TEST FUNCTIONS                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 1: UNIMODAL SMOOTH FUNCTIONS (easy, single minimum)
    // ========================================================================================

    // 1. Simple quadratic: f(x) = (x-2)² + 1
    // Minimum at x=2, f(2)=1
    static Real Opt1D_Quadratic(Real x) { return (x - 2) * (x - 2) + 1; }
    static Real Opt1D_Quadratic_deriv(Real x) { return 2 * (x - 2); }
    const static inline Real Opt1D_Quadratic_min_x = 2.0;
    const static inline Real Opt1D_Quadratic_min_f = 1.0;

    // 2. Quartic: f(x) = (x-1)⁴
    // Minimum at x=1, f(1)=0 (very flat near minimum)
    static Real Opt1D_Quartic(Real x) { 
        Real d = x - 1;
        return d * d * d * d; 
    }
    static Real Opt1D_Quartic_deriv(Real x) { 
        Real d = x - 1;
        return 4 * d * d * d; 
    }
    const static inline Real Opt1D_Quartic_min_x = 1.0;
    const static inline Real Opt1D_Quartic_min_f = 0.0;

    // 3. Shifted exponential: f(x) = exp(x-3) + exp(-(x-3)) (cosh shifted)
    // Minimum at x=3, f(3)=2
    static Real Opt1D_Cosh(Real x) { 
        return std::exp(x - 3) + std::exp(-(x - 3)); 
    }
    static Real Opt1D_Cosh_deriv(Real x) { 
        return std::exp(x - 3) - std::exp(-(x - 3)); 
    }
    const static inline Real Opt1D_Cosh_min_x = 3.0;
    const static inline Real Opt1D_Cosh_min_f = 2.0;

    // 4. Absolute value (non-smooth): f(x) = |x - 1.5| + 2
    // Minimum at x=1.5, f(1.5)=2
    static Real Opt1D_AbsValue(Real x) { return std::abs(x - 1.5) + 2; }
    const static inline Real Opt1D_AbsValue_min_x = 1.5;
    const static inline Real Opt1D_AbsValue_min_f = 2.0;

    // 5. Fourth-degree polynomial: f(x) = x⁴ - 4x³ + 4x² + 1
    // = (x² - 2x)² + 1, minimum at x=0 and x=2, f=1
    static Real Opt1D_PolyDouble(Real x) { 
        Real t = x * x - 2 * x;
        return t * t + 1; 
    }
    static Real Opt1D_PolyDouble_deriv(Real x) { 
        return 4 * (x * x - 2 * x) * (2 * x - 2); 
    }
    const static inline Real Opt1D_PolyDouble_min_x1 = 0.0;
    const static inline Real Opt1D_PolyDouble_min_x2 = 2.0;
    const static inline Real Opt1D_PolyDouble_min_f = 1.0;

    // ========================================================================================
    // CATEGORY 2: MULTIMODAL FUNCTIONS (multiple local minima)
    // ========================================================================================

    // 6. Sinusoidal: f(x) = sin(x) + sin(3x)/3 + 0.5
    // Multiple local minima, global minimum near x ≈ -1.9
    static Real Opt1D_Sinusoidal(Real x) { 
        return std::sin(x) + std::sin(3 * x) / 3 + 0.5; 
    }
    static Real Opt1D_Sinusoidal_deriv(Real x) { 
        return std::cos(x) + std::cos(3 * x); 
    }

    // 7. Rastrigin 1D: f(x) = 10 + x² - 10*cos(2πx)
    // Global minimum at x=0, f(0)=0
    static Real Opt1D_Rastrigin(Real x) { 
        return 10 + x * x - 10 * std::cos(2 * Constants::PI * x); 
    }
    static Real Opt1D_Rastrigin_deriv(Real x) { 
        return 2 * x + 20 * Constants::PI * std::sin(2 * Constants::PI * x); 
    }
    const static inline Real Opt1D_Rastrigin_min_x = 0.0;
    const static inline Real Opt1D_Rastrigin_min_f = 0.0;

    // 8. Ackley 1D: f(x) = -20*exp(-0.2*|x|) - exp(cos(2πx)) + 20 + e
    // Global minimum at x=0, f(0)=0
    static Real Opt1D_Ackley(Real x) { 
        return -20 * std::exp(-0.2 * std::abs(x)) 
               - std::exp(std::cos(2 * Constants::PI * x)) 
               + 20 + Constants::E; 
    }
    const static inline Real Opt1D_Ackley_min_x = 0.0;
    const static inline Real Opt1D_Ackley_min_f = 0.0;

    // 9. Gramacy & Lee (2012): f(x) = sin(10πx)/(2x) + (x-1)⁴
    // Tricky function with local minima, on [0.5, 2.5]
    static Real Opt1D_GramacyLee(Real x) { 
        return std::sin(10 * Constants::PI * x) / (2 * x) + (x - 1) * (x - 1) * (x - 1) * (x - 1); 
    }

    // ========================================================================================
    // CATEGORY 3: CHALLENGING CASES
    // ========================================================================================

    // 10. Very flat: f(x) = exp(-1/(x²+0.01))
    // Nearly flat minimum region around x=0
    static Real Opt1D_VeryFlat(Real x) { 
        return std::exp(-1.0 / (x * x + 0.01)); 
    }
    const static inline Real Opt1D_VeryFlat_min_x = 0.0;

    // 11. Steep valley: f(x) = (x-2)² + 0.1/((x-2)²+0.001)
    // Sharp minimum at x=2
    static Real Opt1D_SteepValley(Real x) { 
        Real d = (x - 2) * (x - 2);
        return d + 0.1 / (d + 0.001); 
    }

    // 12. Log barrier: f(x) = x² - log(x), x > 0
    // Minimum at x = 1/√2 ≈ 0.707, f ≈ 0.847
    static Real Opt1D_LogBarrier(Real x) { 
        return x * x - std::log(x); 
    }
    static Real Opt1D_LogBarrier_deriv(Real x) { 
        return 2 * x - 1.0 / x; 
    }
    const static inline Real Opt1D_LogBarrier_min_x = 0.7071067811865476;  // 1/sqrt(2)
    const static inline Real Opt1D_LogBarrier_min_f = 0.8465735902799727;  // 0.5 + 0.5*ln(2)

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    N-DIMENSIONAL OPTIMIZATION TEST FUNCTIONS                           //
    ///////////////////////////////////////////////////////////////////////////////////////////
    // All classic benchmark functions for optimization algorithm testing
    // Convention: Functions take VectorN<Real,N> and return Real
    // Minima documented with high precision

    // ========================================================================================
    // CATEGORY 1: CONVEX/UNIMODAL (single global minimum)
    // ========================================================================================

    // 13. SPHERE: f(x) = Σ xᵢ²
    // Global minimum: f(0,...,0) = 0
    // Simplest test function, convex, separable
    template<int N>
    static Real OptND_Sphere(const VectorN<Real, N>& x) {
        Real sum = 0;
        for (int i = 0; i < N; ++i)
            sum += x[i] * x[i];
        return sum;
    }
    template<int N>
    static VectorN<Real, N> OptND_Sphere_grad(const VectorN<Real, N>& x) {
        VectorN<Real, N> g;
        for (int i = 0; i < N; ++i)
            g[i] = 2 * x[i];
        return g;
    }

    // 14. SUM OF SQUARES: f(x) = Σ i*xᵢ²
    // Global minimum: f(0,...,0) = 0
    // Weighted version of sphere
    template<int N>
    static Real OptND_SumSquares(const VectorN<Real, N>& x) {
        Real sum = 0;
        for (int i = 0; i < N; ++i)
            sum += (i + 1) * x[i] * x[i];
        return sum;
    }

    // 15. ROTATED HYPER-ELLIPSOID: f(x) = Σᵢ(Σⱼ≤ᵢ xⱼ)²
    // Global minimum: f(0,...,0) = 0
    // Non-separable
    template<int N>
    static Real OptND_RotatedHyperEllipsoid(const VectorN<Real, N>& x) {
        Real sum = 0;
        for (int i = 0; i < N; ++i) {
            Real inner = 0;
            for (int j = 0; j <= i; ++j)
                inner += x[j];
            sum += inner * inner;
        }
        return sum;
    }

    // 16. TRID: f(x) = Σ(xᵢ-1)² - Σxᵢxᵢ₋₁
    // Global minimum depends on dimension
    // For N=2: min at (2,2), f=-2
    // For N=6: min at (6,10,12,12,10,6), f=-50
    template<int N>
    static Real OptND_Trid(const VectorN<Real, N>& x) {
        Real sum1 = 0, sum2 = 0;
        for (int i = 0; i < N; ++i)
            sum1 += (x[i] - 1) * (x[i] - 1);
        for (int i = 1; i < N; ++i)
            sum2 += x[i] * x[i-1];
        return sum1 - sum2;
    }

    // ========================================================================================
    // CATEGORY 2: ROSENBROCK VARIANTS (curved valleys)
    // ========================================================================================

    // 17. ROSENBROCK (Banana function): f(x,y) = 100(y-x²)² + (1-x)²
    // Global minimum: f(1,1) = 0
    // Classic benchmark - long curved valley
    static Real OptND_Rosenbrock2D(const VectorN<Real, 2>& x) {
        Real t1 = x[1] - x[0] * x[0];
        Real t2 = 1 - x[0];
        return 100 * t1 * t1 + t2 * t2;
    }
    static VectorN<Real, 2> OptND_Rosenbrock2D_grad(const VectorN<Real, 2>& x) {
        VectorN<Real, 2> g;
        g[0] = -400 * x[0] * (x[1] - x[0] * x[0]) - 2 * (1 - x[0]);
        g[1] = 200 * (x[1] - x[0] * x[0]);
        return g;
    }
    const static inline VectorN<Real, 2> OptND_Rosenbrock2D_min = {1.0, 1.0};
    const static inline Real OptND_Rosenbrock2D_min_f = 0.0;

    // 18. N-dimensional Rosenbrock: f(x) = Σ[100(xᵢ₊₁-xᵢ²)² + (1-xᵢ)²]
    template<int N>
    static Real OptND_RosenbrockND(const VectorN<Real, N>& x) {
        Real sum = 0;
        for (int i = 0; i < N - 1; ++i) {
            Real t1 = x[i+1] - x[i] * x[i];
            Real t2 = 1 - x[i];
            sum += 100 * t1 * t1 + t2 * t2;
        }
        return sum;
    }
    template<int N>
    static VectorN<Real, N> OptND_RosenbrockND_grad(const VectorN<Real, N>& x) {
        VectorN<Real, N> g;
        for (int i = 0; i < N; ++i) g[i] = 0;
        
        for (int i = 0; i < N - 1; ++i) {
            g[i] += -400 * x[i] * (x[i+1] - x[i] * x[i]) - 2 * (1 - x[i]);
            g[i+1] += 200 * (x[i+1] - x[i] * x[i]);
        }
        return g;
    }

    // ========================================================================================
    // CATEGORY 3: MULTIMODAL FUNCTIONS (many local minima)
    // ========================================================================================

    // 19. RASTRIGIN: f(x) = 10n + Σ[xᵢ² - 10cos(2πxᵢ)]
    // Global minimum: f(0,...,0) = 0
    // Highly multimodal - standard global optimization benchmark
    template<int N>
    static Real OptND_Rastrigin(const VectorN<Real, N>& x) {
        Real sum = 10.0 * N;
        for (int i = 0; i < N; ++i)
            sum += x[i] * x[i] - 10 * std::cos(2 * Constants::PI * x[i]);
        return sum;
    }

    // 20. ACKLEY: f(x) = -20exp(-0.2√(Σxᵢ²/n)) - exp(Σcos(2πxᵢ)/n) + 20 + e
    // Global minimum: f(0,...,0) = 0
    // Multimodal with a global minimum in a large flat region
    template<int N>
    static Real OptND_Ackley(const VectorN<Real, N>& x) {
        Real sum_sq = 0, sum_cos = 0;
        for (int i = 0; i < N; ++i) {
            sum_sq += x[i] * x[i];
            sum_cos += std::cos(2 * Constants::PI * x[i]);
        }
        return -20 * std::exp(-0.2 * std::sqrt(sum_sq / N)) 
               - std::exp(sum_cos / N) 
               + 20 + Constants::E;
    }
    const static inline Real OptND_Ackley_min_f = 0.0;

    // 21. GRIEWANK: f(x) = 1 + Σxᵢ²/4000 - Πcos(xᵢ/√i)
    // Global minimum: f(0,...,0) = 0
    // Product term creates many local minima
    template<int N>
    static Real OptND_Griewank(const VectorN<Real, N>& x) {
        Real sum = 0, prod = 1;
        for (int i = 0; i < N; ++i) {
            sum += x[i] * x[i];
            prod *= std::cos(x[i] / std::sqrt(static_cast<Real>(i + 1)));
        }
        return 1 + sum / 4000 - prod;
    }

    // 22. SCHWEFEL: f(x) = 418.9829n - Σxᵢsin(√|xᵢ|)
    // Global minimum: f(420.9687,...,420.9687) ≈ 0
    // Deceptive - global minimum far from local minima
    template<int N>
    static Real OptND_Schwefel(const VectorN<Real, N>& x) {
        Real sum = 0;
        for (int i = 0; i < N; ++i)
            sum += x[i] * std::sin(std::sqrt(std::abs(x[i])));
        return 418.9829 * N - sum;
    }
    const static inline Real OptND_Schwefel_min_coord = 420.9687;

    // 23. LEVY: Complex multimodal function
    // Global minimum: f(1,...,1) = 0
    template<int N>
    static Real OptND_Levy(const VectorN<Real, N>& x) {
        auto w = [](Real xi) { return 1 + (xi - 1) / 4; };
        
        Real sum = std::pow(std::sin(Constants::PI * w(x[0])), 2);
        
        for (int i = 0; i < N - 1; ++i) {
            Real wi = w(x[i]);
            sum += (wi - 1) * (wi - 1) * (1 + 10 * std::pow(std::sin(Constants::PI * wi + 1), 2));
        }
        
        Real wn = w(x[N-1]);
        sum += (wn - 1) * (wn - 1) * (1 + std::pow(std::sin(2 * Constants::PI * wn), 2));
        
        return sum;
    }

    // ========================================================================================
    // CATEGORY 4: 2D BENCHMARK FUNCTIONS (for visualization)
    // ========================================================================================

    // 24. BEALE: f(x,y) = (1.5-x+xy)² + (2.25-x+xy²)² + (2.625-x+xy³)²
    // Global minimum: f(3, 0.5) = 0
    static Real OptND_Beale(const VectorN<Real, 2>& x) {
        Real t1 = 1.5 - x[0] + x[0] * x[1];
        Real t2 = 2.25 - x[0] + x[0] * x[1] * x[1];
        Real t3 = 2.625 - x[0] + x[0] * x[1] * x[1] * x[1];
        return t1 * t1 + t2 * t2 + t3 * t3;
    }
    const static inline VectorN<Real, 2> OptND_Beale_min = {3.0, 0.5};
    const static inline Real OptND_Beale_min_f = 0.0;

    // 25. BOOTH: f(x,y) = (x+2y-7)² + (2x+y-5)²
    // Global minimum: f(1, 3) = 0
    static Real OptND_Booth(const VectorN<Real, 2>& x) {
        Real t1 = x[0] + 2 * x[1] - 7;
        Real t2 = 2 * x[0] + x[1] - 5;
        return t1 * t1 + t2 * t2;
    }
    const static inline VectorN<Real, 2> OptND_Booth_min = {1.0, 3.0};
    const static inline Real OptND_Booth_min_f = 0.0;

    // 26. MATYAS: f(x,y) = 0.26(x²+y²) - 0.48xy
    // Global minimum: f(0, 0) = 0
    static Real OptND_Matyas(const VectorN<Real, 2>& x) {
        return 0.26 * (x[0] * x[0] + x[1] * x[1]) - 0.48 * x[0] * x[1];
    }
    const static inline VectorN<Real, 2> OptND_Matyas_min = {0.0, 0.0};
    const static inline Real OptND_Matyas_min_f = 0.0;

    // 27. HIMMELBLAU: f(x,y) = (x²+y-11)² + (x+y²-7)²
    // Four identical global minima: f ≈ 0 at
    // (3, 2), (-2.805, 3.131), (-3.779, -3.283), (3.584, -1.848)
    static Real OptND_Himmelblau(const VectorN<Real, 2>& x) {
        Real t1 = x[0] * x[0] + x[1] - 11;
        Real t2 = x[0] + x[1] * x[1] - 7;
        return t1 * t1 + t2 * t2;
    }
    const static inline VectorN<Real, 2> OptND_Himmelblau_min1 = {3.0, 2.0};
    const static inline VectorN<Real, 2> OptND_Himmelblau_min2 = {-2.805118, 3.131312};
    const static inline VectorN<Real, 2> OptND_Himmelblau_min3 = {-3.779310, -3.283186};
    const static inline VectorN<Real, 2> OptND_Himmelblau_min4 = {3.584428, -1.848126};
    const static inline Real OptND_Himmelblau_min_f = 0.0;

    // 28. GOLDSTEIN-PRICE: 
    // f(x,y) = [1+(x+y+1)²(19-14x+3x²-14y+6xy+3y²)] × 
    //          [30+(2x-3y)²(18-32x+12x²+48y-36xy+27y²)]
    // Global minimum: f(0, -1) = 3
    static Real OptND_GoldsteinPrice(const VectorN<Real, 2>& x) {
        Real a = 1 + std::pow(x[0] + x[1] + 1, 2) * 
                 (19 - 14*x[0] + 3*x[0]*x[0] - 14*x[1] + 6*x[0]*x[1] + 3*x[1]*x[1]);
        Real b = 30 + std::pow(2*x[0] - 3*x[1], 2) *
                 (18 - 32*x[0] + 12*x[0]*x[0] + 48*x[1] - 36*x[0]*x[1] + 27*x[1]*x[1]);
        return a * b;
    }
    const static inline VectorN<Real, 2> OptND_GoldsteinPrice_min = {0.0, -1.0};
    const static inline Real OptND_GoldsteinPrice_min_f = 3.0;

    // 29. THREE-HUMP CAMEL: f(x,y) = 2x² - 1.05x⁴ + x⁶/6 + xy + y²
    // Global minimum: f(0, 0) = 0
    static Real OptND_ThreeHumpCamel(const VectorN<Real, 2>& x) {
        Real x2 = x[0] * x[0];
        Real x4 = x2 * x2;
        Real x6 = x4 * x2;
        return 2*x2 - 1.05*x4 + x6/6 + x[0]*x[1] + x[1]*x[1];
    }
    const static inline VectorN<Real, 2> OptND_ThreeHumpCamel_min = {0.0, 0.0};
    const static inline Real OptND_ThreeHumpCamel_min_f = 0.0;

    // 30. SIX-HUMP CAMEL: f(x,y) = (4-2.1x²+x⁴/3)x² + xy + (-4+4y²)y²
    // Two global minima at f ≈ -1.0316: (0.0898, -0.7126), (-0.0898, 0.7126)
    static Real OptND_SixHumpCamel(const VectorN<Real, 2>& x) {
        Real x2 = x[0] * x[0];
        Real y2 = x[1] * x[1];
        return (4 - 2.1*x2 + x2*x2/3)*x2 + x[0]*x[1] + (-4 + 4*y2)*y2;
    }
    const static inline VectorN<Real, 2> OptND_SixHumpCamel_min1 = {0.0898, -0.7126};
    const static inline VectorN<Real, 2> OptND_SixHumpCamel_min2 = {-0.0898, 0.7126};
    const static inline Real OptND_SixHumpCamel_min_f = -1.0316284535;

    // 31. EASOM: f(x,y) = -cos(x)cos(y)exp(-[(x-π)²+(y-π)²])
    // Global minimum: f(π, π) = -1
    // Unimodal but with flat regions
    static Real OptND_Easom(const VectorN<Real, 2>& x) {
        Real dx = x[0] - Constants::PI;
        Real dy = x[1] - Constants::PI;
        return -std::cos(x[0]) * std::cos(x[1]) * std::exp(-(dx*dx + dy*dy));
    }
    const static inline VectorN<Real, 2> OptND_Easom_min = {Constants::PI, Constants::PI};
    const static inline Real OptND_Easom_min_f = -1.0;

    // 32. STYBLINSKI-TANG: f(x) = 0.5 * Σ(xᵢ⁴ - 16xᵢ² + 5xᵢ)
    // Global minimum: f(-2.903534,...,-2.903534) = -39.16599*n
    template<int N>
    static Real OptND_StyblinskiTang(const VectorN<Real, N>& x) {
        Real sum = 0;
        for (int i = 0; i < N; ++i) {
            Real xi = x[i];
            Real xi2 = xi * xi;
            sum += xi2 * xi2 - 16 * xi2 + 5 * xi;
        }
        return 0.5 * sum;
    }
    const static inline Real OptND_StyblinskiTang_min_coord = -2.903534;
    const static inline Real OptND_StyblinskiTang_min_f_per_dim = -39.16599;

    // ========================================================================================
    // CATEGORY 5: ILL-CONDITIONED FUNCTIONS
    // ========================================================================================

    // 33. DIXON-PRICE: f(x) = (x₁-1)² + Σi(2xᵢ²-xᵢ₋₁)²
    // Global minimum: xᵢ = 2^(-(2^i-2)/2^i) for i=1,...,n
    template<int N>
    static Real OptND_DixonPrice(const VectorN<Real, N>& x) {
        Real sum = (x[0] - 1) * (x[0] - 1);
        for (int i = 1; i < N; ++i) {
            Real t = 2 * x[i] * x[i] - x[i-1];
            sum += (i + 1) * t * t;
        }
        return sum;
    }

    // 34. ZAKHAROV: f(x) = Σxᵢ² + (Σ0.5i·xᵢ)² + (Σ0.5i·xᵢ)⁴
    // Global minimum: f(0,...,0) = 0
    template<int N>
    static Real OptND_Zakharov(const VectorN<Real, N>& x) {
        Real sum1 = 0, sum2 = 0;
        for (int i = 0; i < N; ++i) {
            sum1 += x[i] * x[i];
            sum2 += 0.5 * (i + 1) * x[i];
        }
        return sum1 + sum2 * sum2 + sum2 * sum2 * sum2 * sum2;
    }

    // 35. POWELL: f(x) = Σ[(x₄ᵢ₋₃+10x₄ᵢ₋₂)² + 5(x₄ᵢ₋₁-x₄ᵢ)² + (x₄ᵢ₋₂-2x₄ᵢ₋₁)⁴ + 10(x₄ᵢ₋₃-x₄ᵢ)⁴]
    // Works for N divisible by 4. Global minimum: f(0,...,0) = 0
    static Real OptND_Powell4D(const VectorN<Real, 4>& x) {
        Real t1 = x[0] + 10 * x[1];
        Real t2 = x[2] - x[3];
        Real t3 = x[1] - 2 * x[2];
        Real t4 = x[0] - x[3];
        return t1*t1 + 5*t2*t2 + t3*t3*t3*t3 + 10*t4*t4*t4*t4;
    }
    const static inline VectorN<Real, 4> OptND_Powell4D_min = {0.0, 0.0, 0.0, 0.0};
    const static inline Real OptND_Powell4D_min_f = 0.0;

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                              SEARCH DOMAIN BOUNDS                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Typical search domains for benchmark functions
    const static inline Real Domain_Sphere_low = -5.12;
    const static inline Real Domain_Sphere_high = 5.12;

    const static inline Real Domain_Rosenbrock_low = -5.0;
    const static inline Real Domain_Rosenbrock_high = 10.0;

    const static inline Real Domain_Rastrigin_low = -5.12;
    const static inline Real Domain_Rastrigin_high = 5.12;

    const static inline Real Domain_Ackley_low = -32.768;
    const static inline Real Domain_Ackley_high = 32.768;

    const static inline Real Domain_Griewank_low = -600.0;
    const static inline Real Domain_Griewank_high = 600.0;

    const static inline Real Domain_Schwefel_low = -500.0;
    const static inline Real Domain_Schwefel_high = 500.0;

    const static inline Real Domain_Beale_low = -4.5;
    const static inline Real Domain_Beale_high = 4.5;

    const static inline Real Domain_GoldsteinPrice_low = -2.0;
    const static inline Real Domain_GoldsteinPrice_high = 2.0;

    const static inline Real Domain_Himmelblau_low = -5.0;
    const static inline Real Domain_Himmelblau_high = 5.0;

} // namespace MML::TestBeds

#endif // __MML_OPTIMIZATION_DEFS_H
