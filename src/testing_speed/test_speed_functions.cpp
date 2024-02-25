#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

//#include <chrono>


/*
FUNKCIJA
- exp, log, log10, log2,
- pow, sqrt, cbrt
- sin, cos, tan, asin, acos, atan
- sinh, cosh, tanh, asinh, acosh, atanh
- Error and gamma func - erf, tgamma, lgamma
- Special func - assoc_laguerre, assoc_legendre, beta, comp_ellint_1(2,3), cyl_bessel_I(j,k)
                - hermite, legendre, lagurre, riemann_zeta, sph_bessel, sph_legendre
*/

inline double eval_hermite(int n, double x) { return std::hermite(n, x); }

inline double eval_hermite_1(double x) { return std::hermite(1, x); }
inline double eval_hermite_3(double x) { return std::hermite(3, x); }
inline double eval_hermite_5(double x) { return std::hermite(5, x); }
inline double eval_hermite_10(double x) { return std::hermite(10, x); }
inline double eval_hermite_20(double x) { return std::hermite(20, x); }
inline double eval_hermite_30(double x) { return std::hermite(30, x); }

inline double eval_legandre_1(double x) { return std::legendre(1, x); }

inline double eval_laguerre_1(double x) { return std::laguerre(1, x); }

inline double eval_sph_bessel_3(double x) { return std::sph_bessel(3, x); }

inline double eval_sph_legandre_3_1(double x) { return std::sph_legendre(3, 1, x); }

// esencijalno DVOPARAMETARSKE FUNKCIJE???
// pow, beta, ellint

struct FuncToEval
{
    double (*_func)(double);
    std::string _name;
    double _low;
    double _up;
};

std::vector<FuncToEval> vec_func = {
    {sqrt, "sqrt", 0.0, 1e9},
    {cbrt, "cbrt", -1e9, 1e9},
    {exp, "exp", -10, 10},
    {log, "log", 0.0, 1e6},
    {log10, "log10", 0.0, 1e6},
    {sin, "sin", -3.14159, 3.14159},
    {cos, "cos", -3.14159, 3.14159},
    {tan, "tan", -3.14159, 3.14159},
    {asin, "asin", -1, 1},
    {acos, "acos", -1, 1},
    {atan, "atan", -1e3, 1e3},
    {sinh, "sinh", -3.14159, 3.14159},
    {cosh, "cosh", -3.14159, 3.14159},
    {tanh, "tanh", -3.14159, 3.14159},
    {asinh, "asinh", -1, 1},
    {acosh, "acosh", -1, 1},
    {atanh, "atanh", -1e3, 1e3},
    {std::erf, "erf", -1e6, 1e6},
    {std::erfc, "erfc", -1e6, 1e6},
    {std::tgamma, "tgamma", -1e3, 1e3},
    {std::lgamma, "lgamma", -1e3, 1e3},
    {std::riemann_zeta, "riemann_zeta", -1e3, 1e3},
    {std::comp_ellint_1, "comp_ellint_1", -1e3, 1e3},
    {std::comp_ellint_2, "comp_ellint_2", -1e3, 1e3},
    {eval_hermite_1, "Hermite 1", -1e3, 1e3},
    {eval_hermite_3, "Hermite 3", -1e3, 1e3},
    {eval_hermite_5, "Hermite 5", -1e3, 1e3},
    {eval_hermite_10, "Hermite 10", -1e3, 1e3},
    {eval_hermite_20, "Hermite 20", -1e3, 1e3},
    {eval_hermite_30, "Hermite 30", -1e3, 1e3},
    {eval_legandre_1, "Legandre 1", -1e3, 1e3},
    {eval_laguerre_1, "Laguerre 1", -1e3, 1e3},
    {eval_sph_bessel_3, "Sph Bessel 3", -1e3, 1e3},
    {eval_sph_legandre_3_1, "Sph Legendre 3 1", -1e3, 1e3}
};

static const int num_evals = 1000000;

void Test_Speed_Functions()
{
    std::cout << "TESTING SPEED OF FUNCTION EVALUATION\n";

    // Messing up clang build, so commented out for now
    //using std::chrono::duration;
    //using std::chrono::duration_cast;
    //using std::chrono::high_resolution_clock;
    //using std::chrono::milliseconds;
    //double x = 0.0;
    //double y = 0.0;

    //for (auto f : vec_func)
    //{
    //    auto t1 = high_resolution_clock::now();
    //    for (int i = 0; i < num_evals; i++)
    //    {
    //        x = rand() % 1000 * (f._up - f._low) / 1000.0 + f._low;
    //        y = f._func(x);
    //    }
    //    auto t2 = high_resolution_clock::now();

    //    duration<double, std::milli> rand_time = t2 - t1;

    //    int num_repeats = 1000;
    //    double vals[10000];
    //    for (int i = 0; i < 10000; i++)
    //    {
    //        vals[i] = rand() % 1000 * (f._up - f._low) / 1000.0 + f._low;
    //    }

    //    auto t3 = high_resolution_clock::now();
    //    for (int j = 0; j < num_repeats; j++)
    //    {
    //        for (int i = 0; i < 10000; i++)
    //        {
    //            y = f._func(vals[i]);
    //        }
    //    }
    //    auto t4 = high_resolution_clock::now();

    //    duration<double, std::milli> array_time = t4 - t3;

    //    double diff = rand_time.count() - array_time.count(); 
    //    std::cout << std::setw(20) << f._name << " - " << rand_time.count() << " ; " << array_time.count() << " ; " << diff << "\n";
    //}
}
