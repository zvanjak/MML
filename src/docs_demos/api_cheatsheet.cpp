///////////////////////////////////////////////////////////////////////////////////////////
// api_cheatsheet.cpp - Verification of all examples from API_CHEATSHEET.md
//
// This file ensures all code examples in API_CHEATSHEET.md compile and run correctly.
// Each section matches the corresponding section in the markdown file.
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Vector/VectorTypes.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixSym.h"
#include "base/Function.h"
#include "base/Polynom.h"
#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"
#include "core/LinAlgEqSolvers.h"
#include "core/Derivation.h"
#include "core/Integration.h"
#include "algorithms/RootFinding.h"
#include "algorithms/EigenSystemSolvers.h"
#include "mml/algorithms/RootFinding/RootFindingPolynoms.h"
#include "base/InterpolatedFunction.h"
#include "mml/algorithms/ODESolvers/ODEFixedStepIntegrators.h"
#include "mml/algorithms/ODESolvers/ODESystemStepCalculators.h"
#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#include "algorithms/FunctionsAnalyzer.h"
#include "core/Integration/GaussianQuadrature.h"
#include "tools/ConsolePrinter.h"
#endif

using namespace MML;
#include <iostream>
#include <cmath>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 1: VECTORS
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_vectors()
{
    std::cout << "\n=== 1. VECTORS ===" << std::endl;
    
    // Creation
    Vector<Real> v{1.0, 2.0, 3.0};               // Initialize with values
    Vector<Real> zeros(5);                        // 5 zeros
    Vector3Cartesian v3{1.0, 2.0, 3.0};          // 3D Cartesian vector
    
    // Operations
    Vector<Real> w{4.0, 5.0, 6.0};
    Vector<Real> sum = v + w;
    Vector<Real> diff = v - w;
    Vector<Real> scaled = 2.0 * v;
    Real dot = ScalarProduct(v3, Vector3Cartesian{4.0, 5.0, 6.0});  // Dot product
    Vector3Cartesian cross = VectorProduct(v3, Vector3Cartesian{4.0, 5.0, 6.0});  // Cross product
    
    // Norms
    Real l2 = v.NormL2();                         // Euclidean norm
    Real l1 = v.NormL1();                         // Manhattan norm
    Real linf = v.NormLInf();                     // Max norm
    
    std::cout << "v = " << v << std::endl;
    std::cout << "||v||_2 = " << l2 << std::endl;
    std::cout << "v . (4,5,6) = " << dot << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 2: MATRICES
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_matrices()
{
    std::cout << "\n=== 2. MATRICES ===" << std::endl;
    
    // Creation
    Matrix<Real> A{3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9}};  // Row-major
    Matrix<Real> I = Matrix<Real>::Identity(3);     // Identity
    Matrix<Real> Z(3, 3);                                // Zero matrix (default init)
    
    // Operations
    Matrix<Real> B{3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1}};
    Matrix<Real> sum_mat = A + B;
    Matrix<Real> prod = A * B;
    Matrix<Real> transposed = A.transpose();
    
    // Properties
    Real trace = A.trace();
    int rows = A.rows();
    int cols = A.cols();
    
    std::cout << "Matrix A:" << std::endl;
    A.Print(std::cout, 8, 3);
    std::cout << "Trace(A) = " << trace << std::endl;
    std::cout << "Size: " << rows << " x " << cols << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 3: LINEAR ALGEBRA
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_linear_algebra()
{
    std::cout << "\n=== 3. LINEAR ALGEBRA ===" << std::endl;
    
    // LU Solver
    Matrix<Real> A{3, 3, {4, 1, 2, 1, 5, 1, 2, 1, 6}};
    Vector<Real> b{8, 9, 12};
    LUSolver<Real> lu(A);
    Vector<Real> x = lu.Solve(b);
    std::cout << "LU Solution: " << x << std::endl;
    
    // Cholesky (for symmetric positive-definite)
    Matrix<Real> sympos{3, 3, {4, 1, 1, 1, 5, 1, 1, 1, 6}};
    Vector<Real> b2{8, 9, 12};
    CholeskySolver<Real> chol(sympos);
    Vector<Real> x2 = chol.Solve(b2);
    std::cout << "Cholesky Solution: " << x2 << std::endl;
    
    // QR Solver
    QRSolver<Real> qr(A);
    Vector<Real> x3 = qr.Solve(b);
    std::cout << "QR Solution: " << x3 << std::endl;
    
    // Eigenvalues (symmetric matrix)
    MatrixSym<Real> sym(3);
    sym(0,0) = 2; sym(0,1) = 1; sym(0,2) = 0;
                  sym(1,1) = 2; sym(1,2) = 1;
                                sym(2,2) = 2;
    auto result = SymmMatEigenSolverJacobi::Solve(sym);
    std::cout << "Eigenvalues: " << result.eigenvalues << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 4: POLYNOMIALS
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_polynomials()
{
    std::cout << "\n=== 4. POLYNOMIALS ===" << std::endl;
    
    // Create: p(x) = 1 + 2x + 3x²
    PolynomReal p{1.0, 2.0, 3.0};
    
    // Evaluate
    Real y = p(2.0);                              // p(2) = 1 + 4 + 12 = 17
    std::cout << "p(2) = " << y << std::endl;
    
    // Operations
    PolynomReal q{1.0, 1.0};                      // q(x) = 1 + x
    PolynomReal sum = p + q;
    PolynomReal product = p * q;
    
    // Polynomial roots: x² - 5x + 6 = 0  -> roots: 2, 3
    PolynomReal quad{6.0, -5.0, 1.0};             // c + bx + ax²
    auto roots = RootFinding::LaguerreRoots(quad);
    std::cout << "Roots of x^2 - 5x + 6: ";
    for (size_t i = 0; i < roots.size(); i++) {
        if (std::abs(roots[i].imag()) < 1e-10)
            std::cout << roots[i].real() << " ";
        else
            std::cout << roots[i] << " ";
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 5: DERIVATIVES
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_derivatives()
{
    std::cout << "\n=== 5. DERIVATIVES ===" << std::endl;
    
    RealFunction f{ [](Real x) { return x*x*x; } };  // f(x) = x³
    
    // First derivative (different accuracy orders)
    Real x0 = 2.0;
    Real df1 = Derivation::NDer1(f, x0);   // f'(x), 1st order accuracy
    Real df2 = Derivation::NDer2(f, x0);   // f'(x), 2nd order accuracy
    Real df4 = Derivation::NDer4(f, x0);   // f'(x), 4th order accuracy
    
    std::cout << "f(x) = x^3, f'(2) should be 12:" << std::endl;
    std::cout << "  NDer1: " << df1 << std::endl;
    std::cout << "  NDer2: " << df2 << std::endl;
    std::cout << "  NDer4: " << df4 << std::endl;
    
    // Second derivative
    Real d2f = Derivation::NSecDer2(f, x0);  // f''(x), 2nd order accuracy
    std::cout << "f''(2) = " << d2f << " (expected: 12)" << std::endl;
    
    // Third derivative  
    Real d3f = Derivation::NThirdDer2(f, x0);  // f'''(x), 2nd order accuracy
    std::cout << "f'''(2) = " << d3f << " (expected: 6)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 6: INTEGRATION
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_integration()
{
    std::cout << "\n=== 6. INTEGRATION ===" << std::endl;
    
    RealFunction f{ [](Real x) { return std::sin(x); } };
    Real a = 0.0, b = Constants::PI;
    
    // Romberg integration (adaptive, usually best choice)
    IntegrationResult res1 = IntegrateRomberg(f, a, b);
    std::cout << "Romberg: " << res1.value << " (converged: " << res1.converged << ")" << std::endl;
    
    // Trapezoidal
    IntegrationResult res2 = IntegrateTrap(f, a, b);
    std::cout << "Trapezoidal: " << res2.value << std::endl;
    
    // Simpson
    IntegrationResult res3 = IntegrateSimpson(f, a, b);
    std::cout << "Simpson: " << res3.value << std::endl;
    
    // Gauss-Legendre quadrature (adaptive order)
    IntegrationResult gl_result = IntegrateGaussLegendre(f, a, b, 10);
    std::cout << "Gauss-Legendre 10: " << gl_result.value << std::endl;
    
    std::cout << "Expected: 2.0" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 7: ROOT FINDING
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_root_finding()
{
    std::cout << "\n=== 7. ROOT FINDING ===" << std::endl;
    
    RealFunction f{ [](Real x) { return x*x - 2; } };  // Root is sqrt(2)
    
    // Brent's method (recommended)
    Real r1 = RootFinding::FindRootBrent(f, 0.0, 2.0, 1e-10);
    std::cout << "Brent: " << r1 << std::endl;
    
    // Bisection
    Real r2 = RootFinding::FindRootBisection(f, 0.0, 2.0, 1e-10);
    std::cout << "Bisection: " << r2 << std::endl;
    
    // Secant
    Real r3 = RootFinding::FindRootSecant(f, 1.0, 2.0, 1e-10);
    std::cout << "Secant: " << r3 << std::endl;
    
    // Newton's method (uses numerical derivatives)
    Real r4 = RootFinding::FindRootNewton(f, 0.0, 2.0, 1e-10);
    std::cout << "Newton: " << r4 << std::endl;
    
    std::cout << "Expected sqrt(2): " << std::sqrt(2.0) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 8: INTERPOLATION
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_interpolation()
{
    std::cout << "\n=== 8. INTERPOLATION ===" << std::endl;
    
    // Data points
    Vector<Real> x{0.0, 1.0, 2.0, 3.0};
    Vector<Real> y{0.0, 1.0, 4.0, 9.0};  // y = x²
    
    // Linear interpolation
    LinearInterpRealFunc lin_interp(x, y);
    Real y_lin = lin_interp(1.5);
    std::cout << "Linear interp at 1.5: " << y_lin << std::endl;
    
    // Polynomial interpolation (order = number of points)
    PolynomInterpRealFunc poly_interp(x, y, 4);  // 4-point polynomial
    Real y_poly = poly_interp(1.5);
    std::cout << "Polynomial interp at 1.5: " << y_poly << " (expected 2.25)" << std::endl;
    
    // Spline interpolation (natural cubic)
    SplineInterpRealFunc spline_interp(x, y);
    Real y_spline = spline_interp(1.5);
    std::cout << "Spline interp at 1.5: " << y_spline << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 9: ODEs
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_odes()
{
    std::cout << "\n=== 9. ODEs ===" << std::endl;
    
    // Define ODE system: dy/dt = -y (exponential decay)
    class DecayODE : public IODESystem {
    public:
        int getDim() const override { return 1; }
        void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
            dydt[0] = -y[0];
        }
    };
    
    DecayODE ode;
    Vector<Real> y0{1.0};
    Real t0 = 0.0, tend = 2.0;
    Real outputInterval = 0.1;
    
    // Adaptive integrator (Dormand-Prince) - recommended
    ODEAdaptiveIntegrator<DormandPrince5_Stepper> adaptive_solver(ode);
    ODESystemSolution sol1 = adaptive_solver.integrate(y0, t0, tend, outputInterval);
    std::cout << "Dormand-Prince: y(2) = " << sol1.getXValuesAtEnd()[0] << std::endl;
    
    // Fixed-step solver (RK4)
    RungeKutta4_StepCalculator rk4_calc;
    ODESystemFixedStepSolver fixed_solver(ode, rk4_calc);
    ODESystemSolution sol2 = fixed_solver.integrate(y0, t0, tend, 100);  // 100 steps
    std::cout << "RK4 fixed: y(2) = " << sol2.getXValuesAtEnd()[0] << std::endl;
    
    // Access solution data correctly
    std::cout << "Number of steps: " << sol1.getNumSteps() << std::endl;
    std::cout << "t values: " << sol1.getTValue(0) << " to " << sol1.getTValue(sol1.size()-1) << std::endl;
    
    std::cout << "Expected e^(-2): " << std::exp(-2.0) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// SECTION 10: FUNCTION ANALYSIS
///////////////////////////////////////////////////////////////////////////////////////////
void cheatsheet_function_analysis()
{
    std::cout << "\n=== 10. FUNCTION ANALYSIS ===" << std::endl;
    
    // Define function
    RealFunction f{ [](Real x) { return std::sin(x); } };
    Real x1 = 0.0, x2 = 2.0 * Constants::PI;
    
    // Create analyzer with the function
    RealFunctionAnalyzer analyzer(f);
    
    // Find roots (zeros)
    std::vector<Real> zeros = analyzer.GetRoots(x1, x2, 1e-6);
    std::cout << "Zeros of sin(x) in [0, 2*pi]:" << std::endl;
    for (Real z : zeros) {
        std::cout << "  x = " << z << std::endl;
    }
    
    // Find local optima (classified as min/max)
    std::vector<CriticalPoint> optima = analyzer.GetLocalOptimumsClassified(x1, x2, 1e-6);
    std::cout << "Local optima:" << std::endl;
    for (const auto& pt : optima) {
        const char* type_str = (pt.type == CriticalPointType::LOCAL_MAXIMUM) ? "max" : "min";
        std::cout << "  x = " << pt.x << " (" << type_str << "), value = " << pt.value << std::endl;
    }
    
    // Find inflection points
    Vector<Real> inflections = analyzer.GetInflectionPoints(x1, x2, 1e-6);
    std::cout << "Inflection points:" << std::endl;
    for (int i = 0; i < inflections.size(); i++) {
        std::cout << "  x = " << inflections[i] << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// MAIN VERIFICATION FUNCTION
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_ApiCheatsheet()
{
    std::cout << "=====================================================" << std::endl;
    std::cout << "  API CHEATSHEET - Example Verification" << std::endl;
    std::cout << "=====================================================" << std::endl;
    
    cheatsheet_vectors();
    cheatsheet_matrices();
    cheatsheet_linear_algebra();
    cheatsheet_polynomials();
    cheatsheet_derivatives();
    cheatsheet_integration();
    cheatsheet_root_finding();
    cheatsheet_interpolation();
    cheatsheet_odes();
    cheatsheet_function_analysis();
    
    std::cout << "\n=====================================================" << std::endl;
    std::cout << "  ALL API CHEATSHEET EXAMPLES VERIFIED SUCCESSFULLY!" << std::endl;
    std::cout << "=====================================================" << std::endl;
}

int main()
{
    Docs_Demo_ApiCheatsheet();
    return 0;
}
