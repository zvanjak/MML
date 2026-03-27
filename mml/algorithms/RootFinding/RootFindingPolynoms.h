///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFindingPolynoms.h                                               ///
///  Description: Polynomial root finding (Laguerre, Durand-Kerner, companion matrix) ///
///               Real and complex roots for polynomials of any degree                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_POLYNOMS_H
#define MML_ROOTFINDING_POLYNOMS_H

#include "mml/MMLBase.h"

#include "mml/interfaces/IFunction.h"

#include "mml/base/Vector/Vector.h"
#include "mml/base/Matrix/Matrix.h"
#include "mml/base/Polynom.h"

#include "mml/core/Derivation.h"

#include "mml/algorithms/EigenSystemSolvers.h"

namespace MML {
  /***************************************************************************************************
   * LOW-DEGREE POLYNOMIAL SOLVERS (Quadratic, Cubic, Quartic)
   * 
   * Closed-form solutions for polynomials of degree 2, 3, and 4.
   * These use the classical formulas: quadratic formula, Cardano's formula, and Ferrari's method.
   * Complexity: O(1) for all three — direct computation, no iteration.
   ***************************************************************************************************/

  /// Solve quadratic equation: a*x^2 + b*x + c = 0
  /// @return Number of real roots (0 or 2); complex roots returned via x1, x2
  inline int SolveQuadratic(Real a, Real b, Real c, Complex &x1, Complex &x2) {
    // Degenerate: a=0 reduces to linear equation b*x + c = 0
    if (std::abs(a) < Constants::Eps) {
      if (std::abs(b) < Constants::Eps) {
        x1 = x2 = Complex(0);
        return 0;
      }
      x1 = Complex(-c / b);
      x2 = Complex(0);
      return 1;
    }

    Real D = b * b - 4 * a * c;
    if (D >= 0) {
      // Numerically stable formula: avoids catastrophic cancellation when |b| >> |sqrt(D)|
      // Reference: Numerical Recipes, Press et al., "Quadratic and Cubic Equations"
      Real sqrtD = std::sqrt(D);
      Real q = Real(-0.5) * (b + std::copysign(sqrtD, b));
      x1 = q / a;
      x2 = c / q;
      return 2;
    } else {
      Complex sqrtD = std::sqrt(Complex(D));
      x1 = (-b + sqrtD) / (2 * a);
      x2 = (-b - sqrtD) / (2 * a);
      return 0;
    }
  }

  /// Solve quadratic equation with complex coefficients: a*x^2 + b*x + c = 0
  inline void SolveQuadratic(const Complex &a, const Complex &b, const Complex &c,
                             Complex &x1, Complex &x2) {
    // Degenerate: a=0 reduces to linear equation b*x + c = 0
    if (std::abs(a) < Constants::Eps) {
      x1 = (std::abs(b) < Constants::Eps) ? Complex(0) : -c / b;
      x2 = Complex(0);
      return;
    }

    Complex D = b * b - Real(4.0) * a * c;
    Complex sqrtD = std::sqrt(D);
    x1 = (-b + sqrtD) / (Real(2.0) * a);
    x2 = (-b - sqrtD) / (Real(2.0) * a);
  }

  /// Solve cubic equation: a*x^3 + b*x^2 + c*x + d = 0 (Cardano's formula)
  /// @return Number of real roots (1 or 3); all roots returned via x1, x2, x3
  inline int SolveCubic(Real a, Real b, Real c, Real d, Complex &x1, Complex &x2,
                        Complex &x3) {
    // Degenerate: a=0 reduces to quadratic b*x^2 + c*x + d = 0
    if (std::abs(a) < Constants::Eps) {
      int n = SolveQuadratic(b, c, d, x1, x2);
      x3 = Complex(0);
      return n;
    }

    // Normalize the coefficients
    Real A = b / a;
    Real B = c / a;
    Real C = d / a;

    // Calculate the discriminant
    Real Q = (Real(3.0) * B - POW2(A)) / Real(9.0);
    Real R =
        (Real(9.0) * A * B - Real(27.0) * C - Real(2.0) * POW3(A)) / Real(54.0);
    Real D = POW3(Q) + POW2(R); // Discriminant

    if (D >= 0) // Complex or duplicate roots
    {
      Real S = std::cbrt(R + std::sqrt(D));
      Real T = std::cbrt(R - std::sqrt(D));

      x1 = -A / Real(3.0) + (S + T); // Real root
      x2 = -A / Real(3.0) - (S + T) / Real(2.0) +
           Complex(0, std::sqrt(Real(3.0)) * (S - T) / Real(2.0)); // Complex root
      x3 = -A / Real(3.0) - (S + T) / Real(2.0) -
           Complex(0, std::sqrt(Real(3.0)) * (S - T) / Real(2.0)); // Complex root

      return 1;
    } else // Three real roots
    {
      Real theta = std::acos(R / std::sqrt(-POW3(Q)));
      x1 =
          Real(2.0) * std::sqrt(-Q) * std::cos(theta / Real(3.0)) - A / Real(3.0);
      x2 = Real(2.0) * std::sqrt(-Q) *
               std::cos((theta + Real(2.0) * Constants::PI) / Real(3.0)) -
           A / Real(3.0);
      x3 = Real(2.0) * std::sqrt(-Q) *
               std::cos((theta + Real(4.0) * Constants::PI) / Real(3.0)) -
           A / Real(3.0);

      return 3;
    }
  }

  /// Solve quartic equation: a*x^4 + b*x^3 + c*x^2 + d*x + e = 0 (Ferrari's method)
  /// All roots returned via x1, x2, x3, x4
  inline void SolveQuartic(Real a, Real b, Real c, Real d, Real e, Complex &x1,
                           Complex &x2, Complex &x3, Complex &x4) {
    // Degenerate: reduce to cubic
    if (std::abs(a) < Constants::Eps) {
      SolveCubic(b, c, d, e, x1, x2, x3);
      x4 = Complex(0);
      return;
    }

    // Normalize coefficients
    Real A = b / a;
    Real B = c / a;
    Real C = d / a;
    Real D = e / a;

    // Depressed quartic y = x + A/4: y^4 + P y^2 + Q y + R = 0
    Real AA = A * A;
    Real P = B - Real(3.0) * AA / Real(8.0);
    Real Q = C - Real(0.5) * A * B + AA * A / Real(8.0);
    Real R = D - Real(0.25) * A * C + AA * B / Real(16.0) -
             Real(3.0) * AA * AA / Real(256.0);

    // Special case: biquadratic (Q ≈ 0) -> solve t^2 + P t + R = 0 where t = y^2
    if (std::abs(Q) <= Constants::Eps) {
      Complex t1, t2;
      SolveQuadratic(Complex(Real(1.0)), Complex(P), Complex(R), t1, t2);

      x1 = std::sqrt(t1) - A / Real(4.0);
      x2 = -std::sqrt(t1) - A / Real(4.0);
      x3 = std::sqrt(t2) - A / Real(4.0);
      x4 = -std::sqrt(t2) - A / Real(4.0);
      return;
    }

    // General case (Ferrari)
    Complex z1, z2, z3;
    SolveCubic(Real(1.0), -P / Real(2.0), -R,
               R * P / Real(2.0) - Q * Q / Real(8.0), z1, z2, z3);

    auto U_from = [P](const Complex &z) {
      return std::sqrt(Complex(Real(2.0)) * z - P);
    };

    // Choose z to maximize |U| to avoid division by small numbers
    Complex candidates[3] = {z1, z2, z3};
    Complex z = candidates[0];
    Complex U = U_from(z);
    for (int i = 1; i < 3; ++i) {
      Complex Ui = U_from(candidates[i]);
      if (std::abs(Ui) > std::abs(U)) {
        z = candidates[i];
        U = Ui;
      }
    }

    Complex W = Q / (Complex(Real(2.0)) * U);

    // Solve two quadratics in y
    Complex y1, y2, y3, y4;
    SolveQuadratic(Complex(Real(1.0)), U, z - W, y1, y2);
    SolveQuadratic(Complex(Real(1.0)), -U, z + W, y3, y4);

    // Back-substitute x = y - A/4
    x1 = y1 - A / Real(4.0);
    x2 = y2 - A / Real(4.0);
    x3 = y3 - A / Real(4.0);
    x4 = y4 - A / Real(4.0);
  }

	namespace RootFinding {
		/***************************************************************************************************
    * POLYNOMIAL ROOT FINDING METHODS (General Degree)
    ***************************************************************************************************/

		/// Find all roots of a polynomial using Laguerre's method.
		///
		/// Laguerre's method is a globally convergent iterative method for polynomial root finding.
		/// Complexity: O(n² × maxIter) — n roots found by deflation, each requiring up to
		///            maxIter iterations with O(n) polynomial evaluation per iteration.
		///
		/// @param poly      Polynomial with real coefficients
		/// @param tol       Convergence tolerance (default: 1e-10)
		/// @param maxIter   Maximum iterations per root (default: 100)
		/// @return Vector of all complex roots
		static Vector<Complex> LaguerreRoots(const PolynomReal& poly, Real tol = 1e-10, int maxIter = 100) {
			Vector<Complex> roots;
			int degree = poly.degree();

			if (degree <= 0)
				return roots;

			// Work with coefficient array (not Polynom object due to operator[] template issues)
			Vector<Complex> workCoeffs(degree + 1);
			for (int i = 0; i <= degree; i++)
				workCoeffs[i] = std::complex<Real>(poly[i], 0.0);

			// Find roots one by one with deflation
			for (int currentDegree = degree; currentDegree >= 1; currentDegree--) {
				// Initial guess
				Complex x = std::complex<Real>(0.5, 0.5);

				// Laguerre iteration
				bool converged = false;
				for (int iter = 0; iter < maxIter; iter++) {
					// Evaluate polynomial and derivatives using Horner's method
					Complex p = workCoeffs[currentDegree];
					Complex dp = std::complex<Real>(0.0);
					Complex d2p = std::complex<Real>(0.0);

					for (int j = currentDegree - 1; j >= 0; j--) {
						d2p = d2p * x + dp;
						dp = dp * x + p;
						p = p * x + workCoeffs[j];
					}
					d2p *= Real(2.0);

					if (std::abs(p) < tol) {
						converged = true;
						break;
					}

					// Laguerre's formula
					Real n = Real(currentDegree);
					Complex G = dp / p;
					Complex H = G * G - d2p / p;
					Complex sq = std::sqrt((n - Real(1.0)) * (n * H - G * G));

					Complex denom1 = G + sq;
					Complex denom2 = G - sq;
					Complex denom = (std::abs(denom1) > std::abs(denom2)) ? denom1 : denom2;

					Complex dx = n / denom;
					x -= dx;

					if (std::abs(dx) < tol) {
						converged = true;
						break;
					}
				}

				if (!converged)
					throw RootFindingError("Laguerre's method failed to converge");

				// Polish with Newton iterations
				for (int polish = 0; polish < 3; polish++) {
					Complex p = workCoeffs[currentDegree];
					Complex dp = std::complex<Real>(0.0);

					for (int j = currentDegree - 1; j >= 0; j--) {
						dp = dp * x + p;
						p = p * x + workCoeffs[j];
					}

					if (std::abs(dp) > tol)
						x -= p / dp;
				}

				roots.push_back(x);

				// Deflate polynomial
				if (currentDegree > 1) {
					Complex b = workCoeffs[currentDegree];
					for (int j = currentDegree - 1; j >= 0; j--) {
						Complex temp = workCoeffs[j];
						workCoeffs[j] = b;
						b = temp + x * b;
					}
				}
			}

			return roots;
		}

		/// Find all roots using eigenvalue method (companion matrix).
		/// Complexity: O(n³) — constructs n×n companion matrix then solves eigenvalue problem.
		///            Most robust but most expensive method for large polynomials.
		///
		/// @param poly  Polynomial with real coefficients
		/// @return Vector of all complex roots
		static Vector<Complex> EigenvalueRoots(const PolynomReal& poly) {
			Vector<Complex> roots;
			int degree = poly.degree();

			if (degree <= 0)
				return roots;

			if (degree == 1) {
				roots.push_back(std::complex<Real>(-poly[0] / poly[1], 0.0));
				return roots;
			}

			// Construct companion matrix
			Matrix<Real> companion(degree, degree);
			for (int i = 0; i < degree; i++)
				for (int j = 0; j < degree; j++)
					companion[i][j] = 0.0;

			// Superdiagonal of ones
			for (int i = 0; i < degree - 1; i++)
				companion[i][i + 1] = 1.0;

			// Last row: -a_0/a_n, -a_1/a_n, ..., -a_{n-1}/a_n
			Real leadingCoef = poly[degree];
			if (std::abs(leadingCoef) < PrecisionValues<Real>::PolynomialCoeffZeroThreshold)
				throw RootFindingError("Polynomial leading coefficient is zero");

			for (int i = 0; i < degree; i++)
				companion[degree - 1][i] = -poly[i] / leadingCoef;

			// Solve eigenvalue problem
			auto eigenResult = EigenSolver::Solve(companion, 1e-10, 1000);

			if (!eigenResult.converged)
				throw RootFindingError("Eigenvalue solver failed to converge");

			// Extract eigenvalues as roots
			for (const auto& eval : eigenResult.eigenvalues)
				roots.push_back(std::complex<Real>(eval.real, eval.imag));

			return roots;
		}

		/// Find/refine polynomial roots using Bairstow's method.
		/// Complexity: O(n² × maxIter) — extracts n/2 quadratic factors, each requiring
		///            up to maxIter iterations with O(n) synthetic division per iteration.
		///
		/// @param poly          Polynomial with real coefficients
		/// @param initialRoots  Initial root estimates (optional)
		/// @param tol           Convergence tolerance (default: 1e-10)
		/// @param maxIter       Maximum iterations per quadratic (default: 100)
		/// @return Vector of complex roots
		static Vector<Complex> BairstowRoots(const PolynomReal& poly, const Vector<Complex>& initialRoots = Vector<Complex>(),
											 Real tol = 1e-10, int maxIter = 100) {
			Vector<Complex> roots;
			int n = poly.degree();

			if (n <= 0)
				return roots;

			// Work with coefficient array a[0..n] where poly = a[0] + a[1]x + ... + a[n]x^n
			Vector<Real> a(n + 1);
			for (int i = 0; i <= n; i++)
				a[i] = poly[i];

			while (n >= 2) {
				// Initial guess for quadratic factor x² + u*x + v
				// Start with a simple guess - these values tend to work for many polynomials
				Real u = 0.0;
				Real v = 1.0;

				// Use initial roots if provided
				if (initialRoots.size() >= 2 && n == poly.degree()) {
					Complex r1 = initialRoots[roots.size()];
					Complex r2 = (roots.size() + 1 < initialRoots.size()) ? initialRoots[roots.size() + 1] : std::conj(r1);
					u = -(r1 + r2).real();
					v = (r1 * r2).real();
				}

				// Bairstow's method iteration
				// Division: P(x) = (x² + ux + v) * Q(x) + (cx + d)
				// where Q(x) = sum_{i=0}^{n-2} b[i] x^i

				Vector<Real> b(n + 1, 0.0); // Quotient coefficients + work space
				Vector<Real> f(n + 1, 0.0); // Second division for Jacobian

				bool converged = false;
				for (int iter = 0; iter < maxIter; iter++) {
					// First synthetic division: compute b[i]
					// b[n] = b[n-1] = 0 (implied by initialization)
					// b[i] = a[i+2] - u*b[i+1] - v*b[i+2] for i = n-2, ..., 0
					b[n] = 0.0;
					b[n - 1] = 0.0;
					for (int i = n - 2; i >= 0; i--)
						b[i] = a[i + 2] - u * b[i + 1] - v * b[i + 2];

					// Remainder: c*x + d
					Real c = a[1] - u * b[0] - v * b[1];
					Real d = a[0] - v * b[0];

					// Check convergence
					if (std::abs(c) < tol && std::abs(d) < tol) {
						converged = true;
						break;
					}

					// Second synthetic division on Q(x) to get partial derivatives
					// f[n] = f[n-1] = 0 (implied)
					// f[i] = b[i+2] - u*f[i+1] - v*f[i+2] for i = n-2, ..., 0
					f[n] = 0.0;
					f[n - 1] = 0.0;
					for (int i = n - 2; i >= 0; i--)
						f[i] = b[i + 2] - u * f[i + 1] - v * f[i + 2];

					// g and h for Jacobian: g*x + h is remainder of Q(x) / (x² + ux + v)
					Real g = b[1] - u * f[0] - v * f[1];
					Real h = b[0] - v * f[0];

					// Newton update: [u,v]_new = [u,v] - J^{-1} * [c,d]
					// det(J) = v*g² + h*(h - u*g)
					Real det = v * g * g + h * (h - u * g);

					if (std::abs(det) < PrecisionValues<Real>::DeterminantZeroThreshold) {
						// Singular Jacobian - perturb and continue
						u += 0.1;
						v += 0.1;
						continue;
					}

					// Update from Wikipedia formula: [u,v] = [u,v] - J^{-1} * [c,d]
					Real du = (-h * c + g * d) / det;
					Real dv = (-g * v * c + (g * u - h) * d) / det;

					u -= du;
					v -= dv;

					if (std::abs(du) < tol && std::abs(dv) < tol) {
						converged = true;
						break;
					}
				}

				if (!converged)
					throw std::runtime_error("BairstowRoots: failed to converge for quadratic factor");

				// Extract roots from quadratic x² + ux + v = 0
				// x = (-u ± sqrt(u² - 4v)) / 2
				Real discriminant = u * u - Real(4.0) * v;
				if (discriminant >= 0) {
					Real sqrtDisc = std::sqrt(discriminant);
					roots.push_back(std::complex<Real>((-u + sqrtDisc) / Real(2.0), 0.0));
					roots.push_back(std::complex<Real>((-u - sqrtDisc) / Real(2.0), 0.0));
				} else {
					Real sqrtDisc = std::sqrt(-discriminant);
					roots.push_back(std::complex<Real>(-u / Real(2.0), sqrtDisc / Real(2.0)));
					roots.push_back(std::complex<Real>(-u / Real(2.0), -sqrtDisc / Real(2.0)));
				}

				// Deflate: new polynomial is Q(x) with coefficients b[0], b[1], ..., b[n-2]
				if (n > 2) {
					Vector<Real> newA(n - 1);
					for (int i = 0; i <= n - 2; i++)
						newA[i] = b[i];
					a = newA;
				}
				n -= 2;
			}

			// Handle remaining linear term: a[0] + a[1]*x = 0
			if (n == 1) {
				if (std::abs(a[1]) > PrecisionValues<Real>::PolynomialCoeffZeroThreshold)
					roots.push_back(std::complex<Real>(-a[0] / a[1], 0.0));
			}

			return roots;
		} // BairstowRoots
	} // namespace RootFinding
} // namespace MML

#endif // MML_ROOTFINDING_POLYNOMS_H
