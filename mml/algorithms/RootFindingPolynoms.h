///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFindingPolynoms.h                                               ///
///  Description: Polynomial root finding (Laguerre, Durand-Kerner, companion matrix) ///
///               Real and complex roots for polynomials of any degree                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_POLYNOMS_H
#define MML_ROOTFINDING_POLYNOMS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/Polynom.h"

#include "core/Derivation.h"

#include "algorithms/EigenSystemSolvers.h"

namespace MML
{
	namespace RootFinding
	{
    /***************************************************************************************************
    * POLYNOMIAL ROOT FINDING METHODS
    ***************************************************************************************************/

    /**
    * Find all roots of a polynomial using Laguerre's method.
    * 
    * Laguerre's method is a globally convergent iterative method for polynomial root finding.
    * 
    * @param poly      Polynomial with real coefficients
    * @param tol       Convergence tolerance (default: 1e-10)
    * @param maxIter   Maximum iterations per root (default: 100)
    * @return Vector of all complex roots
    */
    static Vector<Complex> LaguerreRoots(const PolynomReal& poly, Real tol = 1e-10, int maxIter = 100)
    {
      Vector<Complex> roots;
      int degree = poly.GetDegree();
      
      if (degree <= 0)
        return roots;
      
      // Work with coefficient array (not Polynom object due to operator[] template issues)
      Vector<Complex> workCoeffs(degree + 1);
      for (int i = 0; i <= degree; i++)
        workCoeffs[i] = std::complex<Real>(poly[i], 0.0);
      
      // Find roots one by one with deflation
      for (int currentDegree = degree; currentDegree >= 1; currentDegree--)
      {
        // Initial guess
        Complex x = std::complex<Real>(0.5, 0.5);
        
        // Laguerre iteration
        bool converged = false;
        for (int iter = 0; iter < maxIter; iter++)
        {
          // Evaluate polynomial and derivatives using Horner's method
          Complex p = workCoeffs[currentDegree];
          Complex dp = std::complex<Real>(0.0);
          Complex d2p = std::complex<Real>(0.0);
          
          for (int j = currentDegree - 1; j >= 0; j--)
          {
            d2p = d2p * x + dp;
            dp = dp * x + p;
            p = p * x + workCoeffs[j];
          }
          d2p *= Real(2.0);
          
          if (std::abs(p) < tol)
          {
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
          
          if (std::abs(dx) < tol)
          {
            converged = true;
            break;
          }
        }
        
        if (!converged)
          throw RootFindingError("Laguerre's method failed to converge");
        
        // Polish with Newton iterations
        for (int polish = 0; polish < 3; polish++)
        {
          Complex p = workCoeffs[currentDegree];
          Complex dp = std::complex<Real>(0.0);
          
          for (int j = currentDegree - 1; j >= 0; j--)
          {
            dp = dp * x + p;
            p = p * x + workCoeffs[j];
          }
          
          if (std::abs(dp) > tol)
            x -= p / dp;
        }
        
        roots.push_back(x);
        
        // Deflate polynomial
        if (currentDegree > 1)
        {
          Complex b = workCoeffs[currentDegree];
          for (int j = currentDegree - 1; j >= 0; j--)
          {
            Complex temp = workCoeffs[j];
            workCoeffs[j] = b;
            b = temp + x * b;
          }
        }
      }
      
      return roots;
    }

    /**
    * Find all roots using eigenvalue method (companion matrix).
    * 
    * @param poly  Polynomial with real coefficients
    * @return Vector of all complex roots
    */
    static Vector<Complex> EigenvalueRoots(const PolynomReal& poly)
    {
      Vector<Complex> roots;
      int degree = poly.GetDegree();
      
      if (degree <= 0)
        return roots;
      
      if (degree == 1)
      {
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

    /**
    * Find/refine polynomial roots using Bairstow's method.
    * 
    * @param poly          Polynomial with real coefficients
    * @param initialRoots  Initial root estimates (optional)
    * @param tol           Convergence tolerance (default: 1e-10)
    * @param maxIter       Maximum iterations per quadratic (default: 100)
    * @return Vector of complex roots
    */
    static Vector<Complex> BairstowRoots(const PolynomReal& poly, 
                                         const Vector<Complex>& initialRoots = Vector<Complex>(),
                                         Real tol = 1e-10, int maxIter = 100)
    {
      Vector<Complex> roots;
      int n = poly.GetDegree();
      
      if (n <= 0)
        return roots;
      
      // Work with coefficient array a[0..n] where poly = a[0] + a[1]x + ... + a[n]x^n
      Vector<Real> a(n + 1);
      for (int i = 0; i <= n; i++)
        a[i] = poly[i];
      
      while (n >= 2)
      {
        // Initial guess for quadratic factor x² + u*x + v
        // Start with a simple guess - these values tend to work for many polynomials
        Real u = 0.0;
        Real v = 1.0;
        
        // Use initial roots if provided
        if (initialRoots.size() >= 2 && n == poly.GetDegree())
        {
          Complex r1 = initialRoots[roots.size()];
          Complex r2 = (roots.size() + 1 < initialRoots.size()) ? 
                      initialRoots[roots.size() + 1] : std::conj(r1);
          u = -(r1 + r2).real();
          v = (r1 * r2).real();
        }
        
        // Bairstow's method iteration
        // Division: P(x) = (x² + ux + v) * Q(x) + (cx + d)
        // where Q(x) = sum_{i=0}^{n-2} b[i] x^i
        
        Vector<Real> b(n + 1, 0.0);  // Quotient coefficients + work space
        Vector<Real> f(n + 1, 0.0);  // Second division for Jacobian
        
        bool converged = false;
        for (int iter = 0; iter < maxIter; iter++)
        {
          // First synthetic division: compute b[i]
          // b[n] = b[n-1] = 0 (implied by initialization)
          // b[i] = a[i+2] - u*b[i+1] - v*b[i+2] for i = n-2, ..., 0
          b[n] = 0.0;
          b[n-1] = 0.0;
          for (int i = n - 2; i >= 0; i--)
            b[i] = a[i + 2] - u * b[i + 1] - v * b[i + 2];
          
          // Remainder: c*x + d
          Real c = a[1] - u * b[0] - v * b[1];
          Real d = a[0] - v * b[0];
          
          // Check convergence
          if (std::abs(c) < tol && std::abs(d) < tol)
          {
            converged = true;
            break;
          }
          
          // Second synthetic division on Q(x) to get partial derivatives
          // f[n] = f[n-1] = 0 (implied)
          // f[i] = b[i+2] - u*f[i+1] - v*f[i+2] for i = n-2, ..., 0
          f[n] = 0.0;
          f[n-1] = 0.0;
          for (int i = n - 2; i >= 0; i--)
            f[i] = b[i + 2] - u * f[i + 1] - v * f[i + 2];
          
          // g and h for Jacobian: g*x + h is remainder of Q(x) / (x² + ux + v)
          Real g = b[1] - u * f[0] - v * f[1];
          Real h = b[0] - v * f[0];
          
          // Newton update: [u,v]_new = [u,v] - J^{-1} * [c,d]
          // det(J) = v*g² + h*(h - u*g)
          Real det = v * g * g + h * (h - u * g);
          
          if (std::abs(det) < PrecisionValues<Real>::DeterminantZeroThreshold)
          {
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
          
          if (std::abs(du) < tol && std::abs(dv) < tol)
          {
            converged = true;
            break;
          }
        }
        
        // Extract roots from quadratic x² + ux + v = 0
        // x = (-u ± sqrt(u² - 4v)) / 2
        Real discriminant = u * u - Real(4.0) * v;
        if (discriminant >= 0)
        {
          Real sqrtDisc = std::sqrt(discriminant);
          roots.push_back(std::complex<Real>((-u + sqrtDisc) / Real(2.0), 0.0));
          roots.push_back(std::complex<Real>((-u - sqrtDisc) / Real(2.0), 0.0));
        }
        else
        {
          Real sqrtDisc = std::sqrt(-discriminant);
          roots.push_back(std::complex<Real>(-u / Real(2.0), sqrtDisc / Real(2.0)));
          roots.push_back(std::complex<Real>(-u / Real(2.0), -sqrtDisc / Real(2.0)));
        }
        
        // Deflate: new polynomial is Q(x) with coefficients b[0], b[1], ..., b[n-2]
        if (n > 2)
        {
          Vector<Real> newA(n - 1);
          for (int i = 0; i <= n - 2; i++)
            newA[i] = b[i];
          a = newA;
        }
        n -= 2;
      }
      
      // Handle remaining linear term: a[0] + a[1]*x = 0
      if (n == 1)
      {
        if (std::abs(a[1]) > PrecisionValues<Real>::PolynomialCoeffZeroThreshold)
          roots.push_back(std::complex<Real>(-a[0] / a[1], 0.0));
      }
      
      return roots;
    } // BairstowRoots
  } // namespace RootFinding
} // namespace MML

#endif // MML_ROOTFINDING_POLYNOMS_H


