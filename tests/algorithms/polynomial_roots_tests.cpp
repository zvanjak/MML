#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/RootFindingPolynoms.h"
#include "base/Polynom.h"
#endif

#include <algorithm>
#include <cmath>

using namespace MML;
using namespace MML::Testing;

// Helper function to check if a complex number is close to expected
bool isCloseComplex(const Complex& actual, const Complex& expected, Real tol = 1e-8)
{
  return std::abs(actual - expected) < tol;
}

// Helper function to sort roots by magnitude for easier comparison
void sortRootsByMagnitude(Vector<Complex>& roots)
{
  std::sort(roots.begin(), roots.end(), [](const Complex& a, const Complex& b) {
    Real magA = std::abs(a);
    Real magB = std::abs(b);
    if (std::abs(magA - magB) < 1e-10)
      return a.real() < b.real();
    return magA < magB;
  });
}

// Helper to verify roots actually satisfy polynomial
bool verifyRoot(const PolynomReal& poly, const Complex& root, Real tol = 1e-6)
{
  // Evaluate polynomial at complex root using Horner's method
  int degree = poly.GetDegree();
  Complex result = Complex(poly[degree], REAL(0.0));
  
  for (int i = degree - 1; i >= 0; i--)
    result = result * root + Complex(poly[i], REAL(0.0));
  
  return std::abs(result) < tol;
}

TEST_CASE("PolynomialRoots::Laguerre - Simple polynomials", "[polynomial][roots][laguerre]")
{
  SECTION("x^2 - 1 (two real roots)")
  {
    PolynomReal poly({-1, 0, 1});  // -1 + 0*x + 1*x^2
    auto roots = RootFinding::LaguerreRoots(poly);
    
    REQUIRE(roots.size() == 2);
    sortRootsByMagnitude(roots);
    
    CHECK(isCloseComplex(roots[0], Complex(-1, 0)));
    CHECK(isCloseComplex(roots[1], Complex(1, 0)));
    
    // Verify roots
    for (const auto& root : roots)
      CHECK(verifyRoot(poly, root));
  }
  
  SECTION("x^2 + 1 (two complex roots)")
  {
    PolynomReal poly({1, 0, 1});  // 1 + 0*x + 1*x^2
    auto roots = RootFinding::LaguerreRoots(poly);
    
    REQUIRE(roots.size() == 2);
    
    // Roots should be ±i
    bool foundPlusI = false, foundMinusI = false;
    for (const auto& root : roots)
    {
      if (isCloseComplex(root, Complex(0, 1)))
        foundPlusI = true;
      if (isCloseComplex(root, Complex(0, -1)))
        foundMinusI = true;
    }
    CHECK(foundPlusI);
    CHECK(foundMinusI);
    
    // Verify roots
    for (const auto& root : roots)
      CHECK(verifyRoot(poly, root));
  }
  
  SECTION("x^3 - 1 (cube roots of unity)")
  {
    PolynomReal poly({-1, 0, 0, 1});  // -1 + 0*x + 0*x^2 + 1*x^3
    auto roots = RootFinding::LaguerreRoots(poly);
    
    REQUIRE(roots.size() == 3);
    
    // One real root at 1, two complex roots at e^(±2πi/3)
    bool foundReal = false, foundComplex1 = false, foundComplex2 = false;
    for (const auto& root : roots)
    {
      if (isCloseComplex(root, Complex(1, 0)))
        foundReal = true;
      if (isCloseComplex(root, Complex(-REAL(0.5), std::sqrt(REAL(3.0))/REAL(2.0))))
        foundComplex1 = true;
      if (isCloseComplex(root, Complex(-REAL(0.5), -std::sqrt(REAL(3.0))/REAL(2.0))))
        foundComplex2 = true;
    }
    CHECK(foundReal);
    CHECK(foundComplex1);
    CHECK(foundComplex2);
  }
}

TEST_CASE("PolynomialRoots::Eigenvalue - Simple polynomials", "[polynomial][roots][eigenvalue]")
{
  SECTION("x^2 - 1")
  {
    PolynomReal poly({-1, 0, 1});
    auto roots = RootFinding::EigenvalueRoots(poly);
    
    REQUIRE(roots.size() == 2);
    sortRootsByMagnitude(roots);
    
    CHECK(isCloseComplex(roots[0], Complex(-1, 0)));
    CHECK(isCloseComplex(roots[1], Complex(1, 0)));
  }
  
  SECTION("x^2 + 1")
  {
    PolynomReal poly({1, 0, 1});
    auto roots = RootFinding::EigenvalueRoots(poly);
    
    REQUIRE(roots.size() == 2);
    
    bool foundPlusI = false, foundMinusI = false;
    for (const auto& root : roots)
    {
      if (isCloseComplex(root, Complex(0, 1)))
        foundPlusI = true;
      if (isCloseComplex(root, Complex(0, -1)))
        foundMinusI = true;
    }
    CHECK(foundPlusI);
    CHECK(foundMinusI);
  }
  
  SECTION("x^3 - 8 (cube roots of 8)")
  {
    PolynomReal poly({-8, 0, 0, 1});  // -8 + 0*x + 0*x^2 + 1*x^3
    auto roots = RootFinding::EigenvalueRoots(poly);
    
    REQUIRE(roots.size() == 3);
    
    // One real root at 2, two complex roots
    bool foundReal = false;
    for (const auto& root : roots)
    {
      if (isCloseComplex(root, Complex(2, 0)))
        foundReal = true;
      
      // Verify the root
      CHECK(verifyRoot(poly, root));
    }
    CHECK(foundReal);
  }
  
  SECTION("Linear polynomial")
  {
    PolynomReal poly({-5, 1});  // -5 + x (root at x=5)
    auto roots = RootFinding::EigenvalueRoots(poly);
    
    REQUIRE(roots.size() == 1);
    CHECK(isCloseComplex(roots[0], Complex(5, 0)));
  }
}

TEST_CASE("PolynomialRoots::Bairstow - Simple polynomials", "[polynomial][roots][bairstow]")
{
  SECTION("x^2 - 1")
  {
    PolynomReal poly({-1, 0, 1});
    auto roots = RootFinding::BairstowRoots(poly);
    
    REQUIRE(roots.size() == 2);
    sortRootsByMagnitude(roots);
    
    CHECK(isCloseComplex(roots[0], Complex(-1, 0)));
    CHECK(isCloseComplex(roots[1], Complex(1, 0)));
  }
  
  SECTION("x^2 + 1")
  {
    PolynomReal poly({1, 0, 1});
    auto roots = RootFinding::BairstowRoots(poly);
    
    REQUIRE(roots.size() == 2);
    
    bool foundPlusI = false, foundMinusI = false;
    for (const auto& root : roots)
    {
      if (isCloseComplex(root, Complex(0, 1)))
        foundPlusI = true;
      if (isCloseComplex(root, Complex(0, -1)))
        foundMinusI = true;
    }
    CHECK(foundPlusI);
    CHECK(foundMinusI);
  }
  
  SECTION("(x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6")
  {
    PolynomReal poly({-6, 11, -6, 1});
    auto roots = RootFinding::BairstowRoots(poly);
    
    REQUIRE(roots.size() == 3);
    
    bool found1 = false, found2 = false, found3 = false;
    for (const auto& root : roots)
    {
      if (isCloseComplex(root, Complex(1, 0)))
        found1 = true;
      if (isCloseComplex(root, Complex(2, 0)))
        found2 = true;
      if (isCloseComplex(root, Complex(3, 0)))
        found3 = true;
    }
    CHECK(found1);
    CHECK(found2);
    CHECK(found3);
  }
}

TEST_CASE("PolynomialRoots::Method comparison - Higher degree", "[polynomial][roots][comparison]")
{
  SECTION("Quartic: (x-1)(x-2)(x+1)(x+2) = x^4 - 5x^2 + 4")
  {
    PolynomReal poly({4, 0, -5, 0, 1});
    
    auto rootsLaguerre = RootFinding::LaguerreRoots(poly);
    auto rootsEigenvalue = RootFinding::EigenvalueRoots(poly);
    auto rootsBairstow = RootFinding::BairstowRoots(poly);
    
    REQUIRE(rootsLaguerre.size() == 4);
    REQUIRE(rootsEigenvalue.size() == 4);
    REQUIRE(rootsBairstow.size() == 4);
    
    // All methods should find roots at -2, -1, 1, 2
    Vector<Complex> expected(4);
    expected[0] = Complex(-2, 0);
    expected[1] = Complex(-1, 0);
    expected[2] = Complex(1, 0);
    expected[3] = Complex(2, 0);
    
    for (const auto& expectedRoot : expected)
    {
      // Check Laguerre
      bool foundInLaguerre = false;
      for (const auto& root : rootsLaguerre)
        if (isCloseComplex(root, expectedRoot))
          foundInLaguerre = true;
      CHECK(foundInLaguerre);
      
      // Check Eigenvalue
      bool foundInEigen = false;
      for (const auto& root : rootsEigenvalue)
        if (isCloseComplex(root, expectedRoot))
          foundInEigen = true;
      CHECK(foundInEigen);
      
      // Check Bairstow
      bool foundInBairstow = false;
      for (const auto& root : rootsBairstow)
        if (isCloseComplex(root, expectedRoot))
          foundInBairstow = true;
      CHECK(foundInBairstow);
    }
  }
  
  SECTION("Quintic with complex roots: x^5 - 1")
  {
    PolynomReal poly({-1, 0, 0, 0, 0, 1});
    
    auto rootsLaguerre = RootFinding::LaguerreRoots(poly);
    auto rootsEigenvalue = RootFinding::EigenvalueRoots(poly);
    
    REQUIRE(rootsLaguerre.size() == 5);
    REQUIRE(rootsEigenvalue.size() == 5);
    
    // Verify roots (5th roots of unity)
    for (const auto& root : rootsLaguerre)
    {
      Complex powered = root;
      for (int i = 1; i < 5; i++)
        powered *= root;
      CHECK(isCloseComplex(powered, Complex(1, 0)));
    }
    
    for (const auto& root : rootsEigenvalue)
    {
      Complex powered = root;
      for (int i = 1; i < 5; i++)
        powered *= root;
      CHECK(isCloseComplex(powered, Complex(1, 0)));
    }
  }
}

TEST_CASE("PolynomialRoots::Chebyshev polynomial T_4", "[polynomial][roots][special]")
{
  // T_4(x) = 8x^4 - 8x^2 + 1
  // Has 4 real roots at cos(π/8), cos(3π/8), cos(5π/8), cos(7π/8)
  PolynomReal poly({1, 0, -8, 0, 8});
  
  auto roots = RootFinding::LaguerreRoots(poly);
  REQUIRE(roots.size() == 4);
  
  // All roots should be real and in [-1, 1]
  for (const auto& root : roots)
  {
    CHECK(std::abs(root.imag()) < 1e-6);  // Nearly real
    CHECK(std::abs(root.real()) <= REAL(1.0));   // In [-1, 1]
    CHECK(verifyRoot(poly, root));
  }
}

TEST_CASE("PolynomialRoots::Edge cases", "[polynomial][roots][edge]")
{
  SECTION("Constant polynomial (degree 0)")
  {
    PolynomReal poly = PolynomReal::Constant(5.0);  // Constant polynomial p(x) = 5
    auto roots = RootFinding::LaguerreRoots(poly);
    CHECK(roots.isEmpty());  // A constant polynomial has no roots
  }
  
  SECTION("Linear polynomial")
  {
    PolynomReal poly({-3, 1});  // x - 3
    auto roots = RootFinding::LaguerreRoots(poly);
    REQUIRE(roots.size() == 1);
    CHECK(isCloseComplex(roots[0], Complex(3, 0)));
  }
  
  SECTION("Double root: (x-2)^2 = x^2 - 4x + 4")
  {
    PolynomReal poly({4, -4, 1});
    auto roots = RootFinding::LaguerreRoots(poly);
    
    REQUIRE(roots.size() == 2);
    // Double roots are numerically challenging - use relaxed tolerance
    for (const auto& root : roots)
      CHECK(isCloseComplex(root, Complex(2, 0), 1e-4));  // Relaxed tolerance for double root
  }
}

TEST_CASE("PolynomialRoots::Wilkinson polynomial (ill-conditioned)", "[polynomial][roots][wilkinson][.]")
{
  // Wilkinson polynomial: (x-1)(x-2)(x-3)(x-4)(x-5)
  // This is a CLASSIC ill-conditioned problem demonstrating numerical instability.
  // The polynomial coefficients are: -120, 274, -225, 85, -15, 1
  // Even tiny coefficient perturbations cause large root changes.
  // 
  // This test is SKIPPED by default (tagged with [.]) because:
  // - It demonstrates a known numerical analysis limitation
  // - Results vary by platform, compiler, and floating-point settings
  // - The test is educational, not a regression test
  //
  // Run with: MML_Tests.exe [wilkinson] to execute this test explicitly
  
  PolynomReal poly({-120, 274, -225, 85, -15, 1});  // Correct expanded form with alternating signs
  
  auto roots = RootFinding::EigenvalueRoots(poly);
  
  REQUIRE(roots.size() == 5);
  
  // Print roots for educational purposes
  INFO("Wilkinson polynomial roots (should be 1,2,3,4,5):");
  for (size_t i = 0; i < roots.size(); ++i) {
    INFO("  root[" << i << "] = " << roots[i].real() << " + " << roots[i].imag() << "i");
  }
  
  // With correct coefficients and eigenvalue method, this should work
  Vector<Complex> expected(5);
  expected[0] = Complex(1, 0);
  expected[1] = Complex(2, 0);
  expected[2] = Complex(3, 0);
  expected[3] = Complex(4, 0);
  expected[4] = Complex(5, 0);
  
  for (const auto& expectedRoot : expected)
  {
    bool found = false;
    for (const auto& root : roots)
      if (isCloseComplex(root, expectedRoot, 0.1))
        found = true;
    CHECK(found);
  }
}
