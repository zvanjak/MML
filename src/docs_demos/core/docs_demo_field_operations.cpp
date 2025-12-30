#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "core/Fields.h"
#include "core/FieldOperations.h"
#include "core/MetricTensor.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////
///                           PREDEFINED FIELDS DEMOS
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Predefined_Potential_Fields()
{
    std::cout << "=== PREDEFINED POTENTIAL FIELDS (Fields.h) ===" << std::endl;
    
    using namespace Fields;
    
    // Inverse radial potential: φ(r) = k/r (gravitational or electrostatic)
    std::cout << "\n1. Inverse Radial Potential Field (Cartesian)" << std::endl;
    
    // Free function approach
    Vec3Cart pos1(3.0, 4.0, 0.0);  // Distance 5 from origin
    Real phi1 = InverseRadialPotentialFieldCart(pos1);  // constant = 1
    Real phi2 = InverseRadialPotentialFieldCart(-10.0, pos1);  // constant = -10
    
    std::cout << "  Position: (3, 4, 0), |r| = " << pos1.NormL2() << std::endl;
    std::cout << "  Potential (k=1):  " << phi1 << " (expected: 0.2 = 1/5)" << std::endl;
    std::cout << "  Potential (k=-10): " << phi2 << " (expected: -2 = -10/5)" << std::endl;
    
    // Class-based approach (implements IScalarFunction<3>)
    std::cout << "\n2. InverseRadialFieldCart class" << std::endl;
    Real GM = 3.986e14;  // Earth GM in m³/s² (approximately)
    InverseRadialFieldCart gravPotential(-GM);
    
    Vec3Cart earthSurface(6.371e6, 0.0, 0.0);  // Earth surface on x-axis
    Real phi_surface = gravPotential(earthSurface);
    
    std::cout << "  Earth gravitational potential at surface:" << std::endl;
    std::cout << "  φ = " << phi_surface << " J/kg" << std::endl;
    std::cout << "  (Expected: ~-6.26×10⁷ J/kg)" << std::endl;
    
    // Spherical coordinates
    std::cout << "\n3. Inverse Radial Potential (Spherical)" << std::endl;
    Vec3Sph posSph(5.0, Constants::PI/4, Constants::PI/2);  // r=5, theta=pi/4, phi=pi/2
    Real phiSph = InverseRadialPotentialFieldSpher(posSph);
    
    std::cout << "  Position: (r=5, θ=π/4, φ=π/2)" << std::endl;
    std::cout << "  Potential: " << phiSph << " (expected: 0.2 = 1/r)" << std::endl;
}

void Docs_Demo_Predefined_Force_Fields()
{
    std::cout << "\n=== PREDEFINED FORCE FIELDS (Fields.h) ===" << std::endl;
    
    using namespace Fields;
    
    // Inverse-square force: F = -k·r/|r|³ (gravitational or Coulomb)
    std::cout << "\n1. Inverse Radial Force Field (Cartesian)" << std::endl;
    
    Vec3Cart pos(3.0, 4.0, 0.0);  // Distance 5 from origin
    Vec3Cart force1 = InverseRadialPotentialForceFieldCart(pos);
    Vec3Cart force2 = InverseRadialPotentialForceFieldCart(100.0, pos);
    
    std::cout << "  Position: (3, 4, 0), |r| = " << pos.NormL2() << std::endl;
    std::cout << "  Force (k=1): "; force1.Print(std::cout, 8, 4); std::cout << std::endl;
    std::cout << "    |F| = " << force1.NormL2() << " (expected: 1/25 = 0.04)" << std::endl;
    std::cout << "  Force (k=100): "; force2.Print(std::cout, 8, 4); std::cout << std::endl;
    std::cout << "    |F| = " << force2.NormL2() << " (expected: 100/25 = 4)" << std::endl;
    
    // Class-based approach
    std::cout << "\n2. InverseRadialForceFieldCart class" << std::endl;
    Real k_coulomb = 8.99e9;  // Coulomb constant
    Real q1_q2 = 1.0e-18;     // Product of two elementary charges
    InverseRadialForceFieldCart coulombForce(k_coulomb * q1_q2);
    
    Vec3Cart separation(1.0e-10, 0.0, 0.0);  // 0.1 nm separation
    Vec3Cart F_coulomb = coulombForce(separation);
    
    std::cout << "  Coulomb force between two electrons at 0.1 nm:" << std::endl;
    std::cout << "  F = "; F_coulomb.Print(std::cout, 8, 4); std::cout << std::endl;
    std::cout << "  |F| = " << F_coulomb.NormL2() << " N (repulsive)" << std::endl;
    
    // Spherical coordinates
    std::cout << "\n3. Inverse Radial Force (Spherical)" << std::endl;
    Vec3Sph posSph(5.0, Constants::PI/4, 0.0);
    Vec3Sph forceSph = InverseRadialPotentialForceFieldSph(posSph);
    
    std::cout << "  Position: (r=5, θ=π/4, φ=0)" << std::endl;
    std::cout << "  Force: (F_r=" << forceSph[0] << ", F_θ=" << forceSph[1] << ", F_φ=" << forceSph[2] << ")" << std::endl;
    std::cout << "  Expected: F_r = -1/r² = -0.04, others = 0" << std::endl;
}

void Docs_Demo_Predefined_Fields_With_Operations()
{
    std::cout << "\n=== PREDEFINED FIELDS + OPERATIONS ===" << std::endl;
    
    using namespace Fields;
    using namespace ScalarFieldOperations;
    using namespace VectorFieldOperations;
    
    // Demonstrate that F = -∇φ for inverse radial fields
    std::cout << "\n1. Verify F = -∇φ (gradient of potential)" << std::endl;
    
    Real k = 10.0;
    InverseRadialFieldCart potential(k);
    InverseRadialForceFieldCart force(k);
    
    Vec3Cart pos(2.0, 1.0, 2.0);
    
    // Compute gradient of potential
    Vec3Cart gradPhi = GradientCart<3>(potential, pos);
    Vec3Cart negGradPhi = -1.0 * gradPhi;
    
    // Get force from predefined field
    Vec3Cart F = force(pos);
    
    std::cout << "  Position: "; pos.Print(std::cout, 8, 4); std::cout << std::endl;
    std::cout << "  -∇φ = "; negGradPhi.Print(std::cout, 10, 6); std::cout << std::endl;
    std::cout << "  F   = "; F.Print(std::cout, 10, 6); std::cout << std::endl;
    std::cout << "  Match: " << ((negGradPhi - F).NormL2() < 1e-6 ? "YES" : "NO") << std::endl;
    
    // Verify ∇²φ = 0 for r > 0 (Laplace equation)
    std::cout << "\n2. Verify ∇²φ = 0 (harmonic away from origin)" << std::endl;
    
    Real laplacian = LaplacianCart<3>(potential, pos);
    std::cout << "  Laplacian ∇²φ = " << laplacian << " (expected: ~0)" << std::endl;
    
    // Verify ∇·F for 1/r² force
    std::cout << "\n3. Verify ∇·F = 0 for inverse-square force" << std::endl;
    
    Real divF = DivCart<3>(force, pos);
    std::cout << "  Divergence ∇·F = " << divF << " (expected: ~0 away from source)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                           GRADIENT DEMOS
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_GradientCartesian()
{
    std::cout << "=== GRADIENT - CARTESIAN ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Example 1: Scalar field f(x,y,z) = x² + y² + z²
    // Gradient should be 2*(x, y, z)
    std::cout << "\n1. Radial distance squared: f(x,y,z) = x² + y² + z²" << std::endl;
    
    auto field1 = [](const Vec3Cart& pos) -> Real {
        return pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
    };
    ScalarFunctionFromStdFunc<3> scalarField1(field1);
    
    Vec3Cart pos1(1.0, 2.0, 3.0);
    Vec3Cart grad1 = GradientCart(scalarField1, pos1);
    
    std::cout << "  Position: (" << pos1[0] << ", " << pos1[1] << ", " << pos1[2] << ")" << std::endl;
    std::cout << "  Gradient: (" << grad1[0] << ", " << grad1[1] << ", " << grad1[2] << ")" << std::endl;
    std::cout << "  Expected: (2, 4, 6)" << std::endl;
    
    // Example 2: Different derivative orders
    std::cout << "\n2. Derivative accuracy comparison:" << std::endl;
    for (int order : {1, 2, 4, 6, 8}) {
        Vec3Cart grad = GradientCart(scalarField1, pos1, order);
        std::cout << "  Order " << order << ": (" << grad[0] << ", " << grad[1] << ", " << grad[2] << ")" << std::endl;
    }
    
    // Example 3: 2D case
    std::cout << "\n3. 2D scalar field: f(x,y) = x*y" << std::endl;
    auto field2D = [](const VectorN<Real,2>& pos) -> Real {
        return pos[0] * pos[1];
    };
    ScalarFunctionFromStdFunc<2> scalarField2D(field2D);
    
    VectorN<Real,2> pos2D{3.0, 4.0};
    VectorN<Real,2> grad2D = GradientCart(scalarField2D, pos2D);
    
    std::cout << "  Position: (" << pos2D[0] << ", " << pos2D[1] << ")" << std::endl;
    std::cout << "  Gradient: (" << grad2D[0] << ", " << grad2D[1] << ")" << std::endl;
    std::cout << "  Expected: (4, 3) -- ∂(xy)/∂x = y, ∂(xy)/∂y = x" << std::endl;
}

void Docs_Demo_GradientSpherical()
{
    std::cout << "\n=== GRADIENT - SPHERICAL ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Example 1: Radial field f(r,θ,φ) = 1/r
    // Gradient: (-1/r², 0, 0)
    std::cout << "\n1. Coulomb-like potential: f(r,θ,φ) = 1/r" << std::endl;
    
    auto radialField = [](const Vec3Sph& pos) -> Real {
        return 1.0 / pos[0];  // pos[0] = r
    };
    ScalarFunctionFromStdFunc<3> phi(radialField);
    
    Vec3Sph pos(2.0, Constants::PI/4, 0.0);  // r=2, θ=π/4, φ=0
    Vec3Sph grad = GradientSpher(phi, pos);
    
    std::cout << "  Position: (r=" << pos[0] << ", θ=" << pos[1] << ", φ=" << pos[2] << ")" << std::endl;
    std::cout << "  Gradient: (∂f/∂r=" << grad[0] << ", grad_θ=" << grad[1] << ", grad_φ=" << grad[2] << ")" << std::endl;
    std::cout << "  Expected: (-0.25, 0, 0) -- derivative of 1/r is -1/r²" << std::endl;
    
    // Example 2: Angular dependence f(r,θ,φ) = r*cos(θ) = z
    std::cout << "\n2. f(r,θ,φ) = r*cos(θ) (same as z in Cartesian)" << std::endl;
    
    auto zField = [](const Vec3Sph& pos) -> Real {
        return pos[0] * cos(pos[1]);  // r*cos(θ)
    };
    ScalarFunctionFromStdFunc<3> z(zField);
    
    Vec3Sph pos2(1.0, Constants::PI/4, 0.0);
    Vec3Sph grad2 = GradientSpher(z, pos2);
    
    std::cout << "  Position: (r=" << pos2[0] << ", θ=π/4, φ=0)" << std::endl;
    std::cout << "  Gradient: (" << grad2[0] << ", " << grad2[1] << ", " << grad2[2] << ")" << std::endl;
    std::cout << "  Expected: (cos(π/4), -sin(π/4), 0) ≈ (0.707, -0.707, 0)" << std::endl;
}

void Docs_Demo_GradientCylindrical()
{
    std::cout << "\n=== GRADIENT - CYLINDRICAL ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Example: f(ρ,φ,z) = ρ² + z²
    // Gradient: (2ρ, 0, 2z)
    std::cout << "\n1. f(ρ,φ,z) = ρ² + z²" << std::endl;
    
    auto field = [](const Vec3Cyl& pos) -> Real {
        return pos[0]*pos[0] + pos[2]*pos[2];  // ρ² + z²
    };
    ScalarFunctionFromStdFunc<3> scalarField(field);
    
    Vec3Cyl pos(2.0, Constants::PI/2, 3.0);  // ρ=2, φ=π/2, z=3
    Vec3Cyl grad = GradientCyl(scalarField, pos);
    
    std::cout << "  Position: (ρ=" << pos[0] << ", φ=" << pos[1] << ", z=" << pos[2] << ")" << std::endl;
    std::cout << "  Gradient: (" << grad[0] << ", " << grad[1] << ", " << grad[2] << ")" << std::endl;
    std::cout << "  Expected: (4, 0, 6) -- (2ρ, 0, 2z)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                           LAPLACIAN DEMOS
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LaplacianCartesian()
{
    std::cout << "\n=== LAPLACIAN - CARTESIAN ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Example 1: Harmonic function f(x,y) = x² - y²
    // Laplacian: 2 - 2 = 0 (harmonic!)
    std::cout << "\n1. Harmonic function: f(x,y) = x² - y²" << std::endl;
    
    auto harmonic = [](const VectorN<Real,2>& pos) -> Real {
        return pos[0]*pos[0] - pos[1]*pos[1];
    };
    ScalarFunctionFromStdFunc<2> f(harmonic);
    
    VectorN<Real,2> pos{1.0, 2.0};
    Real lapl = LaplacianCart(f, pos);
    
    std::cout << "  Position: (" << pos[0] << ", " << pos[1] << ")" << std::endl;
    std::cout << "  Laplacian: " << lapl << std::endl;
    std::cout << "  Expected: 0 (harmonic function!)" << std::endl;
    
    // Example 2: Paraboloid f(x,y,z) = x² + y² + z²
    // Laplacian: 2 + 2 + 2 = 6
    std::cout << "\n2. Paraboloid: f(x,y,z) = x² + y² + z²" << std::endl;
    
    auto paraboloid = [](const Vec3Cart& pos) -> Real {
        return pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
    };
    ScalarFunctionFromStdFunc<3> g(paraboloid);
    
    Vec3Cart pos3D(1.0, 2.0, 3.0);
    Real lapl3D = LaplacianCart(g, pos3D);
    
    std::cout << "  Position: (" << pos3D[0] << ", " << pos3D[1] << ", " << pos3D[2] << ")" << std::endl;
    std::cout << "  Laplacian: " << lapl3D << std::endl;
    std::cout << "  Expected: 6" << std::endl;
}

void Docs_Demo_LaplacianSpherical()
{
    std::cout << "\n=== LAPLACIAN - SPHERICAL ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Example: Radial function f(r,θ,φ) = 1/r
    // Laplacian: 0 for r ≠ 0 (fundamental solution of Laplace equation)
    std::cout << "\n1. Coulomb potential: f(r,θ,φ) = 1/r" << std::endl;
    
    auto coulomb = [](const Vec3Sph& pos) -> Real {
        return 1.0 / pos[0];
    };
    ScalarFunctionFromStdFunc<3> phi(coulomb);
    
    Vec3Sph pos(2.0, Constants::PI/4, 0.0);
    Real lapl = LaplacianSpher(phi, pos);
    
    std::cout << "  Position: (r=" << pos[0] << ", θ=π/4, φ=0)" << std::endl;
    std::cout << "  Laplacian: " << lapl << std::endl;
    std::cout << "  Expected: ~0 (harmonic for r > 0)" << std::endl;
}

void Docs_Demo_LaplacianCylindrical()
{
    std::cout << "\n=== LAPLACIAN - CYLINDRICAL ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Example: f(ρ,φ,z) = ln(ρ) -- 2D Coulomb potential
    // Laplacian: 0 for ρ > 0
    std::cout << "\n1. 2D Coulomb: f(ρ,φ,z) = ln(ρ)" << std::endl;
    
    auto coulomb2D = [](const Vec3Cyl& pos) -> Real {
        return log(pos[0]);
    };
    ScalarFunctionFromStdFunc<3> phi(coulomb2D);
    
    Vec3Cyl pos(2.0, 0.0, 1.0);
    Real lapl = LaplacianCyl(phi, pos);
    
    std::cout << "  Position: (ρ=" << pos[0] << ", φ=" << pos[1] << ", z=" << pos[2] << ")" << std::endl;
    std::cout << "  Laplacian: " << lapl << std::endl;
    std::cout << "  Expected: ~0 (harmonic for ρ > 0)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                           DIVERGENCE DEMOS
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DivergenceCartesian()
{
    std::cout << "\n=== DIVERGENCE - CARTESIAN ===" << std::endl;
    
    using namespace VectorFieldOperations;
    
    // Example 1: Radial field F(x,y,z) = (x, y, z)
    // Divergence: 1 + 1 + 1 = 3 (expanding source)
    std::cout << "\n1. Expanding radial field: F(x,y,z) = (x, y, z)" << std::endl;
    
    auto radialField = [](const Vec3Cart& pos) -> Vec3Cart {
        return pos;
    };
    VectorFunctionFromStdFunc<3> F(radialField);
    
    Vec3Cart pos(1.0, 2.0, 3.0);
    Real div = DivCart(F, pos);
    
    std::cout << "  Position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
    std::cout << "  Divergence: " << div << std::endl;
    std::cout << "  Expected: 3" << std::endl;
    
    // Example 2: Incompressible flow F(x,y,z) = (-y, x, 0)
    // Divergence: 0 (rotation, no sources/sinks)
    std::cout << "\n2. Incompressible rotation: F(x,y,z) = (-y, x, 0)" << std::endl;
    
    auto rotationField = [](const Vec3Cart& pos) -> Vec3Cart {
        return Vec3Cart(-pos[1], pos[0], 0.0);
    };
    VectorFunctionFromStdFunc<3> G(rotationField);
    
    Real divG = DivCart(G, pos);
    
    std::cout << "  Divergence: " << divG << std::endl;
    std::cout << "  Expected: 0 (incompressible flow)" << std::endl;
}

void Docs_Demo_DivergenceSpherical()
{
    std::cout << "\n=== DIVERGENCE - SPHERICAL ===" << std::endl;
    
    using namespace VectorFieldOperations;
    
    // Example: Radial field in spherical coords F = (r, 0, 0)
    // Divergence: (1/r²)∂(r²·r)/∂r = 3
    std::cout << "\n1. Radial field: F = (r, 0, 0) in spherical" << std::endl;
    
    auto radialSph = [](const Vec3Sph& pos) -> Vec3Sph {
        return Vec3Sph(pos[0], 0.0, 0.0);  // F_r = r, F_θ = 0, F_φ = 0
    };
    VectorFunctionFromStdFunc<3> F(radialSph);
    
    Vec3Sph pos(2.0, Constants::PI/4, 0.0);
    Real div = DivSpher(F, pos);
    
    std::cout << "  Position: (r=" << pos[0] << ", θ=π/4, φ=0)" << std::endl;
    std::cout << "  Divergence: " << div << std::endl;
    std::cout << "  Expected: 3" << std::endl;
}

void Docs_Demo_DivergenceCylindrical()
{
    std::cout << "\n=== DIVERGENCE - CYLINDRICAL ===" << std::endl;
    
    using namespace VectorFieldOperations;
    
    // Example: Cylindrical radial field F = (ρ, 0, 0)
    // Divergence: (1/ρ)∂(ρ·ρ)/∂ρ = 2
    std::cout << "\n1. Radial field: F = (ρ, 0, 0) in cylindrical" << std::endl;
    
    auto radialCyl = [](const Vec3Cyl& pos) -> Vec3Cyl {
        return Vec3Cyl(pos[0], 0.0, 0.0);  // F_ρ = ρ
    };
    VectorFunctionFromStdFunc<3> F(radialCyl);
    
    Vec3Cyl pos(2.0, 0.0, 1.0);
    Real div = DivCyl(F, pos);
    
    std::cout << "  Position: (ρ=" << pos[0] << ", φ=" << pos[1] << ", z=" << pos[2] << ")" << std::endl;
    std::cout << "  Divergence: " << div << std::endl;
    std::cout << "  Expected: 2" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                           CURL DEMOS
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_CurlCartesian()
{
    std::cout << "\n=== CURL - CARTESIAN ===" << std::endl;
    
    using namespace VectorFieldOperations;
    
    // Example 1: Rotation field F(x,y,z) = (-y, x, 0)
    // Curl: (0, 0, 2) - constant rotation around z-axis
    std::cout << "\n1. Rotation field: F(x,y,z) = (-y, x, 0)" << std::endl;
    
    auto rotationField = [](const Vec3Cart& pos) -> Vec3Cart {
        return Vec3Cart(-pos[1], pos[0], 0.0);
    };
    VectorFunctionFromStdFunc<3> F(rotationField);
    
    Vec3Cart pos(1.0, 2.0, 0.0);
    Vec3Cart curl = CurlCart(F, pos);
    
    std::cout << "  Position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
    std::cout << "  Curl: (" << curl[0] << ", " << curl[1] << ", " << curl[2] << ")" << std::endl;
    std::cout << "  Expected: (0, 0, 2)" << std::endl;
    
    // Example 2: Conservative field (gradient of potential)
    // F = ∇(x²y + yz²) = (2xy, x² + z², 2yz)
    // Curl of gradient = 0
    std::cout << "\n2. Conservative field (gradient): F = ∇(x²y + yz²)" << std::endl;
    
    auto conservativeField = [](const Vec3Cart& pos) -> Vec3Cart {
        Real x = pos[0], y = pos[1], z = pos[2];
        return Vec3Cart(2*x*y, x*x + z*z, 2*y*z);
    };
    VectorFunctionFromStdFunc<3> G(conservativeField);
    
    Vec3Cart pos2(1.0, 2.0, 3.0);
    Vec3Cart curlG = CurlCart(G, pos2);
    
    std::cout << "  Curl: (" << curlG[0] << ", " << curlG[1] << ", " << curlG[2] << ")" << std::endl;
    std::cout << "  Expected: (0, 0, 0) -- curl of gradient is zero" << std::endl;
    
    Real norm = curlG.NormL2();
    std::cout << "  |Curl|: " << norm << (norm < 1e-10 ? " -- CONSERVATIVE FIELD!" : "") << std::endl;
}

void Docs_Demo_CurlSpherical()
{
    std::cout << "\n=== CURL - SPHERICAL ===" << std::endl;
    
    using namespace VectorFieldOperations;
    
    // Example: Rotational field around z-axis in spherical
    std::cout << "\n1. Azimuthal field: F = (0, 0, r*sin(θ))" << std::endl;
    
    auto azimuthalField = [](const Vec3Sph& pos) -> Vec3Sph {
        return Vec3Sph(0.0, 0.0, pos[0] * sin(pos[1]));
    };
    VectorFunctionFromStdFunc<3> F(azimuthalField);
    
    Vec3Sph pos(1.0, Constants::PI/2, 0.0);
    Vec3Sph curl = CurlSpher(F, pos);
    
    std::cout << "  Position: (r=" << pos[0] << ", θ=π/2, φ=0)" << std::endl;
    std::cout << "  Curl: (r=" << curl[0] << ", θ=" << curl[1] << ", φ=" << curl[2] << ")" << std::endl;
}

void Docs_Demo_CurlCylindrical()
{
    std::cout << "\n=== CURL - CYLINDRICAL ===" << std::endl;
    
    using namespace VectorFieldOperations;
    
    // Example: Cylindrical vortex F = (0, ρ, 0)
    // Curl: (0, 0, 2) - vorticity
    std::cout << "\n1. Solid-body rotation: F = (0, ρ, 0) in cylindrical" << std::endl;
    
    auto vortexField = [](const Vec3Cyl& pos) -> Vec3Cyl {
        return Vec3Cyl(0.0, pos[0], 0.0);  // F_φ = ρ
    };
    VectorFunctionFromStdFunc<3> F(vortexField);
    
    Vec3Cyl pos(2.0, 0.0, 0.0);
    Vec3Cyl curl = CurlCyl(F, pos);
    
    std::cout << "  Position: (ρ=" << pos[0] << ", φ=" << pos[1] << ", z=" << pos[2] << ")" << std::endl;
    std::cout << "  Curl: (ρ=" << curl[0] << ", φ=" << curl[1] << ", z=" << curl[2] << ")" << std::endl;
    std::cout << "  Expected z-component: 2 (solid body rotation vorticity)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                           PHYSICS EXAMPLES
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Electric_Field()
{
    std::cout << "\n=== PHYSICS: ELECTRIC FIELD FROM POTENTIAL ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Point charge potential: φ(r) = k·Q/r
    Real k = 8.99e9;   // Coulomb constant
    Real Q = 1.0e-9;   // 1 nC charge
    
    auto potential = [k, Q](const Vec3Cart& pos) -> Real {
        Real r = pos.NormL2();
        return (r > 1e-10) ? k * Q / r : 0.0;
    };
    ScalarFunctionFromStdFunc<3> phi(potential);
    
    // Electric field: E = -∇φ
    Vec3Cart pos(1.0, 0.0, 0.0);  // 1 meter from charge on x-axis
    Vec3Cart gradPhi = GradientCart(phi, pos);
    Vec3Cart E = -1.0 * gradPhi;
    
    std::cout << "Point charge: Q = " << Q*1e9 << " nC" << std::endl;
    std::cout << "Position: 1 m from charge on x-axis" << std::endl;
    std::cout << "Potential: " << phi(pos) << " V" << std::endl;
    std::cout << "Electric field E = -∇φ = (" << E[0] << ", " << E[1] << ", " << E[2] << ") N/C" << std::endl;
    std::cout << "Expected: E ≈ " << k*Q << " N/C in x-direction" << std::endl;
    
    // Verify: Laplacian away from charge should be ~0
    Real lapl = LaplacianCart(phi, pos);
    std::cout << "Laplacian ∇²φ = " << lapl << " (should be ~0 away from charge)" << std::endl;
}

void Docs_Demo_Incompressible_Flow()
{
    std::cout << "\n=== PHYSICS: INCOMPRESSIBLE FLUID FLOW ===" << std::endl;
    
    using namespace VectorFieldOperations;
    
    // 3D vortex flow: v = (-y, x, 0) / (x² + y²)^0.5
    auto velocity = [](const Vec3Cart& pos) -> Vec3Cart {
        Real r2 = pos[0]*pos[0] + pos[1]*pos[1];
        if (r2 < 1e-10) return Vec3Cart(0, 0, 0);
        Real r = sqrt(r2);
        return Vec3Cart(-pos[1]/r, pos[0]/r, 0.0);
    };
    VectorFunctionFromStdFunc<3> v(velocity);
    
    Vec3Cart pos(2.0, 0.0, 0.0);
    
    // Check incompressibility: ∇·v = 0
    Real div = DivCart(v, pos);
    std::cout << "Vortex flow at (2, 0, 0):" << std::endl;
    std::cout << "  Divergence ∇·v = " << div << std::endl;
    std::cout << "  " << (std::abs(div) < 1e-6 ? "Flow is INCOMPRESSIBLE!" : "Flow has sources/sinks") << std::endl;
    
    // Compute vorticity: ω = ∇×v
    Vec3Cart vorticity = CurlCart(v, pos);
    std::cout << "  Vorticity ∇×v = (" << vorticity[0] << ", " << vorticity[1] << ", " << vorticity[2] << ")" << std::endl;
}

void Docs_Demo_Heat_Diffusion()
{
    std::cout << "\n=== PHYSICS: HEAT DIFFUSION ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Temperature field: T(x,y) = T₀ * sin(πx) * sin(πy)
    Real T0 = 100.0;  // Peak temperature
    
    auto temperature = [T0](const VectorN<Real,2>& pos) -> Real {
        return T0 * sin(Constants::PI * pos[0]) * sin(Constants::PI * pos[1]);
    };
    ScalarFunctionFromStdFunc<2> T(temperature);
    
    // Heat equation: ∂T/∂t = κ·∇²T
    VectorN<Real,2> pos{0.5, 0.5};  // Center of domain
    Real laplacian = LaplacianCart(T, pos);
    
    std::cout << "Temperature at center (0.5, 0.5): " << T(pos) << " K" << std::endl;
    std::cout << "Laplacian ∇²T: " << laplacian << " K/m²" << std::endl;
    
    // Rate of change with thermal diffusivity κ = 1
    Real kappa = 1.0;
    Real dT_dt = kappa * laplacian;
    std::cout << "Rate of cooling: dT/dt = κ∇²T = " << dT_dt << " K/s" << std::endl;
    std::cout << "  (Negative = cooling, positive = heating)" << std::endl;
}

void Docs_Demo_Conservative_Field_Check()
{
    std::cout << "\n=== PHYSICS: CONSERVATIVE FIELD CHECK ===" << std::endl;
    
    using namespace VectorFieldOperations;
    
    // A field is conservative if ∇×F = 0
    // This is true when F = ∇φ for some potential φ
    
    // Test 1: Known conservative field F = ∇(x²y + yz²)
    std::cout << "\n1. Testing F = (2xy, x² + z², 2yz) - gradient of x²y + yz²" << std::endl;
    
    auto conservative = [](const Vec3Cart& pos) -> Vec3Cart {
        Real x = pos[0], y = pos[1], z = pos[2];
        return Vec3Cart(2*x*y, x*x + z*z, 2*y*z);
    };
    VectorFunctionFromStdFunc<3> F_cons(conservative);
    
    Vec3Cart pos(1.0, 2.0, 3.0);
    Vec3Cart curl_cons = CurlCart(F_cons, pos);
    Real norm_cons = curl_cons.NormL2();
    
    std::cout << "  Curl: (" << curl_cons[0] << ", " << curl_cons[1] << ", " << curl_cons[2] << ")" << std::endl;
    std::cout << "  |Curl| = " << norm_cons << std::endl;
    std::cout << "  " << (norm_cons < 1e-8 ? "CONSERVATIVE (curl ≈ 0)" : "NOT conservative") << std::endl;
    
    // Test 2: Non-conservative field F = (-y, x, 0)
    std::cout << "\n2. Testing F = (-y, x, 0) - rotation field" << std::endl;
    
    auto nonconservative = [](const Vec3Cart& pos) -> Vec3Cart {
        return Vec3Cart(-pos[1], pos[0], 0.0);
    };
    VectorFunctionFromStdFunc<3> F_noncons(nonconservative);
    
    Vec3Cart curl_noncons = CurlCart(F_noncons, pos);
    Real norm_noncons = curl_noncons.NormL2();
    
    std::cout << "  Curl: (" << curl_noncons[0] << ", " << curl_noncons[1] << ", " << curl_noncons[2] << ")" << std::endl;
    std::cout << "  |Curl| = " << norm_noncons << std::endl;
    std::cout << "  " << (norm_noncons < 1e-8 ? "CONSERVATIVE (curl ≈ 0)" : "NOT CONSERVATIVE (has rotation)") << std::endl;
}

void Docs_Demo_Spherical_Harmonics()
{
    std::cout << "\n=== PHYSICS: SPHERICAL HARMONICS ANALYSIS ===" << std::endl;
    
    using namespace ScalarFieldOperations;
    
    // Spherical harmonic Y₁₀ ∝ cos(θ)
    // For spherical harmonics: ∇²Yₗₘ = -l(l+1)/r² Yₗₘ on the unit sphere
    
    Real normalization = sqrt(3.0 / (4.0 * Constants::PI));
    
    auto Y10 = [normalization](const Vec3Sph& pos) -> Real {
        return normalization * cos(pos[1]);  // cos(θ)
    };
    ScalarFunctionFromStdFunc<3> harmonic(Y10);
    
    Vec3Sph pos(1.0, Constants::PI/4, 0.0);  // Unit sphere, θ = π/4
    Real value = harmonic(pos);
    Real lapl = LaplacianSpher(harmonic, pos);
    
    // For l=1: eigenvalue is -l(l+1) = -2
    Real expected = -2.0 * value;
    
    std::cout << "Y₁₀(θ=π/4) = " << value << std::endl;
    std::cout << "Laplacian ∇²Y₁₀ = " << lapl << std::endl;
    std::cout << "Expected: -l(l+1)·Y₁₀ = -2·" << value << " = " << expected << std::endl;
    std::cout << "Match: " << (std::abs(lapl - expected) < 0.1 ? "YES" : "APPROXIMATE") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                           MAIN ENTRY POINT
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Field_operations()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*****        SCALAR AND VECTOR FIELD OPERATIONS     *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;

    // Predefined physical fields from Fields.h
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****         PREDEFINED PHYSICAL FIELDS            *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    
    Docs_Demo_Predefined_Potential_Fields();
    Docs_Demo_Predefined_Force_Fields();
    Docs_Demo_Predefined_Fields_With_Operations();

    // Gradient demonstrations
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****              GRADIENT OPERATIONS              *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    
    Docs_Demo_GradientCartesian();
    Docs_Demo_GradientSpherical();
    Docs_Demo_GradientCylindrical();
    
    // Laplacian demonstrations
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****             LAPLACIAN OPERATIONS              *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    
    Docs_Demo_LaplacianCartesian();
    Docs_Demo_LaplacianSpherical();
    Docs_Demo_LaplacianCylindrical();
    
    // Divergence demonstrations
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****            DIVERGENCE OPERATIONS              *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    
    Docs_Demo_DivergenceCartesian();
    Docs_Demo_DivergenceSpherical();
    Docs_Demo_DivergenceCylindrical();
    
    // Curl demonstrations
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****               CURL OPERATIONS                 *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    
    Docs_Demo_CurlCartesian();
    Docs_Demo_CurlSpherical();
    Docs_Demo_CurlCylindrical();
    
    // Physics application examples
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****         PHYSICS APPLICATION EXAMPLES          *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    
    Docs_Demo_Electric_Field();
    Docs_Demo_Incompressible_Flow();
    Docs_Demo_Heat_Diffusion();
    Docs_Demo_Conservative_Field_Check();
    Docs_Demo_Spherical_Harmonics();
    
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****          FIELD OPERATIONS DEMO COMPLETE       *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;
}
