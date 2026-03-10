///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme11_coord_transforms.cpp                                       ///
///  Description: README example - Coordinate Transformations                         ///
///               Demonstrates Cartesian/Spherical/Cylindrical conversions            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/VectorN.h"
#include "base/Function.h"
#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/FieldOperations.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Readme_CoordTransforms()
{
    std::cout << std::endl;
    std::cout << "=== Coordinate Transformations ===" << std::endl;

    // Cartesian point
    Vector3Cartesian cart_pos{1.0, 1.0, 1.0};
    std::cout << "Cartesian point: " << cart_pos << std::endl;

    // Convert to Spherical (r, θ, φ) - Math/ISO convention
    // θ = polar angle from z-axis, φ = azimuthal angle in xy-plane
    Vector3Spherical sph_pos = CoordTransfCartToSpher.transf(cart_pos);
    std::cout << std::endl << "Convert to Spherical (r, θ, φ):" << std::endl;
    std::cout << "  r     = " << sph_pos[0] << " (should be √3 ≈ 1.732)" << std::endl;
    std::cout << "  theta = " << sph_pos[1] << " rad (" << sph_pos[1]*180/Constants::PI << "°)" << std::endl;
    std::cout << "  phi   = " << sph_pos[2] << " rad (" << sph_pos[2]*180/Constants::PI << "°)" << std::endl;

    // Convert to Cylindrical (r, φ, z)
    Vector3Cylindrical cyl_pos = CoordTransfCartToCyl.transf(cart_pos);
    std::cout << std::endl << "Convert to Cylindrical (ρ, φ, z):" << std::endl;
    std::cout << "  rho = " << cyl_pos[0] << " (should be √2 ≈ 1.414)" << std::endl;
    std::cout << "  phi = " << cyl_pos[1] << " rad (45°)" << std::endl;
    std::cout << "  z   = " << cyl_pos[2] << std::endl;

    // Convert back to Cartesian
    Vector3Cartesian back = CoordTransfSpherToCart.transf(sph_pos);
    std::cout << std::endl << "Convert back to Cartesian: " << back << std::endl;
    std::cout << "Roundtrip error: " << (cart_pos - back).NormL2() << std::endl;

    // Field operations in spherical coordinates
    // Inverse-square potential: φ = -1/r
    ScalarFunction<3> pot_spher([](const VectorN<Real, 3>& x) { return -1.0/x[0]; });
    
    // Use a point away from origin for gradient
    Vector3Spherical test_sph{2.0, Constants::PI/4, Constants::PI/4};
    auto grad_spher = ScalarFieldOperations::GradientSpher(pot_spher, test_sph);
    std::cout << std::endl << "Gradient of -1/r in spherical at r=2:" << std::endl;
    std::cout << "  ∂φ/∂r     = " << grad_spher[0] << " (analytical: 1/r² = 0.25)" << std::endl;
    std::cout << "  ∂φ/∂theta = " << grad_spher[1] << " (should be 0)" << std::endl;
    std::cout << "  ∂φ/∂phi   = " << grad_spher[2] << " (should be 0)" << std::endl;

    // Covariant vector transformation (transform gradient to Cartesian coords)
    Vector3Cartesian test_cart = CoordTransfSpherToCart.transf(test_sph);
    auto force_cart = CoordTransfSpherToCart.transfVecCovariant(grad_spher, test_cart);
    std::cout << std::endl << "Force vector in Cartesian: " << force_cart << std::endl;

    std::cout << std::endl;
}
