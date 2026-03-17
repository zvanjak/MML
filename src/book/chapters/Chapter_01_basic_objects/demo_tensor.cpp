///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_tensor.cpp
///  Description: Comprehensive Tensor demonstrations for Chapter 01
///               Index variance, contractions, and physics applications
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/Vector.h"
#include "mml/base/Vector/VectorN.h"
#include "mml/base/Tensor.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        Index Variance - Covariant vs Contravariant                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Tensor_IndexVariance()
{
    std::cout << "\n=== Index Variance ===\n\n";

    std::cout << "In tensor notation, indices can be:\n";
    std::cout << "  - Contravariant (upper): v^i - transforms inversely to basis\n";
    std::cout << "  - Covariant (lower):     v_i - transforms same as basis\n\n";

    // Different index configurations for rank-2 tensors
    std::cout << "--- Rank-2 Tensor Types ---\n";
    
    Tensor2<3> T_contra(0, 2);  // T^{ij} - both indices up
    std::cout << "T^{ij} (2 contravariant): NumCovar=" << T_contra.NumCovar() 
              << ", NumContravar=" << T_contra.NumContravar() << "\n";
    
    Tensor2<3> T_covar(2, 0);   // T_{ij} - both indices down
    std::cout << "T_{ij} (2 covariant):     NumCovar=" << T_covar.NumCovar() 
              << ", NumContravar=" << T_covar.NumContravar() << "\n";
    
    Tensor2<3> T_mixed(1, 1);   // T^i_j or T_i^j - one of each
    std::cout << "T^i_j (mixed):            NumCovar=" << T_mixed.NumCovar() 
              << ", NumContravar=" << T_mixed.NumContravar() << "\n";

    // Query individual index variance
    std::cout << "\nFor mixed tensor T^i_j:\n";
    std::cout << "  Is index 0 contravariant? " << (T_mixed.IsContravar(0) ? "yes" : "no") << "\n";
    std::cout << "  Is index 1 covariant?     " << (T_mixed.IsCovar(1) ? "yes" : "no") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Metric Tensors - Fundamental in Physics                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Metric_Tensors()
{
    std::cout << "\n=== Metric Tensors ===\n\n";

    std::cout << "The metric tensor g_{ij} defines distances and angles in a space.\n";
    std::cout << "ds² = g_{ij} dx^i dx^j\n\n";

    // 2D Euclidean (Cartesian) metric
    std::cout << "--- 2D Euclidean Metric (Cartesian) ---\n";
    Tensor2<2> g_euclidean(2, 0, {
        1, 0,
        0, 1
    });
    std::cout << "g_ij = I (identity):\n" << g_euclidean;
    std::cout << "This gives ds² = dx² + dy² (Pythagorean theorem)\n\n";

    // 2D Polar coordinates metric
    std::cout << "--- 2D Polar Coordinates Metric ---\n";
    double r = 2.0;  // At radius r = 2
    Tensor2<2> g_polar(2, 0, {
        1, 0,
        0, r*r
    });
    std::cout << "At r = " << r << ":\n";
    std::cout << "g_ij = diag(1, r²) =\n" << g_polar;
    std::cout << "This gives ds² = dr² + r²dθ²\n\n";

    // Minkowski metric (Special Relativity)
    std::cout << "--- Minkowski Metric (Special Relativity) ---\n";
    Tensor2<4> eta(2, 0, {
        -1,  0,  0,  0,
         0,  1,  0,  0,
         0,  0,  1,  0,
         0,  0,  0,  1
    });
    std::cout << "η_μν (signature -+++):\n" << eta;
    std::cout << "This gives ds² = -c²dt² + dx² + dy² + dz²\n";
    std::cout << "(Timelike intervals have ds² < 0)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Tensor Construction and Access                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Tensor_Construction()
{
    std::cout << "\n=== Tensor Construction ===\n\n";

    // Rank-2 tensor (3x3)
    std::cout << "--- Rank-2 Tensor (3×3) ---\n";
    Tensor2<3> A(1, 1, { 
        1.0, 2.0, 2.0,
        3.0, 4.0, 0.5,
        5.0, 6.0, 1.5 
    });
    std::cout << "A^i_j =\n" << A;

    // Element access
    std::cout << "Element access:\n";
    std::cout << "  A(0, 0) = " << A(0, 0) << "\n";
    std::cout << "  A(1, 1) = " << A(1, 1) << "\n";
    std::cout << "  A.at(2, 2) = " << A.at(2, 2) << " (bounds-checked)\n";

    // Modification
    A(0, 1) = 10.0;
    std::cout << "\nAfter A(0, 1) = 10:\n" << A;

    // Rank-3 tensor (3×3×3 = 27 elements)
    std::cout << "--- Rank-3 Tensor ---\n";
    Tensor3<3> B(2, 1);  // 2 covariant, 1 contravariant indices
    B(0, 0, 0) = 1.0;
    B(1, 1, 1) = 2.0;
    B(2, 2, 2) = 3.0;
    
    std::cout << "Rank-3 tensor B^i_{jk}:\n";
    std::cout << "  NumCovar = " << B.NumCovar() << ", NumContravar = " << B.NumContravar() << "\n";
    std::cout << "  B(0,0,0) = " << B(0, 0, 0) << "\n";
    std::cout << "  B(1,1,1) = " << B(1, 1, 1) << "\n";
    std::cout << "  B(2,2,2) = " << B(2, 2, 2) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Tensor Arithmetic                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Tensor_Arithmetic()
{
    std::cout << "\n=== Tensor Arithmetic ===\n\n";

    Tensor2<2> A(1, 1, {2, 0, 0, 3});
    Tensor2<2> B(1, 1, {1, 1, 1, 1});

    std::cout << "A =\n" << A;
    std::cout << "B =\n" << B;

    // Addition
    auto sum = A + B;
    std::cout << "A + B =\n" << sum;

    // Subtraction
    auto diff = A - B;
    std::cout << "A - B =\n" << diff;

    // Scalar multiplication
    auto scaled = A * 2.0;
    std::cout << "A * 2.0 =\n" << scaled;

    auto scaled2 = 0.5 * A;
    std::cout << "0.5 * A =\n" << scaled2;

    // Scalar division
    auto divided = A / 2.0;
    std::cout << "A / 2.0 =\n" << divided;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Contraction - The Core Tensor Operation                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Tensor_Contraction()
{
    std::cout << "\n=== Tensor Contraction ===\n\n";

    std::cout << "Contraction: summing over a repeated index (one up, one down)\n";
    std::cout << "Example: A^i_i = A^0_0 + A^1_1 + ... (trace)\n\n";

    // Rank-2 contraction (trace)
    std::cout << "--- Rank-2 Contraction (Trace) ---\n";
    Tensor2<3> A(1, 1, { 
        1.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 3.0 
    });
    std::cout << "A^i_j =\n" << A;
    
    Real trace = A.Contract();
    std::cout << "A^i_i = " << trace << " (1 + 2 + 3 = 6)\n\n";

    // Rank-3 contraction yields a vector
    std::cout << "--- Rank-3 Contraction (→ Vector) ---\n";
    Tensor3<3> T(1, 2);  // T^i_{jk}
    // Fill with some values
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                T(i, j, k) = static_cast<Real>(i + 3*j + 9*k + 1);

    VectorN<Real, 3> v1 = T.Contract(0, 1);  // Sum over first two indices
    std::cout << "T^i_{jk} contracted on indices (0,1): " << v1 << "\n";

    VectorN<Real, 3> v2 = T.Contract(0, 2);  // Sum over indices 0 and 2
    std::cout << "T^i_{jk} contracted on indices (0,2): " << v2 << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Tensor Acting on Vectors                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Tensor_On_Vectors()
{
    std::cout << "\n=== Tensors Acting on Vectors ===\n\n";

    std::cout << "A rank-2 tensor T^i_j acts on vectors:\n";
    std::cout << "  T(v, w) = T^i_j v_i w^j (contraction with both)\n\n";

    // Metric tensor example
    Tensor2<3> g(2, 0, {
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    });
    std::cout << "Euclidean metric g_ij:\n" << g;

    VectorN<Real, 3> v{1.0, 2.0, 3.0};
    VectorN<Real, 3> w{1.0, 0.0, 0.0};

    std::cout << "v = " << v << "\n";
    std::cout << "w = " << w << "\n";

    Real inner_product = g(v, w);
    std::cout << "\ng(v, w) = " << inner_product << " (v·w using metric)\n";

    // Self inner product = squared length
    Real v_squared = g(v, v);
    std::cout << "g(v, v) = " << v_squared << " = |v|² = 1² + 2² + 3² = 14\n";

    // 4D example: Minkowski metric
    std::cout << "\n--- Minkowski Metric ---\n";
    Tensor2<4> eta(2, 0, {
        -1,  0,  0,  0,
         0,  1,  0,  0,
         0,  0,  1,  0,
         0,  0,  0,  1
    });

    // A 4-velocity example (c = 1 units)
    VectorN<Real, 4> u{1.0, 0.0, 0.0, 0.0};  // At rest
    VectorN<Real, 4> p{1.0, 0.5, 0.0, 0.0};  // Moving in x

    std::cout << "4-vector u = " << u << " (at rest)\n";
    std::cout << "4-vector p = " << p << " (moving)\n";
    std::cout << "η(u, u) = " << eta(u, u) << " (should be -1 for proper time)\n";
    std::cout << "η(u, p) = " << eta(u, p) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Physics Applications                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Tensor_Physics()
{
    std::cout << "\n=== Physics Applications ===\n\n";

    // Christoffel symbols (connection coefficients)
    std::cout << "--- Christoffel Symbols (Rank-3) ---\n";
    std::cout << "Γ^i_{jk} - defines parallel transport and geodesics\n\n";
    
    Tensor3<3> Gamma(2, 1);  // 2 covariant, 1 contravariant
    // In curved space, these would be computed from the metric
    // For flat space in Cartesian coords, all are zero
    std::cout << "In Cartesian coordinates (flat space):\n";
    std::cout << "  All Christoffel symbols = 0\n";
    std::cout << "  Γ^0_{11} = " << Gamma(0, 1, 1) << "\n\n";

    // Electromagnetic field tensor
    std::cout << "--- Electromagnetic Field Tensor F^{μν} ---\n";
    // F^{μν} is antisymmetric, encodes E and B fields
    double Ex = 1.0, Ey = 0.0, Ez = 0.0;  // Electric field
    double Bx = 0.0, By = 0.0, Bz = 0.5;  // Magnetic field
    
    Tensor2<4> F(0, 2, {
         0,   Ex,  Ey,  Ez,
        -Ex,  0,   Bz, -By,
        -Ey, -Bz,  0,   Bx,
        -Ez,  By, -Bx,  0
    });
    
    std::cout << "For E = (" << Ex << ", " << Ey << ", " << Ez << "), ";
    std::cout << "B = (" << Bx << ", " << By << ", " << Bz << "):\n";
    std::cout << "F^{μν} =\n" << F;
    std::cout << "(Antisymmetric: F^{μν} = -F^{νμ})\n";

    // Stress-energy tensor
    std::cout << "\n--- Stress-Energy Tensor T^{μν} ---\n";
    std::cout << "Sources gravity in General Relativity:\n";
    std::cout << "  T^{00} = energy density\n";
    std::cout << "  T^{0i} = momentum density\n";
    std::cout << "  T^{ij} = stress (pressure, shear)\n";
    
    double rho = 1.0;  // Energy density
    double p = 0.5;    // Pressure (isotropic)
    Tensor2<4> T_perfect_fluid(0, 2, {
        rho, 0,  0,  0,
         0,  p,  0,  0,
         0,  0,  p,  0,
         0,  0,  0,  p
    });
    std::cout << "\nPerfect fluid (rest frame):\n" << T_perfect_fluid;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Tensor()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                      TENSORS IN MML                           ****\n";
    std::cout << "****           Index Variance, Contraction, and Physics            ****\n";
    std::cout << "***********************************************************************\n";

    Demo_Tensor_IndexVariance();
    Demo_Metric_Tensors();
    Demo_Tensor_Construction();
    Demo_Tensor_Arithmetic();
    Demo_Tensor_Contraction();
    Demo_Tensor_On_Vectors();
    Demo_Tensor_Physics();
}