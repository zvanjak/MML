#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_CLASSIC_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_CLASSIC_DEFS_H

// Classic test matrices from numerical linear algebra literature
// These are standard benchmarks used to test solver robustness and accuracy

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#endif

namespace MML::TestBeds
{
    /***********************************************************************************************/
    /**************                     Hilbert Matrices                              **************/
    /***********************************************************************************************/
    // Hilbert matrix: H[i,j] = 1/(i+j+1)
    // Source: David Hilbert (1894), "Über die Darstellung definiter Formen als Summe von Formenquadraten"
    // Properties: Symmetric, positive definite, notoriously ill-conditioned
    // Condition number grows as O(exp(3.5n))
    // Purpose: Classic test for ill-conditioned systems
    // Note: Exact inverse exists but numerical computation is challenging

    const static inline MML::Matrix<Real> hilbert_3x3{3, 3, { 
        1.0,        0.5,        1.0/3.0,
        0.5,        1.0/3.0,    0.25,
        1.0/3.0,    0.25,       0.2
    }};
    // Condition number ≈ 524.06
    
    const static inline MML::Vector<Real> hilbert_3x3_rhs0{1.0, 0.0, 0.0};
    const static inline MML::Vector<Real> hilbert_3x3_rhs0_sol{ 9.000000000000028422,
                                                             -36.000000000000142109,
                                                             30.000000000000131450 };

    const static inline VectorComplex hilbert_3x3_eigen_val{
        Complex(1.40831893063717, 0),
        Complex(0.122327059972063, 0),
        Complex(0.00268734272409464, 0)
    };

    // Matrix properties for hilbert_3x3
    const static inline Real              hilbert_3x3_det = 0.00046296296296296285;
    const static inline MML::Vector<Real> hilbert_3x3_singular_values{1.4083189306371614, 0.12232705997206326, 0.0026873427240946386};
    const static inline Real              hilbert_3x3_cond_2 = 524.0567775860644;
    const static inline int               hilbert_3x3_rank = 3;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> hilbert_4x4{4, 4, { 
        1.0,        0.5,        1.0/3.0,    0.25,
        0.5,        1.0/3.0,    0.25,       0.2,
        1.0/3.0,    0.25,       0.2,        1.0/6.0,
        0.25,       0.2,        1.0/6.0,    1.0/7.0
    }};
    // Condition number ≈ 15,513.74
    
    const static inline MML::Vector<Real> hilbert_4x4_rhs0{1.0, 1.0, 1.0, 1.0};
    const static inline MML::Vector<Real> hilbert_4x4_rhs0_sol{ -3.99999999999889511,
                                                             59.9999999999874873,
                                                             -179.999999999969873,
                                                             139.999999999980417 };

    const static inline VectorComplex hilbert_4x4_eigen_val{
        Complex(1.50021428498598, 0),
        Complex(0.169141194528054, 0),
        Complex(0.00673882617794846, 0),
        Complex(0.0000967080422693108, 0)
    };

    // Matrix properties for hilbert_4x4
    const static inline Real              hilbert_4x4_det = 1.653439153439153e-07;
    const static inline MML::Vector<Real> hilbert_4x4_singular_values{1.5002142849859763, 0.16914119452805447, 0.006738826177948463, 9.670804226931078e-05};
    const static inline Real              hilbert_4x4_cond_2 = 15513.738738929038;
    const static inline int               hilbert_4x4_rank = 4;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> hilbert_5x5{5, 5, { 
        1.0,        0.5,        1.0/3.0,    0.25,       0.2,
        0.5,        1.0/3.0,    0.25,       0.2,        1.0/6.0,
        1.0/3.0,    0.25,       0.2,        1.0/6.0,    1.0/7.0,
        0.25,       0.2,        1.0/6.0,    1.0/7.0,    0.125,
        0.2,        1.0/6.0,    1.0/7.0,    0.125,      1.0/9.0
    }};
    // Condition number ≈ 476,607.25
    
    const static inline MML::Vector<Real> hilbert_5x5_rhs0{1.0, 0.0, 0.0, 0.0, 0.0};
    const static inline MML::Vector<Real> hilbert_5x5_rhs0_sol{ 24.9999999999885745,
                                                             -299.999999999816851,
                                                             1049.99999999928991,
                                                             -1399.99999999900933,
                                                             629.999999999544343 };

    const static inline VectorComplex hilbert_5x5_eigen_val{
        Complex(1.56705104127683, 0),
        Complex(0.208534218611057, 0),
        Complex(0.0114075076537421, 0),
        Complex(0.000305898135857558, 0),
        Complex(3.28793161393089e-06, 0)
    };

    // Matrix properties for hilbert_5x5
    const static inline Real              hilbert_5x5_det = 3.7492951918427907e-12;
    const static inline MML::Vector<Real> hilbert_5x5_singular_values{1.5670510412768291, 0.20853421861105685, 0.011407507653742109, 0.00030589813585755764, 3.2879316139308927e-06};
    const static inline Real              hilbert_5x5_cond_2 = 476607.25024259434;
    const static inline int               hilbert_5x5_rank = 5;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> hilbert_8x8{8, 8, { 
        1.0,        0.5,        1.0/3.0,    0.25,       0.2,        1.0/6.0,    1.0/7.0,    0.125,
        0.5,        1.0/3.0,    0.25,       0.2,        1.0/6.0,    1.0/7.0,    0.125,      1.0/9.0,
        1.0/3.0,    0.25,       0.2,        1.0/6.0,    1.0/7.0,    0.125,      1.0/9.0,    0.1,
        0.25,       0.2,        1.0/6.0,    1.0/7.0,    0.125,      1.0/9.0,    0.1,        1.0/11.0,
        0.2,        1.0/6.0,    1.0/7.0,    0.125,      1.0/9.0,    0.1,        1.0/11.0,   1.0/12.0,
        1.0/6.0,    1.0/7.0,    0.125,      1.0/9.0,    0.1,        1.0/11.0,   1.0/12.0,   1.0/13.0,
        1.0/7.0,    0.125,      1.0/9.0,    0.1,        1.0/11.0,   1.0/12.0,   1.0/13.0,   1.0/14.0,
        0.125,      1.0/9.0,    0.1,        1.0/11.0,   1.0/12.0,   1.0/13.0,   1.0/14.0,   1.0/15.0
    }};
    // Condition number ≈ 1.53e10
    
    const static inline MML::Vector<Real> hilbert_8x8_rhs0{1.0, 0.5, 1.0/3.0, 0.25, 0.2, 1.0/6.0, 1.0/7.0, 0.125};
    const static inline MML::Vector<Real> hilbert_8x8_rhs0_sol{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    const static inline VectorComplex hilbert_8x8_eigen_val{
        Complex(1.61916364657347, 0),
        Complex(0.258198877582623, 0),
        Complex(0.0201263409990808, 0),
        Complex(0.000893442051546472, 0),
        Complex(0.0000232704524858034, 0),
        Complex(3.50923176262494e-07, 0),
        Complex(2.82226770346168e-09, 0),
        Complex(9.56466020563318e-12, 0)
    };

    // Matrix properties for hilbert_8x8
    const static inline Real              hilbert_8x8_det = 2.737050310006999e-33;
    const static inline MML::Vector<Real> hilbert_8x8_singular_values{1.6191636465734691, 0.2581988775826232, 0.020126340999080806, 0.000893442051546472, 2.327045248580339e-05, 3.509231762624942e-07, 2.822267703461677e-09, 9.564660205633179e-12};
    const static inline Real              hilbert_8x8_cond_2 = 169361458044.76382;
    const static inline int               hilbert_8x8_rank = 8;

    /***********************************************************************************************/
    /**************                     Pascal Matrices                               **************/
    /***********************************************************************************************/
    // Pascal matrix: P[i,j] = C(i+j, i) = binomial coefficient
    // Source: Blaise Pascal, triangle arithmetic (1653)
    // Properties: Symmetric, positive definite, well-conditioned (for small n)
    // Purpose: Tests exact integer arithmetic (all entries and solutions are integers)
    // Note: Can be used to verify solver accuracy - no rounding errors in exact arithmetic

    const static inline MML::Matrix<Real> pascal_3x3{3, 3, {
        1.0,    1.0,    1.0,
        1.0,    2.0,    3.0,
        1.0,    3.0,    6.0
    }};
    // Condition number ≈ 29.93
    
    const static inline MML::Vector<Real> pascal_3x3_rhs0{3.0, 6.0, 10.0};
    const static inline MML::Vector<Real> pascal_3x3_rhs0_sol{1.0, 1.0, 1.0};

    const static inline VectorComplex pascal_3x3_eigen_val{
        Complex(7.61803398874989, 0),
        Complex(1.0, 0),
        Complex(0.381966011250105, 0)
    };

    // Matrix properties for pascal_3x3
    const static inline Real              pascal_3x3_det = 1.0;
    const static inline MML::Vector<Real> pascal_3x3_singular_values{7.872762035082037, 0.8507709622346417, 0.1494315782430217};
    const static inline Real              pascal_3x3_cond_2 = 52.68486290776987;
    const static inline int               pascal_3x3_rank = 3;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> pascal_4x4{4, 4, {
        1.0,    1.0,    1.0,    1.0,
        1.0,    2.0,    3.0,    4.0,
        1.0,    3.0,    6.0,    10.0,
        1.0,    4.0,    10.0,   20.0
    }};
    // Condition number ≈ 296.42
    
    const static inline MML::Vector<Real> pascal_4x4_rhs0{4.0, 10.0, 20.0, 35.0};
    const static inline MML::Vector<Real> pascal_4x4_rhs0_sol{1.0, 1.0, 1.0, 1.0};

    const static inline VectorComplex pascal_4x4_eigen_val{
        Complex(24.5397807763306, 0),
        Complex(3.48698392376551, 0),
        Complex(0.973235297903924, 0),
        Complex(3.21154315992916e-14, 0)
    };

    // Matrix properties for pascal_4x4
    const static inline Real              pascal_4x4_det = 1.0;
    const static inline MML::Vector<Real> pascal_4x4_singular_values{26.089727103276627, 2.4668466068098166, 0.3831706927996587, 0.039877906206478744};
    const static inline Real              pascal_4x4_cond_2 = 654.285632428704;
    const static inline int               pascal_4x4_rank = 4;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> pascal_5x5{5, 5, {
        1.0,    1.0,    1.0,    1.0,    1.0,
        1.0,    2.0,    3.0,    4.0,    5.0,
        1.0,    3.0,    6.0,    10.0,   15.0,
        1.0,    4.0,    10.0,   20.0,   35.0,
        1.0,    5.0,    15.0,   35.0,   70.0
    }};
    // Condition number ≈ 2,934.15
    
    const static inline MML::Vector<Real> pascal_5x5_rhs0{5.0, 15.0, 35.0, 70.0, 126.0};
    const static inline MML::Vector<Real> pascal_5x5_rhs0_sol{1.0, 1.0, 1.0, 1.0, 1.0};

    const static inline VectorComplex pascal_5x5_eigen_val{
        Complex(88.5877825190949, 0),
        Complex(10.3176453932877, 0),
        Complex(0.927623858889912, 0),
        Complex(0.166948228727329, 0),
        Complex(6.14164823655194e-14, 0)
    };

    // Matrix properties for pascal_5x5
    const static inline Real              pascal_5x5_det = 1.0;
    const static inline MML::Vector<Real> pascal_5x5_singular_values{94.54195398399082, 7.256851312892555, 0.8048379619730833, 0.10569063399135155, 0.011289445168814684};
    const static inline Real              pascal_5x5_cond_2 = 8375.934668920896;
    const static inline int               pascal_5x5_rank = 5;

    /***********************************************************************************************/
    /**************                     Vandermonde Matrices                          **************/
    /***********************************************************************************************/
    // Vandermonde matrix: V[i,j] = x[i]^j where x is a set of nodes
    // Source: Alexandre-Théophile Vandermonde (1772)
    // Properties: Often ill-conditioned, related to polynomial interpolation
    // Purpose: Tests solver behavior for polynomial fitting and interpolation problems
    // Note: Condition number grows exponentially with n and depends on node spacing

    const static inline MML::Matrix<Real> vandermonde_3x3{3, 3, {
        1.0,    1.0,    1.0,
        1.0,    2.0,    4.0,
        1.0,    3.0,    9.0
    }};
    // Nodes: x = [1, 2, 3]
    // Condition number ≈ 13.93
    
    const static inline MML::Vector<Real> vandermonde_3x3_rhs0{1.0, 0.0, 0.0};
    const static inline MML::Vector<Real> vandermonde_3x3_rhs0_sol{ 3.0,
                                                             -2.5,
                                                             0.5 };

    const static inline VectorComplex vandermonde_3x3_eigen_val{
        Complex(10.3848189354378, 0),
        Complex(1.0, 0),
        Complex(0.615181064562221, 0)
    };

    // Matrix properties for vandermonde_3x3
    const static inline Real              vandermonde_3x3_det = 2.0;
    const static inline MML::Vector<Real> vandermonde_3x3_singular_values{10.638266291116246, 1.1312948855215682, 0.16617047629946406};
    const static inline Real              vandermonde_3x3_cond_2 = 64.02236697259854;
    const static inline int               vandermonde_3x3_rank = 3;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> vandermonde_4x4{4, 4, {
        1.0,    1.0,    1.0,     1.0,
        1.0,    2.0,    4.0,     8.0,
        1.0,    3.0,    9.0,     27.0,
        1.0,    4.0,    16.0,    64.0
    }};
    // Nodes: x = [1, 2, 3, 4]
    // Condition number ≈ 312.65
    
    const static inline MML::Vector<Real> vandermonde_4x4_rhs0{1.0, 2.0, 3.0, 4.0};
    const static inline MML::Vector<Real> vandermonde_4x4_rhs0_sol{ 0.0,
                                                             1.0,
                                                             0.0,
                                                             0.0 };

    const static inline VectorComplex vandermonde_4x4_eigen_val{
        Complex(67.6524329833645, 0),
        Complex(2.95912959527421, 0),
        Complex(0.388437421361318, 0),
        Complex(1.40166556476996e-13, 0)
    };

    // Matrix properties for vandermonde_4x4
    const static inline Real              vandermonde_4x4_det = 12.0;
    const static inline MML::Vector<Real> vandermonde_4x4_singular_values{69.5629001766858, 2.7626831649139155, 0.23018082437095295, 0.02262652879688948};
    const static inline Real              vandermonde_4x4_cond_2 = 3074.691795753135;
    const static inline int               vandermonde_4x4_rank = 4;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> vandermonde_5x5{5, 5, {
        1.0,    1.0,    1.0,     1.0,     1.0,
        1.0,    2.0,    4.0,     8.0,     16.0,
        1.0,    3.0,    9.0,     27.0,    81.0,
        1.0,    4.0,    16.0,    64.0,    256.0,
        1.0,    5.0,    25.0,    125.0,   625.0
    }};
    // Nodes: x = [1, 2, 3, 4, 5]
    // Condition number ≈ 14,026.03
    
    const static inline MML::Vector<Real> vandermonde_5x5_rhs0{5.0, 15.0, 35.0, 70.0, 126.0};
    const static inline MML::Vector<Real> vandermonde_5x5_rhs0_sol{ 1.0,
                                                             2.08333333333333304,
                                                             1.45833333333333348,
                                                             0.416666666666666685,
                                                             0.0416666666666666644 };

    const static inline VectorComplex vandermonde_5x5_eigen_val{
        Complex(650.010456765544, 0),
        Complex(6.06264621645773, 0),
        Complex(0.231473876551471, 0),
        Complex(0.0954241414532509, 0),
        Complex(2.48335869750195e-13, 0)
    };

    // Matrix properties for vandermonde_5x5
    const static inline Real              vandermonde_5x5_det = 288.0;
    const static inline MML::Vector<Real> vandermonde_5x5_singular_values{673.3959135063234, 6.068058403816063, 0.3062006587217439, 0.015871700050915376, 0.0010389398037693148};
    const static inline Real              vandermonde_5x5_cond_2 = 648162.2009285765;
    const static inline int               vandermonde_5x5_rank = 5;

    /***********************************************************************************************/
    /**************                     Frank Matrices                                **************/
    /***********************************************************************************************/
    // Frank matrix: Upper Hessenberg with F[i,j] = n+1-j for i <= j, F[i,j] = n+1-i for i = j+1, else 0
    // Source: Werner Frank (1958), "Computing eigenvalues of complex matrices"
    // Properties: Upper Hessenberg, all eigenvalues are real and distinct
    // Purpose: Tests eigenvalue computation, especially for eigenvalues of varying magnitude
    // Note: Excellent test for QR algorithm and related eigenvalue methods

    const static inline MML::Matrix<Real> frank_3x3{3, 3, {
        3.0,    2.0,    1.0,
        2.0,    2.0,    1.0,
        0.0,    1.0,    1.0
    }};
    // Condition number ≈ 36.44
    
    const static inline MML::Vector<Real> frank_3x3_rhs0{6.0, 5.0, 2.0};
    const static inline MML::Vector<Real> frank_3x3_rhs0_sol{1.0, 1.0, 1.0};

    const static inline VectorComplex frank_3x3_eigen_val{
        Complex(4.94948316672145, 0),
        Complex(1.10557280900008, 0),
        Complex(-0.0550560157215357, 0)
    };

    // Matrix properties for frank_3x3
    const static inline Real              frank_3x3_det = -0.30000000000000055;
    const static inline MML::Vector<Real> frank_3x3_singular_values{5.147935685095287, 1.3484320232019912, 0.04313437538952867};
    const static inline Real              frank_3x3_cond_2 = 119.35063706339016;
    const static inline int               frank_3x3_rank = 3;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> frank_4x4{4, 4, {
        4.0,    3.0,    2.0,    1.0,
        3.0,    3.0,    2.0,    1.0,
        0.0,    2.0,    2.0,    1.0,
        0.0,    0.0,    1.0,    1.0
    }};
    // Condition number ≈ 278.24
    
    const static inline MML::Vector<Real> frank_4x4_rhs0{10.0, 9.0, 5.0, 2.0};
    const static inline MML::Vector<Real> frank_4x4_rhs0_sol{1.0, 1.0, 1.0, 1.0};

    const static inline VectorComplex frank_4x4_eigen_val{
        Complex(7.65828103570816, 0),
        Complex(2.16331577368876, 0),
        Complex(0.280113127328669, 0),
        Complex(-0.101710036725597, 0)
    };

    // Matrix properties for frank_4x4
    const static inline Real              frank_4x4_det = -0.049999999999999954;
    const static inline MML::Vector<Real> frank_4x4_singular_values{8.003098817046832, 2.0770306888073584, 0.20261050821755314, 0.029854389282313704};
    const static inline Real              frank_4x4_cond_2 = 268.05660277538714;
    const static inline int               frank_4x4_rank = 4;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> frank_5x5{5, 5, {
        5.0,    4.0,    3.0,    2.0,    1.0,
        4.0,    4.0,    3.0,    2.0,    1.0,
        0.0,    3.0,    3.0,    2.0,    1.0,
        0.0,    0.0,    2.0,    2.0,    1.0,
        0.0,    0.0,    0.0,    1.0,    1.0
    }};
    // Condition number ≈ 2,178.91
    
    const static inline MML::Vector<Real> frank_5x5_rhs0{15.0, 14.0, 9.0, 5.0, 2.0};
    const static inline MML::Vector<Real> frank_5x5_rhs0_sol{1.0, 1.0, 1.0, 1.0, 1.0};

    const static inline VectorComplex frank_5x5_eigen_val{
        Complex(11.4949816267438, 0),
        Complex(3.12758815478237, 0),
        Complex(0.486387764798014, 0),
        Complex(0.0415030413691732, 0),
        Complex(-0.149559386692964, 0)
    };

    // Matrix properties for frank_5x5
    const static inline Real              frank_5x5_det = -0.008333333333333327;
    const static inline MML::Vector<Real> frank_5x5_singular_values{11.760929127988423, 2.8177426117551053, 0.28424303706296697, 0.024867896379831584, 0.003005538820706936};
    const static inline Real              frank_5x5_cond_2 = 3912.5971096430165;
    const static inline int               frank_5x5_rank = 5;

    /***********************************************************************************************/
    /**************           Kahan Matrices (severely ill-conditioned)               **************/
    /***********************************************************************************************/
    // Kahan matrix: Tests numerical stability for triangular systems
    // Source: William Kahan (numerical analysis pioneer)
    // Properties: Upper triangular, parametrized by angle θ
    // Purpose: Tests backward stability of triangular solvers
    // Note: Condition number can be made arbitrarily large by choice of θ

    const static inline MML::Matrix<Real> kahan_3x3{3, 3, {
        1.0,        -0.5,       -0.5,
        0.0,        0.866025,   -0.433013,
        0.0,        0.0,        0.75
    }};
    // θ = π/6, s = sin(θ) = 0.5, c = cos(θ) ≈ 0.866025
    // Condition number ≈ 4.62
    
    const static inline MML::Vector<Real> kahan_3x3_rhs0{0.0, 0.433013, 0.75};
    const static inline MML::Vector<Real> kahan_3x3_rhs0_sol{ 1.000000577350538489,
                                                             1.000001154701076755,
                                                             1.000000000000000000 };

    const static inline VectorComplex kahan_3x3_eigen_val{
        Complex(1.0, 0),
        Complex(0.866025403784439, 0),
        Complex(0.75, 0)
    };

    // Matrix properties for kahan_3x3
    const static inline Real              kahan_3x3_det = 0.6495187500000001;
    const static inline MML::Vector<Real> kahan_3x3_singular_values{1.3108974376855447, 0.8770399269579917, 0.5652434916251133};
    const static inline Real              kahan_3x3_cond_2 = 2.3191997618579485;
    const static inline int               kahan_3x3_rank = 3;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Real> kahan_5x5{5, 5, {
        1.0,        -0.5,       -0.5,       -0.5,       -0.5,
        0.0,        0.866025,   -0.433013,  -0.433013,  -0.433013,
        0.0,        0.0,        0.75,       -0.375,     -0.375,
        0.0,        0.0,        0.0,        0.649519,   -0.32476,
        0.0,        0.0,        0.0,        0.0,        0.5625
    }};
    // θ = π/6 for all rows
    // Condition number ≈ 59.94
    
    const static inline MML::Vector<Real> kahan_5x5_rhs0{-2.0, -0.866025, -0.75, -0.649519, 0.5625};
    const static inline MML::Vector<Real> kahan_5x5_rhs0_sol{ -2.687499206143009634,
                                                             -1.124999566986651933,
                                                             -0.749999615099789185,
                                                             -0.499999230199578426,
                                                             1.000000000000000000 };

    const static inline VectorComplex kahan_5x5_eigen_val{
        Complex(1.0, 0),
        Complex(0.866025403784439, 0),
        Complex(0.75, 0),
        Complex(0.649519052838329, 0),
        Complex(0.5625, 0)
    };

    // Matrix properties for kahan_5x5
    const static inline Real              kahan_5x5_det = 0.2371621878051758;
    const static inline MML::Vector<Real> kahan_5x5_singular_values{1.5889817050568896, 0.8773413135200024, 0.6303116648310096, 0.5005118188605684, 0.3367626116437704};
    const static inline Real              kahan_5x5_cond_2 = 4.718811649912687;
    const static inline int               kahan_5x5_rank = 5;

} // namespace MML::TestBeds

#endif // MML_LINEAR_ALG_EQ_SYSTEMS_CLASSIC_DEFS_H
