#if !defined __MML_LINEAR_ALG_EQ_SOLVERS_TEST_BED_H
#define __MML_LINEAR_ALG_EQ_SOLVERS_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#endif

#include <vector>
#include <functional>

#include "linear_alg_eq_systems_real_defs.h"
#include "linear_alg_eq_systems_sym_defs.h"
#include "linear_alg_eq_systems_complex_defs.h"
#include "linear_alg_eq_systems_classic_defs.h"
#include "linear_alg_eq_systems_solverspecific_defs.h"

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * MATRIX CATEGORIES - for filtering and classification
     *******************************************************************************************************************/
    enum class MatrixCategory {
        Random,           // Randomly generated matrices (mat_3x3, mat_5x5, etc.)
        Hilbert,          // Hilbert matrices (ill-conditioned)
        Pascal,           // Pascal matrices (exact arithmetic, integer entries)
        Vandermonde,      // Vandermonde matrices (polynomial interpolation)
        Frank,            // Frank matrices (eigenvalue computation)
        Kahan,            // Kahan matrices (backward stability testing)
        DiagDominant,     // Diagonally dominant matrices (iterative solver friendly)
        Overdetermined,   // m x n matrices where m > n (least squares)
        SPD,              // Symmetric positive definite
        Symmetric         // General symmetric matrices
    };

    /*******************************************************************************************************************
     * TEST SYSTEM STRUCTS
     *******************************************************************************************************************/
    class TestLinearSystem
    {
    public:
        int _n;
        Matrix<Real> _mat;
        Vector<Real> _rhs;
        Vector<Real> _sol;
        Vector<Complex> _eigen_values;
        std::vector<Vector<Complex>> _eigen_vectors;
        
        // Matrix properties for algorithm verification
        Real _determinant = 0.0;           // det(A)
        Vector<Real> _singular_values;     // σ₁ ≥ σ₂ ≥ ... ≥ σₙ (descending)
        Real _cond_2 = 0.0;                // Condition number (2-norm): σ_max/σ_min
        int _rank = 0;                     // Numerical rank

        // Original constructors (backward compatible)
        TestLinearSystem(int n, const MML::Matrix<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol, const MML::Vector<Complex> &eigen_values) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values), _rank(n)
        {}
        TestLinearSystem(int n, const MML::Matrix<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol, const MML::Vector<Complex> &eigen_values, const std::vector<MML::Vector<Complex>> &eigen_vectors) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values), _eigen_vectors(eigen_vectors), _rank(n)
        {}
        
        // Extended constructor with matrix properties
        TestLinearSystem(int n, const MML::Matrix<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol, 
                        const MML::Vector<Complex> &eigen_values, Real determinant, const MML::Vector<Real> &singular_values, Real cond_2, int rank) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values),
            _determinant(determinant), _singular_values(singular_values), _cond_2(cond_2), _rank(rank)
        {}
        TestLinearSystem(int n, const MML::Matrix<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol, 
                        const MML::Vector<Complex> &eigen_values, const std::vector<MML::Vector<Complex>> &eigen_vectors,
                        Real determinant, const MML::Vector<Real> &singular_values, Real cond_2, int rank) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values), _eigen_vectors(eigen_vectors),
            _determinant(determinant), _singular_values(singular_values), _cond_2(cond_2), _rank(rank)
        {}
    };

    class TestLinearSystemMultiRHS
    {
    public:
        int _n;
        Matrix<Real> _mat;
        Matrix<Real> _rhs;
        Matrix<Real> _sol;

        TestLinearSystemMultiRHS(int n, const MML::Matrix<Real> &mat, const MML::Matrix<Real> &rhs, const MML::Matrix<Real> &sol) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol)
        {}
    };    

    class TestLinearSystemSymmetric
    {
    public:
        int _n;
        MatrixSym<Real> _mat;
        Vector<Real> _rhs;
        Vector<Real> _sol;
        Vector<Real> _eigen_values;
        std::vector<Vector<Real>> _eigen_vectors;
        
        // Matrix properties for algorithm verification
        Real _determinant = 0.0;           // det(A)
        Vector<Real> _singular_values;     // For symmetric: |λ₁| ≥ |λ₂| ≥ ... (eigenvalues sorted by magnitude)
        Real _cond_2 = 0.0;                // Condition number (2-norm): |λ_max|/|λ_min|
        int _rank = 0;                     // Numerical rank

        // Original constructors (backward compatible)
        TestLinearSystemSymmetric(int n, const MML::MatrixSym<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol, const MML::Vector<Real> &eigen_values) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values), _rank(n)
        {}
        TestLinearSystemSymmetric(int n, const MML::MatrixSym<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol, const MML::Vector<Real> &eigen_values, const std::vector<MML::Vector<Real>> &eigen_vectors) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values), _eigen_vectors(eigen_vectors), _rank(n)
        {}
        
        // Extended constructor with matrix properties
        TestLinearSystemSymmetric(int n, const MML::MatrixSym<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol, 
                                  const MML::Vector<Real> &eigen_values, Real determinant, const MML::Vector<Real> &singular_values, Real cond_2, int rank) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values),
            _determinant(determinant), _singular_values(singular_values), _cond_2(cond_2), _rank(rank)
        {}
    };

    // Complex linear system
    class TestLinearSystemComplex
    {
    public:
        int _n;
        int _m;  // For overdetermined: m rows, n columns (m >= n)
        Matrix<Complex> _mat;
        Vector<Complex> _rhs;
        Vector<Complex> _sol;
        
        // Matrix properties
        Real _det_abs = 0.0;               // |det(A)| - magnitude of complex determinant
        Vector<Real> _singular_values;     // σ₁ ≥ σ₂ ≥ ... ≥ σₙ (descending)
        Real _cond_2 = 0.0;                // Condition number (2-norm): σ_max/σ_min
        int _rank = 0;                     // Numerical rank

        TestLinearSystemComplex(int n, const MML::Matrix<Complex> &mat, const MML::Vector<Complex> &rhs, const MML::Vector<Complex> &sol) : 
            _n(n), _m(n), _mat(mat), _rhs(rhs), _sol(sol), _rank(n)
        {}
        
        TestLinearSystemComplex(int n, const MML::Matrix<Complex> &mat, const MML::Vector<Complex> &rhs, const MML::Vector<Complex> &sol,
                               Real det_abs, const MML::Vector<Real> &singular_values, Real cond_2, int rank) : 
            _n(n), _m(n), _mat(mat), _rhs(rhs), _sol(sol),
            _det_abs(det_abs), _singular_values(singular_values), _cond_2(cond_2), _rank(rank)
        {}
    };

    // Overdetermined system (m x n, m > n) - for least squares
    class TestLinearSystemOverdetermined
    {
    public:
        int _m;  // rows
        int _n;  // columns (m > n)
        Matrix<Real> _mat;
        Vector<Real> _rhs;
        Vector<Real> _sol;  // Least squares solution
        
        // Matrix properties
        Vector<Real> _singular_values;     // σ₁ ≥ σ₂ ≥ ... ≥ σₙ (descending)
        Real _cond_2 = 0.0;                // Condition number (2-norm)
        int _rank = 0;                     // Numerical rank

        TestLinearSystemOverdetermined(int m, int n, const MML::Matrix<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol) : 
            _m(m), _n(n), _mat(mat), _rhs(rhs), _sol(sol), _rank(n)
        {}
        
        TestLinearSystemOverdetermined(int m, int n, const MML::Matrix<Real> &mat, const MML::Vector<Real> &rhs, const MML::Vector<Real> &sol,
                                       const MML::Vector<Real> &singular_values, Real cond_2, int rank) : 
            _m(m), _n(n), _mat(mat), _rhs(rhs), _sol(sol),
            _singular_values(singular_values), _cond_2(cond_2), _rank(rank)
        {}
    };

    /*******************************************************************************************************************
     * HELPER FUNCTIONS FOR SOLVER-SPECIFIC MATRICES (kept for backward compatibility)
     *******************************************************************************************************************/

    inline TestLinearSystem overdetermined_4x3() {
        return TestLinearSystem(3, overdetermined_4x3_mat, overdetermined_4x3_rhs, overdetermined_4x3_sol, overdetermined_4x3_eigen);
    }

    inline TestLinearSystem overdetermined_5x3_regression() {
        return TestLinearSystem(3, overdetermined_5x3_regression_mat, overdetermined_5x3_regression_rhs, overdetermined_5x3_regression_sol, overdetermined_5x3_regression_eigen);
    }

    inline TestLinearSystem overdetermined_6x4_polynomial() {
        return TestLinearSystem(4, overdetermined_6x4_polynomial_mat, overdetermined_6x4_polynomial_rhs, overdetermined_6x4_polynomial_sol, overdetermined_6x4_polynomial_eigen);
    }

    inline TestLinearSystem diag_dominant_4x4() {
        return TestLinearSystem(4, diag_dominant_4x4_mat, diag_dominant_4x4_rhs, diag_dominant_4x4_sol, diag_dominant_4x4_eigen);
    }

    inline TestLinearSystem diag_dominant_5x5_tridiag() {
        return TestLinearSystem(5, diag_dominant_5x5_tridiag_mat, diag_dominant_5x5_tridiag_rhs, diag_dominant_5x5_tridiag_sol, diag_dominant_5x5_tridiag_eigen);
    }

    inline TestLinearSystem diag_dominant_6x6_poisson2d() {
        return TestLinearSystem(6, diag_dominant_6x6_poisson2d_mat, diag_dominant_6x6_poisson2d_rhs, diag_dominant_6x6_poisson2d_sol, diag_dominant_6x6_poisson2d_eigen);
    }

    inline TestLinearSystemSymmetric spd_3x3() {
        return TestLinearSystemSymmetric(3, spd_3x3_mat, spd_3x3_rhs, spd_3x3_sol, spd_3x3_eigen);
    }

    inline TestLinearSystemSymmetric spd_4x4_correlation() {
        return TestLinearSystemSymmetric(4, spd_4x4_correlation_mat, spd_4x4_correlation_rhs, spd_4x4_correlation_sol, spd_4x4_correlation_eigen);
    }

    inline TestLinearSystemSymmetric spd_5x5_mass_matrix() {
        return TestLinearSystemSymmetric(5, spd_5x5_mass_matrix_mat, spd_5x5_mass_matrix_rhs, spd_5x5_mass_matrix_sol, spd_5x5_mass_matrix_eigen);
    }

    inline TestLinearSystemSymmetric spd_6x6_graph_laplacian() {
        return TestLinearSystemSymmetric(6, spd_6x6_graph_laplacian_mat, spd_6x6_graph_laplacian_rhs, spd_6x6_graph_laplacian_sol, spd_6x6_graph_laplacian_eigen);
    }

    /*******************************************************************************************************************
     * LINEAR ALGEBRA EQUATION TEST BED
     * 
     * Provides access to test matrices organized by:
     * - Type: Real, Symmetric, Complex, Overdetermined
     * - Category: Random, Hilbert, Pascal, Vandermonde, Frank, Kahan, DiagDominant, SPD
     * - Size: getBySize(n), getSmall(), getLarge()
     * - Conditioning: getWellConditioned(), getIllConditioned()
     *******************************************************************************************************************/
    class LinearAlgEqTestBed
    {
    private:
        // Helper to initialize linear systems array - uses construct-on-first-use idiom
        // to avoid static initialization order fiasco
        static const std::vector<std::pair<std::string, TestLinearSystem>>& getLinearSystemsArray() {
            static const std::vector<std::pair<std::string, TestLinearSystem>> systems = {
                { "mat_3x3", { 3, mat_3x3, mat_3x3_rhs0, mat_3x3_rhs0_sol, mat_3x3_eigen_val, mat_3x3_eigen_vecs } },
                { "mat_3x3_1", { 3, mat_3x3_1, mat_3x3_1_rhs0, mat_3x3_1_rhs0_sol, mat_3x3_1_eigen_val, mat_3x3_1_eigen_vecs } },
                { "mat_3x3_2", { 3, mat_3x3_2, mat_3x3_2_rhs0, mat_3x3_2_rhs0_sol, mat_3x3_2_eigen_val, mat_3x3_2_eigen_vecs } },
                { "mat_3x3_3", { 3, mat_3x3_3, mat_3x3_3_rhs0, mat_3x3_3_rhs0_sol, mat_3x3_3_eigen_val, mat_3x3_3_eigen_vecs } },
                { "mat_3x3_4", { 3, mat_3x3_4, mat_3x3_4_rhs0, mat_3x3_4_rhs0_sol, mat_3x3_4_eigen_val, mat_3x3_4_eigen_vecs } },
                { "mat_5x5", { 5, mat_5x5, mat_5x5_rhs0, mat_5x5_rhs0_sol, mat_5x5_eigen_val, mat_5x5_eigen_vecs } },
                { "mat_8x8", { 8, mat_8x8, mat_8x8_rhs0, mat_8x8_rhs0_sol, mat_8x8_eigen_val, mat_8x8_eigen_vecs } },
                { "mat_10x10", { 10, mat_10x10, mat_10x10_rhs0, mat_10x10_rhs0_sol, mat_10x10_eigen_val, mat_10x10_eigen_vecs } },
                { "mat_20x20", { 20, mat_20x20, mat_20x20_rhs0, mat_20x20_rhs0_sol, mat_20x20_eigen_val, mat_20x20_eigen_vecs } },
                { "hilbert_3x3", { 3, hilbert_3x3, hilbert_3x3_rhs0, hilbert_3x3_rhs0_sol, hilbert_3x3_eigen_val } },
                { "hilbert_4x4", { 4, hilbert_4x4, hilbert_4x4_rhs0, hilbert_4x4_rhs0_sol, hilbert_4x4_eigen_val } },
                { "hilbert_5x5", { 5, hilbert_5x5, hilbert_5x5_rhs0, hilbert_5x5_rhs0_sol, hilbert_5x5_eigen_val } },
                { "hilbert_8x8", { 8, hilbert_8x8, hilbert_8x8_rhs0, hilbert_8x8_rhs0_sol, hilbert_8x8_eigen_val } },
                { "pascal_3x3", { 3, pascal_3x3, pascal_3x3_rhs0, pascal_3x3_rhs0_sol, pascal_3x3_eigen_val } },
                { "pascal_4x4", { 4, pascal_4x4, pascal_4x4_rhs0, pascal_4x4_rhs0_sol, pascal_4x4_eigen_val } },
                { "pascal_5x5", { 5, pascal_5x5, pascal_5x5_rhs0, pascal_5x5_rhs0_sol, pascal_5x5_eigen_val } },
                { "vandermonde_3x3", { 3, vandermonde_3x3, vandermonde_3x3_rhs0, vandermonde_3x3_rhs0_sol, vandermonde_3x3_eigen_val } },
                { "vandermonde_4x4", { 4, vandermonde_4x4, vandermonde_4x4_rhs0, vandermonde_4x4_rhs0_sol, vandermonde_4x4_eigen_val } },
                { "vandermonde_5x5", { 5, vandermonde_5x5, vandermonde_5x5_rhs0, vandermonde_5x5_rhs0_sol, vandermonde_5x5_eigen_val } },
                { "frank_3x3", { 3, frank_3x3, frank_3x3_rhs0, frank_3x3_rhs0_sol, frank_3x3_eigen_val } },
                { "frank_4x4", { 4, frank_4x4, frank_4x4_rhs0, frank_4x4_rhs0_sol, frank_4x4_eigen_val } },
                { "frank_5x5", { 5, frank_5x5, frank_5x5_rhs0, frank_5x5_rhs0_sol, frank_5x5_eigen_val } },
                { "kahan_3x3", { 3, kahan_3x3, kahan_3x3_rhs0, kahan_3x3_rhs0_sol, kahan_3x3_eigen_val } },
                { "kahan_5x5", { 5, kahan_5x5, kahan_5x5_rhs0, kahan_5x5_rhs0_sol, kahan_5x5_eigen_val } },
                { "diag_dominant_4x4", { 4, diag_dominant_4x4_mat, diag_dominant_4x4_rhs, diag_dominant_4x4_sol, diag_dominant_4x4_eigen } },
                { "diag_dominant_5x5_tridiag", { 5, diag_dominant_5x5_tridiag_mat, diag_dominant_5x5_tridiag_rhs, diag_dominant_5x5_tridiag_sol, diag_dominant_5x5_tridiag_eigen } },
                { "diag_dominant_6x6_poisson2d", { 6, diag_dominant_6x6_poisson2d_mat, diag_dominant_6x6_poisson2d_rhs, diag_dominant_6x6_poisson2d_sol, diag_dominant_6x6_poisson2d_eigen } }
            };
            return systems;
        }

        static const std::vector<std::pair<std::string, TestLinearSystemSymmetric>>& getSymmetricSystemsArray() {
            static const std::vector<std::pair<std::string, TestLinearSystemSymmetric>> systems = {
                { "symm_mat_3x3", { 3, symm_mat_3x3, symm_mat_3x3_rhs0, symm_mat_3x3_rhs0_sol, symm_mat_3x3_eigen_val, symm_mat_3x3_eigen_vecs } },
                { "symm_mat_5x5", { 5, symm_mat_5x5, symm_mat_5x5_rhs0, symm_mat_5x5_rhs0_sol, symm_mat_5x5_eigen_val, symm_mat_5x5_eigen_vecs } },
                { "symm_mat_10x10", { 10, symm_mat_10x10, symm_mat_10x10_rhs0, symm_mat_10x10_rhs0_sol, symm_mat_10x10_eigen_val, symm_mat_10x10_eigen_vecs } },
                { "spd_3x3", { 3, spd_3x3_mat, spd_3x3_rhs, spd_3x3_sol, spd_3x3_eigen } },
                { "spd_4x4_correlation", { 4, spd_4x4_correlation_mat, spd_4x4_correlation_rhs, spd_4x4_correlation_sol, spd_4x4_correlation_eigen } },
                { "spd_5x5_mass_matrix", { 5, spd_5x5_mass_matrix_mat, spd_5x5_mass_matrix_rhs, spd_5x5_mass_matrix_sol, spd_5x5_mass_matrix_eigen } },
                { "spd_6x6_graph_laplacian", { 6, spd_6x6_graph_laplacian_mat, spd_6x6_graph_laplacian_rhs, spd_6x6_graph_laplacian_sol, spd_6x6_graph_laplacian_eigen } }
            };
            return systems;
        }

        static const std::vector<std::pair<std::string, TestLinearSystemMultiRHS>>& getMultiRHSSystemsArray() {
            static const std::vector<std::pair<std::string, TestLinearSystemMultiRHS>> systems = {
                {
                    "mat_5x5_multi_rhs1",
                    {
                        5, 
                        Matrix<Real>{5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                        1.6, 1.5, 1.1, 0.7, 5.0,
                                        3.8, 8.0, 9.6, 5.4, 8.8,
                                        4.6, 8.2, 8.4, 0.4, 8.0,
                                        2.6, 2.9, 0.1, 9.6, 7.7}},
                        Matrix<Real>{5, 2, {1.1, 1.6, 
                                        4.7, 9.1, 
                                        0.1, 4.0, 
                                        9.3, 8.4, 
                                        0.4, 4.1}},
                        Matrix<Real>{5, 2, {-3.9032710424808688,  15.643114174796667, 
                                        5.2353433160849479, -11.587503332831671, 
                                        -3.2920957702478550,   4.4111268480786325, 
                                        -1.7183300108528281,   0.21432757972725644, 
                                        1.5832710097423177,  -0.70999930382454612}}
                    }
                }
            };
            return systems;
        }

        static const std::vector<std::pair<std::string, TestLinearSystemComplex>>& getComplexSystemsArray() {
            static const std::vector<std::pair<std::string, TestLinearSystemComplex>> systems = {
                { "mat_cmplx_3x3", { 3, mat_cmplx_3x3, mat_cmplx_3x3_rhs0, mat_cmplx_3x3_rhs0_sol, mat_cmplx_3x3_det_abs, mat_cmplx_3x3_singular_values, mat_cmplx_3x3_cond_2, mat_cmplx_3x3_rank } },
                { "mat_cmplx_1_3x3", { 3, mat_cmplx_1_3x3, mat_cmplx_1_3x3_rhs0, mat_cmplx_1_3x3_rhs0_sol, mat_cmplx_1_3x3_det_abs, mat_cmplx_1_3x3_singular_values, mat_cmplx_1_3x3_cond_2, mat_cmplx_1_3x3_rank } },
                { "mat_cmplx_1_5x5", { 5, mat_cmplx_1_5x5, mat_cmplx_1_5x5_rhs0, mat_cmplx_1_5x5_rhs0_sol, mat_cmplx_1_5x5_det_abs, mat_cmplx_1_5x5_singular_values, mat_cmplx_1_5x5_cond_2, mat_cmplx_1_5x5_rank } },
                { "mat_cmplx_2_5x5", { 5, mat_cmplx_2_5x5, mat_cmplx_2_5x5_rhs0, mat_cmplx_2_5x5_rhs0_sol, mat_cmplx_2_5x5_det_abs, mat_cmplx_2_5x5_singular_values, mat_cmplx_2_5x5_cond_2, mat_cmplx_2_5x5_rank } },
                { "mat_cmplx_3_5x5", { 5, mat_cmplx_3_5x5, mat_cmplx_3_5x5_rhs0, mat_cmplx_3_5x5_rhs0_sol, mat_cmplx_3_5x5_det_abs, mat_cmplx_3_5x5_singular_values, mat_cmplx_3_5x5_cond_2, mat_cmplx_3_5x5_rank } },
                { "mat_cmplx_4_5x5", { 5, mat_cmplx_4_5x5, mat_cmplx_4_5x5_rhs0, mat_cmplx_4_5x5_rhs0_sol, mat_cmplx_4_5x5_det_abs, mat_cmplx_4_5x5_singular_values, mat_cmplx_4_5x5_cond_2, mat_cmplx_4_5x5_rank } },
                { "mat_cmplx_1_8x8", { 8, mat_cmplx_1_8x8, mat_cmplx_1_8x8_rhs0, mat_cmplx_1_8x8_rhs0_sol, mat_cmplx_1_8x8_det_abs, mat_cmplx_1_8x8_singular_values, mat_cmplx_1_8x8_cond_2, mat_cmplx_1_8x8_rank } },
                { "mat_cmplx_2_8x8", { 8, mat_cmplx_2_8x8, mat_cmplx_2_8x8_rhs0, mat_cmplx_2_8x8_rhs0_sol, mat_cmplx_2_8x8_det_abs, mat_cmplx_2_8x8_singular_values, mat_cmplx_2_8x8_cond_2, mat_cmplx_2_8x8_rank } },
                { "mat_cmplx_3_8x8", { 8, mat_cmplx_3_8x8, mat_cmplx_3_8x8_rhs0, mat_cmplx_3_8x8_rhs0_sol, mat_cmplx_3_8x8_det_abs, mat_cmplx_3_8x8_singular_values, mat_cmplx_3_8x8_cond_2, mat_cmplx_3_8x8_rank } },
                { "mat_cmplx_1_10x10", { 10, mat_cmplx_1_10x10, mat_cmplx_1_10x10_rhs0, mat_cmplx_1_10x10_rhs0_sol, mat_cmplx_1_10x10_det_abs, mat_cmplx_1_10x10_singular_values, mat_cmplx_1_10x10_cond_2, mat_cmplx_1_10x10_rank } },
                { "mat_cmplx_2_10x10", { 10, mat_cmplx_2_10x10, mat_cmplx_2_10x10_rhs0, mat_cmplx_2_10x10_rhs0_sol, mat_cmplx_2_10x10_det_abs, mat_cmplx_2_10x10_singular_values, mat_cmplx_2_10x10_cond_2, mat_cmplx_2_10x10_rank } },
                { "mat_cmplx_3_10x10", { 10, mat_cmplx_3_10x10, mat_cmplx_3_10x10_rhs0, mat_cmplx_3_10x10_rhs0_sol, mat_cmplx_3_10x10_det_abs, mat_cmplx_3_10x10_singular_values, mat_cmplx_3_10x10_cond_2, mat_cmplx_3_10x10_rank } }
            };
            return systems;
        }

        static const std::vector<std::pair<std::string, TestLinearSystemOverdetermined>>& getOverdeterminedSystemsArray() {
            static const std::vector<std::pair<std::string, TestLinearSystemOverdetermined>> systems = {
                { "overdetermined_4x3", { 4, 3, overdetermined_4x3_mat, overdetermined_4x3_rhs, overdetermined_4x3_sol, overdetermined_4x3_singular_values, overdetermined_4x3_cond_2, overdetermined_4x3_rank } },
                { "overdetermined_5x3_regression", { 5, 3, overdetermined_5x3_regression_mat, overdetermined_5x3_regression_rhs, overdetermined_5x3_regression_sol, overdetermined_5x3_regression_singular_values, overdetermined_5x3_regression_cond_2, overdetermined_5x3_regression_rank } },
                { "overdetermined_6x4_polynomial", { 6, 4, overdetermined_6x4_polynomial_mat, overdetermined_6x4_polynomial_rhs, overdetermined_6x4_polynomial_sol, overdetermined_6x4_polynomial_singular_values, overdetermined_6x4_polynomial_cond_2, overdetermined_6x4_polynomial_rank } }
            };
            return systems;
        }

    public:
        //==============================================================================================================
        // ORIGINAL INTERFACE (backward compatible)
        //==============================================================================================================
        static int numLinAlgEqSystems()     { return static_cast<int>(getLinearSystemsArray().size()); }
        
        static const TestLinearSystem& getLinAlgEqSystem(int index) { return getLinearSystemsArray()[index].second; }
        static const TestLinearSystem& getLinAlgEqSystem(std::string sysName)
        {
            for (const auto &sys : getLinearSystemsArray())
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("LinearAlgEqTestBed::getLinAlgEqSystem: system with name " + sysName + " not found");
        }
        
        static int numLinAlgEqSystemsSymmetric() { return static_cast<int>(getSymmetricSystemsArray().size()); }
        static const TestLinearSystemSymmetric& getLinAlgEqSystemSymmetric(int index) { return getSymmetricSystemsArray()[index].second; }
        static const TestLinearSystemSymmetric& getLinAlgEqSystemSymmetric(std::string sysName)
        {
            for (const auto &sys : getSymmetricSystemsArray())
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("LinearAlgEqTestBed::getLinAlgEqSystemSymmetric: system with name " + sysName + " not found");
        }
        
        static int numLinAlgEqSystemsMultiRHS() { return static_cast<int>(getMultiRHSSystemsArray().size()); }
        static const TestLinearSystemMultiRHS& getLinAlgEqSystemMultiRHS(int index) { return getMultiRHSSystemsArray()[index].second; }
        static const TestLinearSystemMultiRHS& getLinAlgEqSystemMultiRHS(std::string sysName)
        {
            for (const auto &sys : getMultiRHSSystemsArray())
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("LinearAlgEqTestBed::getLinAlgEqSystemMultiRHS: system with name " + sysName + " not found");
        }

        //==============================================================================================================
        // COMPLEX MATRICES
        //==============================================================================================================
        static int numComplexSystems() { return static_cast<int>(getComplexSystemsArray().size()); }
        static const TestLinearSystemComplex& getComplexSystem(int index) { return getComplexSystemsArray()[index].second; }
        static const TestLinearSystemComplex& getComplexSystem(std::string sysName)
        {
            for (const auto &sys : getComplexSystemsArray())
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("LinearAlgEqTestBed::getComplexSystem: system with name " + sysName + " not found");
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystemComplex*>> getComplexSystemsBySize(int n)
        {
            std::vector<std::pair<std::string, const TestLinearSystemComplex*>> result;
            for (const auto &sys : getComplexSystemsArray())
            {
                if (sys.second._n == n)
                    result.push_back({sys.first, &sys.second});
            }
            return result;
        }

        //==============================================================================================================
        // OVERDETERMINED MATRICES (for least squares)
        //==============================================================================================================
        static int numOverdeterminedSystems() { return static_cast<int>(getOverdeterminedSystemsArray().size()); }
        static const TestLinearSystemOverdetermined& getOverdeterminedSystem(int index) { return getOverdeterminedSystemsArray()[index].second; }
        static const TestLinearSystemOverdetermined& getOverdeterminedSystem(std::string sysName)
        {
            for (const auto &sys : getOverdeterminedSystemsArray())
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("LinearAlgEqTestBed::getOverdeterminedSystem: system with name " + sysName + " not found");
        }

        //==============================================================================================================
        // BY CATEGORY
        //==============================================================================================================
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getRandomMatrices()
        {
            return filterSystems([](const std::string& name) { 
                return name.find("mat_") == 0 && name.find("hilbert") == std::string::npos; 
            });
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getHilbertMatrices()
        {
            return filterSystems([](const std::string& name) { return name.find("hilbert") != std::string::npos; });
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getPascalMatrices()
        {
            return filterSystems([](const std::string& name) { return name.find("pascal") != std::string::npos; });
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getVandermondeMatrices()
        {
            return filterSystems([](const std::string& name) { return name.find("vandermonde") != std::string::npos; });
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getFrankMatrices()
        {
            return filterSystems([](const std::string& name) { return name.find("frank") != std::string::npos; });
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getKahanMatrices()
        {
            return filterSystems([](const std::string& name) { return name.find("kahan") != std::string::npos; });
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getDiagDominantMatrices()
        {
            return filterSystems([](const std::string& name) { return name.find("diag_dominant") != std::string::npos; });
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystemSymmetric*>> getSPDMatrices()
        {
            return filterSymmetricSystems([](const std::string& name) { return name.find("spd_") != std::string::npos; });
        }
        
        // All classic test matrices (Hilbert, Pascal, Vandermonde, Frank, Kahan)
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getClassicMatrices()
        {
            return filterSystems([](const std::string& name) { 
                return name.find("hilbert") != std::string::npos ||
                       name.find("pascal") != std::string::npos ||
                       name.find("vandermonde") != std::string::npos ||
                       name.find("frank") != std::string::npos ||
                       name.find("kahan") != std::string::npos;
            });
        }

        //==============================================================================================================
        // BY SIZE
        //==============================================================================================================
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getSystemsBySize(int n)
        {
            return filterSystems([n](const std::string&, const TestLinearSystem& sys) { return sys._n == n; });
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystemSymmetric*>> getSymmetricSystemsBySize(int n)
        {
            return filterSymmetricSystems([n](const std::string&, const TestLinearSystemSymmetric& sys) { return sys._n == n; });
        }
        
        // Small matrices (n <= 4) - good for debugging, visual inspection
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getSmallMatrices()
        {
            return filterSystems([](const std::string&, const TestLinearSystem& sys) { return sys._n <= 4; });
        }
        
        // Medium matrices (5 <= n <= 8)
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getMediumMatrices()
        {
            return filterSystems([](const std::string&, const TestLinearSystem& sys) { return sys._n >= 5 && sys._n <= 8; });
        }
        
        // Large matrices (n >= 10) - stress testing
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getLargeMatrices()
        {
            return filterSystems([](const std::string&, const TestLinearSystem& sys) { return sys._n >= 10; });
        }

        //==============================================================================================================
        // BY CONDITIONING (requires _cond_2 to be set)
        //==============================================================================================================
        // Well-conditioned: κ₂ < 100
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getWellConditioned()
        {
            return filterSystems([](const std::string&, const TestLinearSystem& sys) { 
                return sys._cond_2 > 0 && sys._cond_2 < 100; 
            });
        }
        
        // Moderately conditioned: 100 <= κ₂ < 1000
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getModeratelyConditioned()
        {
            return filterSystems([](const std::string&, const TestLinearSystem& sys) { 
                return sys._cond_2 >= 100 && sys._cond_2 < 1000; 
            });
        }
        
        // Ill-conditioned: κ₂ >= 1000
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getIllConditioned()
        {
            return filterSystems([](const std::string&, const TestLinearSystem& sys) { 
                return sys._cond_2 >= 1000; 
            });
        }

        //==============================================================================================================
        // UTILITY - Get all systems as vectors (easier for range-based for loops)
        //==============================================================================================================
        static std::vector<std::pair<std::string, const TestLinearSystem*>> getAllSystems()
        {
            std::vector<std::pair<std::string, const TestLinearSystem*>> result;
            for (const auto &sys : getLinearSystemsArray())
                result.push_back({sys.first, &sys.second});
            return result;
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystemSymmetric*>> getAllSymmetricSystems()
        {
            std::vector<std::pair<std::string, const TestLinearSystemSymmetric*>> result;
            for (const auto &sys : getSymmetricSystemsArray())
                result.push_back({sys.first, &sys.second});
            return result;
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystemComplex*>> getAllComplexSystems()
        {
            std::vector<std::pair<std::string, const TestLinearSystemComplex*>> result;
            for (const auto &sys : getComplexSystemsArray())
                result.push_back({sys.first, &sys.second});
            return result;
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystemOverdetermined*>> getAllOverdeterminedSystems()
        {
            std::vector<std::pair<std::string, const TestLinearSystemOverdetermined*>> result;
            for (const auto &sys : getOverdeterminedSystemsArray())
                result.push_back({sys.first, &sys.second});
            return result;
        }

    private:
        //==============================================================================================================
        // FILTER HELPERS
        //==============================================================================================================
        static std::vector<std::pair<std::string, const TestLinearSystem*>> filterSystems(
            std::function<bool(const std::string&)> predicate)
        {
            std::vector<std::pair<std::string, const TestLinearSystem*>> result;
            for (const auto &sys : getLinearSystemsArray())
            {
                if (predicate(sys.first))
                    result.push_back({sys.first, &sys.second});
            }
            return result;
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystem*>> filterSystems(
            std::function<bool(const std::string&, const TestLinearSystem&)> predicate)
        {
            std::vector<std::pair<std::string, const TestLinearSystem*>> result;
            for (const auto &sys : getLinearSystemsArray())
            {
                if (predicate(sys.first, sys.second))
                    result.push_back({sys.first, &sys.second});
            }
            return result;
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystemSymmetric*>> filterSymmetricSystems(
            std::function<bool(const std::string&)> predicate)
        {
            std::vector<std::pair<std::string, const TestLinearSystemSymmetric*>> result;
            for (const auto &sys : getSymmetricSystemsArray())
            {
                if (predicate(sys.first))
                    result.push_back({sys.first, &sys.second});
            }
            return result;
        }
        
        static std::vector<std::pair<std::string, const TestLinearSystemSymmetric*>> filterSymmetricSystems(
            std::function<bool(const std::string&, const TestLinearSystemSymmetric&)> predicate)
        {
            std::vector<std::pair<std::string, const TestLinearSystemSymmetric*>> result;
            for (const auto &sys : getSymmetricSystemsArray())
            {
                if (predicate(sys.first, sys.second))
                    result.push_back({sys.first, &sys.second});
            }
            return result;
        }
    };
}

#endif