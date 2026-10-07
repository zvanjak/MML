///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        sparse_matrix_tests.cpp                                             ///
///  Description: Tests for sparse matrix classes (COO, CSR, CSC)                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <mml/base/SparseMatrix/SparseMatrix.h>
#include <mml/base/Vector/Vector.h>

#include <cmath>

using namespace MML;
using namespace MML::SparseMatrix;
using Catch::Approx;

///////////////////////////////////////////////////////////////////////////
///                         COO Format Tests                            ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("SparseMatrixCOO - Construction", "[pde][sparse][coo]") {
    SparseMatrixCOO<double> coo(5, 5);
    
    REQUIRE(coo.rows() == 5);
    REQUIRE(coo.cols() == 5);
    REQUIRE(coo.nnz() == 0);
}

TEST_CASE("SparseMatrixCOO - Add entries", "[pde][sparse][coo]") {
    SparseMatrixCOO<double> coo(3, 3);
    
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(1, 1, 2.0);
    coo.addEntry(2, 2, 3.0);
    coo.addEntry(0, 1, 0.5);
    
    REQUIRE(coo.nnz() == 4);
}

TEST_CASE("SparseMatrixCOO - Skip zeros", "[pde][sparse][coo]") {
    SparseMatrixCOO<double> coo(3, 3);
    
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(1, 1, 0.0);  // Should be skipped
    coo.addEntry(2, 2, 3.0);
    
    REQUIRE(coo.nnz() == 2);
}

TEST_CASE("SparseMatrixCOO - Duplicate handling", "[pde][sparse][coo]") {
    SparseMatrixCOO<double> coo(3, 3);
    
    coo.addEntry(1, 1, 2.0);
    coo.addEntry(1, 1, 3.0);  // Duplicate - should be summed
    
    REQUIRE(coo.nnz() == 2);  // Before consolidation
    
    coo.consolidate();
    
    REQUIRE(coo.nnz() == 1);  // After consolidation
    REQUIRE(coo.triplets()[0].value == Approx(5.0));
}

TEST_CASE("SparseMatrixCOO - Consolidate and sort", "[pde][sparse][coo]") {
    SparseMatrixCOO<double> coo(3, 3);
    
    // Add entries out of order
    coo.addEntry(2, 2, 3.0);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(1, 1, 2.0);
    
    coo.consolidate();
    
    REQUIRE(coo.isSorted());
    const auto& triplets = coo.triplets();
    REQUIRE(triplets[0].row == 0);
    REQUIRE(triplets[1].row == 1);
    REQUIRE(triplets[2].row == 2);
}

TEST_CASE("SparseMatrixCOO - Stencil helpers", "[pde][sparse][coo]") {
    SparseMatrixCOO<double> coo(5, 5);
    
    // Add 1D Laplacian for interior node i=2
    coo.addLaplacian1D(2, 2, 1.0);
    
    REQUIRE(coo.nnz() == 3);
    
    coo.consolidate();
    const auto& triplets = coo.triplets();
    
    // Should have entries at (2,1), (2,2), (2,3)
    REQUIRE(triplets[0].row == 2);
    REQUIRE(triplets[0].col == 1);
    REQUIRE(triplets[0].value == Approx(1.0));
    
    REQUIRE(triplets[1].row == 2);
    REQUIRE(triplets[1].col == 2);
    REQUIRE(triplets[1].value == Approx(-2.0));
    
    REQUIRE(triplets[2].row == 2);
    REQUIRE(triplets[2].col == 3);
    REQUIRE(triplets[2].value == Approx(1.0));
}

TEST_CASE("Sparse matrices - Structural validation", "[pde][sparse][validation]") {
    REQUIRE_THROWS_AS(SparseMatrixCOO<double>(-1, 2), MatrixDimensionError);
    REQUIRE_THROWS_AS(SparseMatrixCSR<double>(-1, 2), MatrixDimensionError);
    REQUIRE_THROWS_AS(SparseMatrixCSC<double>(2, -1), MatrixDimensionError);

    SparseMatrixCOO<double> coo(2, 2);
    coo.addEntry(1, 0, 1.0);
    coo.addEntry(0, 1, 2.0);
    REQUIRE(coo.validate());

    REQUIRE_THROWS_AS((SparseMatrixCSR<double>(2, 2,
        {1.0, 2.0}, {0, 1}, {0, 2, 1})), ArgumentError);
    REQUIRE_THROWS_AS((SparseMatrixCSR<double>(2, 2,
        {1.0}, {2}, {0, 1, 1})), ArgumentError);
    REQUIRE_THROWS_AS((SparseMatrixCSR<double>(2, 2,
        {1.0, 2.0}, {1, 0}, {0, 2, 2})), ArgumentError);

    REQUIRE_THROWS_AS((SparseMatrixCSC<double>(2, 2,
        {1.0, 2.0}, {0, 1}, {0, 2, 1})), ArgumentError);
    REQUIRE_THROWS_AS((SparseMatrixCSC<double>(2, 2,
        {1.0}, {2}, {0, 1, 1})), ArgumentError);
    REQUIRE_THROWS_AS((SparseMatrixCSC<double>(2, 2,
        {1.0, 2.0}, {1, 0}, {0, 2, 2})), ArgumentError);

    SparseMatrixCSR<double> mutableCsr(2, 2, {1.0}, {0}, {0, 1, 1});
    REQUIRE(mutableCsr.validate());
    mutableCsr.colIndices()[0] = 2;
    REQUIRE_FALSE(mutableCsr.validate());
}

TEST_CASE("SparseMatrixCOO - 2D Laplacian respects grid boundaries", "[pde][sparse][coo][laplacian]") {
    SparseMatrixCOO<double> coo(9, 9);
    coo.addLaplacian2D_5pt(0, 0, 0, 3, 3, 1.0);
    coo.consolidate();

    REQUIRE(coo.nnz() == 3);
    REQUIRE(coo.validate());
    REQUIRE(coo.triplets()[0].col == 0);
    REQUIRE(coo.triplets()[1].col == 1);
    REQUIRE(coo.triplets()[2].col == 3);
}

///////////////////////////////////////////////////////////////////////////
///                         CSR Format Tests                            ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("SparseMatrixCSR - Construction from COO", "[pde][sparse][csr]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSR<double> csr(coo);
    
    REQUIRE(csr.rows() == 3);
    REQUIRE(csr.cols() == 3);
    REQUIRE(csr.nnz() == 5);
}

TEST_CASE("SparseMatrixCSR - Element access", "[pde][sparse][csr]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSR<double> csr(coo);
    
    REQUIRE(csr(0, 0) == Approx(1.0));
    REQUIRE(csr(0, 1) == Approx(0.0));  // Zero element
    REQUIRE(csr(0, 2) == Approx(2.0));
    REQUIRE(csr(1, 1) == Approx(3.0));
    REQUIRE(csr(2, 0) == Approx(4.0));
    REQUIRE(csr(2, 2) == Approx(5.0));
}

TEST_CASE("SparseMatrixCSR - Row pointers", "[pde][sparse][csr]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 2, 2.0);  // Row 0: 2 entries
    coo.addEntry(1, 1, 3.0);  // Row 1: 1 entry
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(2, 2, 5.0);  // Row 2: 2 entries
    
    SparseMatrixCSR<double> csr(coo);
    
    const auto& rp = csr.rowPointers();
    REQUIRE(rp[0] == 0);
    REQUIRE(rp[1] == 2);
    REQUIRE(rp[2] == 3);
    REQUIRE(rp[3] == 5);
    
    REQUIRE(csr.rowNnz(0) == 2);
    REQUIRE(csr.rowNnz(1) == 1);
    REQUIRE(csr.rowNnz(2) == 2);
}

TEST_CASE("SparseMatrixCSR - Matrix-vector product", "[pde][sparse][csr]") {
    // Create matrix:
    // [1 0 2]
    // [0 3 0]
    // [4 0 5]
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSR<double> csr(coo);
    
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y;
    
    csr.multiply(x, y);
    
    // y = A*x = [1*1 + 2*3, 3*2, 4*1 + 5*3] = [7, 6, 19]
    REQUIRE(y[0] == Approx(7.0));
    REQUIRE(y[1] == Approx(6.0));
    REQUIRE(y[2] == Approx(19.0));
}

TEST_CASE("SparseMatrixCSR - MML Vector operations", "[pde][sparse][csr][vector]") {
    SparseMatrixCSR<double> matrix(2, 2, {2.0, 3.0}, {0, 1}, {0, 1, 2});
    Vector<double> input{4.0, 5.0};
    Vector<double> output;

    matrix.multiply(input, output);
    REQUIRE(output[0] == Approx(8.0));
    REQUIRE(output[1] == Approx(15.0));

    output = matrix * input;
    REQUIRE(output[0] == Approx(8.0));
    REQUIRE(output[1] == Approx(15.0));

    matrix.gemv(2.0, input, 1.0, output);
    REQUIRE(output[0] == Approx(24.0));
    REQUIRE(output[1] == Approx(45.0));

    Vector<double> transposeOutput;
    matrix.multiplyTranspose(input, transposeOutput);
    REQUIRE(transposeOutput[0] == Approx(8.0));
    REQUIRE(transposeOutput[1] == Approx(15.0));

    Vector<double> solution;
    Vector<double> rightHandSide{8.0, 15.0};
    matrix.solveLower(rightHandSide, solution);
    REQUIRE(solution[0] == Approx(4.0));
    REQUIRE(solution[1] == Approx(5.0));
    matrix.solveUpper(rightHandSide, solution);
    REQUIRE(solution[0] == Approx(4.0));
    REQUIRE(solution[1] == Approx(5.0));
}

TEST_CASE("SparseMatrixCSR - Transpose-vector product", "[pde][sparse][csr]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSR<double> csr(coo);
    
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y;
    
    csr.multiplyTranspose(x, y);
    
    // y = A^T*x = [1*1 + 4*3, 3*2, 2*1 + 5*3] = [13, 6, 17]
    REQUIRE(y[0] == Approx(13.0));
    REQUIRE(y[1] == Approx(6.0));
    REQUIRE(y[2] == Approx(17.0));
}

TEST_CASE("SparseMatrixCSR - Diagonal extraction", "[pde][sparse][csr]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 1, 0.5);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 1, 0.5);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSR<double> csr(coo);
    
    auto diag = csr.getDiagonal();
    
    REQUIRE(diag[0] == Approx(1.0));
    REQUIRE(diag[1] == Approx(3.0));
    REQUIRE(diag[2] == Approx(5.0));
}

TEST_CASE("SparseMatrixCSR - Norms", "[pde][sparse][csr]") {
    SparseMatrixCOO<double> coo(2, 2);
    coo.addEntry(0, 0, 3.0);
    coo.addEntry(0, 1, 4.0);
    coo.addEntry(1, 0, 0.0);
    coo.addEntry(1, 1, 5.0);
    
    SparseMatrixCSR<double> csr(coo);
    
    // Frobenius norm: sqrt(9 + 16 + 25) = sqrt(50)
    REQUIRE(csr.normFrobenius() == Approx(std::sqrt(50.0)));
    
    // Infinity norm: max row sum = max(7, 5) = 7
    REQUIRE(csr.normInf() == Approx(7.0));
    
    // 1-norm: max column sum = max(3, 9) = 9
    REQUIRE(csr.norm1() == Approx(9.0));
}

TEST_CASE("SparseMatrixCSR - Identity matrix", "[pde][sparse][csr]") {
    auto I = identityCSR<double>(4);
    
    REQUIRE(I.rows() == 4);
    REQUIRE(I.cols() == 4);
    REQUIRE(I.nnz() == 4);
    
    for (int i = 0; i < 4; ++i) {
        REQUIRE(I(i, i) == Approx(1.0));
        if (i > 0) REQUIRE(I(i, i-1) == Approx(0.0));
    }
    
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
    auto y = I * x;
    
    for (int i = 0; i < 4; ++i) {
        REQUIRE(y[i] == Approx(x[i]));
    }
}

TEST_CASE("SparseMatrixCSR - Symmetry check", "[pde][sparse][csr]") {
    // Symmetric matrix
    SparseMatrixCOO<double> coo_sym(3, 3);
    coo_sym.addEntry(0, 0, 2.0);
    coo_sym.addEntry(0, 1, 1.0);
    coo_sym.addEntry(1, 0, 1.0);
    coo_sym.addEntry(1, 1, 3.0);
    coo_sym.addEntry(1, 2, 1.0);
    coo_sym.addEntry(2, 1, 1.0);
    coo_sym.addEntry(2, 2, 2.0);
    
    SparseMatrixCSR<double> csr_sym(coo_sym);
    REQUIRE(csr_sym.isSymmetric());
    
    // Non-symmetric matrix
    SparseMatrixCOO<double> coo_nonsym(3, 3);
    coo_nonsym.addEntry(0, 0, 2.0);
    coo_nonsym.addEntry(0, 1, 1.0);
    coo_nonsym.addEntry(1, 0, 2.0);  // Different from (0,1)
    coo_nonsym.addEntry(1, 1, 3.0);
    
    SparseMatrixCSR<double> csr_nonsym(coo_nonsym);
    REQUIRE_FALSE(csr_nonsym.isSymmetric());
}

TEST_CASE("SparseMatrixCSR - Lower triangular solve", "[pde][sparse][csr]") {
    // L = [2 0 0]
    //     [1 3 0]
    //     [1 1 4]
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 2.0);
    coo.addEntry(1, 0, 1.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 0, 1.0);
    coo.addEntry(2, 1, 1.0);
    coo.addEntry(2, 2, 4.0);
    
    SparseMatrixCSR<double> L(coo);
    
    std::vector<double> b = {2.0, 4.0, 6.0};
    std::vector<double> x;
    
    L.solveLower(b, x);
    
    // Verify L * x = b
    auto Lx = L * x;
    for (int i = 0; i < 3; ++i) {
        REQUIRE(Lx[i] == Approx(b[i]).margin(1e-10));
    }
}

TEST_CASE("SparseMatrixCSR - Triangular solves reject invalid diagonals", "[pde][sparse][csr]") {
    const std::vector<double> b = {1.0, 2.0};
    std::vector<double> x;

    SparseMatrixCSR<double> missingDiagonal(2, 2,
        {1.0, 1.0}, {0, 0}, {0, 1, 2});
    REQUIRE_THROWS_AS(missingDiagonal.solveLower(b, x), SingularMatrixError);
    REQUIRE_THROWS_AS(missingDiagonal.solveUpper(b, x), SingularMatrixError);

    SparseMatrixCSR<double> zeroDiagonal(2, 2,
        {1.0, 0.0}, {0, 1}, {0, 1, 2});
    REQUIRE_THROWS_AS(zeroDiagonal.solveLower(b, x), SingularMatrixError);
    REQUIRE_THROWS_AS(zeroDiagonal.solveUpper(b, x), SingularMatrixError);
}

TEST_CASE("SparseMatrixCSR - gemv operation", "[pde][sparse][csr]") {
    SparseMatrixCOO<double> coo(2, 2);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 1, 2.0);
    coo.addEntry(1, 0, 3.0);
    coo.addEntry(1, 1, 4.0);
    
    SparseMatrixCSR<double> A(coo);
    
    std::vector<double> x = {1.0, 1.0};
    std::vector<double> y = {10.0, 20.0};
    
    // y = 2*A*x + 3*y = 2*[3, 7] + 3*[10, 20] = [6, 14] + [30, 60] = [36, 74]
    A.gemv(2.0, x, 3.0, y);
    
    REQUIRE(y[0] == Approx(36.0));
    REQUIRE(y[1] == Approx(74.0));
}

///////////////////////////////////////////////////////////////////////////
///                         CSC Format Tests                            ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("SparseMatrixCSC - Construction from COO", "[pde][sparse][csc]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSC<double> csc(coo);
    
    REQUIRE(csc.rows() == 3);
    REQUIRE(csc.cols() == 3);
    REQUIRE(csc.nnz() == 5);
    
    // Check column pointers
    const auto& cp = csc.colPointers();
    REQUIRE(cp[0] == 0);
    REQUIRE(cp[1] == 2);  // Column 0 has 2 entries
    REQUIRE(cp[2] == 3);  // Column 1 has 1 entry
    REQUIRE(cp[3] == 5);  // Column 2 has 2 entries
}

TEST_CASE("SparseMatrixCSC - Element access", "[pde][sparse][csc]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSC<double> csc(coo);
    
    REQUIRE(csc(0, 0) == Approx(1.0));
    REQUIRE(csc(2, 0) == Approx(4.0));
    REQUIRE(csc(1, 1) == Approx(3.0));
    REQUIRE(csc(0, 2) == Approx(2.0));
    REQUIRE(csc(2, 2) == Approx(5.0));
    REQUIRE(csc(0, 1) == Approx(0.0));  // Zero entry
}

TEST_CASE("SparseMatrixCSC - Matrix-vector product", "[pde][sparse][csc]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSC<double> csc(coo);
    
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y;
    
    csc.multiply(x, y);
    
    // Same result as CSR
    REQUIRE(y[0] == Approx(7.0));
    REQUIRE(y[1] == Approx(6.0));
    REQUIRE(y[2] == Approx(19.0));
}

TEST_CASE("SparseMatrixCSC - MML Vector operations", "[pde][sparse][csc][vector]") {
    SparseMatrixCSC<double> matrix(2, 2, {2.0, 3.0}, {0, 1}, {0, 1, 2});
    Vector<double> input{4.0, 5.0};
    Vector<double> output;

    matrix.multiply(input, output);
    REQUIRE(output[0] == Approx(8.0));
    REQUIRE(output[1] == Approx(15.0));

    output = matrix * input;
    REQUIRE(output[0] == Approx(8.0));
    REQUIRE(output[1] == Approx(15.0));

    matrix.multiplyTranspose(input, output);
    REQUIRE(output[0] == Approx(8.0));
    REQUIRE(output[1] == Approx(15.0));
}

TEST_CASE("SparseMatrixCSC - CSR conversion round-trip", "[pde][sparse][csc]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSR<double> csr_orig(coo);
    SparseMatrixCSC<double> csc(csr_orig);
    SparseMatrixCSR<double> csr_back = csc.toCSR();
    
    // Check same elements
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            REQUIRE(csr_back(i, j) == Approx(csr_orig(i, j)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////
///                         Factory Function Tests                      ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("Laplacian1D - Tridiagonal structure", "[pde][sparse][laplacian]") {
    auto L = laplacian1D<double>(5, 1.0);
    
    REQUIRE(L.rows() == 5);
    REQUIRE(L.cols() == 5);
    REQUIRE(L.nnz() == 13);  // 5 diagonal + 4 upper + 4 lower
    
    // Check diagonal
    for (int i = 0; i < 5; ++i) {
        REQUIRE(L(i, i) == Approx(2.0));
    }
    
    // Check off-diagonals
    for (int i = 0; i < 4; ++i) {
        REQUIRE(L(i, i+1) == Approx(-1.0));
        REQUIRE(L(i+1, i) == Approx(-1.0));
    }
}

TEST_CASE("Laplacian1D - Eigenvalue test", "[pde][sparse][laplacian]") {
    // For 1D Laplacian with h=1, eigenvalues are 4*sin²(k*π/(2(n+1))) for k=1..n
    // For n=10: max eigenvalue ≈ 4*sin²(10π/22) ≈ 3.83
    
    auto L = laplacian1D<double>(10);
    
    // Power iteration to estimate largest eigenvalue
    std::vector<double> v(10, 1.0);
    for (int iter = 0; iter < 100; ++iter) {
        auto Lv = L * v;
        double norm = norm2(Lv);
        for (size_t i = 0; i < v.size(); ++i) v[i] = Lv[i] / norm;
    }
    
    auto Lv = L * v;
    double lambda_max = dot(v, Lv) / dot(v, v);
    
    // Largest eigenvalue should be in range [3.5, 4.0] for this discretization
    REQUIRE(lambda_max > 3.5);
    REQUIRE(lambda_max <= 4.0);
}

TEST_CASE("Laplacian2D - 5-point stencil structure", "[pde][sparse][laplacian]") {
    auto L = laplacian2D<double>(3, 3);
    
    REQUIRE(L.rows() == 9);
    REQUIRE(L.cols() == 9);
    
    // Center node (1,1) -> index 4
    // Should have 5 non-zeros: itself and 4 neighbors
    REQUIRE(L(4, 4) == Approx(4.0));
    REQUIRE(L(4, 3) == Approx(-1.0));  // Left
    REQUIRE(L(4, 5) == Approx(-1.0));  // Right
    REQUIRE(L(4, 1) == Approx(-1.0));  // Bottom
    REQUIRE(L(4, 7) == Approx(-1.0));  // Top
    
    // Corner node (0,0) -> index 0
    // Should have 3 non-zeros
    REQUIRE(L(0, 0) == Approx(4.0));
    REQUIRE(L(0, 1) == Approx(-1.0));
    REQUIRE(L(0, 3) == Approx(-1.0));
}

TEST_CASE("Laplacian2D - Symmetry", "[pde][sparse][laplacian]") {
    auto L = laplacian2D<double>(5, 5);
    
    REQUIRE(L.isSquare());
    REQUIRE(L.isSymmetric());
}

TEST_CASE("Laplacian3D - 7-point stencil", "[pde][sparse][laplacian]") {
    auto L = laplacian3D<double>(3, 3, 3);
    
    REQUIRE(L.rows() == 27);
    REQUIRE(L.cols() == 27);
    
    // Center node (1,1,1) -> index 13
    REQUIRE(L(13, 13) == Approx(6.0));
    REQUIRE(L(13, 12) == Approx(-1.0));  // x-1
    REQUIRE(L(13, 14) == Approx(-1.0));  // x+1
    REQUIRE(L(13, 10) == Approx(-1.0));  // y-1
    REQUIRE(L(13, 16) == Approx(-1.0));  // y+1
    REQUIRE(L(13, 4) == Approx(-1.0));   // z-1
    REQUIRE(L(13, 22) == Approx(-1.0));  // z+1
}

///////////////////////////////////////////////////////////////////////////
///                         Vector Utility Tests                        ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector utilities - norm2", "[pde][sparse][util]") {
    std::vector<double> v = {3.0, 4.0};
    REQUIRE(norm2(v) == Approx(5.0));
}

TEST_CASE("Vector utilities - dot", "[pde][sparse][util]") {
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {4.0, 5.0, 6.0};
    REQUIRE(dot(a, b) == Approx(32.0));
}

TEST_CASE("Vector utilities - axpy", "[pde][sparse][util]") {
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y = {10.0, 20.0, 30.0};
    
    axpy(2.0, x, y);  // y = 2*x + y
    
    REQUIRE(y[0] == Approx(12.0));
    REQUIRE(y[1] == Approx(24.0));
    REQUIRE(y[2] == Approx(36.0));
}

TEST_CASE("Vector utilities - residual", "[pde][sparse][util]") {
    auto I = identityCSR<double>(3);
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> b = {1.0, 2.0, 3.0};
    
    auto r = residual(I, x, b);  // r = b - I*x = 0
    
    REQUIRE(norm2(r) == Approx(0.0).margin(1e-14));
}

///////////////////////////////////////////////////////////////////////////
///                         Memory and Performance                      ///
///////////////////////////////////////////////////////////////////////////

TEST_CASE("SparseMatrixCSR - Memory usage", "[pde][sparse][memory]") {
    auto L = laplacian2D<double>(100, 100);
    
    // 10000 nodes, ~5 nnz per node (except boundaries)
    REQUIRE(L.nnz() < 50000);
    
    size_t mem = L.memoryUsage();
    // Should be much less than dense: 10000 * 10000 * 8 = 800 MB
    REQUIRE(mem < 1000000);  // Less than 1 MB
}

TEST_CASE("SparseMatrixCSR - Dense conversion", "[pde][sparse][dense]") {
    SparseMatrixCOO<double> coo(3, 3);
    coo.addEntry(0, 0, 1.0);
    coo.addEntry(0, 2, 2.0);
    coo.addEntry(1, 1, 3.0);
    coo.addEntry(2, 0, 4.0);
    coo.addEntry(2, 2, 5.0);
    
    SparseMatrixCSR<double> csr(coo);
    auto dense = csr.toDense();
    
    REQUIRE(dense[0][0] == Approx(1.0));
    REQUIRE(dense[0][1] == Approx(0.0));
    REQUIRE(dense[0][2] == Approx(2.0));
    REQUIRE(dense[1][1] == Approx(3.0));
    REQUIRE(dense[2][0] == Approx(4.0));
    REQUIRE(dense[2][2] == Approx(5.0));
}
