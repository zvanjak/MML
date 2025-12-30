#if !defined __MML_EIGENVALUE_TEST_BED_H
#define __MML_EIGENVALUE_TEST_BED_H

#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include <algorithm>

#include "eigenvalue_defs.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                        EIGENVALUE TEST DATA STRUCTURE                                  //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Represents a test case for eigenvalue solvers
     * 
     * Contains the matrix, known eigenvalues/eigenvectors, and metadata
     * for categorization and difficulty assessment.
     */
    struct TestEigenSystem
    {
        std::string name;                               ///< Descriptive name
        std::function<Matrix<Real>()> getMatrix;       ///< Factory to get test matrix
        
        // Eigenvalue data
        Vector<Real> realEigenvalues;                  ///< Real eigenvalues (if all real)
        std::vector<std::complex<Real>> complexEigenvalues;  ///< Complex eigenvalues (general)
        std::vector<Vector<Real>> eigenvectors;        ///< Eigenvectors (if known/computed)
        
        // Matrix properties
        int dimension;                                  ///< Matrix dimension n x n
        bool isSymmetric = false;                      ///< Symmetric matrices have real eigenvalues
        bool isPositiveDefinite = false;               ///< All eigenvalues positive
        bool hasRealEigenvalues = true;                ///< All eigenvalues are real
        bool isDefective = false;                      ///< Fewer eigenvectors than eigenvalues
        bool hasRepeatedEigenvalues = false;           ///< Multiple eigenvalues exist
        bool hasClusteredEigenvalues = false;          ///< Eigenvalues very close together
        
        // Numerical properties
        Real conditionNumber = 1.0;                    ///< Condition number (if known)
        Real eigenvalueSpread = 1.0;                   ///< Ratio of largest/smallest |λ|
        int multiplicityMax = 1;                       ///< Maximum algebraic multiplicity
        
        // Metadata
        std::string category;                          ///< standard, repeated, defective, etc.
        std::string description;                       ///< Detailed description
        int difficulty = 1;                            ///< 1=easy, 2=medium, 3=hard, 4=extreme
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                        CLASSIC/STANDARD TEST GENERATORS                                //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestEigenSystem getPascal5x5Test()
    {
        TestEigenSystem test;
        test.name = "Pascal 5x5";
        test.getMatrix = []() { return classic_pascal_5x5; };
        test.realEigenvalues = classic_pascal_5x5_eigenvalues;
        test.dimension = 5;
        test.isSymmetric = true;
        test.isPositiveDefinite = true;
        test.hasRealEigenvalues = true;
        test.conditionNumber = classic_pascal_5x5_eigenvalues[4] / classic_pascal_5x5_eigenvalues[0];
        test.eigenvalueSpread = test.conditionNumber;
        test.category = "standard";
        test.description = "Pascal matrix - symmetric positive definite, moderate conditioning.";
        test.difficulty = 1;
        return test;
    }

    inline TestEigenSystem getHilbert10x10Test()
    {
        TestEigenSystem test;
        test.name = "Hilbert 10x10";
        test.getMatrix = getHilbert10x10;
        test.realEigenvalues = classic_hilbert_10x10_eigenvalues;
        test.dimension = 10;
        test.isSymmetric = true;
        test.isPositiveDefinite = true;
        test.hasRealEigenvalues = true;
        test.conditionNumber = 1.6e13;  // Extremely ill-conditioned
        test.eigenvalueSpread = 1e10;
        test.category = "standard";
        test.description = "Hilbert matrix - extremely ill-conditioned, tests numerical precision.";
        test.difficulty = 3;
        return test;
    }

    inline TestEigenSystem getToeplitz20x20Test()
    {
        TestEigenSystem test;
        test.name = "Toeplitz Tridiag 20x20";
        test.getMatrix = getToeplitz20x20;
        test.realEigenvalues = getToeplitz20x20_eigenvalues();
        test.dimension = 20;
        test.isSymmetric = true;
        test.isPositiveDefinite = true;
        test.hasRealEigenvalues = true;
        test.conditionNumber = (2 + 2*std::cos(Constants::PI/21)) / (2 + 2*std::cos(20*Constants::PI/21));
        test.eigenvalueSpread = test.conditionNumber;
        test.category = "standard";
        test.description = "Tridiagonal Toeplitz - analytically known eigenvalues, good reference.";
        test.difficulty = 1;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                      REPEATED EIGENVALUE TEST GENERATORS                               //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestEigenSystem getRepeatedDouble4x4Test()
    {
        TestEigenSystem test;
        test.name = "Repeated Double 4x4";
        test.getMatrix = []() { return repeated_double_4x4; };
        test.realEigenvalues = repeated_double_4x4_eigenvalues;
        test.dimension = 4;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 2;
        test.category = "repeated";
        test.description = "Diagonal with double eigenvalue at 2. Tests multiplicity handling.";
        test.difficulty = 1;
        return test;
    }

    inline TestEigenSystem getRepeatedTriple5x5Test()
    {
        TestEigenSystem test;
        test.name = "Repeated Triple 5x5";
        test.getMatrix = []() { return repeated_triple_5x5; };
        test.realEigenvalues = repeated_triple_5x5_eigenvalues;
        test.dimension = 5;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 3;
        test.category = "repeated";
        test.description = "Diagonal with triple eigenvalue. Tests high multiplicity.";
        test.difficulty = 2;
        return test;
    }

    inline TestEigenSystem getRepeatedQuadruple4x4Test()
    {
        TestEigenSystem test;
        test.name = "Repeated Quadruple 4x4";
        test.getMatrix = []() { return repeated_quadruple_4x4; };
        test.realEigenvalues = repeated_quadruple_4x4_eigenvalues;
        test.dimension = 4;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 2;
        test.category = "repeated";
        test.description = "Non-diagonal symmetric with repeated eigenvalues. Tests orthogonalization.";
        test.difficulty = 2;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                       DEFECTIVE MATRIX TEST GENERATORS                                 //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestEigenSystem getDefectiveJordan2x2Test()
    {
        TestEigenSystem test;
        test.name = "Jordan Block 2x2";
        test.getMatrix = []() { return defective_jordan_2x2; };
        test.realEigenvalues = defective_jordan_2x2_eigenvalues;
        test.dimension = 2;
        test.isSymmetric = false;
        test.hasRealEigenvalues = true;
        test.isDefective = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 2;
        test.category = "defective";
        test.description = "Classic Jordan block - defective, only one eigenvector for double eigenvalue.";
        test.difficulty = 3;
        return test;
    }

    inline TestEigenSystem getDefectiveJordan3x3Test()
    {
        TestEigenSystem test;
        test.name = "Jordan Block 3x3";
        test.getMatrix = []() { return defective_jordan_3x3; };
        test.realEigenvalues = defective_jordan_3x3_eigenvalues;
        test.dimension = 3;
        test.isSymmetric = false;
        test.hasRealEigenvalues = true;
        test.isDefective = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 3;
        test.category = "defective";
        test.description = "3x3 Jordan block - maximally defective, single eigenvector.";
        test.difficulty = 3;
        return test;
    }

    inline TestEigenSystem getDefectiveMixed3x3Test()
    {
        TestEigenSystem test;
        test.name = "Mixed Jordan 3x3";
        test.getMatrix = []() { return defective_mixed_3x3; };
        test.realEigenvalues = defective_mixed_3x3_eigenvalues;
        test.dimension = 3;
        test.isSymmetric = false;
        test.hasRealEigenvalues = true;
        test.isDefective = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 2;
        test.category = "defective";
        test.description = "J_2(2) ⊕ J_1(3) structure. Mixed Jordan blocks.";
        test.difficulty = 3;
        return test;
    }

    inline TestEigenSystem getDefectiveTwoBlocks4x4Test()
    {
        TestEigenSystem test;
        test.name = "Two Jordan Blocks 4x4";
        test.getMatrix = []() { return defective_two_blocks_4x4; };
        test.realEigenvalues = defective_two_blocks_4x4_eigenvalues;
        test.dimension = 4;
        test.isSymmetric = false;
        test.hasRealEigenvalues = true;
        test.isDefective = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 2;
        test.category = "defective";
        test.description = "Two separate 2x2 Jordan blocks. Tests block detection.";
        test.difficulty = 3;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                     NEARLY DEFECTIVE TEST GENERATORS                                   //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestEigenSystem getNearlyDefective3x3Test(Real epsilon = 1e-8)
    {
        TestEigenSystem test;
        test.name = "Nearly Defective 3x3 (ε=" + std::to_string(epsilon) + ")";
        test.getMatrix = [epsilon]() { return getNearlyDefective3x3(epsilon); };
        test.dimension = 3;
        test.isSymmetric = false;
        test.hasRealEigenvalues = true;
        test.hasClusteredEigenvalues = true;
        test.category = "nearly_defective";
        test.description = "Perturbed Jordan block - eigenvalues very close, ill-conditioned eigenvectors.";
        test.difficulty = 4;
        return test;
    }

    inline TestEigenSystem getWilkinsonTest(int n = 21)
    {
        TestEigenSystem test;
        test.name = "Wilkinson W⁻ " + std::to_string(n) + "x" + std::to_string(n);
        test.getMatrix = [n]() { return getWilkinsonMinus(n); };
        test.dimension = n;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.hasClusteredEigenvalues = true;
        test.category = "nearly_defective";
        test.description = "Wilkinson matrix - famous for extremely close eigenvalue pairs.";
        test.difficulty = 4;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                      CLUSTERED EIGENVALUE TEST GENERATORS                              //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestEigenSystem getClusteredTight4x4Test()
    {
        TestEigenSystem test;
        test.name = "Clustered Tight 4x4";
        test.getMatrix = []() { return clustered_tight_4x4; };
        test.realEigenvalues = clustered_tight_4x4_eigenvalues;
        test.dimension = 4;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.hasClusteredEigenvalues = true;
        test.category = "clustered";
        test.description = "Three eigenvalues within 0.002 of each other. Tests discrimination.";
        test.difficulty = 2;
        return test;
    }

    inline TestEigenSystem getClusteredVeryTight5x5Test()
    {
        TestEigenSystem test;
        test.name = "Clustered Very Tight 5x5";
        test.getMatrix = []() { return clustered_very_tight_5x5; };
        test.realEigenvalues = clustered_very_tight_5x5_eigenvalues;
        test.dimension = 5;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.hasClusteredEigenvalues = true;
        test.category = "clustered";
        test.description = "Four eigenvalues within 0.0003. Extreme clustering.";
        test.difficulty = 3;
        return test;
    }

    inline TestEigenSystem getGradedClusters6x6Test()
    {
        TestEigenSystem test;
        test.name = "Graded Clusters 6x6";
        test.getMatrix = getGradedClusters6x6;
        test.realEigenvalues = clustered_graded_6x6_eigenvalues;
        test.dimension = 6;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.hasClusteredEigenvalues = true;
        test.eigenvalueSpread = 100.0;
        test.category = "clustered";
        test.description = "Three clusters at scales 0.1, 1, 10. Tests multi-scale handling.";
        test.difficulty = 2;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    COMPLEX EIGENVALUE TEST GENERATORS                                  //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestEigenSystem getRotation2x2Test(Real theta = Constants::PI / 4)
    {
        TestEigenSystem test;
        test.name = "Rotation 2x2 (θ=" + std::to_string(theta) + ")";
        test.getMatrix = [theta]() { return getRotation2x2(theta); };
        test.complexEigenvalues = {
            {std::cos(theta), std::sin(theta)},
            {std::cos(theta), -std::sin(theta)}
        };
        test.dimension = 2;
        test.isSymmetric = false;
        test.hasRealEigenvalues = false;
        test.category = "complex";
        test.description = "Pure rotation matrix - complex conjugate eigenvalues on unit circle.";
        test.difficulty = 2;
        return test;
    }

    inline TestEigenSystem getComplexPair3x3Test()
    {
        TestEigenSystem test;
        test.name = "Complex Pair 3x3";
        test.getMatrix = []() { return complex_pair_3x3; };
        test.complexEigenvalues = complex_pair_3x3_eigenvalues;
        test.dimension = 3;
        test.isSymmetric = false;
        test.hasRealEigenvalues = false;
        test.category = "complex";
        test.description = "One real eigenvalue (3) and one complex pair (1±i).";
        test.difficulty = 2;
        return test;
    }

    inline TestEigenSystem getComplexTwoPairs4x4Test()
    {
        TestEigenSystem test;
        test.name = "Two Complex Pairs 4x4";
        test.getMatrix = []() { return complex_two_pairs_4x4; };
        test.complexEigenvalues = complex_two_pairs_4x4_eigenvalues;
        test.dimension = 4;
        test.isSymmetric = false;
        test.hasRealEigenvalues = false;
        test.category = "complex";
        test.description = "All complex eigenvalues: ±i and ±2i (purely imaginary).";
        test.difficulty = 2;
        return test;
    }

    inline TestEigenSystem getCompanionTest()
    {
        TestEigenSystem test;
        test.name = "Companion x⁴-1";
        test.getMatrix = []() { return companion_x4_minus_1; };
        test.complexEigenvalues = companion_x4_minus_1_eigenvalues;
        test.dimension = 4;
        test.isSymmetric = false;
        test.hasRealEigenvalues = false;
        test.category = "complex";
        test.description = "Companion matrix for x⁴-1. Eigenvalues are 4th roots of unity.";
        test.difficulty = 2;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                   SPECIAL STRUCTURE TEST GENERATORS                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestEigenSystem getCirculant4x4Test()
    {
        TestEigenSystem test;
        test.name = "Circulant 4x4";
        test.getMatrix = getCirculant4x4;
        test.complexEigenvalues = circulant_4x4_eigenvalues;
        test.dimension = 4;
        test.isSymmetric = false;
        test.hasRealEigenvalues = false;
        test.category = "special";
        test.description = "Circulant matrix - eigenvalues are DFT of first row.";
        test.difficulty = 2;
        return test;
    }

    inline TestEigenSystem getBlockDiagonal6x6Test()
    {
        TestEigenSystem test;
        test.name = "Block Diagonal 6x6";
        test.getMatrix = getBlockDiagonal6x6;
        test.realEigenvalues = block_diagonal_6x6_eigenvalues;
        test.dimension = 6;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.category = "special";
        test.description = "Block diagonal - eigenvalues are union of block eigenvalues.";
        test.difficulty = 1;
        return test;
    }

    inline TestEigenSystem getStochastic4x4Test()
    {
        TestEigenSystem test;
        test.name = "Stochastic 4x4";
        test.getMatrix = []() { return stochastic_4x4; };
        test.dimension = 4;
        test.isSymmetric = false;
        test.hasRealEigenvalues = true;  // Perron-Frobenius applies
        test.category = "special";
        test.description = "Row stochastic - dominant eigenvalue = 1 (Perron-Frobenius).";
        test.difficulty = 2;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                      PATHOLOGICAL TEST GENERATORS                                      //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestEigenSystem getNilpotent3x3Test()
    {
        TestEigenSystem test;
        test.name = "Nilpotent 3x3";
        test.getMatrix = []() { return nilpotent_3x3; };
        test.realEigenvalues = Vector<Real>{0.0, 0.0, 0.0};
        test.dimension = 3;
        test.isSymmetric = false;
        test.hasRealEigenvalues = true;
        test.isDefective = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 3;
        test.category = "pathological";
        test.description = "Nilpotent - all eigenvalues zero but matrix is not zero.";
        test.difficulty = 3;
        return test;
    }

    inline TestEigenSystem getTinyEigenvalue3x3Test(Real eps = 1e-15)
    {
        TestEigenSystem test;
        test.name = "Tiny Eigenvalue 3x3";
        test.getMatrix = [eps]() { return getTinyEigenvalue3x3(eps); };
        test.realEigenvalues = Vector<Real>{eps, 1.0, 2.0};
        test.dimension = 3;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.eigenvalueSpread = 2.0 / eps;
        test.category = "pathological";
        test.description = "One eigenvalue at machine epsilon scale. Tests underflow handling.";
        test.difficulty = 3;
        return test;
    }

    inline TestEigenSystem getHugeSpread3x3Test(Real scale = 1e15)
    {
        TestEigenSystem test;
        test.name = "Huge Spread 3x3";
        test.getMatrix = [scale]() { return getHugeSpread3x3(scale); };
        test.realEigenvalues = Vector<Real>{1.0/scale, 1.0, scale};
        test.dimension = 3;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.eigenvalueSpread = scale * scale;
        test.conditionNumber = scale * scale;
        test.category = "pathological";
        test.description = "Eigenvalue spread of 10^30. Extreme dynamic range.";
        test.difficulty = 4;
        return test;
    }

    inline TestEigenSystem getAlmostSingular4x4Test(Real eps = 1e-12)
    {
        TestEigenSystem test;
        test.name = "Almost Singular 4x4";
        test.getMatrix = [eps]() { return getAlmostSingular4x4(eps); };
        test.dimension = 4;
        test.isSymmetric = true;
        test.hasRealEigenvalues = true;
        test.hasRepeatedEigenvalues = true;
        test.multiplicityMax = 3;
        test.hasClusteredEigenvalues = true;
        test.category = "pathological";
        test.description = "Nearly rank-deficient - three tiny eigenvalues, one large.";
        test.difficulty = 4;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                        CATEGORY RETRIEVAL FUNCTIONS                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Get all eigenvalue test cases
     */
    inline std::vector<TestEigenSystem> getAllEigenTests()
    {
        return {
            // Standard/Classic
            getPascal5x5Test(),
            getHilbert10x10Test(),
            getToeplitz20x20Test(),
            // Repeated eigenvalues
            getRepeatedDouble4x4Test(),
            getRepeatedTriple5x5Test(),
            getRepeatedQuadruple4x4Test(),
            // Defective
            getDefectiveJordan2x2Test(),
            getDefectiveJordan3x3Test(),
            getDefectiveMixed3x3Test(),
            getDefectiveTwoBlocks4x4Test(),
            // Nearly defective
            getNearlyDefective3x3Test(),
            getWilkinsonTest(),
            // Clustered
            getClusteredTight4x4Test(),
            getClusteredVeryTight5x5Test(),
            getGradedClusters6x6Test(),
            // Complex
            getRotation2x2Test(),
            getComplexPair3x3Test(),
            getComplexTwoPairs4x4Test(),
            getCompanionTest(),
            // Special structure
            getCirculant4x4Test(),
            getBlockDiagonal6x6Test(),
            getStochastic4x4Test(),
            // Pathological
            getNilpotent3x3Test(),
            getTinyEigenvalue3x3Test(),
            getHugeSpread3x3Test(),
            getAlmostSingular4x4Test()
        };
    }

    /**
     * @brief Get symmetric matrix tests only
     */
    inline std::vector<TestEigenSystem> getSymmetricEigenTests()
    {
        std::vector<TestEigenSystem> result;
        for (const auto& test : getAllEigenTests())
            if (test.isSymmetric) result.push_back(test);
        return result;
    }

    /**
     * @brief Get defective matrix tests
     */
    inline std::vector<TestEigenSystem> getDefectiveEigenTests()
    {
        return {
            getDefectiveJordan2x2Test(),
            getDefectiveJordan3x3Test(),
            getDefectiveMixed3x3Test(),
            getDefectiveTwoBlocks4x4Test(),
            getNilpotent3x3Test()
        };
    }

    /**
     * @brief Get complex eigenvalue tests
     */
    inline std::vector<TestEigenSystem> getComplexEigenTests()
    {
        return {
            getRotation2x2Test(),
            getComplexPair3x3Test(),
            getComplexTwoPairs4x4Test(),
            getCompanionTest(),
            getCirculant4x4Test()
        };
    }

    /**
     * @brief Get challenging tests (difficulty >= 3)
     */
    inline std::vector<TestEigenSystem> getChallengingEigenTests()
    {
        std::vector<TestEigenSystem> result;
        for (const auto& test : getAllEigenTests())
            if (test.difficulty >= 3) result.push_back(test);
        return result;
    }

    /**
     * @brief Get easy tests for basic validation
     */
    inline std::vector<TestEigenSystem> getEasyEigenTests()
    {
        return {
            getPascal5x5Test(),
            getToeplitz20x20Test(),
            getRepeatedDouble4x4Test(),
            getBlockDiagonal6x6Test()
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                            VERIFICATION UTILITIES                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Verify eigenvalues match expected (real case)
     * @param computed Computed eigenvalues (may be unsorted)
     * @param expected Expected eigenvalues
     * @param tolerance Acceptable error
     * @return True if all eigenvalues match within tolerance
     */
    inline bool verifyRealEigenvalues(Vector<Real> computed, Vector<Real> expected, Real tolerance = 1e-10)
    {
        if (computed.size() != expected.size()) return false;
        
        // Sort both for comparison
        std::vector<Real> comp(computed.size()), exp(expected.size());
        for (size_t i = 0; i < computed.size(); ++i) {
            comp[i] = computed[i];
            exp[i] = expected[i];
        }
        std::sort(comp.begin(), comp.end());
        std::sort(exp.begin(), exp.end());
        
        for (size_t i = 0; i < comp.size(); ++i)
            if (std::abs(comp[i] - exp[i]) > tolerance)
                return false;
        return true;
    }

    /**
     * @brief Compute eigenvalue error (max absolute difference)
     */
    inline Real computeEigenvalueError(Vector<Real> computed, Vector<Real> expected)
    {
        std::vector<Real> comp(computed.size()), exp(expected.size());
        for (size_t i = 0; i < computed.size(); ++i) {
            comp[i] = computed[i];
            exp[i] = expected[i];
        }
        std::sort(comp.begin(), comp.end());
        std::sort(exp.begin(), exp.end());
        
        Real maxErr = 0;
        for (size_t i = 0; i < comp.size(); ++i)
            maxErr = std::max(maxErr, std::abs(comp[i] - exp[i]));
        return maxErr;
    }

    /**
     * @brief Verify eigenvector orthogonality (for symmetric matrices)
     * @param V Matrix with eigenvectors as columns
     * @param tolerance Acceptable deviation from orthonormality
     */
    inline bool verifyOrthonormality(const Matrix<Real>& V, Real tolerance = 1e-10)
    {
        int n = V.RowNum();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Real dot = 0;
                for (int k = 0; k < n; ++k)
                    dot += V(k, i) * V(k, j);
                Real expected = (i == j) ? 1.0 : 0.0;
                if (std::abs(dot - expected) > tolerance)
                    return false;
            }
        }
        return true;
    }

    /**
     * @brief Verify A*v = λ*v for eigenpair
     */
    inline Real computeResidual(const Matrix<Real>& A, const Vector<Real>& v, Real lambda)
    {
        int n = A.RowNum();
        Real residual = 0;
        for (int i = 0; i < n; ++i) {
            Real Av_i = 0;
            for (int j = 0; j < n; ++j)
                Av_i += A(i, j) * v[j];
            residual += (Av_i - lambda * v[i]) * (Av_i - lambda * v[i]);
        }
        return std::sqrt(residual);
    }

} // namespace MML::TestBeds

#endif // __MML_EIGENVALUE_TEST_BED_H
