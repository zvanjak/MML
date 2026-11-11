#include <mml/MMLBase.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/BaseUtils.h>

#include <mml/core/LinAlgEqSolvers.h>

#include <mml/algorithms/EigenSystemSolvers.h>
#include <mml/tools/persistence/VectorBinary.h>
#include <mml/tools/persistence/MatrixIO.h>

#include "test_data/linear_alg_eq_systems_real_defs.h"
#include "test_data/linear_alg_eq_systems_complex_defs.h"
#include "test_data/linear_alg_eq_systems_classic_defs.h"

#include <filesystem>
#include <chrono>

using namespace MML;
using namespace MML::Persistence;

// TODO - LOW, EASY - generate random complex matrix, i ispis da se lako copy paste u complex_defs.h
Matrix<Real> GenerateRandomMatrix(int n, int m, Real min, Real max)
{
    Matrix<Real> ret(n, m);

    srand((unsigned int)time(NULL));

    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
        {
            Real rnd = min + rand() / (Real)RAND_MAX * (max - min);
            ret(i,j) = rnd;
        }

    return ret;
}

void MatrixAnalyzer(std::string matName, const Matrix<Real> &origMat, const Vector<Real> &rhs)
{
    //std::cout << "*******************************************************************\n";
    //std::cout << "*****   MATRIX ANALYZER FOR MATRIX   " << origMat.rows() << "x" << origMat.cols() << "   *****\n";

    //int width = 15;
    //int prec  = 10;

    //int dim = origMat.rows();
    //Matrix<Real>  matcopy(origMat);
    //Vector<Real>  rhscopy(rhs);

    //std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    //std::cout << "Right side:\n";        rhscopy.Print(std::cout,10,3);

    //Vector<Real>     vecSol(rhs.size());

    //LUDecompositionSolver<Real> luSolver(matcopy);
    //luSolver.Solve(rhscopy, vecSol);

    //std::cout << "\nSolution:\n" << vecSol << std::endl;
    //std::cout << "Multiplying solution with matrix: " << origMat * vecSol << std::endl;

    //Matrix<Real> inv(rhs.size(),rhs.size());
    //luSolver.inverse(inv);

    //std::cout << "Inverse " << inv << std::endl;
    //std::cout << "Orig * inv = " << origMat * inv << std::endl;    

    //EigenSolver eigen_solver(matcopy, true, false);

    //std::cout << "Num real eigenvalues    : " << eigen_solver.getNumReal() << std::endl;
    //std::cout << "Num complex eigenvalues : " << eigen_solver.getNumComplex() << "\n\n";
    //
    //std::cout << "Eigenvalues         : "; eigen_solver.getEigenvalues().Print(std::cout,width,prec); std::cout << std::endl;
    //std::cout << "Real eigenvalues    : "; eigen_solver.getRealEigenvalues().Print(std::cout,width,prec); std::cout << std::endl;
    //std::cout << "Complex eigenvalues : "; eigen_solver.getComplexEigenvalues().Print(std::cout,width,prec); std::cout << "\n\n";
    // 
    //std::cout << "Eigenvectors matrix:\n"; eigen_solver.zz.Print(std::cout, width, 7); std::cout << std::endl;

    //Matrix<Complex> mat_cmplx = Utils::CmplxMatFromRealMat(matcopy);

    //for(int i=0; i<dim; i++ )
    //{
    //    if( eigen_solver.isRealEigenvalue(i) )
    //    {
    //        double        eigen_val = eigen_solver.getEigenvalues()[i].real();
    //        Vector<Real>  eigen_vec = eigen_solver.getRealPartEigenvector(i);

    //        std::cout << "Eigen value  " << i+1 << " (real)    = " << eigen_val << std::endl;
    //        std::cout << "Eigen vector " << i+1 << " (real)    = "; eigen_vec.Print(std::cout,width,prec); std::cout << std::endl;

    //        std::cout << "Mat * eig_vec            = "; (matcopy * eigen_vec).Print(std::cout,width,prec); std::cout << std::endl;
    //        std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";            
    //    }
    //    else
    //    {
    //        Complex          eigen_val = eigen_solver.getEigenvalues()[i];
    //        Vector<Complex>  eigen_vec = eigen_solver.getEigenvector(i);
    //        
    //        std::cout << "Eigen value  " << i+1 << " (complex) = " << eigen_val << std::endl;
    //        std::cout << "Eigen vector " << i+1 << " (complex) = "; eigen_solver.getEigenvector(i).Print(std::cout, width, prec); std::cout << std::endl;

    //        std::cout << "Mat * eig_vec            = "; (mat_cmplx * eigen_vec).Print(std::cout, width, prec); std::cout << std::endl;
    //        std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";            
    //    }
    //}
    //
    //std::cout << "    const static inline MML::Vector     " << matName << "_rhs0_sol{";
    //for(int i=0; i<dim; i++ )
    //{
    //    std::cout << std::setprecision(15) << vecSol[i];
    //    if( i < dim-1 )
    //        std::cout << ", ";
    //    else
    //        std::cout << "};\n\n";
    //}

    //std::cout << "    const static inline VectorComplex   " << matName << "_eigen_val{ Complex(" << eigen_solver.getEigenvalues()[0].real() << ", " << eigen_solver.getEigenvalues()[0].imag() << "),\n";
    //for(int i=1; i<dim; i++ )
    //{
    //    std::cout << "                                                             Complex(" << eigen_solver.getEigenvalues()[i].real() << ", " << eigen_solver.getEigenvalues()[i].imag() << ")";
    //    if( i < dim-1 )
    //        std::cout << ",\n";
    //    else
    //        std::cout << "};\n";
    //}

    //std::cout << "    const static inline std::vector<VectorComplex> " << matName << "_eigen_vecs\n";
    //std::cout << std::setprecision(10) << "    {\n";
    //for(int i=0; i<dim; i++ )
    //{
    //    std::cout << "        VectorComplex{ Complex(" << eigen_solver.getEigenvector(i)[0].real() << ", " << eigen_solver.getEigenvector(i)[0].imag() << "),";
    //    for(int j=1; j<dim; j++ )
    //    {
    //        std::cout << "Complex(" << eigen_solver.getEigenvector(i)[j].real() << ", " << eigen_solver.getEigenvector(i)[j].imag() << ")";
    //        if( j < dim-1 )
    //            std::cout << ", ";
    //        else
    //            std::cout << " },\n";
    //    }
    //    // if( i < dim-1 )
    //    //     std::cout << ",\n";
    //    // else
    //    //     std::cout << "\n";
    //}
    //std::cout << "    };\n";
}


void SymmMatrixAnalyzer(std::string matName, const MatrixSym<Real> &origMat, const Vector<Real> &rhs)
{
    //std::cout << "*******************************************************************\n";
    //std::cout << "*****   SYMMETRIC MATRIX ANALYZER FOR MATRIX   " << origMat.rows() << "x" << origMat.cols() << "   *****\n";

    //int width = 15;
    //int prec  = 10;

    //int dim = origMat.rows();
    //Vector<Real>  rhscopy(rhs);

    //std::cout << "Initial matrix:\n";    origMat.Print(std::cout,10,3);
    //std::cout << "Right side:\n";        rhscopy.Print(std::cout,10,3);

    //auto matcopy = origMat.GetAsMatrix();
    //Vector<Real>   vecSol(rhs.size());
    //
    //LUDecompositionSolver<Real> luSolver(matcopy);
    //luSolver.Solve(rhscopy, vecSol);

    //std::cout << "\nSolution:\n" << vecSol << std::endl;
    //std::cout << "Multiplying solution with matrix: " << origMat * vecSol << std::endl;

    //Matrix<Real> inv(rhs.size(),rhs.size());
    //luSolver.inverse(inv);

    //std::cout << "Inverse " << inv << std::endl;
    //std::cout << "Orig * inv = " << origMat * inv << std::endl;    

    //SymmMatEigenSolver eigen_solver(origMat, true);

    //std::cout << "Eigenvalues         : "; eigen_solver.getEigenvalues().Print(std::cout,width,prec); std::cout << std::endl;
    //std::cout << "Eigenvectors matrix:\n"; eigen_solver.z.Print(std::cout, width, 7); std::cout << std::endl;

    //auto matcopy2 = origMat.GetAsMatrix();
    //for(int i=0; i<dim; i++ )
    //{
    //    double        eigen_val = eigen_solver.getEigenvalues()[i];
    //    Vector<Real>  eigen_vec = eigen_solver.getEigenvector(i);

    //    std::cout << "Eigen value  " << i+1 << " (real)    = " << eigen_val << std::endl;
    //    std::cout << "Eigen vector " << i+1 << " (real)    = "; eigen_vec.Print(std::cout,width,prec); std::cout << std::endl;

    //    std::cout << "Mat * eig_vec            = "; (matcopy2 * eigen_vec).Print(std::cout,width,prec); std::cout << std::endl;
    //    std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";            
    //}
    //
    //std::cout << "    const static inline MML::Vector     " << matName << "_rhs0_sol{";
    //for(int i=0; i<dim; i++ )
    //{
    //    std::cout << std::setprecision(15) << vecSol[i];
    //    if( i < dim-1 )
    //        std::cout << ", ";
    //    else
    //        std::cout << "};\n\n";
    //}

    //std::cout << "    const static inline Vector<Real> " << matName << "_eigen_val{ " << eigen_solver.getEigenvalues()[0] << ", ";
    //for(int i=1; i<dim; i++ )
    //{
    //    std::cout << eigen_solver.getEigenvalues()[i];
    //    if( i < dim-1 )
    //        std::cout << ", ";
    //    else
    //        std::cout << " };\n";
    //}

    //std::cout << "    const static inline std::vector<Vector<Real>> " << matName << "_eigen_vecs\n";
    //std::cout << std::setprecision(10) << "    {\n";
    //for(int i=0; i<dim; i++ )
    //{
    //    std::cout << "        Vector<Real>{ " << eigen_solver.getEigenvector(i)[0] << ", ";
    //    for(int j=1; j<dim; j++ )
    //    {
    //        std::cout  << eigen_solver.getEigenvector(i)[j];
    //        if( j < dim-1 )
    //            std::cout << ", ";
    //        else
    //            std::cout << " },\n";
    //    }
    //}
    //std::cout << "    };\n";
}


void ComplexMatrixAnalyzer(std::string matName, const Matrix<Complex> &origMat, const Vector<Complex> &rhs)
{
 /*   std::cout << "*******************************************************************\n";
    std::cout << "*****   COMPLEX MATRIX ANALYZER FOR MATRIX   " << origMat.rows() << "x" << origMat.cols() << "   *****\n";

    int width = 15;
    int prec  = 10;

    int dim = origMat.rows();
    Matrix<Complex>  matcopy(origMat);
    Vector<Complex>  rhscopy(rhs);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,15,3);
    std::cout << "Right side:\n";        rhscopy.Print(std::cout,15,3);

    Vector<Complex>     vecSol(rhs.size());

    LUDecompositionSolver<Complex> luSolver(matcopy);
    luSolver.Solve(rhscopy, vecSol);

    std::cout << "\nSolution:\n" << vecSol << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * vecSol << std::endl;

    Matrix<Complex> inv(rhs.size(),rhs.size());
    luSolver.inverse(inv);

    std::cout << "Inverse:\n";
    inv.Print(std::cout, 20, 5);
    std::cout << "\nOrig * inv:\n";
    (origMat * inv).Print(std::cout, 20, 5, 1e-10);
    
    std::cout << "\n    const static inline MML::Vector<Complex>  " << matName << "_rhs0_sol{ ";
    for(int i=0; i<dim; i++ )
    {
        std::cout << "Complex";
        std::cout << std::setprecision(15) << vecSol[i];
        if( i < dim-1 )
            std::cout << ", \n";
        else
            std::cout << "};\n\n";
    }

    std::cout << "    };\n";*/
}


void TestMatrix(const std::string& name, const Matrix<Real>& mat, const Vector<Real>& rhs, const Vector<Real>& expected_sol)
{
    LUSolver<Real> solver(mat);
    Vector<Real> computed_sol = solver.Solve(rhs);
    
    Vector<Real> residual = mat * computed_sol - rhs;
    double residual_norm = residual.NormL2();
    
    bool matches = computed_sol.IsEqualTo(expected_sol, 1e-13);
    
    std::cout << name << ":\n";
    std::cout << "  Residual norm: " << std::scientific << std::setprecision(6) << residual_norm << "\n";
    std::cout << "  Matches expected (tol=1e-13): " << (matches ? "YES" : "NO") << "\n";
    
    if (!matches) {
        std::cout << "  Expected solution:\n";
        for (int i = 0; i < expected_sol.size(); i++)
            std::cout << "    [" << i << "] " << std::setprecision(18) << std::fixed << expected_sol[i] << "\n";
        std::cout << "  Computed solution:\n";
        for (int i = 0; i < computed_sol.size(); i++)
            std::cout << "    [" << i << "] " << std::setprecision(18) << std::fixed << computed_sol[i] << "\n";
    }
    std::cout << "\n";
}

/******************************************************************************
 * INTEGRATION TEST DATA GENERATOR
 * 
 * Generates binary test data files for large linear system integration tests.
 * Creates RANDOM matrices (not diagonally dominant!) with known solutions.
 *****************************************************************************/

/**
 * @brief Generate a random matrix using deterministic seeding
 * Uses a simple LCG for reproducibility across platforms
 */
Matrix<Real> GenerateRandomMatrix(int n, unsigned int seed = 12345) {
    Matrix<Real> A(n, n);
    
    // Simple LCG for reproducible random numbers
    unsigned int state = seed;
    auto nextRandom = [&state]() -> Real {
        state = state * 1103515245 + 12345;
        // Map to [-1, 1] range
        return Real(2.0) * (Real(state & 0x7FFFFFFF) / Real(0x7FFFFFFF)) - Real(1.0);
    };
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A(i, j) = nextRandom();
        }
    }
    
    return A;
}

/**
 * @brief Generate a random solution vector (deterministic)
 */
Vector<Real> GenerateRandomSolution(int n, unsigned int seed = 54321) {
    Vector<Real> x(n);
    
    unsigned int state = seed;
    auto nextRandom = [&state]() -> Real {
        state = state * 1103515245 + 12345;
        return Real(2.0) * (Real(state & 0x7FFFFFFF) / Real(0x7FFFFFFF)) - Real(1.0);
    };
    
    for (int i = 0; i < n; i++) {
        x[i] = nextRandom();
    }
    return x;
}

/**
 * @brief Generate and save a complete linear system test case
 * 
 * Saves three files:
 *   - <basename>_matrix.mmlbin  - coefficient matrix A (RANDOM!)
 *   - <basename>_rhs.mmlbin     - right-hand side b = A*x
 *   - <basename>_solution.mmlbin - known solution x
 * 
 * @param n Matrix dimension
 * @param outputDir Output directory path
 * @param basename Base name for files (e.g., "system_1000")
 * @param matrixSeed Seed for matrix generation
 * @param solutionSeed Seed for solution vector generation
 * @return true if all files saved successfully
 */
bool GenerateLinearSystemTestData(int n, const std::string& outputDir, const std::string& basename,
                                   unsigned int matrixSeed = 12345, unsigned int solutionSeed = 54321) {
    std::cout << "Generating " << n << "x" << n << " RANDOM linear system test data...\n";
    std::cout << "  (Matrix seed: " << matrixSeed << ", Solution seed: " << solutionSeed << ")\n";
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Generate RANDOM matrix (not diagonally dominant!)
    std::cout << "  Creating random matrix..." << std::flush;
    Matrix<Real> A = GenerateRandomMatrix(n, matrixSeed);
    auto matrixTime = std::chrono::high_resolution_clock::now();
    std::cout << " done (" 
              << std::chrono::duration_cast<std::chrono::milliseconds>(matrixTime - startTime).count() 
              << " ms)\n";
    
    // Generate random solution
    std::cout << "  Creating random solution vector..." << std::flush;
    Vector<Real> x = GenerateRandomSolution(n, solutionSeed);
    std::cout << " done\n";
    
    // Compute RHS: b = A * x
    std::cout << "  Computing RHS vector (b = A*x)..." << std::flush;
    Vector<Real> b = A * x;
    auto rhsTime = std::chrono::high_resolution_clock::now();
    std::cout << " done (" 
              << std::chrono::duration_cast<std::chrono::milliseconds>(rhsTime - matrixTime).count() 
              << " ms)\n";
    
    // Create output directory if needed
    std::filesystem::create_directories(outputDir);
    
    // Save files
    std::string matrixFile = outputDir + "/" + basename + "_matrix.mmlbin";
    std::string rhsFile = outputDir + "/" + basename + "_rhs.mmlbin";
    std::string solutionFile = outputDir + "/" + basename + "_solution.mmlbin";
    
    std::cout << "  Saving matrix to " << matrixFile << "..." << std::flush;
    auto res = SaveMatrixToBinary(A, matrixFile);

    if (!res.success) {
        std::cerr << " FAILED!\n";
        return false;
    }
    std::cout << " done (" << std::filesystem::file_size(matrixFile) / (1024*1024) << " MB)\n";
    
    std::cout << "  Saving RHS to " << rhsFile << "..." << std::flush;
    if (!SaveBinary(b, rhsFile).success) {
        std::cerr << " FAILED!\n";
        return false;
    }
    std::cout << " done\n";
    
    std::cout << "  Saving solution to " << solutionFile << "..." << std::flush;
    if (!SaveBinary(x, solutionFile).success) {
        std::cerr << " FAILED!\n";
        return false;
    }
    std::cout << " done\n";
    
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "  Total time: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() 
              << " ms\n\n";
    
    return true;
}

/******************************************************************************
 * EIGENVALUE TEST DATA GENERATOR
 * 
 * Generates symmetric matrices and computes reference eigenvalues using
 * Jacobi solver (high accuracy for symmetric matrices).
 *****************************************************************************/

/**
 * @brief Generate a random symmetric matrix using deterministic seeding
 */
MatrixSym<Real> GenerateRandomSymmetricMatrix(int n, unsigned int seed = 12345) {
    MatrixSym<Real> A(n);
    
    // Simple LCG for reproducible random numbers
    unsigned int state = seed;
    auto nextRandom = [&state]() -> Real {
        state = state * 1103515245 + 12345;
        return Real(2.0) * (Real(state & 0x7FFFFFFF) / Real(0x7FFFFFFF)) - Real(1.0);
    };
    
    // Fill upper triangle (symmetric matrix stores only upper)
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            A(i, j) = nextRandom();
        }
    }
    
    return A;
}

/**
 * @brief Generate and save eigenvalue test case for symmetric matrix
 * 
 * Saves two files:
 *   - <basename>_matrix.mmlbin  - symmetric matrix (stored as full matrix)
 *   - <basename>_eigenvalues.mmlb - complex eigenvalues (real for symmetric)
 * 
 * Uses Jacobi solver as reference (best accuracy for symmetric matrices)
 */
bool GenerateEigenvalueTestData(int n, const std::string& outputDir, const std::string& basename,
                                 unsigned int seed = 12345) {
    std::cout << "Generating " << n << "x" << n << " symmetric eigenvalue test data...\n";
    std::cout << "  (Seed: " << seed << ")\n";
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Generate random symmetric matrix
    std::cout << "  Creating random symmetric matrix..." << std::flush;
    MatrixSym<Real> Asym = GenerateRandomSymmetricMatrix(n, seed);
    auto matrixTime = std::chrono::high_resolution_clock::now();
    std::cout << " done (" 
              << std::chrono::duration_cast<std::chrono::milliseconds>(matrixTime - startTime).count() 
              << " ms)\n";
    
    // Convert to full matrix for saving and for solver
    Matrix<Real> A = Asym.ToFullMatrix();
    
    // Compute reference eigenvalues using Jacobi solver
    // Note: Use looser tolerance for initial convergence
    std::cout << "  Computing eigenvalues (Jacobi)..." << std::flush;
    EigenSolverConfig config;
    config.tolerance = 1e-10;  // Reasonable accuracy for reference
    config.max_iterations = 100;  // Standard Jacobi sweeps
    config.sort_eigenvalues = true;
    
    auto result = SymmMatEigenSolverJacobi::Solve(A, config);  // Use Matrix directly
    auto eigenTime = std::chrono::high_resolution_clock::now();
    
    if (!result.converged) {
        std::cerr << " FAILED (did not converge after " << result.iterations << " sweeps)!\n";
        return false;
    }
    std::cout << " done (" 
              << std::chrono::duration_cast<std::chrono::milliseconds>(eigenTime - matrixTime).count() 
              << " ms, " << result.iterations << " sweeps)\n";
    
    // Create output directory if needed
    std::filesystem::create_directories(outputDir);
    
    // Save matrix
    std::string matrixFile = outputDir + "/" + basename + "_matrix.mmlbin";
    std::cout << "  Saving matrix to " << matrixFile << "..." << std::flush;
    auto res = SaveMatrixToBinary(A, matrixFile);
    if (!res.success) {
        std::cerr << " FAILED!\n";
        return false;
    }
    std::cout << " done (" << std::filesystem::file_size(matrixFile) / 1024 << " KB)\n";
    
    // Save eigenvalues as complex vector (real for symmetric - imag parts are zero)
    std::string eigenFile = outputDir + "/" + basename + "_eigenvalues.mmlb";
    std::cout << "  Saving eigenvalues to " << eigenFile << "..." << std::flush;
    auto resEigen = Persistence::SaveVectorAsComplex(result.eigenvalues, eigenFile);
    if (!resEigen.success) {
        std::cerr << " FAILED!\n";
        return false;
    }
    std::cout << " done (" << n << " eigenvalues, sorted ascending)\n";
    
    // Print some stats
    std::cout << "  Eigenvalue range: [" << result.eigenvalues[0] << ", " 
              << result.eigenvalues[n-1] << "]\n";
    
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "  Total time: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() 
              << " ms\n\n";
    
    return true;
}

/**
 * @brief Generate all integration test linear system data files
 */
void GenerateAllIntegrationTestData(const std::string& baseDir) {
    std::string linearDir = baseDir + "/integration_tests/linear_systems";
    std::string eigenDir = baseDir + "/integration_tests/eigenvalues";
    
    std::cout << "=======================================================\n";
    std::cout << "MML Integration Test Data Generator\n";
    std::cout << "=======================================================\n\n";
    
    // Linear systems
    std::cout << "--- LINEAR SYSTEM TEST DATA ---\n";
    std::cout << "Output directory: " << linearDir << "\n";
    std::cout << "Generating RANDOM matrices (NOT diagonally dominant!)\n\n";
    
    // Generate 1000x1000 system (~8 MB) with unique seeds
    GenerateLinearSystemTestData(1000, linearDir, "system_1000", 1000001, 1000002);
    
    // Generate 2000x2000 system (~32 MB) with different seeds
    GenerateLinearSystemTestData(2000, linearDir, "system_2000", 2000001, 2000002);
    
    // Eigenvalue systems
    std::cout << "--- EIGENVALUE TEST DATA ---\n";
    std::cout << "Output directory: " << eigenDir << "\n";
    std::cout << "Generating SYMMETRIC matrices with Jacobi-computed eigenvalues\n\n";
    
    // Generate symmetric matrices for eigenvalue testing
    GenerateEigenvalueTestData(200, eigenDir, "symm_200", 200001);
    GenerateEigenvalueTestData(300, eigenDir, "symm_300", 300001);
    GenerateEigenvalueTestData(500, eigenDir, "symm_500", 500001);
    
    std::cout << "=======================================================\n";
    std::cout << "All test data generated successfully!\n";
    std::cout << "=======================================================\n";
}

void PrintUsage(const char* programName) {
    std::cout << "MML Utilities Application\n\n";
    std::cout << "Usage:\n";
    std::cout << "  " << programName << " generate-test-data <base_dir>\n";
    std::cout << "      Generate integration test data files\n";
    std::cout << "      <base_dir> is the MML project root directory\n\n";
    std::cout << "  " << programName << " eigen-test\n";
    std::cout << "      Run eigenvalue test on random 100x100 matrix\n\n";
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        PrintUsage(argv[0]);
        return 1;
    }
    
    std::string command = argv[1];
    
    // Command: generate-test-data <base_dir>
    if (command == "generate-test-data") {
        if (argc < 3) {
            std::cerr << "Error: Missing base directory argument\n\n";
            PrintUsage(argv[0]);
            return 1;
        }
        std::string baseDir = argv[2];
        GenerateAllIntegrationTestData(baseDir);
        return 0;
    }
    
    // Command: eigen-test (original functionality)
    if (command == "eigen-test") {
        int n = 100;
        auto mat = GenerateRandomMatrix(n, n, -5, 5);
        
        std::cout << "Testing eigenvalue solver on " << n << "x" << n << " random matrix...\n";
        std::cout << mat;
        
        EigenSolver eigenSolver;
        auto res = eigenSolver.Solve(mat);
        
        std::cout << "\nEigenvalues:\n";
        for(int i = 0; i < n; i++)
            std::cout << res.eigenvalues[i].real << " + i * " << res.eigenvalues[i].imag << std::endl;
        
        return 0;
    }
    
    // Command: test-matrix (original commented-out functionality)
    if (command == "test-matrix") {
        std::cout << "Testing classic matrix solutions...\n\n";
        
        TestMatrix("hilbert_3x3", TestBeds::hilbert_3x3(), TestBeds::hilbert_3x3_rhs0(), TestBeds::hilbert_3x3_rhs0_sol());
        TestMatrix("hilbert_4x4", TestBeds::hilbert_4x4(), TestBeds::hilbert_4x4_rhs0(), TestBeds::hilbert_4x4_rhs0_sol());
        TestMatrix("hilbert_5x5", TestBeds::hilbert_5x5(), TestBeds::hilbert_5x5_rhs0(), TestBeds::hilbert_5x5_rhs0_sol());
        TestMatrix("hilbert_8x8", TestBeds::hilbert_8x8(), TestBeds::hilbert_8x8_rhs0(), TestBeds::hilbert_8x8_rhs0_sol());
        
        TestMatrix("vandermonde_3x3", TestBeds::vandermonde_3x3(), TestBeds::vandermonde_3x3_rhs0(), TestBeds::vandermonde_3x3_rhs0_sol());
        TestMatrix("vandermonde_4x4", TestBeds::vandermonde_4x4(), TestBeds::vandermonde_4x4_rhs0(), TestBeds::vandermonde_4x4_rhs0_sol());
        TestMatrix("vandermonde_5x5", TestBeds::vandermonde_5x5(), TestBeds::vandermonde_5x5_rhs0(), TestBeds::vandermonde_5x5_rhs0_sol());
        
        TestMatrix("kahan_3x3", TestBeds::kahan_3x3(), TestBeds::kahan_3x3_rhs0(), TestBeds::kahan_3x3_rhs0_sol());
        TestMatrix("kahan_5x5", TestBeds::kahan_5x5(), TestBeds::kahan_5x5_rhs0(), TestBeds::kahan_5x5_rhs0_sol());
        
        return 0;
    }
    
    std::cerr << "Unknown command: " << command << "\n\n";
    PrintUsage(argv[0]);
    return 1;
}
