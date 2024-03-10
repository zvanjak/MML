#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/LinAlgEqSolvers.h"

#include "algorithms/MatrixAlg.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

void Docs_Demo_LinAlgSolvers_GaussJordan()
{
	std::cout << "SOLVING SYSTEMS VIA GAUSS-JORDAN ELIMINATION:\n";

	Matrix<Real>     origMat = TestBeds::mat_5x5;
	Matrix<Real>     matcopy(origMat);
	Vector<Real>     rhscopy(TestBeds::mat_5x5_rhs0);

	std::cout << "Initial matrix:\n";    matcopy.Print(std::cout, 10, 3);
	std::cout << "Right side:\n";        rhscopy.Print(std::cout, 10, 3);

	GaussJordanSolver<Real>::Solve(matcopy, rhscopy);

	std::cout << "\nSolution:\n" << rhscopy << std::endl;
	std::cout << "Mat * sol: \n"; (origMat * rhscopy).Print(std::cout, 10, 3);

	// matcopy now containes inverted matrix
	std::cout << "\nInverse:\n" << matcopy << std::endl;
	std::cout << "Orig * inv:\n" << origMat * matcopy << std::endl;

	std::cout << "SOLVING MULTIPLE SYSTEMS VIA GAUSS-JORDAN ELIMINATION:\n";
	Matrix<Real>     origMat2 = TestBeds::mat_5x5;
	Matrix<Real>     matcopy2(origMat2);
	Matrix<Real>     rhscopy2(TestBeds::mat_5x5_rhs_multi);

	std::cout << "Initial matrix:\n";    matcopy2.Print(std::cout, 10, 3);
	std::cout << "Right side:\n";        rhscopy2.Print(std::cout, 10, 3);

	GaussJordanSolver<Real>::Solve(matcopy2, rhscopy2);

	std::cout << "Solution  :\n" << rhscopy2 << std::endl;

	/* OUTPUT
	SOLVING SYSTEMS VIA GAUSS-JORDAN ELIMINATION:
	Initial matrix:
	Rows: 5 Cols: 5
	[        1.4,        2.1,        2.1,        7.4,        9.6,  ]
	[        1.6,        1.5,        1.1,        0.7,          5,  ]
	[        3.8,          8,        9.6,        5.4,        8.8,  ]
	[        4.6,        8.2,        8.4,        0.4,          8,  ]
	[        2.6,        2.9,        0.1,        9.6,        7.7,  ]
	Right side:
	[       1.1,        4.7,        0.1,        9.3,        0.4]
	Solution:
	[   -3.903271042,     5.235343316,     -3.29209577,    -1.718330011,      1.58327101]
	Mat * sol:
	[       1.1,        4.7,        0.1,        9.3,        0.4]
	Inverse:
	Rows: 5 Cols: 5
	[      -1.68,       2.71,       1.47,      -1.62,      0.338,  ]
	[       1.16,      -2.37,      -1.42,       1.64,    0.00331,  ]
	[     -0.517,       1.01,      0.809,     -0.807,    -0.0967,  ]
	[      -0.15,      0.129,      0.226,     -0.238,     0.0936,  ]
	[      0.324,     -0.197,     -0.256,      0.239,     -0.101,  ]

	Orig * inv:
	Rows: 5 Cols: 5
	[          1,   1.11e-15,   4.44e-16,  -4.44e-16,   4.44e-16,  ]
	[  -2.22e-16,          1,          0,          0,   1.11e-16,  ]
	[  -8.88e-16,   2.22e-15,          1,          0,   3.33e-16,  ]
	[  -4.44e-16,  -1.11e-15,   4.44e-16,          1,   3.33e-16,  ]
	[  -4.44e-16,   6.66e-16,  -4.44e-16,  -2.22e-16,          1,  ]

	SOLVING MULTIPLE SYSTEMS VIA GAUSS-JORDAN ELIMINATION:
	Initial matrix:
	Rows: 5 Cols: 5
	[        1.4,        2.1,        2.1,        7.4,        9.6,  ]
	[        1.6,        1.5,        1.1,        0.7,          5,  ]
	[        3.8,          8,        9.6,        5.4,        8.8,  ]
	[        4.6,        8.2,        8.4,        0.4,          8,  ]
	[        2.6,        2.9,        0.1,        9.6,        7.7,  ]
	Right side:
	Rows: 5 Cols: 2
	[        1.1,        1.6,  ]
	[        4.7,        9.1,  ]
	[        0.1,          4,  ]
	[        9.3,        8.4,  ]
	[        0.4,        4.1,  ]
	Solution  :
	Rows: 5 Cols: 2
	[       -3.9,       15.6,  ]
	[       5.24,      -11.6,  ]
	[      -3.29,       4.41,  ]
	[      -1.72,      0.214,  ]
	[       1.58,      -0.71,  ]
	*/
}

void Docs_Demo_LinAlgSolvers_LU_decomposition()
{
	std::cout << "\nSOLVING VIA LU DECOMPOSITION:\n";

	Matrix<Real>  mat = TestBeds::mat_5x5;
	Vector<Real>  rhs(TestBeds::mat_5x5_rhs0);

	std::cout << "Initial matrix:\n";    mat.Print(std::cout, 10, 3);
	std::cout << "Right side:\n";        rhs.Print(std::cout, 10, 3);

	LUDecompositionSolver<Real> luSolver(mat);

	Vector<Real>	vecSol = luSolver.Solve(rhs);

	std::cout << "\nSolution: " << vecSol << std::endl;
	std::cout << "Mat * sol :\n "; (mat * vecSol).Print(std::cout, 8, 4); std::cout << std::endl;

	Matrix<Real> inv(5, 5);
	luSolver.inverse(inv);

	std::cout << "Inverse:\n" << inv << std::endl;
	std::cout << "Orig * inv:\n" << mat * inv << std::endl;

	std::cout << "\nSOLVING MULTIPLE SYSTEMS VIA LU DECOMPOSITION:\n";

	Matrix<Real>     origMat = TestBeds::mat_5x5;
	Matrix<Real>     matcopy(origMat);
	Matrix<Real>     rhscopy(TestBeds::mat_5x5_rhs_multi);

	std::cout << "Initial matrix:\n";    matcopy.Print(std::cout, 10, 3);
	std::cout << "Right side:\n";        rhscopy.Print(std::cout, 10, 3);

	Matrix<Real>     matSol(5, 2);

	LUDecompositionSolver<Real> luSolver2(matcopy);
	luSolver2.Solve(rhscopy, matSol);

	std::cout << "Solution:\n" << matSol << std::endl;

	/* OUTPUT
	SOLVING VIA LU DECOMPOSITION:
	Initial matrix:
	Rows: 5 Cols: 5
	[        1.4,        2.1,        2.1,        7.4,        9.6,  ]
	[        1.6,        1.5,        1.1,        0.7,          5,  ]
	[        3.8,          8,        9.6,        5.4,        8.8,  ]
	[        4.6,        8.2,        8.4,        0.4,          8,  ]
	[        2.6,        2.9,        0.1,        9.6,        7.7,  ]
	Right side:
	[       1.1,        4.7,        0.1,        9.3,        0.4]
	Solution: [   -3.903271042,     5.235343316,     -3.29209577,    -1.718330011,      1.58327101]
	Mat * sol :
	 [     1.1,      4.7,      0.1,      9.3,      0.4]
	Inverse:
	Rows: 5 Cols: 5
	[      -1.68,       2.71,       1.47,      -1.62,      0.338,  ]
	[       1.16,      -2.37,      -1.42,       1.64,    0.00331,  ]
	[     -0.517,       1.01,      0.809,     -0.807,    -0.0967,  ]
	[      -0.15,      0.129,      0.226,     -0.238,     0.0936,  ]
	[      0.324,     -0.197,     -0.256,      0.239,     -0.101,  ]

	Orig * inv:
	Rows: 5 Cols: 5
	[          1,          0,   4.44e-16,  -4.44e-16,   1.11e-16,  ]
	[   2.22e-16,          1,  -2.22e-16,  -2.22e-16,          0,  ]
	[  -4.44e-16,          0,          1,  -4.44e-16,   2.22e-16,  ]
	[   4.44e-16,   4.44e-16,          0,          1,  -1.11e-16,  ]
	[   8.88e-16,  -4.44e-16,  -6.66e-16,  -4.44e-16,          1,  ]


	SOLVING MULTIPLE SYSTEMS VIA LU DECOMPOSITION:
	Initial matrix:
	Rows: 5 Cols: 5
	[        1.4,        2.1,        2.1,        7.4,        9.6,  ]
	[        1.6,        1.5,        1.1,        0.7,          5,  ]
	[        3.8,          8,        9.6,        5.4,        8.8,  ]
	[        4.6,        8.2,        8.4,        0.4,          8,  ]
	[        2.6,        2.9,        0.1,        9.6,        7.7,  ]
	Right side:
	Rows: 5 Cols: 2
	[        1.1,        1.6,  ]
	[        4.7,        9.1,  ]
	[        0.1,          4,  ]
	[        9.3,        8.4,  ]
	[        0.4,        4.1,  ]
	Solution:
	Rows: 5 Cols: 2
	[       -3.9,       15.6,  ]
	[       5.24,      -11.6,  ]
	[      -3.29,       4.41,  ]
	[      -1.72,      0.214,  ]
	[       1.58,      -0.71,  ]
	*/
}

void Docs_Demo_LinAlgSolvers_QR_decomposition()
{
	std::cout << "\nSOLVING VIA QR DECOMPOSITION:\n";

	Matrix<Real>   mat4(5, 5, { 1.4, 2.1, 2.1, 7.4, 9.6,
															1.6, 1.5, 1.1, 0.7, 5.0,
															3.8, 8.0, 9.6, 5.4, 8.8,
															4.6, 8.2, 8.4, 0.4, 8.0,
															2.6, 2.9, 0.1, 9.6, 7.7 });

	Vector<Real>     solVec(5), vecrhs4{ 1.1, 4.7, 0.1, 9.3, 0.4 };

	Matrix<Real>     origMat = mat4;
	Matrix<Real>     matcopy(origMat);

	std::cout << "Initial matrix:\n";    matcopy.Print(std::cout, 10, 3);
	std::cout << "Right side:\n" << vecrhs4 << std::endl;

	QRDecompositionSolver qrSolver(matcopy);
	qrSolver.Solve(vecrhs4, solVec);

	std::cout << "Solution:\n" << solVec << std::endl;
	std::cout << "Mat * sol :\n" << origMat * solVec << std::endl;

	/* OUTPUT
	SOLVING VIA QR DECOMPOSITION:
Initial matrix:
Rows: 5 Cols: 5
[        1.4,        2.1,        2.1,        7.4,        9.6,  ]
[        1.6,        1.5,        1.1,        0.7,          5,  ]
[        3.8,          8,        9.6,        5.4,        8.8,  ]
[        4.6,        8.2,        8.4,        0.4,          8,  ]
[        2.6,        2.9,        0.1,        9.6,        7.7,  ]
Right side:
[            1.1,             4.7,             0.1,             9.3,             0.4]
Solution:
[   -3.903271042,     5.235343316,     -3.29209577,    -1.718330011,      1.58327101]
Mat * sol :
[            1.1,             4.7,             0.1,             9.3,             0.4]
	*/
}

void Docs_Demo_LinAlgSolvers_SVD_decomposition()
{
	std::cout << "\nSOLVING VIA SVD DECOMPOSITION:\n";

	Matrix<Real>     mat4(5, 5, { 1.4, 2.1, 2.1, 7.4, 9.6,
																1.6, 1.5, 1.1, 0.7, 5.0,
																3.8, 8.0, 9.6, 5.4, 8.8,
																4.6, 8.2, 8.4, 0.4, 8.0,
																2.6, 2.9, 0.1, 9.6, 7.7 });


	Vector<Real>     vecrhs4{ 1.1, 4.7, 0.1, 9.3, 0.4 }, solVec(5);

	Matrix<Real>     origMat = mat4;
	Matrix<Real>     matcopy(origMat);

	std::cout << "Initial matrix:\n";    matcopy.Print(std::cout, 10, 3);
	std::cout << "Right side:\n";        std::cout << vecrhs4 << std::endl;

	SVDecompositionSolver svdSolver(matcopy);
	svdSolver.Solve(vecrhs4, solVec);

	std::cout << "Solution:\n" << solVec << std::endl;
	std::cout << "Mat * sol :\n" << origMat * solVec << std::endl;

	/* OUTPUT
	SOLVING VIA SVD DECOMPOSITION:
	Initial matrix:
	Rows: 5 Cols: 5
	[        1.4,        2.1,        2.1,        7.4,        9.6,  ]
	[        1.6,        1.5,        1.1,        0.7,          5,  ]
	[        3.8,          8,        9.6,        5.4,        8.8,  ]
	[        4.6,        8.2,        8.4,        0.4,          8,  ]
	[        2.6,        2.9,        0.1,        9.6,        7.7,  ]
	Right side:
	[            1.1,             4.7,             0.1,             9.3,             0.4]
	Solution:
	[   -3.903271042,     5.235343316,     -3.29209577,    -1.718330011,      1.58327101]
	Mat * sol :
	[            1.1,             4.7,             0.1,             9.3,             0.4]
	*/
}

void Docs_Demo_LinAlgSolvers_Cholesky_decomposition()
{
	std::cout << "\nSOLVING VIA Cholesky DECOMPOSITION:\n";

	// example positive definite matrix
	Matrix<Real>  mat4{ 3, 3, { 2.0, -1.0,  0.0,
														 -1.0,  3.0, -2.0,
															0.0, -2.0,  4.0 } };

	Vector<Real>     solVec(3), vecrhs4{ 1.1, 4.7, 0.1 };

	Matrix<Real>     origMat = mat4;
	Matrix<Real>     matcopy(origMat);

	std::cout << "Initial matrix:\n";    matcopy.Print(std::cout, 10, 3);
	std::cout << "Is positive definite: " << MatrixUtils::IsPositiveDefinite(mat4) << std::endl;

	std::cout << "Right side:\n" << vecrhs4 << std::endl;

	CholeskyDecompositionSolver choleskySolver(matcopy);
	choleskySolver.Solve(vecrhs4, solVec);

	std::cout << "Solution:\n" << solVec << std::endl;
	std::cout << "Multiplying solution with matrix: " << origMat * solVec << std::endl;

	/* OUTPUT
SOLVING VIA Cholesky DECOMPOSITION:
Initial matrix:
Rows: 3 Cols: 3
[          2,         -1,          0,  ]
[         -1,          3,         -2,  ]
[          0,         -2,          4,  ]
Is positive definite: 1
Right side:
[            1.1,             4.7,             0.1]
Solution:
[    2.316666667,     3.533333333,     1.791666667]
Multiplying solution with matrix: [            1.1,             4.7,             0.1]
*/
}

void Docs_Demo_LinAlgSolvers()
{
	//Docs_Demo_LinAlgSolvers_GaussJordan();
	//Docs_Demo_LinAlgSolvers_LU_decomposition();
	Docs_Demo_LinAlgSolvers_QR_decomposition();
	Docs_Demo_LinAlgSolvers_SVD_decomposition();
	Docs_Demo_LinAlgSolvers_Cholesky_decomposition();
}