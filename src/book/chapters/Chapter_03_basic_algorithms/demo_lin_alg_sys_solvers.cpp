#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/core/LinAlgEqSolvers.h"

#endif

using namespace MML;


void Demo_Lin_alg_sys_solvers()
{
	// initialize a 4x4 matrix, defining the system of linear equations
	Matrix<Real> matOrig(4, 4, { -2,  4,  3,  5,
																3,  8, -7, -2,
																1, -3,  7,  4,
																5, -4,  9,  1 });
	Matrix<Real> mat1(matOrig), mat2(matOrig);
	std::cout << "Initial matrix:\n";  matOrig.Print(std::cout, 10, 3);

	// right side vector
	Vector<Real> rhs({ 3, -2, 1, 5 }), sol1(rhs);
	std::cout << "\nRight side: ";     rhs.Print(std::cout, 10, 3);

	// solution is returned in given rhs vector (in place solution)
	GaussJordanSolver<Real>::SolveInPlace(mat1, sol1);

	// solution is returned in a new vector (rhs is not modified)
	Vector<Real> sol2 = GaussJordanSolver<Real>::Solve(mat2, rhs);

	// solution is returned in a new vector, with the original matrix unchanged
	Vector<Real> sol3 = GaussJordanSolver<Real>::SolveConst(matOrig, rhs);

	std::cout << "\nSol1 = " << sol1 << std::endl;
	std::cout << "Sol2 = " << sol2 << std::endl;
	std::cout << "Sol3 = " << sol3 << std::endl;

	std::cout << "mat * sol1 = "; (matOrig * sol1).Print(std::cout, 10, 3);
	std::cout << "\nmat * sol2 = "; (matOrig * sol2).Print(std::cout, 10, 3);
	std::cout << "\nmat * sol3 = "; (matOrig * sol3).Print(std::cout, 10, 3);
}