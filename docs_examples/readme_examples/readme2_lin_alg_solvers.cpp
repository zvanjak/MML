#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/LinAlgEqSolvers.h"

#include "algorithms/EigenSystemSolvers.h"
#endif
#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;


void Readme_linear_system_solvers()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                  README - ling.alg. solvers                   ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Matrix<Real>    mat{ 5, 5, {  0.2,  4.1, -2.1, -7.4,  1.6,
																 1.6,  1.5, -1.1,  0.7,  5.0,
																-3.8, -8.0,  9.6, -5.4, -7.8,
																 4.6, -8.2,  8.4,  0.4,  8.0,
																-2.6,  2.9,  0.1, -9.6, -2.7 } };
	Vector<Real> 	rhs{ 1.1, 4.7, 0.1, 9.3, 0.4 };

	LUDecompositionSolver<Real> luSolver(mat);

	Vector<Real>	vecSol = luSolver.Solve(rhs);

	std::cout << "Solution   = " << vecSol << std::endl;
	std::cout << "Right side = "; rhs.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "Mat * sol  = "; (mat * vecSol).Print(std::cout, 8, 4); std::cout << std::endl;

	Matrix<Real>  matcopy(mat);

	EigenSolver eigenSolver(matcopy, true, false);

	std::cout << "\nNum real eigenvalues    : " << eigenSolver.getNumReal();
	std::cout << "\nNum complex eigenvalues : " << eigenSolver.getNumComplex() << "\n";

	std::cout << "\nEigenvalues : "; eigenSolver.getEigenvalues().Print(std::cout, 10, 5);
	std::cout << "\nReal        : "; eigenSolver.getRealEigenvalues().Print(std::cout, 15, 10);
	std::cout << "\nComplex     : "; eigenSolver.getComplexEigenvalues().Print(std::cout, 15, 10);

	std::cout << "\n";

	/* OUTPUT
			Solution   = [   -5.568500786,    -5.944693206,    -5.007620645,    -1.393638021,     3.598760994]
			Right side = [     1.1,      4.7,      0.1,      9.3,      0.4]
			Mat * sol  = [     1.1,      4.7,      0.1,      9.3,      0.4]

			Num real eigenvalues    : 3
			Num complex eigenvalues : 2

			Eigenvalues : [(12.974,0), (0.99944,0), (-0.033184,0), (-2.4701,12.994), (-2.4701,-12.994)]
			Real        : [    12.97392154,    0.9994371124,  -0.03318390189]
			Complex     : [(-2.470087376,12.99433106), (-2.470087376,-12.99433106)]
	*/
}