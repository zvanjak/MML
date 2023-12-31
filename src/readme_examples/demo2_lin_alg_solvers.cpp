#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "core/LinAlgEqSolvers.h"

#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;


void Readme_linear_system_solvers()
{
 	Matrix<Real>    mat{5, 5, { 1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7 } };
	Vector<Real> 	rhs{1.1, 4.7, 0.1, 9.3, 0.4};
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);

	Vector<Real>    res_rhs = mat * vecSol;
}