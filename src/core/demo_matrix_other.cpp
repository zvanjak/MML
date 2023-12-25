#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/MatrixSym.h"
#include "core/MatrixBandDiag.h"

#include "core/LinAlgEqSolvers.h"
#endif

using namespace MML;

// TODO - BIG, EMPTY!!! - Demo matrix others

void Demo_Matrix_Sym()
{
    MatrixSym<Real> a(3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    
    a.Print(std::cout, 10, 3);
}

void Demo_Matrix_Tridiag()
{
    TridiagonalMatrix<Real> a(4, {0.0, 4.5, 9.0, 10.0}, {4.0, 1.5, 6.0, 7.0}, {1.0, 2.0, 3.0, 0.0}  );
    TridiagonalMatrix<Real> b(4, {      1.0, 2.0, 3.0,
                                   4.0, 1.5, 6.0, 7.0, 
                                        4.5, 9.0, 10.0 });
    // Dim: 4
    // [          4,          1,          0,          0,  ]e
    // [        4.5,        1.5,          2,          0,  ]
    // [          0,          9,          6,          3,  ]
    // [          0,          0,         10,          7,  ]

    Matrix<Real> c(4, 4, {  4.0, 1.0, 0.0, 0.0,
                            4.5, 1.5, 2.0, 0.0, 
                            0.0, 9.0, 6.0, 3.0,
                            0.0, 0.0, 10.0, 7.0 });   
    Vector<Real> rhs{1.0, 2.0, 3.0, 4.0};

    a.Print(std::cout, 10, 3);
    b.Print(std::cout, 10, 3);
    std::cout << "rhs: " << rhs << std::endl;

    Vector<Real> sol_a(4);
    Vector<Real> sol_b(4);

    a.solve(rhs, sol_a);
    b.solve(rhs, sol_b);

    std::cout << "sol_a: " << sol_a << std::endl;
    std::cout << "sol_b: " << sol_b << std::endl;
    
    std::cout << "Tridiag sol_a mul" << c * sol_a << std::endl;
    std::cout << "Tridiag sol_b mul" << c * sol_b << std::endl;
    
    
    c.Print(std::cout, 10, 3);
    Vector<Real> solLU(4);

    LUDecompositionSolver<Real> solver(c); 
    solver.Solve(rhs, solLU);
    std::cout << "LU sol: " << solLU << std::endl;

    std::cout << "Tridiag solLU mul" << c * solLU << std::endl;
}

void Demo_Matrix_BandDiag()
{
}

void Demo_Matrix_Other()
{
    // Demo_Matrix_Sym();
    Demo_Matrix_Tridiag();
    Demo_Matrix_BandDiag();
}
