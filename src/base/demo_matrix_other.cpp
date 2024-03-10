#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/MatrixSym.h"
#include "base/MatrixBandDiag.h"

#include "core/LinAlgEqSolvers.h"
#endif

using namespace MML;

// TODO 0.9 - HIGH, - Demo matrix others

void Demo_Matrix_Sym()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******                  Matrix Sym                 *******" << std::endl;    
    std::cout << "***********************************************************" << std::endl;
    
    MatrixSym<Real> a(3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    
    a.Print(std::cout, 10, 3);

    MatrixSym<Real> b(3, {1.0, 
                          2.0, 3.0, 
                          4.0, 5.0, 6.0}); 

    b.Print(std::cout, 10, 3);

}

void Demo_Matrix_Tridiag()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******               Matrix Tridiag                *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;

    // initializing with 3 vectors
//    TridiagonalMatrix<Real> a(4, {0.0, 4.5, 9.0, 10.0}, {4.0, 1.5, 6.0, 7.0}, {1.0, 2.0, 3.0, 0.0}  );

    // initializing with values in single initializer list
    TridiagonalMatrix<Real> b(4, {  4.0, 1.0,
                                    4.5, 1.5,  2.0, 
                                         9.0,  6.0, 3.0,
                                              10.0, 7.0 });
    // Dim: 4
    // [          4,          1,          0,          0,  ]
    // [        4.5,        1.5,          2,          0,  ]
    // [          0,          9,          6,          3,  ]
    // [          0,          0,         10,          7,  ]

    Matrix<Real> c(4, 4, {  4.0, 1.0, 0.0, 0.0,
                            4.5, 1.5, 2.0, 0.0, 
                            0.0, 9.0, 6.0, 3.0,
                            0.0, 0.0, 10.0, 7.0 });   
    Vector<Real> rhs{1.0, 2.0, 3.0, 4.0};

    //a.Print(std::cout, 10, 3);
    b.Print(std::cout, 10, 3);
    std::cout << "rhs: " << rhs << std::endl;

    Vector<Real> sol_a(4);
    Vector<Real> sol_b(4);

    //a.Solve(rhs, sol_a);
    b.Solve(rhs, sol_b);

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
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******            Matrix Band diagonal             *******" << std::endl;
    std::cout << "***********************************************************" << std::endl;

    // Dim: 4
    // [          4,          1,          0,          0,  ]
    // [        4.5,        1.5,          2,          0,  ]
    // [          0,          9,          6,          3,  ]
    // [          0,          0,         10,          7,  ]
    Matrix<Real> mat(4, 3, { 0.0, 4.0, 1.0, 
                             4.5, 1.5, 2.0,
                             9.0, 6.0, 3.0,
                            10.0, 7.0, 0.0 } );
    Vector<Real> rhs{1.0, 2.0, 3.0, 4.0};
    Vector<Real> sol(4);

    BandDiagLUSolver solver(mat, 1, 1);
    solver.Solve(rhs, sol);
    
    std::cout << "BandDiag sol: " << sol << std::endl;
    //std::cout << "BandDiag sol mul: " << mat * sol << std::endl;

    BandDiagonalMatrix m(4, 1, 1, mat);
    m.Print(std::cout, 10, 3);

    Matrix<Real> mat2(7, 4, { 0.0, 0.0, 3.0, 1.0,
                              0.0, 4.0, 1.0, 5.0,
                              9.0, 2.0, 6.0, 5.0,
                              4.5, 1.5, 2.0, 5.0,
                              4.5, 1.5, 2.0, 4.0,
                              3.0, 8.0, 4.0, 6.0,
                              2.0, 4.0, 4.0, 0.0 } );

    // Dim: 7
    // [      3,      1,      0,      0,      0,      0,      0,  ]
    // [      4,      1,      5,      0,      0,      0,      0,  ]
    // [      9,      2,      6,      5,      0,      0,      0,  ]
    // [      0,      4,      2,      2,      5,      0,      0,  ]
    // [      0,      0,      4,      2,      2,      4,      0,  ]
    // [      0,      0,      0,      3,      8,      4,      6,  ]
    // [      0,      0,      0,      0,      2,      4,      4,  ]

    BandDiagonalMatrix m2(7, 2, 1, mat2);
    m2.Print(std::cout, 6, 1);
}

void Demo_Matrix_Other()
{
    Demo_Matrix_Sym();
    Demo_Matrix_Tridiag();
    Demo_Matrix_BandDiag();
}
