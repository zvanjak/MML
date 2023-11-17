#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/MatrixSym.h"
#include "core/MatrixBandDiag.h"
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
}

void Demo_Matrix_BandDiag()
{
}

void Demo_Matrix_Other()
{
    Demo_Matrix_Sym();
    Demo_Matrix_Tridiag();
    Demo_Matrix_BandDiag();
}
