#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "core/Function.h"
#include "core/InterpolatedFunction.h"

#include "core/Integration.h"
#include "algorithms/FunctionAnalyzers.h"
#endif

using namespace std;
using namespace MML;


void Test_Interpolation_2DFunc()
{
    const int N=5;
    const double PI=3.141592653589793238;
    int i,j;
    Real dy,f,x1,x2,y;
    Vector<Real> x1a(N),x2a(N);
    Matrix<Real> ya(N,N);

    for (i=0;i<N;i++) {
        x1a[i]=(i+1)*PI/N;
        for (j=0;j<N;j++) {
            x2a[j]=1.0*(j+1)/N;
            ya[i][j]=sin(x1a[i])*exp(x2a[j]);
        }
    }

    BilinInterpScalarFunction2D bilin_int(x1a,x2a, ya);
    PolynomInterpScalarFunction2D poly_int(x1a,x2a,ya,3,3);
    SplineInterpScalarFunction2D spline_int(x1a,x2a,ya);

    // test 2-dimensional interpolation
    cout << endl << "Two dimensional interpolation of sin(x1)exp(x2)";
    cout << endl;
    cout << setw(9) << "x1" << setw(13) << "x2" << setw(14) << "f(x)";
    cout << "        polint2       error       Bilin         error        Poly        error        Spline        error " << endl;
    cout << fixed << setprecision(6);
    for (i=0;i<4;i++) {
        x1=(-0.1+(i+1)/5.0)*PI;
        for (j=0;j<4;j++) {
            x2 = -0.1+(j+1)/5.0;
            f=sin(x1)*exp(x2);

            polin2(x1a,x2a,ya,x1,x2,y,dy);
            
            double y_bilin = bilin_int.interp(x1, x2);
            double y_poly = poly_int.interp(x1, x2);
            double y_spline = spline_int.interp(x1, x2);
            
            cout << setw(12) << x1 << setw(13) << x2 << setw(13) << f;
            cout << setw(13) << y  << setw(13) << dy;
            cout << setw(13) << y_bilin  << setw(13) << y_bilin - f;
            cout << setw(13) << y_poly   << setw(13) << y_poly - f;
            cout << setw(13) << y_spline << setw(13) << y_spline - f << endl;
        }
        cout << endl << "***********************************" << endl;
    }
}

void Demo_Interpolators()
{
    std::cout << endl;
    std::cout << "***********************************************************************" << endl;
    std::cout << "****                       INTERPOLATION                           ****" << endl;
    std::cout << "***********************************************************************" << endl;

    Test_Interpolation_2DFunc();   
}