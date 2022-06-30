#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "MMLBase.h"
#include "algorithms/Interpolators.h"
#include "basic_types/Vector.h"
#include "basic_types/Matrix.h"
#endif

using namespace std;
using namespace MML;


void Test_Interpolation_RealFunc_Polynomial()
{
    const double PI=3.141592653589793238;
    int i,n,nfunc;
    double dy,f,x,y;

    cout << "generation of interpolation tables" << endl;
    cout << " ... sin(x)    0<x<PI" << endl;
    cout << " ... exp(x)    0<x<1 " << endl;
    // cout << "how many entries go in these tables?" << endl;
    // cin >> n;
    // cin.get();
    n = 10;
    if (n < 1) 
        return;
        
    Vector xa(n);
    Vector ya(n);
    cout << fixed << setprecision(6);
    for (nfunc=0;nfunc<2;nfunc++) {
        if (nfunc == 0) {
        cout << endl << "sine function from 0 to PI" << endl;
        for (i=0;i<n;i++) {
            xa[i]=(i+1)*PI/n;
            ya[i]=sin(xa[i]);
        }
        } else if (nfunc == 1) {
        cout << endl << "exponential function from 0 to 1" << endl;
        for (i=0;i<n;i++) {
            xa[i]=(i+1)*1.0/n;
            ya[i]=exp(xa[i]);
        }
        } else {
        break;
        }
        cout << setw(9) << "x" << setw(14) << "f(x)";
        cout << setw(17) << "interpolated" << setw(14) << "error" << endl;
        for (i=0;i<10;i++) {
        if (nfunc == 0) {
            x=(-0.05+(i+1)/10.0)*PI;
            f=sin(x);
        } else if (nfunc == 1) {
            x=(-0.05+(i+1)/10.0);
            f=exp(x);
        }
        MML::polint(xa,ya,x,y,dy);
        cout << setw(12) << x << setw(13) << f << setw(13) << y;
        cout << "     " << setw(12) << dy << endl;
        }
        cout << endl << "***********************************" << endl;
        cout << "press RETURN" << endl;
        //cin.get();
    }    
}

double test_func(const double x, const double eps)
{
        return x*exp(-x)/(SQR(x-1.0)+eps*eps);
}

void Test_Interpolation_RealFunc_Rational()
{
    const int NPT=6;
    const double EPS=1.0;
    int i;
    double dyy,xx,yexp,yy;
    Vector x(NPT),y(NPT);

    for (i=0;i<NPT;i++) {
        x[i]=(i+1)*2.0/NPT;
        y[i]=test_func(x[i],EPS);
    }
    cout << endl << "Diagonal rational function interpolation" << endl;
    cout << endl << setw(5) << "x" << setw(15) << "interp.";
    cout << setw(15) << "accuracy" << setw(13) << "actual" << endl;
    cout << fixed << setprecision(6);
    for (i=0;i<10;i++) {
        xx=0.2*(i+1);
        
        MML::ratint(x,y,xx,yy,dyy);
        
        yexp=test_func(xx,EPS);
        cout << setw(8) << xx << setw(13) << yy;
        cout << setw(14) << dyy << setw(14) << yexp << endl;
    } 
}

void Test_Interpolation_2DFunc()
{
    const int N=5;
    const double PI=3.141592653589793238;
    int i,j;
    double dy,f,x1,x2,y;
    Vector x1a(N),x2a(N);
    Matrix ya(N,N);

    for (i=0;i<N;i++) {
        x1a[i]=(i+1)*PI/N;
        for (j=0;j<N;j++) {
            x2a[j]=1.0*(j+1)/N;
            ya[i][j]=sin(x1a[i])*exp(x2a[j]);
        }
    }
    // test 2-dimensional interpolation
    cout << endl << "Two dimensional interpolation of sin(x1)exp(x2)";
    cout << endl;
    cout << setw(9) << "x1" << setw(13) << "x2" << setw(14) << "f(x)";
    cout << setw(17) << "interpolated" << setw(12) << "error" << endl;
    cout << fixed << setprecision(6);
    for (i=0;i<4;i++) {
        x1=(-0.1+(i+1)/5.0)*PI;
        for (j=0;j<4;j++) {
            x2 = -0.1+(j+1)/5.0;
            f=sin(x1)*exp(x2);

            MML::polin2(x1a,x2a,ya,x1,x2,y,dy);
            
            cout << setw(12) << x1 << setw(13) << x2 << setw(13) << f;
            cout << setw(13) << y << setw(13) << dy << endl;
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

    Test_Interpolation_RealFunc_Polynomial();   
    Test_Interpolation_2DFunc();   
    Test_Interpolation_RealFunc_Rational();
}