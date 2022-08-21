#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "algorithms/Integration.h"
#endif

using namespace std;

int Test_qtrap(void)
{
    // Test function
    MML::RealFunction func_qtrap( [](double x) { return x*x*(x*x-2.0)*sin(x); } );
    // Integral of test function
    MML::RealFunction fint_qtap( [](double x) { return 4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x); } );

    const double PIO2=1.570796326794896619;
    double a=0.0,b=PIO2,s;

    cout << "Integral of func computed with QTRAP" << endl << endl;
    cout << fixed << setprecision(6);
    cout << "Actual value of integral is ";
    cout << setw(12) << (fint_qtap(b)-fint_qtap(a)) << endl;

    s=MML::Integration::IntegrateTrap(func_qtrap,a,b);

    cout << "Result from routine QTRAP is " << setw(12) << s << endl;
    return 0;
}

// Driver for routine qsimp
int Test_qsimp(void)
{
    // Test function
    MML::RealFunction func_qsimp( [](double x) { return x*x*(x*x-2.0)*sin(x); } );
    // Integral of test function
    MML::RealFunction fint_qsimp( [](double x) { return 4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x); } );

    const double PIO2=1.570796326794896619;
    double a=0.0,b=PIO2,s;

    cout << "Integral of func computed with QSIMP" << endl << endl;
    cout << "Actual value of integral is "
        << setw(12) << (fint_qsimp(b)-fint_qsimp(a)) << endl;
    
    s=MML::Integration::IntegrateSimpson(func_qsimp,a,b);
    
    cout << "Result from routine QSIMP is " << setw(12) << s << endl;
    return 0;
}

// Driver for routine qromb
int Test_qromb(void)
{
    // Test function
    MML::RealFunction func_qromb( [] (double x) { return x*x*(x*x-2.0)*sin(x); } );
    // Integral of test function func
    MML::RealFunction fint_qromb( [] (double x) { return 4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x); } );
    
    const double PIO2=1.570796326794896619;
    double a=0.0,b=PIO2,s;

    cout << "Integral of func computed with QROMB" << endl << endl;
    cout << fixed << setprecision(6);
    cout << "Actual value of integral is ";
    cout << setw(12) << (fint_qromb(b)-fint_qromb(a)) << endl;

    s=MML::Integration::IntegrateRomberg(func_qromb,a,b);

    cout << "Result from routine QROMB is " << setw(12) << s << endl;
    return 0;
}

void Demo_Integration()
{
    std::cout << endl;
    std::cout << "***********************************************************************" << endl;
    std::cout << "****                         INTEGRATION                           ****" << endl;
    std::cout << "***********************************************************************" << endl;

    Test_qtrap();
    Test_qsimp();
    Test_qromb();
}