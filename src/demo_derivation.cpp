#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "algorithms/Derivation.h"
#endif
class DemoDerivationTestMemberFunc
{
    private:
        double _param;
    public:
        DemoDerivationTestMemberFunc(double param) : _param(param) { }
    
        double operator()(double x ) { return _param * sin(x); }
        double derivation(double x ) { return _param * cos(x); }
};

void Demo_Derivation()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         DERIVATION                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    DemoDerivationTestMemberFunc      funcObj(3.0);
    
    std::function<double(double)> func{funcObj};

    MML::RealFunctionFromStdFunc g{func};

    for( double h=1e-6; h>=1e-10; h/=10.0 )
    {
        std::cout << "h = " << h << std::endl;
        for( double x=0; x<=3.0; x+=0.25)
        {
            double num_der = MML::Derivation::NDer1(g, x, h);
            std::cout << std::setw(10) << std::setprecision(7) << x << "  " 
                    << std::setw(13) << std::setprecision(8) << num_der << "   " 
                    << std::setw(13) << std::setprecision(8) << funcObj.derivation(x) << "  "
                    << std::setw(13) << std::setprecision(8) << num_der - funcObj.derivation(x) << std::endl;

        }
    }
}