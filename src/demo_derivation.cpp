#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "algorithms/Derivation.h"
#include "../tests/test_data/functions_test_bed.h"

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

void Demo_Derivation_func_ptr()
{

}

void Demo_Derivation_func_interpolated()
{

}

void Demo_Derivation_member_fun()
{
    DemoDerivationTestMemberFunc      funcObj(3.0);
    
    std::function<double(double)> func{funcObj};

    MML::RealFunctionFromStdFunc g{func};

    for( double h=1e-5; h>=1e-10; h/=10.0 )
    {
        double err_sum = 0.0;
        std::cout << "h = " << h << std::endl;
        for( double x=0; x<=3.0; x+=0.5)
        {
            double num_der = MML::Derivation::NDer4(g, x, h);
            double err = num_der - funcObj.derivation(x);
            std::cout << std::setw(10) << std::setprecision(7) << x << "  " 
                    << std::setw(13) << std::setprecision(8) << num_der << "   " 
                    << std::setw(13) << std::setprecision(8) << funcObj.derivation(x) << "  "
                    << std::setw(13) << std::setprecision(8) << err << std::endl;
            
            err_sum += std::abs(err);
        }
        std::cout << "Total error = " << err_sum << std::endl;
    }
}

void Demo_Partial_derivation()
{

}

// Comparing precision of derivation routines of different order, for simple function
void Demo_Derivation_precision_comparison_single()
{
    MML::RealFunction g     = MML::Tests::FunctionsTestBed::_listFuncReal[0]._func;
    MML::RealFunction g_der = MML::Tests::FunctionsTestBed::_listFuncReal[0]._funcDerived;

    for(int i=0; i<5; i++ )
    {
        double num_der;
        double err_sum = 0.0;
        for( double x=0; x<=3.0; x+=0.5)
        {
            switch(i) {
                case 0:  num_der = MML::Derivation::NDer1(g, x); break;
                case 1:  num_der = MML::Derivation::NDer2(g, x); break;
                case 2:  num_der = MML::Derivation::NDer4(g, x); break;
                case 3:  num_der = MML::Derivation::NDer6(g, x); break;
                case 4:  num_der = MML::Derivation::NDer8(g, x); break;
            }
            double err = num_der - g_der(x);
            std::cout << std::setw(10) << std::setprecision(7) << x << "  " 
                    << std::setw(13) << std::setprecision(8) << num_der << "   " 
                    << std::setw(13) << std::setprecision(8) << g_der(x) << "  "
                    << std::setw(13) << std::setprecision(8) << err << std::endl;
            
            err_sum += std::abs(err);
        }
        std::cout << "Total error = " << err_sum << std::endl;
    }
}

void Demo_Derivation_precision_comparison()
{
    double start_x = -10.0;
    double end_x = 10.0;
    double step = 0.5;
    int count = (int) ((end_x - start_x) / step);

    std::cout << "AVERAGE DERIVATION ERROR FOR DIFFERENT ORDERS" << std::endl;
    std::cout << "Function                  Nder1           NDer2           NDer4           NDer6           NDer8" << std::endl;

    for(auto&& test_func : MML::Tests::FunctionsTestBed::_listFuncReal)
    {
        MML::RealFunction g     = test_func._func;
        MML::RealFunction g_der = test_func._funcDerived;


        
        std::cout << std::setw(15) << test_func._funcExpr << " - ";

        for(int i=0; i<5; i++ )
        {
            double num_der;
            double err_sum = 0.0;
            for( double x=start_x; x<=end_x; x+=step)
            {
                switch(i) {
                    case 0:  num_der = MML::Derivation::NDer1(g, x); break;
                    case 1:  num_der = MML::Derivation::NDer2(g, x); break;
                    case 2:  num_der = MML::Derivation::NDer4(g, x); break;
                    case 3:  num_der = MML::Derivation::NDer6(g, x); break;
                    case 4:  num_der = MML::Derivation::NDer8(g, x); break;
                }
                double err = num_der - g_der(x);
                // std::cout << std::setw(10) << std::setprecision(7) << x << "  " 
                //         << std::setw(13) << std::setprecision(8) << num_der << "   " 
                //         << std::setw(13) << std::setprecision(8) << g_der(x) << "  "
                //         << std::setw(13) << std::setprecision(8) << err << std::endl;
                
                err_sum += std::abs(err);
            }
            std::cout << " " << std::setw(15) << err_sum/count ;
        }
        std::cout << std::endl;
    }
}

void Demo_Derivation()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         DERIVATION                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    Demo_Derivation_precision_comparison();
}