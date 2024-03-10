#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Function.h"
#include "core/Derivation.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"

using namespace MML;
using namespace MML::TestBeds;

// For given RealFunction
// and given interval of test points
// produces table of first derivation error order for different derivation orders
void NDer_Error_Order_Diff_Der_Orders_Single_Func(std::string funcName, const TestFunctionReal  &inFunc, std::vector<Real> intervalPoints)
{
    std::vector<int> nder_orders{1, 2, 4, 6, 8};

    const RealFunction  &f     = inFunc._func;
    const RealFunction  &f_der = inFunc._funcDerived;

    double err, err_sum, pr_err;

    std::cout << "\nORDER OF DERIVATION ERROR FOR DIFFERENT DERIVATION ORDERS - function - " << funcName << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Ord \\ X  |";
    for(auto x : intervalPoints) {
        std::cout << std::setw(7) << x << " ";
    }
    std::cout << "  |  Avg. abs. error" << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

    for(auto ord : nder_orders)
    {
        err_sum = 0.0;
        std::cout << "NDer" << ord << "    |  ";
        
        for(auto x : intervalPoints) 
        {
            double exact_der = f_der(x);
            double num_der = 0.0;

            switch(ord)
            {
                case 1: num_der = MML::Derivation::NDer1(f, x); break;
                case 2: num_der = MML::Derivation::NDer2(f, x); break;
                case 4: num_der = MML::Derivation::NDer4(f, x); break;
                case 6: num_der = MML::Derivation::NDer6(f, x); break;
                case 8: num_der = MML::Derivation::NDer8(f, x); break;
                default: std::cout << "Wrong order" << std::endl;
            }
            
            err      = num_der - exact_der;
            err_sum += std::abs(err);
            pr_err   = pow(10,ceil(log10(std::abs(err))));

            std::cout << std::scientific << std::setw(8) << std::setprecision(0) << std::abs(pr_err) << "   ";
        }
        std::cout << "|   " << std::setw(11) << std::setprecision(6) << err_sum / intervalPoints.size() << std::fixed << std::endl;
    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
}


// For given RealFunction
// and given interval of test points
// produces table of second derivation error order for different derivation orders
void NSecDer_Error_Order_Diff_Der_Orders_Single_Func(std::string funcName, const TestFunctionReal  &inFunc, std::vector<Real> intervalPoints)
{
    std::vector<int> nder_orders{1, 2, 4, 6, 8};

    const RealFunction  &f     = inFunc._func;
    const RealFunction  &f_sec_der = inFunc._funcSecDer;

    double err, err_sum, pr_err;

    std::cout << "\nORDER OF SECOND DERIVATION ERROR FOR DIFFERENT DERIVATION ORDERS - function - " << funcName << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Ord \\ X  |";
    for(auto x : intervalPoints) {
        std::cout << std::setw(10) << x << " ";
    }
    std::cout << "  |  Avg. abs. error" << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

    for(auto ord : nder_orders)
    {
        err_sum = 0.0;
        std::cout << "NSecDer" << ord << " |  ";
        
        for(auto x : intervalPoints) 
        {
            double exact_der = f_sec_der(x);
            double num_der = 0.0;

            switch(ord)
            {
                case 1: num_der = MML::Derivation::NSecDer1(f, x); break;
                case 2: num_der = MML::Derivation::NSecDer2(f, x); break;
                case 4: num_der = MML::Derivation::NSecDer4(f, x); break;
                case 6: num_der = MML::Derivation::NSecDer6(f, x); break;
                case 8: num_der = MML::Derivation::NSecDer8(f, x); break;
                default: std::cout << "Wrong order" << std::endl;
            }
            
            err      = num_der - exact_der;
            err_sum += std::abs(err);
            pr_err   = pow(10,ceil(log10(std::abs(err))));

            std::cout << std::scientific << std::setw(8) << std::setprecision(0) << std::abs(pr_err) << "   ";
        }
        std::cout << "|   " << std::setw(11) << std::setprecision(6) << err_sum / intervalPoints.size() << std::fixed << std::endl;
    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
}
void NDer_Average_Error_Diff_Orders_Single_Func(std::string funcName, const TestFunctionReal  &inFunc)
{
    const RealFunction  &f     = inFunc._func;
    const RealFunction  &f_der = inFunc._funcDerived;
    double x1 = inFunc._intervalTest->getLowerBound();
    double x2 = inFunc._intervalTest->getUpperBound();

    const int numPntForEval = 20;

    double err_sum1 = 0.0;
    double err_sum2 = 0.0;
    double err_sum4 = 0.0;
    double err_sum6 = 0.0;
    double err_sum8 = 0.0;

    std::cout << "\nAVERAGE DERIVATION ERROR FOR DIFFERENT ORDERS" << std::endl;
    std::cout << "    X    Exact der.       Nder1        NDer1 err.          Nder2         NDer2 err.          Nder4         NDer4 err.          Nder6         NDer6 err.          Nder8         NDer8 err.         " << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

    for(int i=0; i<numPntForEval; i++ ){
        double x = x1 + (x2-x1)*i/(numPntForEval-1);

        double exact_der = f_der(x);

        double num_der1 = MML::Derivation::NDer1(f, x);
        double num_der2 = MML::Derivation::NDer2(f, x);
        double num_der4 = MML::Derivation::NDer4(f, x);
        double num_der6 = MML::Derivation::NDer6(f, x);
        double num_der8 = MML::Derivation::NDer8(f, x);

        double err1 = num_der1 - exact_der;
        double err2 = num_der2 - exact_der;
        double err4 = num_der4 - exact_der;
        double err6 = num_der6 - exact_der;
        double err8 = num_der8 - exact_der;

        std::cout << std::setw(6)  << std::setprecision(3) << x 
                << std::setw(13)   << std::setprecision(8) << exact_der << " "
                << std::setw(13)   << std::setprecision(8) << num_der1 << "   " 
                << std::scientific << std::setw(15)        << err1 << "   " << std::fixed
                << std::setw(13)   << std::setprecision(8) << num_der2 << "   " 
                << std::scientific << std::setw(15)        << err2 << "   " << std::fixed
                << std::setw(13)   << std::setprecision(8) << num_der4 << "   " 
                << std::scientific << std::setw(15)        << err4 << "   " << std::fixed
                << std::setw(13)   << std::setprecision(8) << num_der6 << "   " 
                << std::scientific << std::setw(15)        << err6 << "   " << std::fixed
                << std::setw(13)   << std::setprecision(8) << num_der8 << "   " 
                << std::scientific << std::setw(15)        << err8 << "   " << std::fixed
                << std::endl;

        err_sum1 += std::abs(err1);
        err_sum2 += std::abs(err2);
        err_sum4 += std::abs(err4);
        err_sum6 += std::abs(err6);
        err_sum8 += std::abs(err8);
    }    

    std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Total abs error =                    " << std::scientific << err_sum1 << "                    " 
                                                        << err_sum2 << "                    " 
                                                        << err_sum4 << "                    " 
                                                        << err_sum6 << "                    " 
                                                        << err_sum8 << std::endl;
}

void NDer_Error_Diff_h_Single_Func(int order, std::string funcName, const TestFunctionReal  &inFunc)
{
    const RealFunction  &f     = inFunc._func;
    const RealFunction  &f_der = inFunc._funcDerived;
    double x1 = inFunc._intervalTest->getLowerBound();
    double x2 = inFunc._intervalTest->getUpperBound();

    const int numPntForEval = 20;

    Vector<Real> h_list{0.1, 0.01, 0.001, 0.0001, 0.00001};
    Vector<Real> h_list2{0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001};
    Vector<Real> err_list(h_list.size(), 0.0);

    std::cout << "\nNDer" << order << " DERIVATION ERROR FOR DIFFERENT STEP SIZES" << std::endl;
    std::cout << "   X     Exact der.";
    for(int j=0; j<h_list.size(); j++ )
        std::cout << "   h=" << std::setw(10) << std::setprecision(8) << h_list[j] << "        err.       ";
    std::cout << "\n---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

    for(int i=0; i<numPntForEval; i++ ){
        double x = x1 + (x2-x1)*i/(numPntForEval-1);

        double exact_der = f_der(x);

        std::cout << std::setw(6)  << std::setprecision(3) << x << std::setw(13) << std::setprecision(8) << exact_der << " ";

        for(int j=0; j<h_list.size(); j++ )
        {
            double h       = h_list[j];
            double num_der = 0.0;

            switch(order)
            {
                case 1: num_der = MML::Derivation::NDer1(f, x, h); break;
                case 2: num_der = MML::Derivation::NDer2(f, x, h); break;
                case 4: num_der = MML::Derivation::NDer4(f, x, h); break;
                case 6: num_der = MML::Derivation::NDer6(f, x, h); break;
                case 8: num_der = MML::Derivation::NDer8(f, x, h); break;
                default: std::cout << "Wrong order" << std::endl;
            }
            double err = num_der - exact_der;

            std::cout << std::setw(13)   << std::setprecision(8) << num_der << "   " 
                      << std::scientific << std::setw(15) << err << "   " << std::fixed;

            err_list[j] += std::abs(err);

        }
        std::cout << std::endl;
    }    

    std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Total abs error =                   ";
    for(int j=0; j<h_list.size(); j++ )
        std::cout << std::scientific << std::setw(15) << err_list[j] << "                   " << std::fixed;
    std::cout << std::endl;
}

// TODO - finish NDer_Error_Diff_Orders_Multi_Func_comparison
void NDer_Error_Diff_Orders_Multi_Func_comparison()
{
    // 10 funkcija po redovima
    const int numPntEval = 50;
    // ime, interval, točna vrijednost, za svaki red absolutna sumarna greška i average greška 
}

void Demo_Second_derivation()
{
    MML::RealFunction g{[](Real x) { return sin(x); } };
    MML::RealFunction g_sec_der{[](Real x) { return -sin(x); } };

    for( double h=1e-4; h>=1e-10; h/=10.0 )
    {
        double err_sum = 0.0;
        std::cout << "h = " << h << std::endl;
        std::cout << "         x     num.sec.der.   exact.sec.der.     error     " << std::endl;
        std::cout << "-----------------------------------------------------------" << std::endl;
        for( double x=0; x<=3.0; x+=0.5)
        {
            double num_der = MML::Derivation::NSecDer1(g, x, h);
            double err = num_der - g_sec_der(x);
            std::cout << std::setw(10) << std::setprecision(7) << x << "  " 
                      << std::setw(15) << std::setprecision(10) << num_der << "   " 
                      << std::setw(15) << std::setprecision(10) << g_sec_der(x) << "  "
                      << std::setw(15) << std::setprecision(10) << err << std::endl;
            
            err_sum += std::abs(err);
        }
        std::cout << "Total error = " << err_sum << std::endl;
    }
}

void Demo_Third_derivation()
{
    MML::RealFunction g{[](Real x) { return sin(x); } };
    MML::RealFunction g_third_der{[](Real x) { return -cos(x); } };

    for( double h=1e-4; h>=1e-6; h/=10.0 )
    {
        double err_sum = 0.0;
        std::cout << "h = " << h << std::endl;
        std::cout << "         x     num.3rd.der.   exact.3rd.der.     error     " << std::endl;
        std::cout << "-----------------------------------------------------------" << std::endl;
        for( double x=0; x<=3.0; x+=0.5)
        {
            double num_der = MML::Derivation::NThirdDer1(g, x, h);
            double err = num_der - g_third_der(x);
            std::cout << std::setw(10) << std::setprecision(7) << x << "  " 
                      << std::setw(15) << std::setprecision(10) << num_der << "   " 
                      << std::setw(15) << std::setprecision(10) << g_third_der(x) << "  "
                      << std::setw(15) << std::setprecision(10) << err << std::endl;
            
            err_sum += std::abs(err);
        }
        std::cout << "Total error = " << err_sum << std::endl;
    }
}

void Test_Precision_Derivation()
{
    //std::vector<Real> intervalPoints{-10.0, -5.0, -3.0, -1.0, -0.5, -0.1, -0.025, 0.0, 0.025, 0.1, 0.5, 1.0, 3.0, 5.0, 10.0};
    std::vector<Real> intervalPoints{0.0, 0.5, 1.0, 2.0, 3.0, 5.0 };

    NDer_Error_Order_Diff_Der_Orders_Single_Func("Sin", RealFunctionsTestBed::getTestFunctionReal("Sin"), intervalPoints);
    // NSecDer_Error_Order_Diff_Der_Orders_Single_Func("Sin", RealFunctionsTestBed::getTestFunctionReal("Sin"), intervalPoints);

    NDer_Average_Error_Diff_Orders_Single_Func("Sin", RealFunctionsTestBed::getTestFunctionReal(0));

    // NDer_Error_Diff_h_Single_Func(1, "Sin", RealFunctionsTestBed::getTestFunctionReal("Sin"));
    // NDer_Error_Diff_h_Single_Func(2, "Sin", RealFunctionsTestBed::getTestFunctionReal("Sin"));
    // NDer_Error_Diff_h_Single_Func(4, "Sin", RealFunctionsTestBed::getTestFunctionReal("Sin"));
    // NDer_Error_Diff_h_Single_Func(6, "Sin", RealFunctionsTestBed::getTestFunctionReal("Sin"));
    // NDer_Error_Diff_h_Single_Func(8, "Sin", RealFunctionsTestBed::getTestFunctionReal("Sin"));

    // Demo_Second_derivation();
    // Demo_Third_derivation();
}