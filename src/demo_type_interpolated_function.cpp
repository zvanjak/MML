#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "utilities/DataContainers.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "basic_types/Function.h"
#include "basic_types/InterpolatedFunction.h"

#include "algorithms/Integration.h"
#include "algorithms/FunctionAnalyzer.h"
#endif

using namespace MML;
using namespace std;

double TestFunc2(double x)
{
    return sin(x);
}

// TODO - naci skup  kompleksnijih funkcija, npr. sin(x) + cos(x) + exp(x) + x^2 + x^3, u intervalu -10, 10
// TODO - testirati 2D interpolaciju
// TODO - testirati parametric curve interpolaciju

double test_func(const double x)
{
    const double eps = 1.0;
    return x*exp(-x)/(SQR(x-1.0)+eps*eps);
}

void CreateInterpolatedValues(RealFunction f, Real x1, Real x2, int numPnt, Vector<Real> &outX, Vector<Real> &outY)
{
    outX.Resize(numPnt);
    outY.Resize(numPnt);

    for (int i=0;i<numPnt;i++) {
        outX[i] = x1 + i * (x2 - x1) / (numPnt - 1);
        outY[i] = f(outX[i]);
    } 
}

void Test_RealFunction_interp()
{
    Vector<Real> vec_x, vec_y;
    const double x1 = 0.0;
    const double x2 = 2.0;
    const int    numPnt = 5;

    CreateInterpolatedValues(test_func, x1, x2, numPnt, vec_x, vec_y);
    
    LinearInterpRealFunc    linear_interp(vec_x,vec_y);
    PolynomInterpRealFunc   poly_interp(vec_x,vec_y,3);
    SplineInterpRealFunc    spline_interp(vec_x,vec_y);
    RationalInterpRealFunc  rat_interp(vec_x,vec_y, 4);
    BaryRatInterpRealFunc   barry_interp(vec_x,vec_y, 4);

    std::cout << "\nINTERPOLATION PRECISION:\n";

    // TODO - interval -10, 10 : 10 pnts, 50 pnts, 100 pnts

    DataSeriesMultiRow<double, double> data("x", 8, 3, 
                                                { "exact", 
                                                  "linear", "abs.err", "rel.err", 
                                                  "poly", "abs.err", "rel.err", 
                                                  "spline", "abs.err", "rel.err", 
                                                  "rat", "abs.err", "rel.err", 
                                                  "barry", "abs.err", "rel.err"
                                                },
                                                { {10,5}, 
                                                  {11,5}, {9,5}, {9,3},
                                                  {11,5}, {9,5}, {9,3},
                                                  {11,5}, {9,5}, {9,3},
                                                  {11,5}, {9,5}, {9,3},
                                                  {11,5}, {9,5}, {9,3}
                                                } 
                                            );


    const int printPnt = 20;
    for (int i=0;i<printPnt;i++) 
    {
        double x = i * (x2 - x1) / (printPnt - 1);
        
        double y_exact  = test_func(x);
        
        double y_linear = linear_interp(x);
        double y_poly   = poly_interp(x);
        double y_spline = spline_interp(x);
        double y_rat    = rat_interp(x);
        double y_barry  = barry_interp(x);

        double linear_err       = y_linear - y_exact;
        double poly_err         = y_poly - y_exact;
        double spline_err       = y_spline - y_exact;
        double rat_err          = y_rat - y_exact;
        double barry_err        = y_barry - y_exact;

        double linear_perc_err  = y_exact == 0 ? 0.0 : linear_err / y_exact * 100.0;
        double poly_perc_err    = y_exact == 0 ? 0.0 : poly_err / y_exact * 100.0;
        double spline_perc_err  = y_exact == 0 ? 0.0 : spline_err / y_exact * 100.0;
        double rat_perc_err     = y_exact == 0 ? 0.0 : rat_err / y_exact * 100.0;
        double barry_perc_err   = y_exact == 0 ? 0.0 : barry_err / y_exact * 100.0;

        data.addRow(x, { y_exact, 
                        y_linear, linear_err, linear_perc_err, 
                        y_poly, poly_err, poly_perc_err, 
                        y_spline, spline_err, spline_perc_err, 
                        y_rat, rat_err, rat_perc_err, 
                        y_barry, barry_err, barry_perc_err } );

        // cout << fixed << setw(8) << setprecision(3) << x << setw(10) << setprecision(5) << y_exact;
        // cout << " " << setw(11) << setprecision(5) << y_linear   << "  " << setw(9) << setprecision(5) << right << linear_err << "   " << setw(6) << setprecision(2) << linear_perc_err;
        // cout << " " << setw(11) << setprecision(5) << y_poly     << "  " << setw(9) << setprecision(5) << right << poly_err << "   " << setw(6) << setprecision(2) << poly_perc_err;
        // cout << " " << setw(11) << setprecision(5) << y_spline   << "  " << setw(9) << setprecision(5) << right << spline_err << "   " << setw(6) << setprecision(2) << spline_perc_err;
        // cout << " " << setw(11) << setprecision(5) << y_rat      << "  " << setw(9) << setprecision(5) << right << rat_err << "   " << setw(6) << setprecision(2) << rat_perc_err;
        // cout << " " << setw(11) << setprecision(5) << y_barry    << "  " << setw(9) << setprecision(5) << right << barry_err << "   " << setw(6) << setprecision(2) << barry_perc_err << endl;
    }
    data.Print();

    RealFunction test(test_func);
    
    std::cout << "Linear diff = " << FunctionAnalyzer::FuncDiff(linear_interp, test, 0.0, 2.0) << std::endl;
    std::cout << "Poly diff   = " << FunctionAnalyzer::FuncDiff(poly_interp, test, 0.0, 2.0) << std::endl;
    std::cout << "Spline diff = " << FunctionAnalyzer::FuncDiff(spline_interp, test, 0.0, 2.0) << std::endl;
    std::cout << "Rat diff    = " << FunctionAnalyzer::FuncDiff(rat_interp, test, 0.0, 2.0) << std::endl;
    std::cout << "Barry diff  = " << FunctionAnalyzer::FuncDiff(barry_interp, test, 0.0, 2.0) << std::endl;
}

void Demo_Interpolated_Function()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     INTERPOLATED FUNCTION                     ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Test_RealFunction_interp();
}