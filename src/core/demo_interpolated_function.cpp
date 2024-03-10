#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/Function.h"
#include "core/InterpolatedFunction.h"
#include "core/Integration.h"

#include "core/ConsolePrinter.h"
#include "core/Serializer.h"
#include "core/Visualizer.h"

#include "algorithms/FunctionAnalyzers.h"
#endif

using namespace MML;
using namespace std;

// TODO - naci skup  kompleksnijih funkcija, npr. sin(x) + cos(x) + exp(x) + x^2 + x^3, u intervalu -10, 10
// TODO - testirati 2D interpolaciju
// TODO - testirati parametric curve interpolaciju

Real test_func(const Real x)
{
    const double eps = 1.0;
    return x*exp(-x)/(POW2(x-1.0)+eps*eps);
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
    
    LinearInterpRealFunc    linear_interp(vec_x, vec_y);
    PolynomInterpRealFunc   poly_interp(vec_x, vec_y, 3);
    SplineInterpRealFunc    spline_interp(vec_x, vec_y);
    RationalInterpRealFunc  rat_interp(vec_x, vec_y, 4);
    BaryRatInterpRealFunc   barry_interp(vec_x, vec_y, 4);

    std::cout << "\nINTERPOLATION PRECISION:\n";

    // TODO - interval -10, 10 : 10 pnts, 50 pnts, 100 pnts

    TablePrinter<double, double> data("x", 8, 3, 
                                                { "exact", 
                                                  "linear", "abs.err", "rel.err", 
                                                  "poly", "abs.err", "rel.err", 
                                                  "spline", "abs.err", "rel.err", 
                                                  "rat", "abs.err", "rel.err", 
                                                  "barry", "abs.err", "rel.err"
                                                },
                                                { {10,5,'F'},
                                                  {11,5,'F'}, {9,5,'F'}, {9,3,'F'},
                                                  {11,5,'F'}, {9,5,'F'}, {9,3,'F'},
                                                  {11,5,'F'}, {9,5,'F'}, {9,3,'F'},
                                                  {11,5,'F'}, {9,5,'F'}, {9,3,'F'},
                                                  {11,5,'F'}, {9,5,'F'}, {9,3,'F'}
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

        double linear_err = y_linear - y_exact;
        double poly_err   = y_poly - y_exact;
        double spline_err = y_spline - y_exact;
        double rat_err    = y_rat - y_exact;
        double barry_err  = y_barry - y_exact;

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
    }
    data.Print();

    RealFunction test(test_func);
    
    std::cout << "Linear diff = " << RealFunctionComparer::getIntegratedDiff(linear_interp, test, 0.0, 2.0) << std::endl;
    std::cout << "Poly diff   = " << RealFunctionComparer::getIntegratedDiff(poly_interp, test, 0.0, 2.0) << std::endl;
    std::cout << "Spline diff = " << RealFunctionComparer::getIntegratedDiff(spline_interp, test, 0.0, 2.0) << std::endl;
    std::cout << "Rat diff    = " << RealFunctionComparer::getIntegratedDiff(rat_interp, test, 0.0, 2.0) << std::endl;
    std::cout << "Barry diff  = " << RealFunctionComparer::getIntegratedDiff(barry_interp, test, 0.0, 2.0) << std::endl;
}

void Test_RealFunction_Linear_interp()
{
    RealFunction f{[](Real x) { return sin(x); } };

    Vector<Real> vec_x, vec_y;
    const double x1 = 0.0;
    const double x2 = 2.0 * Constants::PI;
    const int    numPnt = 4;

    CreateInterpolatedValues(f, x1, x2, numPnt, vec_x, vec_y);
    std::cout << "Interpolated values:\n";
    for( int i=0; i<numPnt; i++)
    {
        std::cout << vec_x[i] << "   " << vec_y[i] << std::endl;
    }
    std::cout <<  std::endl;

    
    // PolynomInterpRealFunc    linear_interp(vec_x, vec_y, 2);
    LinearInterpRealFunc    linear_interp(vec_x, vec_y);
    // PolynomInterpRealFunc   poly_interp(vec_x, vec_y, 3);
    // SplineInterpRealFunc    spline_interp(vec_x, vec_y);
    // RationalInterpRealFunc  rat_interp(vec_x, vec_y, 4);
    // BaryRatInterpRealFunc   barry_interp(vec_x, vec_y, 4);

    std::cout << "\nINTERPOLATION PRECISION:\n";

    TablePrinter<double, double> data("x", 8, 3, 
                                                { "exact", 
                                                  "linear", "abs.err", "rel.err", "rel.perc.err"
                                                //   "poly", "abs.err", "rel.err", 
                                                //   "spline", "abs.err", "rel.err", 
                                                //   "rat", "abs.err", "rel.err", 
                                                //   "barry", "abs.err", "rel.err"
                                                },
                                                { {10,5,'F'},
                                                  {11,5,'F'}, {9,5,'F'}, {9,3,'F'}, {15,3,'F'}
                                                //   {11,5}, {9,5}, {9,3},
                                                //   {11,5}, {9,5}, {9,3},
                                                //   {11,5}, {9,5}, {9,3},
                                                //   {11,5}, {9,5}, {9,3}
                                                } 
                                            );

    const int printPnt = 50;
    for (int i=0;i<printPnt;i++) 
    {
        double x = i * (x2 - x1) / (printPnt - 1);
        

        double y_exact  = f(x);
        
        double y_linear = linear_interp(x);
        // double y_poly   = poly_interp(x);
        // double y_spline = spline_interp(x);
        // double y_rat    = rat_interp(x);
        // double y_barry  = barry_interp(x);

        double linear_abs_err = y_exact - y_linear;
        // double poly_err   = y_poly - y_exact;
        // double spline_err = y_spline - y_exact;
        // double rat_err    = y_rat - y_exact;
        // double barry_err  = y_barry - y_exact;

        double linear_rel_err  = y_exact == 0 ? 0.0 : linear_abs_err / y_exact;
        double linear_perc_err  = y_exact == 0 ? 0.0 : linear_abs_err / y_exact * 100.0;
        // double poly_perc_err    = y_exact == 0 ? 0.0 : poly_err / y_exact * 100.0;
        // double spline_perc_err  = y_exact == 0 ? 0.0 : spline_err / y_exact * 100.0;
        // double rat_perc_err     = y_exact == 0 ? 0.0 : rat_err / y_exact * 100.0;
        // double barry_perc_err   = y_exact == 0 ? 0.0 : barry_err / y_exact * 100.0;

        data.addRow(x, { y_exact, 
                        y_linear, linear_abs_err, linear_rel_err, linear_perc_err 
                        // y_poly, poly_err, poly_perc_err, 
                        // y_spline, spline_err, spline_perc_err, 
                        // y_rat, rat_err, rat_perc_err, 
                        // y_barry, barry_err, barry_perc_err 
                        } );
    }
    data.Print();

    RealFunctionComparer      comparer(f, linear_interp);

    double totalAbsDiff = comparer.getAbsDiffSum(0.0, 2.0*Constants::PI, 100);
    double avgAbsDiff = comparer.getAbsDiffAvg(0.0, 2.0*Constants::PI, 100);
    double maxAbsDiff = comparer.getAbsDiffMax(0.0, 2.0*Constants::PI, 100);

    double totalRelDiff = comparer.getRelDiffSum(0.0, 2.0*Constants::PI, 100);
    double avgRelDiff = comparer.getRelDiffAvg(0.0, 2.0*Constants::PI, 100);
    double maxRelDiff = comparer.getRelDiffMax(0.0, 2.0*Constants::PI, 100);  
 
    std::cout << "Total abs diff = " << totalAbsDiff << std::endl;
    std::cout << "Avg abs diff   = " << avgAbsDiff << std::endl;
    std::cout << "Max abs diff   = " << maxAbsDiff << std::endl;
    std::cout << "Total rel diff = " << totalRelDiff << std::endl;
    std::cout << "Avg rel diff   = " << avgRelDiff << std::endl;
    std::cout << "Max rel diff   = " << maxRelDiff << std::endl;
    
    Visualizer::VisualizeRealFunction(linear_interp, "sin(x) linear interpolation", x1, x2, 100, "func_sin_x_lin_interp.txt");
}

void Demo_Interpolated_Function()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     INTERPOLATED FUNCTION                     ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    //Test_RealFunction_interp();
    Test_RealFunction_Linear_interp();
}