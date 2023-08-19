#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/InterpolatedFunction.h"

#include "algorithms/Derivation.h"
#endif

#include "../test_data/real_functions_test_bed.h"

using namespace MML;

// TODO - simple precision comparison for selected function

//////////////////////          REAL FUNCTION DERIVATION           ////////////////////////
double DemoDerRealFunc_TestFunc(double x) 
{ 
    return sin(x)*(1.0 + 0.5*x*x); 
}

void Demo_Derivation_func_ptr()
{
    // creating a function object from a already existing (standalone) function
    RealFunction f1(DemoDerRealFunc_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    // creating a function object from a std::function object
    std::function<double(double)> h{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
    RealFunctionFromStdFunc f3(h);

    double der_f12 = Derivation::NDer2(f1, 0.5);
    double der_f14 = Derivation::NDer4(f1, 0.5);
    double der_f18 = Derivation::NDer8(f1, 0.5);
    double der_f21 = Derivation::NDer1(f2, 0.5);
    double der_f36 = Derivation::NDer6(f3, 0.5);

    // we can use default Derive routine (default set to NDer4)
    double num_der4 = Derivation::Derive(f1, 0.5, nullptr);
}

// If you CAN change the class where your data for calculation is
class ClassProvidingFuncToDerive
{
    private:
        double _param;
    public:
        ClassProvidingFuncToDerive(double param) : _param(param) { }
    
        // ust provide operator() for your class
        double operator()(double x ) { return _param * sin(x); }
};

void Demo_Derivation_member_fun()
{
    ClassProvidingFuncToDerive   funcObj(3.0);
    
    MML::RealFunctionFromStdFunc g(std::function<double(double)>{funcObj});
    double der_g = MML::Derivation::NDer4(g, 0.5);
}

// If you CAN'T change the class where your data for calculation is
class BigComplexClassYouCantChange
{
    // has data for calculating function you want to derive
};

// Create a helper wrapper class that will provide operator() for your class
// and implement calculation of function value
class BigComplexDerivFunc
{
    const BigComplexClassYouCantChange &_ref;
public:
    BigComplexDerivFunc(const BigComplexClassYouCantChange &bigClass) : _ref(bigClass) { }

    double operator()(double x ) 
    {
        double val = 0.0; 
        // complex calculations using _ref
        return val; 
    }
};

void Demo_Derivation_member_fun2(const BigComplexClassYouCantChange &ref)
{
    BigComplexDerivFunc          funcObj(ref);  
    MML::RealFunctionFromStdFunc func_to_derive(std::function<double(double)>{funcObj});
    
    double der_1 = MML::Derivation::NDer4(func_to_derive, 0.5);
    double der_2 = MML::Derivation::Derive(func_to_derive, 0.5, nullptr);
}

void Demo_Derivation_Interpolated_RealFunc()
{
    // creating interpolated function
    Vector<double> x_val(100);
    Vector<double> y_val(100);
    for( int i=0; i<100; i++ )
    {
        x_val[i] = i / 100.0;
        y_val[i] = sin(x_val[i])*(1.0 + 0.5*x_val[i]*x_val[i]);
    }

    LinearInterpRealFunc f_linear(x_val, y_val);

    double der_f1 = Derivation::NDer1(f_linear, 0.5);
    double der_f2 = Derivation::NDer2(f_linear, 0.5);
    double der_f3 = Derivation::NDer4(f_linear, 0.5);
    double der_f4 = Derivation::NDer6(f_linear, 0.5);
    double der_f5 = Derivation::NDer8(f_linear, 0.5);
    double der_f6 = Derivation::Derive(f_linear, 0.5, nullptr);
}

void Demo_Derivation_RealFunc_Second_and_Third()
{
    RealFunction f1(DemoDerRealFunc_TestFunc);
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
    
    std::function<double(double)>   h{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
    RealFunctionFromStdFunc         f3(h);

    ClassProvidingFuncToDerive      funcObj(3.0);
    MML::RealFunctionFromStdFunc    f4(std::function<double(double)>{funcObj});

    double sec_der_f1 = Derivation::NSecDer1(f1, 0.5);
    double sec_der_f2 = Derivation::NSecDer2(f2, 0.5);
    double sec_der_f3 = Derivation::NSecDer4(f3, 0.5);
    double sec_der_f4 = Derivation::NSecDer6(f4, 0.5);
    double sec_der_f5 = Derivation::DeriveSec(f1, 0.5, nullptr);

    double third_der_f1 = Derivation::NThirdDer1(f1, 0.5);
    double third_der_f2 = Derivation::NThirdDer2(f2, 0.5);
    double third_der_f3 = Derivation::NThirdDer4(f3, 0.5);
    double third_der_f4 = Derivation::NThirdDer6(f4, 0.5);
    double third_der_f5 = Derivation::DeriveThird(f1, 0.5, nullptr);
}

// simple example for a single function of various precisions
void Demo_Derivation_RealFunc_precision_comparison()
{

}

/////////////////////          SCALAR FUNCTION DERIVATION           ///////////////////////
void Demo_Derivation_Scalar_func_partial()
{
    ScalarFunction<3> f{[](const VectorN<double, 3> &x) { return sin(x[0])*(1.0 + 0.5*x[1]*x[2]); } };

    VectorN<double, 3> pos{0.5, 0.5, 0.5};

    double der_f1 = Derivation::NDer1Partial(f, 0, pos);
    double der_f2 = Derivation::NDer1Partial(f, 1, pos);
    double der_f3 = Derivation::NDer1Partial(f, 2, pos);

    VectorN<double, 3> der_f_all = Derivation::NDer1PartialByAll<3>(f, pos, nullptr);
}

void Demo_Derivation_Scalar_func_second()
{
    ScalarFunction<3> f{[](const VectorN<double, 3> &x) { return sin(x[0])*(1.0 + 0.5*x[1]*x[2]); } };

    VectorN<double, 3> pos{0.5, 0.5, 0.5};

    double der_d2f_dxdx = Derivation::NSecDer1Partial(f, 0, 0, pos);
    double der_d2f_dydy = Derivation::NSecDer1Partial(f, 1, 1, pos);
    double der_d2f_dzdx = Derivation::NSecDer1Partial(f, 2, 0, pos);
}

/////////////////////          VECTOR FUNCTION DERIVATION           ///////////////////////
void Demo_Derivation_Vector_func()
{
    VectorFunction<3> f{[](const VectorN<double, 3> &x) { return VectorN<double, 3>{sin(x[0]), cos(x[1]), x[2]}; } };

    VectorN<double, 3> pos{0.5, 0.5, 0.5};

    double der_f11 = Derivation::NDer1Partial(f, 0, 0, pos);

    VectorN<double, 3> der_f1 = Derivation::NDer1PartialByAll(f, 0, pos);
    VectorN<double, 3> der_f2 = Derivation::NDer1PartialByAll(f, 1, pos);
    VectorN<double, 3> der_f3 = Derivation::NDer1PartialByAll(f, 2, pos);

    MatrixNM<double, 3, 3> der_f_all = Derivation::NDer1PartialAllByAll<3>(f, pos, nullptr);
}

////////////////////          PARAMETRIC CURVE DERIVATION           ///////////////////////
void Demo_Derivation_Parametric_curve()
{
    // circle in plane
    ParametricCurve<2> f{[](double t) { return VectorN<double, 2>{sin(t), cos(t)}; } };

    VectorN<Real, 2> der_f1 = Derivation::NDer1(f, 0.5);
    VectorN<Real, 2> der_f2 = Derivation::NDer2(f, 0.5);
    VectorN<Real, 2> der_f3 = Derivation::NDer4(f, 0.5);
    VectorN<Real, 2> der_f4 = Derivation::NDer6(f, 0.5);
    VectorN<Real, 2> der_f5 = Derivation::DeriveCurve<2>(f, 0.5, nullptr);
}

void Demo_Derivation_Parametric_curve_second_and_third()
{
    // helix
    ParametricCurve<3> f{[](double t) { return VectorN<double, 3>{sin(t), cos(t), t}; } };

    VectorN<Real, 3> der_f1 = Derivation::NDer1(f, 0.5);
    VectorN<Real, 3> der_f2 = Derivation::NDer2(f, 0.5);
    VectorN<Real, 3> der_f3 = Derivation::NDer6(f, 0.5);
    VectorN<Real, 3> der_f4 = Derivation::NDer8(f, 0.5);
    VectorN<Real, 3> der_f5 = Derivation::DeriveCurve<3>(f, 0.5, nullptr);

    VectorN<Real, 3> sec_der_f1 = Derivation::NSecDer1(f, 0.5);
    VectorN<Real, 3> sec_der_f2 = Derivation::NSecDer2(f, 0.5);
    // VectorN<Real, 3> sec_der_f3 = Derivation::NSecDer6(f, 0.5);
    // VectorN<Real, 3> sec_der_f4 = Derivation::NSecDer8(f, 0.5);
    VectorN<Real, 3> sec_der_f5 = Derivation::DeriveCurveSec<3>(f, 0.5, nullptr);

    VectorN<Real, 3> third_der_f1 = Derivation::NThirdDer1(f, 0.5);
    VectorN<Real, 3> third_der_f2 = Derivation::NThirdDer2(f, 0.5);
    // VectorN<Real, 3> third_der_f3 = Derivation::NThirdDer6(f, 0.5);
    // VectorN<Real, 3> third_der_f4 = Derivation::NThirdDer6(f, 0.5);
    VectorN<Real, 3> third_der_f5 = Derivation::DeriveCurveThird<3>(f, 0.5, nullptr);
}

void Demo_Derivation()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         DERIVATION                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_Derivation_func_ptr();
    Demo_Derivation_member_fun();
    Demo_Derivation_member_fun2(BigComplexClassYouCantChange{});
    Demo_Derivation_Interpolated_RealFunc();
    Demo_Derivation_RealFunc_Second_and_Third();
    Demo_Derivation_RealFunc_precision_comparison();

    Demo_Derivation_Scalar_func_partial();
    Demo_Derivation_Scalar_func_second();

    Demo_Derivation_Vector_func();

    Demo_Derivation_Parametric_curve();
    Demo_Derivation_Parametric_curve_second_and_third();
}