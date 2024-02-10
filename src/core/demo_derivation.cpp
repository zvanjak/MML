#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/InterpolatedFunction.h"
#include "core/Derivation.h"

#include "core/MetricTensor.h"
#endif

#include "../test_data/real_functions_test_bed.h"

using namespace MML;

/*****************************************************************************************/
/*************                REAL FUNCTION DERIVATION                       *************/
/*****************************************************************************************/

////////////////////////////////////////////////////////////
////                Simple func pointer
Real DemoDerRealFunc_TestFunc(Real x) 
{ 
    return sin(x)*(1.0 + 0.5*x*x); 
}

void Demo_Derivation_func_ptr()
{
    // creating a function object from a already existing (standalone) function
    RealFunction f1(DemoDerRealFunc_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](Real x) { return sin(x)*(1 + x*x / 2); } };

    // creating a function object from a std::function object
    std::function<Real(Real)> h{[](Real x) { return sin(x)*(1 + x*x / 2); } };
    RealFunctionFromStdFunc f3(h);

    Real der_f12 = Derivation::NDer2(f1, 0.5);
    Real der_f14 = Derivation::NDer4(f1, 0.5);
    Real der_f18 = Derivation::NDer8(f1, 0.5);
    Real der_f21 = Derivation::NDer1(f2, 0.5);
    Real der_f36 = Derivation::NDer6(f3, 0.5);

    // we can use default Derive routine (default set to NDer4)
    Real num_der4 = Derivation::Derive(f1, 0.5, nullptr);
}

////////////////////////////////////////////////////////////
////                class you CAN change

// Option 1 - define operator() for your class and use RealFunctionFromStdFunc
class ClassProvidingFuncToDerive
{
    private:
        Real _param;
    public:
        ClassProvidingFuncToDerive(Real param) : _param(param) { }
    
        Real operator()(Real x ) const { return _param * sin(x); }
};

void Demo_Derivation_member_fun()
{
    ClassProvidingFuncToDerive   funcObj(3.0);
    
    RealFunctionFromStdFunc g(std::function<Real(Real)>{funcObj});
    
    Real der_g = Derivation::NDer4(g, 0.5);
}

// Option 2 - make your class inherit IRearFunction interface, and use the object itself as RealFunction
class ClassProvidingFuncToDerive2 : public IRealFunction
{
    private:
        Real _param;
    public:
        ClassProvidingFuncToDerive2(Real param) : _param(param) { }
    
        // just provide operator() for your class
        Real operator()(Real x ) const { return _param * sin(x); }
};

void Demo_Derivation_member_fun1()
{
    ClassProvidingFuncToDerive2   funcObj(3.0);
    
    Real der_g = Derivation::NDer4(funcObj, 0.5);
}

////////////////////////////////////////////////////////////
////                Class you CAN'T change
class BigComplexClassYouCantChange
{
    // has data for calculating function you want to derive
};

// Option 1 - Create a helper wrapper class that will provide operator() for your class
class BigComplexDerivFunc
{
    const BigComplexClassYouCantChange &_ref;
public:
    BigComplexDerivFunc(const BigComplexClassYouCantChange &bigClass) : _ref(bigClass) { }

    Real operator()(Real x ) 
    {
        Real val = 0.0; 
        // complex calculations using _ref
        return val; 
    }
};

void Demo_Derivation_member_fun2(const BigComplexClassYouCantChange &ref)
{
    BigComplexDerivFunc         funcObj(ref);  
    RealFunctionFromStdFunc     func_to_derive(std::function<Real(Real)>{funcObj});
    
    Real der_1 = Derivation::NDer4(func_to_derive, 0.5);
    Real der_2 = Derivation::Derive(func_to_derive, 0.5, nullptr);
}

// Option 2 - Create a helper wrapper class, inherit from IRealFunction and use it as RealFunction
class BigComplexDerivFunc2 : public IRealFunction
{
    const BigComplexClassYouCantChange &_ref;
public:
    BigComplexDerivFunc2(const BigComplexClassYouCantChange &bigClass) : _ref(bigClass) { }

    Real operator()(Real x ) const
    {
        Real val = 0.0; 
        // complex calculations using _ref
        return val; 
    }
};

void Demo_Derivation_member_fun3(const BigComplexClassYouCantChange &ref)
{
    BigComplexDerivFunc2    funcObj(ref);  
    
    Real der_1 = Derivation::NDer4(funcObj, 0.5);
    Real der_2 = Derivation::Derive(funcObj, 0.5, nullptr);
}

////////////////////////////////////////////////////////////
////                MULTIPLE functions class
class BigComplexDerivMultipleFunc
{
    const BigComplexClassYouCantChange &_ref;
public:
    BigComplexDerivMultipleFunc(const BigComplexClassYouCantChange &bigClass) : _ref(bigClass) { }

    Real func1(Real x ) 
    {
        Real val = 0.0; 
        // first complex calculations using _ref
        return val; 
    }
    Real func2(Real x ) 
    {
        Real val = 0.0; 
        // second complex calculations using _ref
        return val; 
    }
};

void Demo_Derivation_member_fun4(const BigComplexClassYouCantChange &ref)
{
    BigComplexDerivMultipleFunc     funcObj(ref);  
    RealFunctionFromStdFunc         func_to_derive1(std::function<Real(Real)>{std::bind(&BigComplexDerivMultipleFunc::func1, funcObj, std::placeholders::_1)});
    RealFunctionFromStdFunc         func_to_derive2(std::function<Real(Real)>{std::bind(&BigComplexDerivMultipleFunc::func2, funcObj, std::placeholders::_1)});
    
    Real der_1 = Derivation::NDer4(func_to_derive1, 0.5);
    Real der_2 = Derivation::Derive(func_to_derive2, 0.5, nullptr);
}

////////////////////////////////////////////////////////////
////                Interpolated function
void Demo_Derivation_Interpolated_RealFunc()
{
    // creating interpolated function
    Vector<Real> x_val(100);
    Vector<Real> y_val(100);
    for( int i=0; i<100; i++ )
    {
        x_val[i] = i / 100.0;
        y_val[i] = sin(x_val[i])*(1 + x_val[i]*x_val[i] / 2);
    }

    LinearInterpRealFunc f_linear(x_val, y_val);

    Real der_f1 = Derivation::NDer1(f_linear, 0.5);
    Real der_f2 = Derivation::NDer2(f_linear, 0.5);
    Real der_f3 = Derivation::NDer4(f_linear, 0.5);
    Real der_f4 = Derivation::NDer6(f_linear, 0.5);
    Real der_f5 = Derivation::NDer8(f_linear, 0.5);
    Real der_f6 = Derivation::Derive(f_linear, 0.5, nullptr);
}

void Demo_Derivation_RealFunc_Second_and_Third()
{
    RealFunction f1(DemoDerRealFunc_TestFunc);
    RealFunction f2{[](Real x) { return sin(x)*(1 + x*x / 2); } };
    
    std::function<Real(Real)>   h{[](Real x) { return sin(x)*(1 + x*x / 2); } };
    RealFunctionFromStdFunc     f3(h);

    ClassProvidingFuncToDerive      funcObj(3.0);
    RealFunctionFromStdFunc         f4(std::function<double(double)>{funcObj});

    Real sec_der_f1 = Derivation::NSecDer1(f1, 0.5);
    Real sec_der_f2 = Derivation::NSecDer2(f2, 0.5);
    Real sec_der_f3 = Derivation::NSecDer4(f3, 0.5);
    Real sec_der_f4 = Derivation::NSecDer6(f4, 0.5);
    Real sec_der_f5 = Derivation::DeriveSec(f1, 0.5, nullptr);

    Real third_der_f1 = Derivation::NThirdDer1(f1, 0.5);
    Real third_der_f2 = Derivation::NThirdDer2(f2, 0.5);
    Real third_der_f3 = Derivation::NThirdDer4(f3, 0.5);
    Real third_der_f4 = Derivation::NThirdDer6(f4, 0.5);
    Real third_der_f5 = Derivation::DeriveThird(f1, 0.5, nullptr);
}

// TODO - simple precision comparison for selected function
// simple example for a single function of various precisions
void Demo_Derivation_RealFunc_precision_comparison()
{

}

/*****************************************************************************************/
/*************               SCALAR FUNCTION DERIVATION                      *************/
/*****************************************************************************************/
void Demo_Derivation_Scalar_func_partial()
{
    ScalarFunction<3> f{[](const VectorN<Real, 3> &x) { return sin(x[0])*(1 + x[1]*x[2] / 2); } };

    VectorN<Real, 3> pos{0.5, 0.5, 0.5};

    Real der_f1 = Derivation::NDer1Partial(f, 0, pos);
    Real der_f2 = Derivation::NDer1Partial(f, 1, pos);
    Real der_f3 = Derivation::NDer1Partial(f, 2, pos);

    VectorN<Real, 3> der_f_all = Derivation::NDer1PartialByAll<3>(f, pos, nullptr);
}

void Demo_Derivation_Scalar_func_second()
{
    ScalarFunction<3> f{[](const VectorN<Real, 3> &x) { return sin(x[0])*(1 + x[1]*x[2] / 2); } };

    VectorN<Real, 3> pos{0.5, 0.5, 0.5};

    Real der_d2f_dxdx = Derivation::NSecDer1Partial(f, 0, 0, pos);
    Real der_d2f_dydy = Derivation::NSecDer1Partial(f, 1, 1, pos);
    Real der_d2f_dzdx = Derivation::NSecDer1Partial(f, 2, 0, pos);
}

/*****************************************************************************************/
/*************               VECTOR FUNCTION DERIVATION                      *************/
/*****************************************************************************************/
void Demo_Derivation_Vector_func()
{
    VectorFunction<3> f{[](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{sin(x[0]), cos(x[1]), x[2]}; } };

    VectorN<Real, 3> pos{0.5, 0.5, 0.5};

    double der_f11 = Derivation::NDer1Partial(f, 0, 0, pos);

    VectorN<Real, 3> der_f1 = Derivation::NDer1PartialByAll(f, 0, pos);
    VectorN<Real, 3> der_f2 = Derivation::NDer1PartialByAll(f, 1, pos);
    VectorN<Real, 3> der_f3 = Derivation::NDer1PartialByAll(f, 2, pos);

    MatrixNM<Real, 3, 3> der_f_all = Derivation::NDer1PartialAllByAll<3>(f, pos, nullptr);

    VectorN<Real, 3> der_f1_2 = Derivation::NDer2PartialByAll(f, 0, pos);
    VectorN<Real, 3> der_f2_2 = Derivation::NDer2PartialByAll(f, 1, pos);
    VectorN<Real, 3> der_f3_2 = Derivation::NDer2PartialByAll(f, 2, pos);

    MatrixNM<Real, 3, 3> der_f_all_2 = Derivation::NDer2PartialAllByAll<3>(f, pos, nullptr);    
    
}

/*****************************************************************************************/
/*************               TENSOR FIELD  DERIVATION                      *************/
/*****************************************************************************************/
void Demo_Derivation_Tensor_field()
{
    MetricTensorCartesian<3> metricCart;
    VectorN<Real, 3> pos{0.5, 0.5, 0.5};

    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
            Real g_der_ij_x = Derivation::NDer1Partial(metricCart, i, j, 0, pos);
            Real g_der_ij_y = Derivation::NDer2Partial(metricCart, i, j, 1, pos);
            Real g_der_ij_z = Derivation::NDer4Partial(metricCart, i, j, 2, pos);
        }
}

/*****************************************************************************************/
/*************              PARAMETRIC CURVE DERIVATION                      *************/
/*****************************************************************************************/
void Demo_Derivation_Parametric_curve()
{
    // circle in plane
    ParametricCurve<2> f{[](Real t) { return VectorN<Real, 2>{sin(t), cos(t)}; } };

    VectorN<Real, 2> der_f1 = Derivation::NDer1(f, 0.5);
    VectorN<Real, 2> der_f2 = Derivation::NDer2(f, 0.5);
    VectorN<Real, 2> der_f3 = Derivation::NDer4(f, 0.5);
    VectorN<Real, 2> der_f4 = Derivation::NDer6(f, 0.5);
    VectorN<Real, 2> der_f5 = Derivation::DeriveCurve<2>(f, 0.5, nullptr);
}

void Demo_Derivation_Parametric_curve_second_and_third()
{
    // helix
    ParametricCurve<3> f{[](Real t) { return VectorN<Real, 3>{sin(t), cos(t), t}; } };

    VectorN<Real, 3> der_f1 = Derivation::NDer1(f, 0.5);
    VectorN<Real, 3> der_f2 = Derivation::NDer2(f, 0.5);
    VectorN<Real, 3> der_f3 = Derivation::NDer6(f, 0.5);
    VectorN<Real, 3> der_f4 = Derivation::NDer8(f, 0.5);
    VectorN<Real, 3> der_f5 = Derivation::DeriveCurve<3>(f, 0.5, nullptr);

    VectorN<Real, 3> sec_der_f1 = Derivation::NSecDer1(f, 0.5);
    VectorN<Real, 3> sec_der_f2 = Derivation::NSecDer2(f, 0.5);
    VectorN<Real, 3> sec_der_f3 = Derivation::NSecDer6(f, 0.5);
    VectorN<Real, 3> sec_der_f4 = Derivation::NSecDer8(f, 0.5);
    VectorN<Real, 3> sec_der_f5 = Derivation::DeriveCurveSec<3>(f, 0.5, nullptr);

    VectorN<Real, 3> third_der_f1 = Derivation::NThirdDer1(f, 0.5);
    VectorN<Real, 3> third_der_f2 = Derivation::NThirdDer2(f, 0.5);
    VectorN<Real, 3> third_der_f3 = Derivation::NThirdDer6(f, 0.5);
    VectorN<Real, 3> third_der_f4 = Derivation::NThirdDer6(f, 0.5);
    VectorN<Real, 3> third_der_f5 = Derivation::DeriveCurveThird<3>(f, 0.5, nullptr);
}

/*****************************************************************************************/
/*************                    DEMO  DERIVATION                           *************/
/*****************************************************************************************/
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

    Demo_Derivation_Scalar_func_partial();
    Demo_Derivation_Scalar_func_second();

    Demo_Derivation_Vector_func();
    
    Demo_Derivation_Tensor_field();

    Demo_Derivation_Parametric_curve();
    Demo_Derivation_Parametric_curve_second_and_third();

    Demo_Derivation_RealFunc_precision_comparison();
}
