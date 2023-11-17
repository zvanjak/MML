#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Function.h"
#endif

using namespace MML;

namespace Tests::Core::FunctionTests
{

double FunctionTests_TestFunc(double x) 
{ 
    return sin(x)*(1.0 + 0.5*x*x); 
}

TEST_CASE("Test_Function_func_pointer", "[simple]") {
    // creating a function object from a already existing (standalone) function
    RealFunction f1(FunctionTests_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    REQUIRE(f1(1.0) == FunctionTests_TestFunc(1.0));
    REQUIRE(f2(1.0) == FunctionTests_TestFunc(1.0));
}

class ClassProvidingFuncToDerive
{
    private:
        double _param;
    public:
        ClassProvidingFuncToDerive(double param) : _param(param) { }
        void setParam(double inParam ) { _param = inParam; }
    
        double operator()(double x ) { return _param * sin(x); }
};

TEST_CASE("Test_Function_class_obj_overload_op()", "[simple]") 
{
    ClassProvidingFuncToDerive   funcObj(3.0);
    
    RealFunctionFromStdFunc f1(std::function<double(double)>{funcObj});
    
    REQUIRE(f1(1.0) == funcObj(1.0));
    REQUIRE(f1(5.0) == funcObj(5.0));

    // VERY IMPORTANT!!!!!!!!!!!!!!!!
    funcObj.setParam(6.0);
    REQUIRE(f1(1.0) != funcObj(1.0));
    REQUIRE(f1(5.0) != funcObj(5.0));

    // BUT!!!
    RealFunctionFromStdFunc f2(std::function<double(double)>{funcObj});
    REQUIRE(f2(1.0) == funcObj(1.0));
    REQUIRE(f2(5.0) == funcObj(5.0));
}

class ClassProvidingFuncToDerive2 : public IRealFunction
{
    private:
        double _param;
    public:
        ClassProvidingFuncToDerive2(double param) : _param(param) { }
        void setParam(double inParam ) { _param = inParam; }

        double operator()(double x ) const { return _param * sin(x); }
};

TEST_CASE("Test_Function_class_inheriting_IRealFunction", "[simple]") 
{
    ClassProvidingFuncToDerive2   funcObj(3.0);

    IRealFunction& f1 = funcObj;
    
    REQUIRE(f1(1.0) == funcObj(1.0));
    REQUIRE(f1(5.0) == funcObj(5.0));

    funcObj.setParam(6.0);
    REQUIRE(f1(1.0) == funcObj(1.0));
    REQUIRE(f1(5.0) == funcObj(5.0));
}

class BigComplexClassYouCantChange
{
    // has data for calculating function you want to derive
    public:
        double _param;
        double complexCalc(double x) const { return cos(x); }
};

class BigComplexDerivFunc
{
    const BigComplexClassYouCantChange &_ref;
public:
    BigComplexDerivFunc(const BigComplexClassYouCantChange &bigClass) : _ref(bigClass) { }

    double operator()(double x ) 
    {
        return _ref._param * _ref.complexCalc(x);
    }
};

TEST_CASE("Test_Function_class_unchangeable_overload_op()", "[simple]") 
{
    BigComplexClassYouCantChange    ref;
    ref._param = 3.0;

    BigComplexDerivFunc             funcObj(ref);  
    RealFunctionFromStdFunc         f1(std::function<double(double)>{funcObj});
    
    REQUIRE(f1(1.0) == funcObj(1.0));
    REQUIRE(f1(5.0) == funcObj(5.0));

    // f1 object is valid even if you changed the referenced class object
    ref._param = 6.0;
    REQUIRE(f1(1.0) == funcObj(1.0));
    REQUIRE(f1(5.0) == funcObj(5.0));
}

} // namespace Tests::Core::FunctionTests