#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Function.h"
#endif

using namespace MML;

namespace MML::Tests::Core::FunctionTests
{
	/*********************************************************************/
	/*****                   RealFunction tests                      *****/
	/*********************************************************************/

	Real FunctionTests_TestFunc(Real x)
	{
		return sin(x) * (1.0 + 0.5 * x * x);
	}

	TEST_CASE("Test_Function_func_pointer", "[simple]") {
		// creating a function object from a already existing (standalone) function
		RealFunction f1(FunctionTests_TestFunc);

		// or creating a function object directly
		RealFunction f2{ [](Real x) { return sin(x) * (1 + x * x / 2); } };

		REQUIRE(f1(1.0) == FunctionTests_TestFunc(1.0));
		REQUIRE(f2(1.0) == FunctionTests_TestFunc(1.0));
	}

	class ClassProvidingFuncToDerive
	{
	private:
		double _param;
	public:
		ClassProvidingFuncToDerive(double param) : _param(param) { }
		void setParam(double inParam) { _param = inParam; }

		Real operator()(Real x) { return _param * sin(x); }
	};

	TEST_CASE("Test_Function_class_obj_overload_op()", "[simple]")
	{
		ClassProvidingFuncToDerive   funcObj(3.0);

		RealFunctionFromStdFunc f1(std::function<double(double)>{funcObj});

		REQUIRE(3.0 * sin(1.0) == f1(1.0));
		REQUIRE(3.0 * sin(1.0) == funcObj(1.0));

		REQUIRE(f1(-2.0) == funcObj(-2.0));
		REQUIRE(f1(1.0) == funcObj(1.0));
		REQUIRE(f1(5.0) == funcObj(5.0));

		// VERY IMPORTANT - our f1 has a COPY of funcObj, so if we change funcObj, f1 will not change
		funcObj.setParam(6.0);

		REQUIRE(3.0 * sin(1.0) == f1(1.0));   // WE HAVE AN OLD VALUE OF PARAMETAR!
		REQUIRE(6.0 * sin(1.0) != f1(1.0));
		REQUIRE(6.0 * sin(1.0) == funcObj(1.0));

		REQUIRE(f1(1.0) != funcObj(1.0));
		REQUIRE(f1(5.0) != funcObj(5.0));

		// if we create a new copy f2, it gets the new value of parametar
		RealFunctionFromStdFunc f2(std::function<double(double)>{funcObj});

		REQUIRE(6.0 * sin(1.0) == f2(1.0));
		REQUIRE(f2(1.0) == funcObj(1.0));
		REQUIRE(f2(5.0) == funcObj(5.0));
	}

	class ClassProvidingFuncToDerive2 : public IRealFunction
	{
	private:
		Real _param;
	public:
		ClassProvidingFuncToDerive2(Real param) : _param(param) { }
		void setParam(Real inParam) { _param = inParam; }

		Real operator()(Real x) const { return _param * sin(x); }
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
		// Has data and functionality for calculating function value we want
		// BUT, we have to manually implement exact calculation we need
	public:
		double _param;
		Real complexCalc(Real x) const { return cos(x); }
	};

	class BigComplexDerivFunc
	{
		const BigComplexClassYouCantChange& _ref;
	public:
		BigComplexDerivFunc(const BigComplexClassYouCantChange& bigClass) : _ref(bigClass) { }

		// here we implement function with behavior we want
		Real operator()(Real x)
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

		REQUIRE(3.0 * cos(1.0) == f1(1.0));

		REQUIRE(f1(1.0) == funcObj(1.0));
		REQUIRE(f1(5.0) == funcObj(5.0));

		// f1 object is valid even if you changed the referenced class object
		ref._param = 6.0;
		REQUIRE(f1(1.0) == funcObj(1.0));
		REQUIRE(f1(5.0) == funcObj(5.0));
	}

	/*********************************************************************/
	/*****                 ScalarFunction tests                    *****/
	/*********************************************************************/

	/*********************************************************************/
	/*****                 VectorFunction tests                    *****/
	/*********************************************************************/

	/*********************************************************************/
	/*****                 VectorNMFunction tests                    *****/
	/*********************************************************************/

	/*********************************************************************/
	/*****                 ParametricCurve tests                    *****/
	/*********************************************************************/

	/*********************************************************************/
	/*****                 ParametricSurface tests                    *****/
	/*********************************************************************/

} // namespace Tests::Core::FunctionTests