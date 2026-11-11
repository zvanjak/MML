#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLConcepts.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#endif

#include <string>

using namespace MML;

namespace
{
	template<class VectorType>
	concept HasUnitVector = requires {
		VectorType::UnitVector(3, 1);
	};
}

namespace MML::Tests::Base::ConceptsTests
{
	TEST_CASE("MMLConcepts - scalar concepts classify supported numeric types", "[Concepts][C++20]")
	{
		static_assert(MMLReal<float>);
		static_assert(MMLReal<double>);
		static_assert(MMLReal<long double>);
		static_assert(!MMLReal<int>);

		static_assert(MMLComplex<std::complex<float>>);
		static_assert(MMLComplex<std::complex<double>>);
		static_assert(MMLComplex<std::complex<long double>>);
		static_assert(MMLComplex<Complex>);
		static_assert(!MMLComplex<Real>);

		static_assert(MMLArithmetic<int>);
		static_assert(MMLArithmetic<Real>);
		static_assert(!MMLArithmetic<Complex>);

		static_assert(MMLNumeric<int>);
		static_assert(MMLNumeric<Real>);
		static_assert(MMLNumeric<Complex>);
		static_assert(!MMLNumeric<std::string>);

		static_assert(MMLScalar<Real>);
		static_assert(MMLScalar<Complex>);
		static_assert(!MMLScalar<int>);
		static_assert(RealFunctionCallable<Real(*)(Real)>);
		static_assert(RealFunctionCallable<decltype([](Real x) -> Real { return x; })>);
		static_assert(!RealFunctionCallable<decltype([](const VectorN<Real, 2>& x) -> Real { return x[0]; })>);

		REQUIRE(is_MML_simple_numeric<Real>);
		REQUIRE(is_MML_simple_numeric<Complex>);
		REQUIRE_FALSE(is_MML_simple_numeric<std::string>);
	}

	TEST_CASE("MMLConcepts - shape concepts accept MML vector and matrix storage", "[Concepts][C++20]")
	{
		static_assert(VectorLike<Vector<Real>>);
		static_assert(VectorLike<VectorN<Real, 3>>);
		static_assert(MatrixLike<Matrix<Real>>);

		REQUIRE(Vector<Real>(3).size() == 3);
		REQUIRE(Matrix<Real>(2, 3).rows() == 2);
	}

	TEST_CASE("MMLConcepts - integer scalar conversion preserves the destination type", "[Concepts][precision]")
	{
		static_assert(std::is_same_v<decltype(ScalarDetail::ScalarFromInteger<float>(1)), float>);
		static_assert(std::is_same_v<decltype(ScalarDetail::ScalarFromInteger<double>(1)), double>);
		static_assert(std::is_same_v<decltype(ScalarDetail::ScalarFromInteger<std::complex<float>>(1)), std::complex<float>>);
		CHECK(ScalarDetail::ScalarFromInteger<float>(-3) == -3.0f);
		CHECK(ScalarDetail::ScalarFromInteger<double>(16777217) == 16777217.0);
		CHECK(ScalarDetail::ScalarFromInteger<std::complex<float>>(-3) == std::complex<float>(-3.0f, 0.0f));
		CHECK(ScalarDetail::ScalarFromInteger<std::complex<double>>(16777217) == std::complex<double>(16777217.0, 0.0));
	}

	TEST_CASE("MMLConcepts - requires clauses constrain arithmetic-only vector APIs", "[Concepts][C++20]")
	{
		static_assert(HasUnitVector<Vector<Real>>);
		static_assert(HasUnitVector<Vector<int>>);
		static_assert(!HasUnitVector<Vector<std::string>>);

		Vector<Real> e1 = Vector<Real>::UnitVector(3, 1);
		REQUIRE(e1[0] == static_cast<MML::Real>(0.0));
		REQUIRE(e1[1] == static_cast<MML::Real>(1.0));
		REQUIRE(e1[2] == static_cast<MML::Real>(0.0));
	}
}
