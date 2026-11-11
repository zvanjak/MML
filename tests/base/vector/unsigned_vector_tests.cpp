#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#endif

#include <iomanip>
#include <limits>
#include <sstream>
#include <type_traits>

using namespace MML;

TEMPLATE_TEST_CASE("Unsigned magnitudes preserve the configured Real return contract",
                   "[Abs][UnsignedVector]", unsigned char, unsigned short,
                   unsigned int, unsigned long, unsigned long long)
{
	STATIC_REQUIRE(std::is_same_v<decltype(Abs(TestType{})), Real>);
	const TestType maximum = std::numeric_limits<TestType>::max();
	const TestType highBit = maximum / 2 + 1;
	REQUIRE(Abs(TestType(0)) == Real(0));
	REQUIRE(Abs(TestType(1)) == Real(1));
	REQUIRE(Abs(highBit) == static_cast<Real>(highBit));
	REQUIRE(Abs(maximum) == static_cast<Real>(maximum));

	const Vector<TestType> dynamic{TestType(3), TestType(4), TestType(0)};
	const VectorN<TestType, 3> fixed{TestType(3), TestType(4), TestType(0)};
	REQUIRE(IsNearlyEqual(dynamic.NormL1(), Real(7)));
	REQUIRE(IsNearlyEqual(dynamic.NormL2(), Real(5)));
	REQUIRE(IsNearlyEqual(dynamic.NormLInf(), Real(4)));
	REQUIRE(IsNearlyEqual(fixed.NormL1(), Real(7)));
	REQUIRE(IsNearlyEqual(fixed.NormL2(), Real(5)));
	REQUIRE(IsNearlyEqual(fixed.NormLInf(), Real(4)));
	REQUIRE(Vector<TestType>{TestType(0)}.IsZero());
	REQUIRE_FALSE(dynamic.IsZero());
	REQUIRE_FALSE(dynamic.IsNearZero());
}

TEMPLATE_TEST_CASE("Unsigned vectors print and stringify native integer values without signed narrowing",
                   "[UnsignedVector][Vector][VectorN][Printing]",
                   unsigned short, unsigned int, unsigned long, unsigned long long)
{
	const TestType maximum = std::numeric_limits<TestType>::max();
	const Vector<TestType> dynamic{maximum, TestType(0)};
	const VectorN<TestType, 2> fixed{maximum, TestType(0)};
	const std::string compact = "[" + std::to_string(maximum) + ", 0]";
	std::ostringstream printed;
	REQUIRE(&dynamic.Print(printed, 0, 6, Real(2)) == &printed);
	REQUIRE(printed.str() == compact);
	REQUIRE(dynamic.to_string(0, 6) == compact);
	std::ostringstream fixedPrinted;
	fixed.Print(fixedPrinted, 0, 6, Real(2));
	REQUIRE(fixedPrinted.str() == compact);
	REQUIRE(fixed.to_string(0, 6) == compact);

	std::ostringstream expected;
	expected << '[' << std::setw(Defaults::VectorPrintWidth) << maximum
	         << ", " << std::setw(Defaults::VectorPrintWidth) << Real(0) << ']';
	std::ostringstream streamed;
	streamed << dynamic;
	REQUIRE(streamed.str() == expected.str());
	REQUIRE(dynamic == Vector<TestType>{maximum, TestType(0)});

	const Vector<TestType> small{TestType(3), TestType(4), TestType(0)};
	std::ostringstream thresholded;
	small.Print(thresholded, 0, 6, Real(3));
	REQUIRE(thresholded.str() == "[0, 4, 0]");
}

TEST_CASE("Unsigned magnitude support preserves signed, floating and complex behavior", "[Abs][UnsignedVector]")
{
	REQUIRE(Abs(-7) == Real(7));
	REQUIRE(Abs(-Real(3)) == Real(3));
	REQUIRE(Abs(Complex(3, 4)) == Real(5));
	REQUIRE(std::isinf(Abs(-std::numeric_limits<Real>::infinity())));
	REQUIRE(std::isnan(Abs(std::numeric_limits<Real>::quiet_NaN())));
	REQUIRE_FALSE(std::signbit(Abs(-Real(0))));
}
