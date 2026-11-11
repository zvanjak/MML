struct Real { int applicationValue; };
struct Complex { int applicationValue; };
inline int GetEnv(const char*) { return 7; }

#include <MML.h>
#include <catch2/catch_test_macros.hpp>

#ifdef REAL
#error MML single header must not define a global REAL macro.
#endif

TEST_CASE("Namespaced single-header scalars coexist with application global names", "[single-header][namespace]")
{
	::Real application{1};
	::Complex applicationComplex{2};
	MML::Vector<MML::Real> values{MML::Real(3), MML::Real(4)};
	REQUIRE(application.applicationValue == 1);
	REQUIRE(applicationComplex.applicationValue == 2);
	REQUIRE(::GetEnv("application-owned") == 7);
	REQUIRE(MML::IsNearlyEqual(values.NormL2(), MML::Real(5)));
	STATIC_REQUIRE(std::is_same_v<MML::Complex, std::complex<MML::Real>>);
}
