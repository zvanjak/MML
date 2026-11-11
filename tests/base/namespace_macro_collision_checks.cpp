constexpr int ApplicationRealLiteral(int value) { return value + 1; }
#define REAL(value) ApplicationRealLiteral(value)

#include <mml/MMLBase.h>
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>

static_assert(REAL(3) == 4);
static_assert(std::is_same_v<MML::Complex, std::complex<MML::Real>>);

#undef REAL

#pragma push_macro("INFINITE")
#pragma push_macro("DOMAIN")
#undef INFINITE
#undef DOMAIN
#define INFINITE(value) ((value) + 17)
#define DOMAIN(value) ((value) ^ 8)

#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#include <mml/tools/analyzers/FunctionAnalysisReports.h>

#if !defined(INFINITE) || !defined(DOMAIN)
#error Function analyzer headers must preserve function-like caller macros.
#else
static_assert(INFINITE(3) == 20);
static_assert(DOMAIN(3) == 11);
#endif

#pragma pop_macro("DOMAIN")
#pragma pop_macro("INFINITE")
