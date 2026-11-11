#include <catch2/catch_all.hpp>
#include <cmath>
#include <sstream>

#if defined(_WIN32) && defined(MML_SINGLE_HEADER_INCLUDE)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#pragma push_macro("INFINITE")
#pragma push_macro("DOMAIN")
#undef INFINITE
#undef DOMAIN
#define INFINITE 0x13572468u
#define DOMAIN 731

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/base/Function.h>
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#include <mml/tools/analyzers/FunctionAnalysisReports.h>
#endif

#if !defined(INFINITE) || !defined(DOMAIN)
#error Function analyzer headers must preserve caller INFINITE and DOMAIN macros.
#else
static_assert(INFINITE == 0x13572468u);
static_assert(DOMAIN == 731);
static_assert(static_cast<int>(MML::DiscontinuityType::JUMP) == 0);
static_assert(static_cast<int>(MML::DiscontinuityType::REMOVABLE) == 1);
static_assert(static_cast<int>(MML::DiscontinuityType::Infinite) == 2);
static_assert(static_cast<int>(MML::DiscontinuityType::OSCILLATORY) == 3);
static_assert(static_cast<int>(MML::DiscontinuityType::UNKNOWN) == 4);

#pragma push_macro("INFINITE")
#undef INFINITE
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4996)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
static_assert(MML::DiscontinuityType::INFINITE == MML::DiscontinuityType::Infinite);
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
#pragma pop_macro("INFINITE")
static_assert(INFINITE == 0x13572468u);

TEST_CASE("Function analyzer and reports preserve object-like caller macros",
          "[FunctionAnalyzer][function-analyzer-macro-contract]")
{
	using namespace MML;
	CHECK(INFINITE == 0x13572468u);
	CHECK(DOMAIN == 731);
	RealFunctionFromStdFunc pole(+[](Real x) { return Precision::OverflowDetectionThreshold / x; });
	const auto point = RealFunctionAnalyzer(pole).ClassifyDiscontinuity(0);
	CHECK(point.type == DiscontinuityType::Infinite);
	CHECK(std::isinf(point.leftLimit));
	CHECK(std::isinf(point.rightLimit));
	const ContinuityAnalysisResult analysis{-1, 1, 1, Precision::LooseTolerance, {point}};
	std::ostringstream report;
	AnalysisReports::PrintContinuityAnalysis(report, analysis);
	CHECK(report.str().find("Estimated type: INFINITE") != std::string::npos);
}
#endif

#pragma pop_macro("DOMAIN")
#pragma pop_macro("INFINITE")
