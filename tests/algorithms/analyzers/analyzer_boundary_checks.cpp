#include <cmath>
#pragma push_macro("INFINITE")
#pragma push_macro("DOMAIN")
#undef INFINITE
#undef DOMAIN

#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#include <mml/algorithms/Analyzers/FieldAnalyzers.h>

#include <string>
#include <type_traits>

#if defined(MML_FUNCTION_ANALYSIS_REPORTS_H) || defined(MML_MATRIX_ANALYSIS_REPORTS_H)
#error Numerical analyzers must not include presentation headers.
#endif

#if defined(INFINITE) || defined(DOMAIN)
#error Numerical analyzer headers must not introduce INFINITE or DOMAIN macros.
#endif

#define INFINITE 0x24681357u
#define DOMAIN 317
#include <mml/tools/analyzers/FunctionAnalysisReports.h>
#if !defined(INFINITE) || !defined(DOMAIN)
#error Reports included after numerical analyzers must preserve caller macros.
#else
static_assert(INFINITE == 0x24681357u);
static_assert(DOMAIN == 317);
static_assert(static_cast<int>(MML::DiscontinuityType::Infinite) == 2);
#endif

template<class Analyzer>
concept HasPointPrinter = requires(const Analyzer& analyzer) {
	analyzer.PrintPointAnalysis(MML::Real(0));
};

template<class Analyzer>
concept HasIntervalPrinter = requires(const Analyzer& analyzer) {
	analyzer.PrintIntervalAnalysis(MML::Real(0), MML::Real(1), 10);
};

template<class Analyzer>
concept HasDetailedPrinter = requires(const Analyzer& analyzer) {
	analyzer.PrintDetailedIntervalAnalysis(MML::Real(0), MML::Real(1), 10);
};

template<class Analyzer>
concept HasContinuityPrinter = requires(Analyzer& analyzer) {
	analyzer.PrintContinuityAnalysis(MML::Real(0), MML::Real(1), 10);
};

template<class Result>
concept HasReportText = requires(const Result& result) {
	result.report;
};

static_assert(!HasPointPrinter<MML::RealFunctionAnalyzer>);
static_assert(!HasIntervalPrinter<MML::RealFunctionAnalyzer>);
static_assert(!HasDetailedPrinter<MML::RealFunctionAnalyzer>);
static_assert(!HasContinuityPrinter<MML::RealFunctionAnalyzer>);
static_assert(!HasReportText<MML::MatrixAlg::MatrixAnalysis<MML::Real>>);
static_assert(!std::is_constructible_v<MML::RealFunctionAnalyzer, const MML::IRealFunction&, std::string>);

#pragma pop_macro("DOMAIN")
#pragma pop_macro("INFINITE")
