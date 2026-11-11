#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Function.h>
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#include <mml/tools/analyzers/FunctionAnalysisReports.h>
#include <mml/tools/analyzers/MatrixAnalysisReports.h>
#endif

#include <array>
#include <functional>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <streambuf>

using namespace MML;

namespace {
	struct CommaDecimal : std::numpunct<char> {
		char do_decimal_point() const override { return ','; }
	};

	struct ThrowingOutput : std::streambuf {
		std::streamsize xsputn(const char*, std::streamsize) override {
			throw std::runtime_error("report output failed");
		}
		int_type overflow(int_type) override {
			throw std::runtime_error("report output failed");
		}
	};
}

TEST_CASE("Analysis reports consume snapshots without evaluating functions", "[AnalysisReports][FunctionAnalyzer]")
{
	int evaluations = 0;
	RealFunctionFromStdFunc function([&](Real x) {
		++evaluations;
		return x * x;
	});
	const RealFunctionAnalyzer analyzer(function);
	const auto point = analyzer.AnalyzePoint(Real(1));
	const auto interval = analyzer.AnalyzeInterval(Real(-1), Real(2), 3);
	const auto details = analyzer.AnalyzePoints(Real(-1), Real(2), 3);
	const auto continuity = analyzer.AnalyzeContinuity(Real(-1), Real(2), 3);
	const int before = evaluations;
	REQUIRE(before > 0);

	REQUIRE(interval.minValue == Real(0));
	REQUIRE(interval.maxValue == Real(1));
	REQUIRE(details.size() == 3);
	REQUIRE(details.back().x == Real(1));
	REQUIRE(point.tolerance == Precision::LooseTolerance);
	REQUIRE(interval.tolerance == Precision::LooseTolerance);
	REQUIRE(continuity.numSubintervals == 3);

	std::ostringstream first;
	std::ostringstream repeated;
	for (auto* output : {&first, &repeated}) {
		AnalysisReports::PrintPointAnalysis(*output, point, "x^2");
		AnalysisReports::PrintIntervalAnalysis(*output, interval);
		AnalysisReports::PrintDetailedIntervalAnalysis(*output, details);
		AnalysisReports::PrintContinuityAnalysis(*output, continuity);
	}
	REQUIRE(evaluations == before);
	REQUIRE(first.str() == repeated.str());
	REQUIRE(first.str().find("Sampled minimum: 0") != std::string::npos);
	REQUIRE(first.str().find("Sampled maximum: 1") != std::string::npos);
	REQUIRE(first.str().find("right endpoint excluded") != std::string::npos);
	REQUIRE(first.str().find("not proofs") != std::string::npos);
}

TEST_CASE("Every analyzer reporter restores formatting even on output failure", "[AnalysisReports][StreamState]")
{
	RealFunction function(+[](Real x) { return x * x; });
	const RealFunctionAnalyzer analyzer(function);
	const auto point = analyzer.AnalyzePoint(Real(1));
	const auto interval = analyzer.AnalyzeInterval(Real(0), Real(1), 4);
	const auto details = analyzer.AnalyzePoints(Real(0), Real(1), 4);
	const auto continuity = analyzer.AnalyzeContinuity(Real(0), Real(1), 4);
	const auto matrix = MatrixAnalyzer<Real>(Matrix<Real>{2, 2, {2, 0, 0, 3}}).Analyze();
	const std::array<std::function<void(std::ostream&)>, 5> reporters{{
		[&](std::ostream& out) { AnalysisReports::PrintPointAnalysis(out, point); },
		[&](std::ostream& out) { AnalysisReports::PrintIntervalAnalysis(out, interval); },
		[&](std::ostream& out) { AnalysisReports::PrintDetailedIntervalAnalysis(out, details); },
		[&](std::ostream& out) { AnalysisReports::PrintContinuityAnalysis(out, continuity); },
		[&](std::ostream& out) { AnalysisReports::PrintMatrixAnalysis(out, matrix); }
	}};
	for (const auto& report : reporters) {
		for (bool fail : {false, true}) {
			ThrowingOutput buffer;
			std::ostringstream successful;
			std::ostream failing(&buffer);
			std::ostream& out = fail ? failing : successful;
			out.imbue(std::locale(std::locale::classic(), new CommaDecimal));
			out.setf(std::ios::scientific | std::ios::showpos | std::ios::uppercase);
			out.precision(3);
			out.width(19);
			out.fill('#');
			out.exceptions(std::ios::badbit | std::ios::failbit);
			const auto locale = out.getloc();
			const auto flags = out.flags();
			const auto exceptions = out.exceptions();
			if (fail)
				REQUIRE_THROWS(report(out));
			else {
				report(out);
				std::ostringstream expected;
				report(expected);
				REQUIRE(successful.str() == expected.str());
			}
			REQUIRE(out.getloc() == locale);
			REQUIRE(out.flags() == flags);
			REQUIRE(out.precision() == 3);
			REQUIRE(out.width() == 19);
			REQUIRE(out.fill() == '#');
			REQUIRE(out.exceptions() == exceptions);
		}
	}
}

TEST_CASE("Reports distinguish sampled evidence from mathematical certification", "[AnalysisReports][FunctionAnalyzer]")
{
	RealFunction narrowJump(+[](Real x) {
		return x > Real(1) / 9 && x < Real(1) / 7 ? Real(-100) : x;
	});
	const RealFunctionAnalyzer analyzer(narrowJump);
	REQUIRE(analyzer.IsMonotonicOnSamples(Real(0), Real(1), 4));
	REQUIRE(narrowJump(Real(1) / 8) < narrowJump(Real(0)));
	const auto result = analyzer.AnalyzeContinuity(Real(0), Real(1), 4);
	REQUIRE(result.candidates.empty());
	std::ostringstream out;
	AnalysisReports::PrintContinuityAnalysis(out, result, "narrow jump");
	REQUIRE(out.str().find("No discontinuity candidates were detected") != std::string::npos);
	REQUIRE(out.str().find("not proofs") != std::string::npos);
	REQUIRE(out.str().find("CONTINUOUS throughout") == std::string::npos);
}

TEST_CASE("Continuity reporting covers every candidate classification", "[AnalysisReports][FunctionAnalyzer]")
{
	ContinuityAnalysisResult result{Real(-1), Real(1), 16, Precision::LooseTolerance, {}};
	for (auto type : {DiscontinuityType::JUMP, DiscontinuityType::REMOVABLE,
	                  DiscontinuityType::Infinite, DiscontinuityType::OSCILLATORY,
	                  DiscontinuityType::UNKNOWN})
		result.candidates.push_back({Real(0), type, Real(-1), Real(1),
		                             std::numeric_limits<Real>::quiet_NaN(), Real(2)});
	std::ostringstream out;
	AnalysisReports::PrintContinuityAnalysis(out, result);
	for (const char* type : {"JUMP", "REMOVABLE", "INFINITE", "OSCILLATORY", "UNKNOWN"})
		REQUIRE(out.str().find(type) != std::string::npos);
	REQUIRE(out.str().find("estimated jump: 2") != std::string::npos);
	REQUIRE(out.str().find("Left limit estimate: -1") != std::string::npos);
	REQUIRE(out.str().find("Right limit estimate: 1") != std::string::npos);
}

TEST_CASE("Reports represent absent finite samples and empty detailed results explicitly", "[AnalysisReports][FunctionAnalyzer]")
{
	RealFunction nonFinite(+[](Real) { return std::numeric_limits<Real>::quiet_NaN(); });
	const auto result = RealFunctionAnalyzer(nonFinite).AnalyzeInterval(Real(0), Real(1), 4);
	std::ostringstream out;
	AnalysisReports::PrintIntervalAnalysis(out, result);
	REQUIRE(out.str().find("unavailable (no finite samples)") != std::string::npos);
	REQUIRE(out.str().find("Sampled minimum:") == std::string::npos);
	REQUIRE(out.str().find("Non-finite sample locations:") != std::string::npos);
	REQUIRE(out.str().find("Monotonic on samples: unavailable") != std::string::npos);
	std::ostringstream empty;
	AnalysisReports::PrintDetailedIntervalAnalysis(empty, {});
	REQUIRE(empty.str().find("(0 points)") != std::string::npos);
}

TEST_CASE("Structured continuity search rejects invalid search metadata", "[AnalysisReports][FunctionAnalyzer]")
{
	RealFunction function(+[](Real x) { return x; });
	const RealFunctionAnalyzer analyzer(function);
	REQUIRE_THROWS_AS(analyzer.AnalyzeContinuity(Real(1), Real(0)), std::invalid_argument);
	REQUIRE_THROWS_AS(analyzer.AnalyzeContinuity(Real(0), Real(1), 0), std::invalid_argument);
	REQUIRE_THROWS_AS(analyzer.AnalyzeContinuity(Real(0), Real(1), -1), std::invalid_argument);
	REQUIRE_THROWS_AS(analyzer.AnalyzeContinuity(Real(0), Real(1), 4, Real(0)), std::invalid_argument);
	REQUIRE_THROWS_AS(analyzer.AnalyzeContinuity(Real(0), Real(1), 4, Real(-1)), std::invalid_argument);
	REQUIRE_THROWS_AS(analyzer.AnalyzeContinuity(std::numeric_limits<Real>::quiet_NaN(), Real(1)), std::invalid_argument);
	REQUIRE_THROWS_AS(analyzer.AnalyzeContinuity(Real(0), std::numeric_limits<Real>::infinity()), std::invalid_argument);
	REQUIRE_THROWS_AS(analyzer.AnalyzeContinuity(Real(0), Real(1), 4, std::numeric_limits<Real>::infinity()), std::invalid_argument);
}

TEMPLATE_TEST_CASE("Matrix reports preserve numerical snapshots and optional classifications",
                   "[AnalysisReports][MatrixAnalyzer]", Real, Complex)
{
	const MatrixAnalyzer<TestType> analyzer(Matrix<TestType>{2, 2, {TestType(2), TestType(0), TestType(0), TestType(3)}});
	const auto result = analyzer.Analyze();
	const auto* cached = &analyzer.SVDDecompose();
	std::ostringstream out;
	AnalysisReports::PrintMatrixAnalysis(out, result, "diagonal");
	REQUIRE(out.str().find("Matrix: 2x2") != std::string::npos);
	REQUIRE(out.str().find("Rank: 2, nullity: 0") != std::string::npos);
	REQUIRE(out.str().find("Condition number: 1.5") != std::string::npos);
	REQUIRE(out.str().find("Definiteness: positive definite") != std::string::npos);
	REQUIRE(cached == &analyzer.SVDDecompose());

	MatrixAlg::MatrixAnalysis<TestType> rectangular;
	rectangular.rows = 3;
	rectangular.cols = 2;
	std::ostringstream absent;
	AnalysisReports::PrintMatrixAnalysis(absent, rectangular);
	REQUIRE(absent.str().find("Definiteness:") == std::string::npos);
	const std::array<std::pair<MatrixAlg::Definiteness, const char*>, 6> classifications{{
		{MatrixAlg::Definiteness::PositiveDefinite, "positive definite"},
		{MatrixAlg::Definiteness::PositiveSemidefinite, "positive semidefinite"},
		{MatrixAlg::Definiteness::NegativeDefinite, "negative definite"},
		{MatrixAlg::Definiteness::NegativeSemidefinite, "negative semidefinite"},
		{MatrixAlg::Definiteness::Indefinite, "indefinite"},
		{MatrixAlg::Definiteness::ZeroSemidefinite, "zero semidefinite"}
	}};
	for (const auto& [type, expected] : classifications) {
		rectangular.definiteness = type;
		std::ostringstream classified;
		AnalysisReports::PrintMatrixAnalysis(classified, rectangular);
		REQUIRE(classified.str().find(std::string("Definiteness: ") + expected + '\n') != std::string::npos);
	}
}
