#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/Function.h>
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#include <mml/tools/analyzers/FunctionAnalysisReports.h>
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#include <mml/tools/analyzers/MatrixAnalysisReports.h>
#endif

#include <sstream>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                    BASIC FUNCTION ANALYSIS                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_FunctionAnalyzer_Basic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Basic Function Analysis\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nRealFunctionAnalyzer provides comprehensive function analysis.\n";
	
	// Simple polynomial
	RealFunctionFromStdFunc f([](Real x) { return x * x * x - 3.0 * x + 1.0; });
	RealFunctionAnalyzer analyzer(f);
	
	std::cout << "\n--- f(x) = x³ - 3x + 1 ---\n";
	
	// Point analysis
	std::cout << "\nPoint analysis at x = 0:\n";
	AnalysisReports::PrintPointAnalysis(std::cout, analyzer.AnalyzePoint(0.0));
	
	std::cout << "\nPoint analysis at x = 1:\n";
	AnalysisReports::PrintPointAnalysis(std::cout, analyzer.AnalyzePoint(1.0));
}

void Docs_Demo_FunctionAnalyzer_Interval()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Interval Analysis\n";
	std::cout << "==========================================================================\n";
	
	// Quadratic function
	RealFunctionFromStdFunc f([](Real x) { return (x - 1.0) * (x - 1.0) - 4.0; });
	RealFunctionAnalyzer analyzer(f);
	
	std::cout << "\n--- f(x) = (x-1)² - 4 ---\n";
	AnalysisReports::PrintIntervalAnalysis(std::cout, analyzer.AnalyzeInterval(-2.0, 4.0, 20));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    FINDING ZEROS (ROOTS)                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_FunctionAnalyzer_Zeros()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Finding Function Zeros\n";
	std::cout << "==========================================================================\n";
	
	// Function with known roots
	RealFunctionFromStdFunc f([](Real x) { 
		return (x - 1.0) * (x + 2.0) * (x - 3.0);  // roots at 1, -2, 3
	});
	RealFunctionAnalyzer analyzer(f);
	
	std::cout << "\n--- f(x) = (x-1)(x+2)(x-3)  (roots at -2, 1, 3) ---\n";
	
	// Find zeros in interval
	auto zeros = analyzer.GetRoots(-5.0, 5.0, 1e-8);
	
	std::cout << "\nFound zeros:\n";
	for (const auto& z : zeros)
		std::cout << "  x = " << z << ", f(x) = " << f(z) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    FINDING EXTREMA                                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_FunctionAnalyzer_Extrema()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Finding Extrema (Min/Max)\n";
	std::cout << "==========================================================================\n";
	
	// Function with local extrema
	RealFunctionFromStdFunc f([](Real x) { 
		return std::sin(x) + 0.2 * x;  // oscillating with trend
	});
	RealFunctionAnalyzer analyzer(f);
	
	std::cout << "\n--- f(x) = sin(x) + 0.2x ---\n";
	
	// Find critical points (classified)
	auto critPoints = analyzer.GetLocalOptimumsClassified(-10.0, 10.0);
	
	std::cout << "\nCritical points:\n";
	for (const auto& cp : critPoints) {
		std::string type;
		switch (cp.type) {
			case CriticalPointType::LOCAL_MINIMUM: type = "minimum"; break;
			case CriticalPointType::LOCAL_MAXIMUM: type = "maximum"; break;
			case CriticalPointType::SADDLE_POINT: type = "saddle"; break;
		}
		std::cout << "  x = " << cp.x << ", f(x) = " << cp.value << " (" << type << ")\n";
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    INFLECTION POINTS                                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_FunctionAnalyzer_Inflection()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Finding Inflection Points\n";
	std::cout << "==========================================================================\n";
	
	// Cubic function (inflection at x=0)
	RealFunctionFromStdFunc f([](Real x) { return x * x * x; });
	RealFunctionAnalyzer analyzer(f);
	
	std::cout << "\n--- f(x) = x³  (inflection at x=0) ---\n";
	
	// Find inflection points
	auto inflections = analyzer.GetInflectionPoints(-3.0, 3.0);
	
	std::cout << "\nInflection points:\n";
	for (int i = 0; i < inflections.size(); i++)
		std::cout << "  x = " << inflections[i] << ", f(x) = " << f(inflections[i]) << std::endl;
	
	// Verify: second derivative changes sign
	std::cout << "\nVerification (f''(x) changes sign):\n";
	std::cout << "  f''(-0.1) ≈ " << 6.0 * (-0.1) << " (< 0)\n";
	std::cout << "  f''(+0.1) ≈ " << 6.0 * (0.1) << " (> 0)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    DISCONTINUITY DETECTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_FunctionAnalyzer_Discontinuities()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Detecting Discontinuities\n";
	std::cout << "==========================================================================\n";
	
	// Function with jump discontinuity (floor function)
	RealFunctionFromStdFunc f([](Real x) { return std::floor(x); });
	RealFunctionAnalyzer analyzer(f);
	
	std::cout << "\n--- f(x) = floor(x)  (jump discontinuities at integers) ---\n";
	
	// Find discontinuities
	const auto continuity = analyzer.AnalyzeContinuity(-2.5, 2.5, 100);
	const auto& discs = continuity.candidates;
	
	std::cout << "\nDiscontinuities found:\n";
	for (const auto& d : discs) {
		std::string type;
		switch (d.type) {
			case DiscontinuityType::JUMP: type = "jump"; break;
			case DiscontinuityType::REMOVABLE: type = "removable"; break;
			case DiscontinuityType::INFINITE: type = "infinite"; break;
			case DiscontinuityType::OSCILLATORY: type = "oscillatory"; break;
			case DiscontinuityType::UNKNOWN: type = "unknown"; break;
		}
		std::cout << "  x = " << d.x << " (" << type << "), "
		          << "left = " << d.leftLimit << ", right = " << d.rightLimit << std::endl;
	}
	AnalysisReports::PrintContinuityAnalysis(std::cout, continuity, "floor(x)");
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    MONOTONICITY                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_FunctionAnalyzer_Monotonicity()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Monotonicity Analysis\n";
	std::cout << "==========================================================================\n";
	
	// Increasing function
	RealFunctionFromStdFunc f1([](Real x) { return std::exp(x); });
	RealFunctionAnalyzer analyzer1(f1);
	
	// Oscillating function
	RealFunctionFromStdFunc f2([](Real x) { return std::sin(x); });
	RealFunctionAnalyzer analyzer2(f2);
	
	std::cout << "\n--- f(x) = exp(x) on [0, 5] ---\n";
	std::cout << "Monotonic on samples: " << (analyzer1.IsMonotonicOnSamples(0.0, 5.0, 100) ? "yes" : "no") << std::endl;
	
	std::cout << "\n--- f(x) = sin(x) on [0, 5] ---\n";
	std::cout << "Monotonic on samples: " << (analyzer2.IsMonotonicOnSamples(0.0, 5.0, 100) ? "yes" : "no") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    DETAILED ANALYSIS                                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_FunctionAnalyzer_Detailed()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Detailed Point-by-Point Analysis\n";
	std::cout << "==========================================================================\n";
	
	// Function with interesting features
	RealFunctionFromStdFunc f([](Real x) { 
		return x * x * x - 2.0 * x * x - x + 2.0;  // = (x-1)(x+1)(x-2)
	});
	RealFunctionAnalyzer analyzer(f);
	
	std::cout << "\n--- f(x) = x³ - 2x² - x + 2 = (x-1)(x+1)(x-2) ---\n";
	AnalysisReports::PrintDetailedIntervalAnalysis(std::cout, analyzer.AnalyzePoints(-2.0, 3.0, 10));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_AnalyzerReports()
{
	RealFunction function(+[](Real x) { return x * x * x - 3 * x + 1; });
	const RealFunctionAnalyzer analyzer(function);
	const auto point = analyzer.AnalyzePoint(Real(1));
	const auto interval = analyzer.AnalyzeInterval(Real(-2), Real(2), 100);
	const auto details = analyzer.AnalyzePoints(Real(-2), Real(2), 20);
	const auto continuity = analyzer.AnalyzeContinuity(Real(-2), Real(2), 100);
	AnalysisReports::PrintPointAnalysis(std::cout, point, "x^3 - 3x + 1");
	AnalysisReports::PrintIntervalAnalysis(std::cout, interval);
	AnalysisReports::PrintDetailedIntervalAnalysis(std::cout, details);
	AnalysisReports::PrintContinuityAnalysis(std::cout, continuity);
	std::ostringstream savedReport;
	AnalysisReports::PrintIntervalAnalysis(savedReport, interval);
	std::cout << "Reused snapshot report size: " << savedReport.str().size() << '\n';
	const auto sampled = analyzer.AnalyzeInterval(Real(-3), Real(3), 100, Precision::LooseTolerance);
	AnalysisReports::PrintIntervalAnalysis(std::cout, sampled, "x^3 - 3x + 1");
	RealFunctionFromStdFunc step([](Real x) -> Real {
		if (x < 0) return 0;
		if (x > 0) return 1;
		return Real(1) / 2;
	});
	AnalysisReports::PrintIntervalAnalysis(
		std::cout, RealFunctionAnalyzer(step).AnalyzeInterval(Real(-2), Real(2), 100, Precision::LooseTolerance));

	RealFunction absoluteValue(+[](Real x) { return std::abs(x); });
	AnalysisReports::PrintPointAnalysis(std::cout, RealFunctionAnalyzer(absoluteValue).AnalyzePoint(Real(0)), "|x|");
	RealFunction lorentzian(+[](Real x) { return Real(1) / (Real(1) + x * x); });
	AnalysisReports::PrintIntervalAnalysis(
		std::cout, RealFunctionAnalyzer(lorentzian).AnalyzeInterval(Real(-5), Real(5), 1000), "1/(1+x^2)");

	Matrix<Real> matrix(2, 2, {2, 0, 0, 3});
	MatrixAnalyzer<Real> matrixAnalyzer(matrix);
	const auto matrixAnalysis = matrixAnalyzer.Analyze();
	AnalysisReports::PrintMatrixAnalysis(std::cout, matrixAnalysis, "diagonal");
}

void Docs_Demo_FunctionAnalyzer()
{
	std::cout << "\n##########################################################################\n";
	std::cout << "#                   FUNCTION ANALYZER DEMOS                              #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_FunctionAnalyzer_Basic();
	Docs_Demo_FunctionAnalyzer_Interval();
	Docs_Demo_FunctionAnalyzer_Zeros();
	Docs_Demo_FunctionAnalyzer_Extrema();
	Docs_Demo_FunctionAnalyzer_Inflection();
	Docs_Demo_FunctionAnalyzer_Discontinuities();
	Docs_Demo_FunctionAnalyzer_Monotonicity();
	Docs_Demo_FunctionAnalyzer_Detailed();
	
	std::cout << "\n=== All Function Analyzer Demos Complete ===\n";
}
