///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        tools/analyzers/FunctionAnalysisReports.h                             ///
///  Description: Stream reports for already computed numerical function analysis      ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_ANALYSIS_REPORTS_H
#define MML_FUNCTION_ANALYSIS_REPORTS_H

#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#include <mml/tools/TextIOFormat.h>

#include <iomanip>
#include <ostream>
#include <string_view>
#include <vector>

namespace MML::AnalysisReports {
	inline void PrintPointAnalysis(std::ostream& os, const PointAnalysisResult& result,
	                               std::string_view functionName = {}) {
		TextIO::ScopedTextFormat format(os, 6);
		if (!functionName.empty())
			os << "f(x) = " << functionName << '\n';
		os << "Numerical point analysis at x = " << result.x << '\n'
		   << "  Probe tolerance: " << result.tolerance << '\n'
		   << "  Finite value: " << (result.isDefined ? "yes" : "no") << '\n'
		   << "  Appears continuous: " << (result.isContinuous ? "yes" : "no") << '\n'
		   << "  Numerical derivatives agree: " << (result.isDerivativeDefined ? "yes" : "no") << '\n'
		   << "  Inflection candidate: " << (result.isInflectionPoint ? "yes" : "no") << '\n'
		   << "  Value: " << result.value << '\n'
		   << "  First derivative estimate: " << result.derivative << '\n'
		   << "  Second derivative estimate: " << result.secondDerivative << '\n';
	}

	inline void PrintIntervalAnalysis(std::ostream& os, const IntervalAnalysisResult& result,
	                                  std::string_view functionName = {}) {
		TextIO::ScopedTextFormat format(os, 6);
		if (!functionName.empty())
			os << "f(x) = " << functionName << '\n';
		os << "Sampled interval analysis [" << result.x1 << ", " << result.x2 << ")\n"
		   << "  Sample points: " << result.numPoints << " (right endpoint excluded)\n"
		   << "  Probe tolerance: " << result.tolerance << '\n'
		   << "  Finite at all sampled points: " << (result.isDefinedEverywhere ? "yes" : "no") << '\n'
		   << "  Continuity probes passed on samples: " << (result.isContinuousEverywhere ? "yes" : "no") << '\n';
		if (result.undefinedPoints.empty())
			os << "  Monotonic on samples: " << (result.isMonotonic ? "yes" : "no") << '\n';
		else
			os << "  Monotonic on samples: unavailable (non-finite samples)\n";
		if (result.numPoints > 0 && result.undefinedPoints.size() < static_cast<std::size_t>(result.numPoints))
			os << "  Sampled minimum: " << result.minValue << '\n'
			   << "  Sampled maximum: " << result.maxValue << '\n';
		else
			os << "  Sampled extrema: unavailable (no finite samples)\n";
		if (!result.undefinedPoints.empty()) {
			os << "  Non-finite sample locations:";
			for (Real x : result.undefinedPoints)
				os << ' ' << x;
			os << '\n';
		}
		if (!result.discontinuityPoints.empty()) {
			os << "  Discontinuity candidate locations:";
			for (Real x : result.discontinuityPoints)
				os << ' ' << x;
			os << '\n';
		}
		os << "Sampling and numerical probes are not proofs over the interval.\n";
	}

	inline void PrintDetailedIntervalAnalysis(std::ostream& os, const std::vector<PointAnalysisResult>& results,
	                                          std::string_view functionName = {}) {
		TextIO::ScopedTextFormat format(os, 6);
		if (!functionName.empty())
			os << "f(x) = " << functionName << '\n';
		os << "Sampled point analysis (" << results.size() << " points)\n"
		   << "         x          value             d1             d2   finite  continuity-probes  derivative-agreement  inflection-candidate\n";
		for (const auto& result : results)
			os << std::setw(10) << result.x << ' ' << std::setw(14) << result.value
			   << ' ' << std::setw(14) << result.derivative << ' ' << std::setw(14) << result.secondDerivative
			   << "   " << (result.isDefined ? "yes" : "no")
			   << "   " << (result.isContinuous ? "pass" : "fail")
			   << "   " << (result.isDerivativeDefined ? "yes" : "no")
			   << "   " << (result.isInflectionPoint ? "yes" : "no") << '\n';
	}

	inline void PrintContinuityAnalysis(std::ostream& os, const ContinuityAnalysisResult& result,
	                                    std::string_view functionName = {}) {
		TextIO::ScopedTextFormat format(os, 6);
		if (!functionName.empty())
			os << "f(x) = " << functionName << '\n';
		os << "Sampled continuity search [" << result.x1 << ", " << result.x2 << "]\n"
		   << "  Search subdivisions: " << result.numSubintervals << '\n'
		   << "  Probe tolerance: " << result.tolerance << '\n'
		   << "Candidate classifications and limits are numerical estimates, not proofs.\n";
		if (result.candidates.empty())
			os << "No discontinuity candidates were detected by the sampling/probing procedure.\n";
		for (const auto& candidate : result.candidates) {
			os << "Discontinuity candidate at x = " << candidate.x << "\n  Estimated type: ";
			switch (candidate.type) {
				case DiscontinuityType::JUMP: os << "JUMP (estimated jump: " << candidate.jumpSize << ')'; break;
				case DiscontinuityType::REMOVABLE: os << "REMOVABLE"; break;
				case DiscontinuityType::Infinite: os << "INFINITE"; break;
				case DiscontinuityType::OSCILLATORY: os << "OSCILLATORY"; break;
				case DiscontinuityType::UNKNOWN: os << "UNKNOWN"; break;
			}
			os << "\n  Left limit estimate: " << candidate.leftLimit
			   << "\n  Right limit estimate: " << candidate.rightLimit
			   << "\n  Value at candidate: " << candidate.valueAtPoint << '\n';
		}
	}
}

#endif // MML_FUNCTION_ANALYSIS_REPORTS_H
