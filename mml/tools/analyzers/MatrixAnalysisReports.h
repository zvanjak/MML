///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        tools/analyzers/MatrixAnalysisReports.h                               ///
///  Description: Stream reports for already computed numerical matrix analysis        ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_ANALYSIS_REPORTS_H
#define MML_MATRIX_ANALYSIS_REPORTS_H

#include <mml/algorithms/MatrixAnalysisTypes.h>
#include <mml/tools/TextIOFormat.h>

#include <ostream>
#include <string_view>

namespace MML::AnalysisReports {
	template<MMLScalar Scalar>
	void PrintMatrixAnalysis(std::ostream& os, const MatrixAlg::MatrixAnalysis<Scalar>& analysis,
	                         std::string_view matrixName = {}) {
		TextIO::ScopedTextFormat format(os, 6);
		if (!matrixName.empty())
			os << matrixName << '\n';
		os << "Matrix: " << analysis.rows << 'x' << analysis.cols << '\n'
		   << "Rank: " << analysis.rank << ", nullity: " << analysis.nullity << '\n'
		   << "Condition number: " << analysis.conditionNumber << '\n';
		if (analysis.definiteness.has_value()) {
			os << "Definiteness: ";
			switch (*analysis.definiteness) {
				case MatrixAlg::Definiteness::PositiveDefinite: os << "positive definite"; break;
				case MatrixAlg::Definiteness::PositiveSemidefinite: os << "positive semidefinite"; break;
				case MatrixAlg::Definiteness::NegativeDefinite: os << "negative definite"; break;
				case MatrixAlg::Definiteness::NegativeSemidefinite: os << "negative semidefinite"; break;
				case MatrixAlg::Definiteness::Indefinite: os << "indefinite"; break;
				case MatrixAlg::Definiteness::ZeroSemidefinite: os << "zero semidefinite"; break;
			}
			os << '\n';
		}
	}
}

#endif // MML_MATRIX_ANALYSIS_REPORTS_H
