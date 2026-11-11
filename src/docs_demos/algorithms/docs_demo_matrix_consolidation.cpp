#include <mml/algorithms/MatrixAlg/MatrixInversion.h>
#include <mml/algorithms/MatrixSingularity.h>
#include <mml/core/DifferentialGeometry/Metric.h>
#include <mml/base/DifferentialGeometry/FrameTags.h>
#include <mml/tools/ConsolePrinter.h>

#include <iostream>
#include <sstream>
#include <stdexcept>

void Docs_Demo_MatrixConsolidation()
{
	using namespace MML;
	Matrix<Real> matrix(2, 2, {4, 7, 2, 6});
	const auto inverse = MatrixAlg::Inverse(matrix);
	if (!(matrix * inverse).IsIdentity(Precision::LooseTolerance))
		throw std::runtime_error("LU inversion residual regression");
	MatrixAlg::Invert(matrix);
	MatrixNM<Real, 2, 2> fixed{2, 0, 0, 4};
	const auto fixedInverse = MatrixAlg::Inverse(fixed);
	const Metric<2, Cartesian2> metric(fixed);
	if (!IsNearlyEqual(metric.contravariantComponents()(1, 1), fixedInverse(1, 1)))
		throw std::runtime_error("Typed metric inversion regression");
	MatrixPrintFormat format = MatrixPrintFormat::Compact();
	format.width = 0;
	format.zeroThreshold = Precision::GeometryEpsilon;
	std::ostringstream rendered;
	matrix.Print(rendered, format);
	std::cout << "Algorithm-only inverse with shared formatting:\n" << rendered.str();
	Vector<Real> vector{1, 2, 3};
	vector.Print(std::cout, 0, 3);
	std::cout << "\nTable numeric value: "
	          << ColumnFormat("value").format(FormatType::Fixed).precision(3).formatValue(Real(1) / 4) << '\n';
}
