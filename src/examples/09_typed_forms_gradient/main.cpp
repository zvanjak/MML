/******************************************************************************
 * MML Example: Typed Differential Forms - Gradient and Directional Derivative
 * ============================================================================
 *
 * Demonstrates the typed forms layer:
 *
 *   1. df is a one-form produced by exterior_derivative(f)
 *   2. df(v) is the directional derivative along a tangent vector
 *   3. grad(f) is sharp(df), so it requires a metric
 *
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/DifferentialGeometry/Metric.h>
#include <mml/base/DifferentialGeometry/TypedFields.h>
#include <mml/core/DifferentialGeometry/FieldOperations.h>
#endif

#include <iomanip>
#include <iostream>

using namespace MML;

class TemperatureField : public IScalarField<3, Cartesian3>
{
public:
	Real operator()(const Point<Real, 3, Cartesian3>& point) const override
	{
		return point[0] * point[0] + REAL(2.0) * point[1] - point[2];
	}
};

class PlanarOneFormField : public IFormField<2, 1, Cartesian2>
{
public:
	Form1<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>& point) const override
	{
		Form1<2, Cartesian2> alpha;
		alpha.Component(0) = point[0] * point[1];
		alpha.Component(1) = point[0] * point[0] + point[1];
		return alpha;
	}
};

template<class VectorLike>
void PrintVector(const char* label, const VectorLike& vector)
{
	std::cout << "  " << std::left << std::setw(24) << label << " = (";
	for (int i = 0; i < vector.size(); i++) {
		if (i > 0)
			std::cout << ", ";
		std::cout << std::right << std::setw(9) << vector[i];
	}
	std::cout << ")\n";
}

int main()
{
	std::cout << std::fixed << std::setprecision(6);
	std::cout << "Typed forms: gradient and directional derivative\n\n";

	TemperatureField temperature;
	Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Euclidean();
	Point<Real, 3, Cartesian3> point{ REAL(1.0), REAL(2.0), REAL(3.0) };
	TangentVector<3, Cartesian3> direction{ REAL(0.5), -REAL(1.0), REAL(2.0) };

	Form1<3, Cartesian3> df = exterior_derivative(temperature)(point);
	Covector<3, Cartesian3> dfCovector = ToCovector(df);
	TangentVector<3, Cartesian3> grad = gradient(temperature, metric)(point);
	Real directionalFromForm = df(direction);
	Real directionalFromHelper = DirectionalDerivative(temperature, point, direction);

	PrintVector("point", point.coordinates());
	PrintVector("direction", direction);
	PrintVector("df components", dfCovector);
	PrintVector("grad = sharp(df)", grad);

	std::cout << "\n";
	std::cout << "  df(direction)            = " << directionalFromForm << "\n";
	std::cout << "  DirectionalDerivative()  = " << directionalFromHelper << "\n";
	std::cout << "\nExpected for f=x^2+2y-z at (1,2,3): df=(2,2,-1), grad=(2,2,-1).\n";

	PlanarOneFormField alpha;
	Point<Real, 2, Cartesian2> planarPoint{ REAL(2.0), REAL(3.0) };
	Form2<2, Cartesian2> dAlpha = exterior_derivative(alpha)(planarPoint);

	std::cout << "\nExterior derivative of alpha = xy dx + (x^2+y) dy at (2,3):\n";
	std::cout << "  d_alpha(0,1)          = " << dAlpha.Component(0, 1) << "\n";
	std::cout << "  d_alpha(1,0)          = " << dAlpha.Component(1, 0) << "\n";
	std::cout << "  Expected: d_alpha = 2 dx^dy with antisymmetric components (+2,-2).\n";

	return 0;
}