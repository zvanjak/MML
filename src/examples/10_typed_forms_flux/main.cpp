/******************************************************************************
 * MML Example: Typed Differential Forms - Flux and Surface Normal
 * ============================================================================
 *
 * Demonstrates the typed forms layer in 3D Cartesian geometry:
 *
 *   1. flat(v) turns a velocity tangent vector into a covector
 *   2. hodge_star(flat(v)) is the flux two-form
 *   3. cross(ru, rv) is derived from wedge + Hodge star + sharp
 *
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/base/DifferentialGeometry/Hodge.h>
#include <mml/base/DifferentialGeometry/Metric.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#include <mml/core/DifferentialGeometry/CoordinateMap.h>
#endif

#include <iomanip>
#include <iostream>

using namespace MML;

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
	std::cout << "Typed forms: flux and surface normal\n\n";

	Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Euclidean();
	Orientation orientation = Orientation::Positive;

	TangentVector<3, Cartesian3> velocity{ REAL(2.0), -REAL(1.0), REAL(3.0) };
	TangentVector<3, Cartesian3> ru{ REAL(1.0), REAL(0.0), REAL(0.0) };
	TangentVector<3, Cartesian3> rv{ REAL(0.0), REAL(1.0), REAL(0.0) };

	Form1<3, Cartesian3> velocityFlat = ToForm(Flat(metric, velocity));
	Form2<3, Cartesian3> fluxForm = hodge_star(velocityFlat, metric, orientation);
	TangentVector<3, Cartesian3> normal = cross(ru, rv, metric, orientation);
	Form2<3, Cartesian3> tangentArea = wedge(ToForm(Flat(metric, ru)), ToForm(Flat(metric, rv)));

	Real fluxThroughPatch = fluxForm(ru, rv);
	Real velocityDotNormal = Inner(metric, velocity, normal);

	PrintVector("velocity", velocity);
	PrintVector("ru", ru);
	PrintVector("rv", rv);
	PrintVector("normal = cross(ru,rv)", normal);

	std::cout << "\n";
	std::cout << "  flux_form(ru, rv)        = " << fluxThroughPatch << "\n";
	std::cout << "  inner(velocity, normal)  = " << velocityDotNormal << "\n";
	std::cout << "  area_form(ru, rv)        = " << tangentArea(ru, rv) << "\n";
	std::cout << "\nFor the xy patch, positive orientation gives normal=(0,0,1), so flux=v_z=3.\n";

	CoordTransfCylindricalToCartesian cylToCart;
	CylindricalToCartesian3DMap map = MakeCylindricalToCartesian3DMap(cylToCart);
	Point<Real, 3, Cylindrical3> cylPoint{ REAL(2.0), Constants::PI / REAL(2.0), REAL(5.0) };
	Form3<3, Cartesian3> cartVolume;
	cartVolume.SetAlternatingComponent(REAL(1.0), 0, 1, 2);
	Form3<3, Cylindrical3> cylVolume = pull_back(map, cartVolume, cylPoint);

	std::cout << "\nCoordinate pull-back of Cartesian volume form to cylindrical at r=2:\n";
	std::cout << "  pulled_volume(0,1,2)    = " << cylVolume.Component(0, 1, 2) << "\n";
	std::cout << "  Expected Jacobian factor = r = 2.\n";

	return 0;
}