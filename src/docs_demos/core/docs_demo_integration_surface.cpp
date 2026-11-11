#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/Function.h>
#include <mml/base/Geometry/Geometry3DBodies.h>

#include <mml/core/Integration/SurfaceIntegration.h>
#include <mml/core/Fields/Fields.h>
#endif


using namespace MML;

namespace {
	class BoundedFluxPlane : public IParametricSurfaceRect<3> {
		Real scale_, extent_;
	public:
		BoundedFluxPlane(Real scale, Real extent) : scale_(scale), extent_(extent) {}
		Real getMinU() const override { return 0; }
		Real getMaxU() const override { return extent_; }
		Real getMinW() const override { return 0; }
		Real getMaxW() const override { return extent_; }
		VectorN<Real, 3> operator()(Real u, Real v) const override {
			if (u < 0 || u > extent_ || v < 0 || v > extent_)
				throw ArgumentError("BoundedFluxPlane: parameter outside domain");
			return {scale_ * u, scale_ * v, 0};
		}
	};
}

void Docs_Demo_GeometricSurfaceContracts()
{
	auto check = [](Real scale, Real extent, Real fieldZ) {
		const BoundedFluxPlane plane(scale, extent);
		int calls = 0;
		VectorFunctionFromStdFunc<3> field([&](const VectorN<Real, 3>&) {
			++calls;
			return VectorN<Real, 3>{0, 0, fieldZ};
		});
		const auto flux = SurfaceIntegration::SurfaceIntegralDetailed(field, plane, 0, extent, 0, extent, 2, 3);
		if (std::abs(flux.value - 1) > Precision::SurfaceIntegrationTolerance ||
			flux.status != AlgorithmStatus::Completed || flux.converged || flux.error_estimated ||
			HasEstimate(flux.error_estimate) || flux.function_evaluations != calls || calls != 6)
			throw NumericalMethodError("Geometric surface documentation contract failed");
		std::cout << "Bounded midpoint flux: " << flux.value << "; field calls: " << calls << '\n';
	};
	check(Real(1e-8), 1, Real(1e16));
	check(1, Real(1e-7), Real(1e14));
}

void Docs_Demo_Integration_surface()
{
  std::cout << std::endl;
  std::cout << "***********************************************************************" << std::endl;
  std::cout << "****                  SURFACE INTEGRATION                          ****" << std::endl;
  std::cout << "***********************************************************************" << std::endl;

  Cube3D cube(7);
  static VectorFunction<3> field([](const VectorN<Real, 3>& x_cart) { return Fields::InverseRadialPotentialForceFieldCart(1 / (4 * Constants::PI), x_cart); });

  Real integral = SurfaceIntegration::SurfaceIntegral(field, cube);

  std::cout << "Surface integral of the vector field over the cube: " << integral << std::endl;

}
