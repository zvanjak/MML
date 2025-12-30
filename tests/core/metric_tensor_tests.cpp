#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Tensor.h"

#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/MetricTensor.h"

#include "base/Geometry3D.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::MetricTensorTests
{
	// TODO REAL(0.9) - HIGH, BIG!!! verify that generating tensor from transf works
	TEST_CASE("Test_Metric_Tensors", "[simple]") {
			TEST_PRECISION_INFO();
		MetricTensorCartesian3D metricCart;
		MetricTensorSpherical metricSpher;
		MetricTensorCylindrical metricCyl;

		CoordTransfSphericalToCartesian coordTransfSpherToCart;

		MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> metricSpherFromCart(coordTransfSpherToCart);
		// using static 
		MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> metricSpherFromCart2(CoordTransfSpherToCart);

		Vector3Cartesian pos(REAL(1.0), REAL(2.0), -REAL(1.0));
		Vector3Spherical posSpher = CoordTransfSpherToCart.transf(pos);
		Vector3Cylindrical posCyl = CoordTransfCylToCart.transf(pos);

		auto cart_metric = metricCart(pos);

		auto spher_metric = metricSpher(posSpher);

		auto cyl_metric = metricCyl(posCyl);
	}
} // namespace MML::Tests::Core::MetricTensorTests