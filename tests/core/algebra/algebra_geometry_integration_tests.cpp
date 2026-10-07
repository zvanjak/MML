#include <catch2/catch_test_macros.hpp>

#include "../../TestPrecision.h"
#include <mml/base/Algebra_base.h>
#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/base/Tensor/Tensor2.h>

namespace MML::Tests::Core::AlgebraTests
{
	constexpr Real AlgebraGeometryTolerance = Testing::Tol(REAL(1e-12), REAL(1e-6));

	TEST_CASE("SO3 preserves typed volume forms and isotropic tensors", "[Algebra][Integration]")
	{
		const auto rotation = Algebra::SO3::Exp({REAL(0.2), REAL(-0.3), REAL(0.4)});
		TangentVector<3, Cartesian3> e1(rotation.apply({1, 0, 0}));
		TangentVector<3, Cartesian3> e2(rotation.apply({0, 1, 0}));
		TangentVector<3, Cartesian3> e3(rotation.apply({0, 0, 1}));
		Form3<3, Cartesian3> volume;
		volume.SetAlternatingComponent(REAL(1.0), 0, 1, 2);

		REQUIRE(std::abs(volume(e1, e2, e3) - REAL(1.0)) < AlgebraGeometryTolerance);

		MatrixNM<Real, 3, 3> identity = MatrixNM<Real, 3, 3>::Identity();
		const auto transformed = rotation.matrix() * identity * rotation.matrix().transpose();
		REQUIRE(transformed.IsEqualTo(identity, AlgebraGeometryTolerance));

		Tensor2<3> isotropic(2, 0);
		for (int i = 0; i < 3; i++) isotropic(i, i) = REAL(1.0);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				REQUIRE(std::abs(transformed(i, j) - isotropic(i, j)) < AlgebraGeometryTolerance);
	}
}