#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include <mml/core/CoordTransf/CoordTransf2D.h>
#include <mml/core/CoordTransf/CoordTransf3D.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#include <mml/core/Curves.h>
#include <mml/core/Fields/FieldOperations.h>

using namespace MML;

namespace MML::Tests::Core::MathematicalIdentityTests {
	namespace {
		template<int N, class Transform, class Point>
		void RequireInverseJacobians(const Transform& transform, const Point& point) {
			const auto mapped = transform.transf(point);
			const auto back = transform.transfInverse(mapped);
			const auto forward = transform.jacobian(point);
			const auto inverse = transform.inverseJacobian(mapped);
			for (int i = 0; i < N; ++i) {
				CAPTURE(i);
				REQUIRE(IsWithinAbs(back[i], point[i], Precision::LooseTolerance));
				for (int j = 0; j < N; ++j) {
					CAPTURE(j);
					Real left = 0, right = 0;
					for (int k = 0; k < N; ++k) {
						left += inverse(i, k) * forward(k, j);
						right += forward(i, k) * inverse(k, j);
					}
					const Real expected = i == j ? 1 : 0;
					REQUIRE(IsWithinAbs(left, expected, Precision::LooseTolerance));
					REQUIRE(IsWithinAbs(right, expected, Precision::LooseTolerance));
					const Real numerical = Derivation::NDer4Partial(transform.coordTransfFunc(i), j, point);
					REQUIRE(IsWithinAbs(forward(i, j), numerical, Precision::LooseTolerance));
				}
			}
		}
	} // namespace

	TEST_CASE("Mathematical identities - every core coordinate transform has mutually inverse Jacobians",
						"[identities][coordtransf][wewi.2.6]") {
		for (Real angle : {Real(-0.7), Real(0.4), Real(1.2)}) {
			CAPTURE(angle);
			const Vec2Cart p2{Real(1.25), Real(-0.75)};
			RequireInverseJacobians<2>(CoordTransfPolarToCartesian2D{}, Vector2Polar{2, angle});
			RequireInverseJacobians<2>(CoordTransfCart2DRotation(angle), p2);
			RequireInverseJacobians<2>(CoordTransfObliqueToCartesian2D(angle), p2);
			RequireInverseJacobians<2>(CoordTransfCartesianToOblique2D(angle), p2);

			const Vec3Cart p3{Real(1.25), Real(-0.75), Real(2.5)};
			RequireInverseJacobians<3>(CoordTransfSpherToCart, Vec3Sph{2, Real(1.1), angle});
			RequireInverseJacobians<3>(CoordTransfCartToSpher, p3);
			RequireInverseJacobians<3>(CoordTransfCylToCart, Vec3Cyl{2, angle, -1});
			RequireInverseJacobians<3>(CoordTransfCartToCyl, p3);
			RequireInverseJacobians<3>(CoordTransfCart3DRotationXAxis(angle), p3);
			RequireInverseJacobians<3>(CoordTransfCart3DRotationYAxis(angle), p3);
			RequireInverseJacobians<3>(CoordTransfCart3DRotationZAxis(angle), p3);
			RequireInverseJacobians<3>(CoordTransfCart3DRotationEuler(angle, 2 * angle, -angle), p3);
			RequireInverseJacobians<3>(CoordTransfCart3DRotationQuaternion(Vec3Cart{1, 2, -1}, angle), p3);
			RequireInverseJacobians<3>(CoordTransf3DCartOrthogonal({0, 1, 0}, {-1, 0, 0}, {0, 0, 1}), p3);
			RequireInverseJacobians<3>(CoordTransf3DCartGeneral(MatrixNM<Real, 3, 3>{2, 1, 0, 0, 3, 1, 1, 0, 2}), p3);
			RequireInverseJacobians<3>(CoordTransfCartesianToOblique3D({2, 0, 0}, {1, 3, 0}, {0, 1, 2}), p3);
		}
	}

	TEST_CASE("Mathematical identities - curvilinear Euclidean metrics have symmetric connection and zero curvature",
						"[identities][MetricTensor][wewi.2.6]") {
		MetricTensorCartesian3D cartesian;
		MetricTensorSpherical spherical;
		MetricTensorCylindrical cylindrical;
		for (const MetricTensorField<3>* metric :
				 {static_cast<const MetricTensorField<3>*>(&cartesian), static_cast<const MetricTensorField<3>*>(&spherical),
					static_cast<const MetricTensorField<3>*>(&cylindrical)}) {
			for (const VectorN<Real, 3>& point : {VectorN<Real, 3>{2, Real(0.8), Real(-0.7)}, {Real(3.5), Real(1.4), Real(0.5)}}) {
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						for (int k = 0; k < 3; ++k) {
							CAPTURE(point, i, j, k);
							REQUIRE(IsWithinAbs(metric->GetChristoffelSymbolSecondKind(i, j, k, point),
																	metric->GetChristoffelSymbolSecondKind(i, k, j, point), Precision::EqualityTolerance));
							for (int l = 0; l < 3; ++l) {
								CAPTURE(l);
								REQUIRE(IsWithinAbs(metric->GetRiemannCurvatureTensor(i, j, k, l, point), 0, Precision::LooseTolerance));
							}
						}
			}
		}
	}

	TEST_CASE("Mathematical identities - divergence of both gradient conventions equals the analytic Laplacian",
						"[identities][fields][wewi.2.6]") {
		ScalarFunction<3> cartesian([](const VectorN<Real, 3>& p) { return p[0] * p[0] + 2 * p[1] * p[1] + 3 * p[2] * p[2] + p[0] * p[1]; });
		ScalarFunctionFromStdFunc<3> spherical([&](const VectorN<Real, 3>& p) { return cartesian(CoordTransfSpherToCart.transf(Vec3Sph(p))); });
		ScalarFunctionFromStdFunc<3> cylindrical([&](const VectorN<Real, 3>& p) { return cartesian(CoordTransfCylToCart.transf(Vec3Cyl(p))); });
		MetricTensorCartesian3D cartMetric;
		MetricTensorSpherical sphMetric;
		MetricTensorCylindrical cylMetric;
		for (int coordinates = 0; coordinates < 3; ++coordinates) {
			const IScalarFunction<3>& scalar = coordinates == 0		? static_cast<const IScalarFunction<3>&>(cartesian)
																				 : coordinates == 1 ? static_cast<const IScalarFunction<3>&>(spherical)
																														: cylindrical;
			const MetricTensorField<3>& metric = coordinates == 0		? static_cast<const MetricTensorField<3>&>(cartMetric)
																					 : coordinates == 1 ? static_cast<const MetricTensorField<3>&>(sphMetric)
																															: cylMetric;
			VectorFunctionFromStdFunc<3> coordinateGradient(
				[&](const VectorN<Real, 3>& p) { return ScalarFieldOperations::Gradient(scalar, p, metric); });
			VectorFunctionFromStdFunc<3> physicalGradient([&](const VectorN<Real, 3>& p) -> VectorN<Real, 3> {
				if (coordinates == 0)
					return ScalarFieldOperations::GradientCart(scalar, p);
				if (coordinates == 1)
					return ScalarFieldOperations::GradientSpher(scalar, Vec3Sph(p));
				return ScalarFieldOperations::GradientCyl(scalar, Vec3Cyl(p));
			});
			VectorFunctionFromStdFunc<3> convertedPhysical(
				[&](const VectorN<Real, 3>& p) { return VectorFieldOperations::ToPhysical(coordinateGradient(p), p, metric); });
			VectorFunctionFromStdFunc<3> convertedCoordinate(
				[&](const VectorN<Real, 3>& p) { return VectorFieldOperations::ToContravariant(physicalGradient(p), p, metric); });
			for (const VectorN<Real, 3>& point : {VectorN<Real, 3>{2, Real(0.8), Real(-0.7)}, {Real(3.5), Real(1.4), Real(0.5)}}) {
				CAPTURE(coordinates, point);
				REQUIRE(convertedPhysical(point).IsEqualTo(physicalGradient(point), Precision::LooseTolerance));
				REQUIRE(convertedCoordinate(point).IsEqualTo(coordinateGradient(point), Precision::LooseTolerance));
				const Real laplacian = coordinates == 0		? ScalarFieldOperations::LaplacianCart(scalar, point)
															 : coordinates == 1 ? ScalarFieldOperations::LaplacianSpher(scalar, Vec3Sph(point))
																									: ScalarFieldOperations::LaplacianCyl(scalar, Vec3Cyl(point));
				const Real physicalDiv = coordinates == 0		? VectorFieldOperations::DivCart(physicalGradient, point)
																 : coordinates == 1 ? VectorFieldOperations::DivSpher(physicalGradient, point)
																										: VectorFieldOperations::DivCyl(physicalGradient, point);
				const Real convertedDiv = coordinates == 0	 ? VectorFieldOperations::DivCart(convertedPhysical, point)
																	: coordinates == 1 ? VectorFieldOperations::DivSpher(convertedPhysical, point)
																										 : VectorFieldOperations::DivCyl(convertedPhysical, point);
				for (Real actual : {laplacian, physicalDiv, VectorFieldOperations::Divergence(coordinateGradient, point, metric),
														VectorFieldOperations::DivergenceViaDet(coordinateGradient, point, metric), convertedDiv,
														VectorFieldOperations::Divergence(convertedCoordinate, point, metric)})
					REQUIRE(IsWithinAbs(actual, 12, Precision::CoarseTolerance));
			}
		}
	}

	TEST_CASE("Field components - metric scale factors convert basis components without lowering indices", "[fields][components][wewi.8.1]") {
		class DiagonalMetric : public MetricTensorField<2> {
		public:
			Real component(int i, int j, const VectorN<Real, 2>&) const override { return i != j ? 0 : i == 0 ? 4 : 9; }
		} metric;
		const VectorN<Real, 2> point{2, -1};
		const VectorN<Real, 2> coordinate{3, -4};
		const auto physical = VectorFieldOperations::ToPhysical(coordinate, point, metric);
		REQUIRE(physical == VectorN<Real, 2>{6, -12});
		REQUIRE(VectorFieldOperations::ToContravariant(physical, point, metric).IsEqualTo(coordinate));
		REQUIRE(metric.LowerIndex(coordinate, point) == VectorN<Real, 2>{12, -36});
		REQUIRE(IsNearlyEqual(physical.NormL2() * physical.NormL2(), 180));
		REQUIRE(point == VectorN<Real, 2>{2, -1});
		REQUIRE(coordinate == VectorN<Real, 2>{3, -4});

		MetricTensorCartesian3D cartesian;
		const Vec3Cart vector{2, -3, 4};
		REQUIRE(VectorFieldOperations::ToPhysical(vector, vector, cartesian) == vector);
		REQUIRE(VectorFieldOperations::ToContravariant(vector, vector, cartesian) == vector);
	}

	TEST_CASE("Field components - orthogonality checks use relative basis lengths and tolerate roundoff", "[fields][components][wewi.8.1]") {
		class ScaledMetric : public MetricTensorField<2> {
		public:
			Real firstScale = std::sqrt(Precision::DivisionSafetyThreshold);
			Real secondScale = std::sqrt(std::numeric_limits<Real>::max()) / 4;
			Real correlation = Precision::OrthogonalityTolerance / 4;
			Real component(int i, int j, const VectorN<Real, 2>&) const override {
				return i == j ? (i == 0 ? firstScale * firstScale : secondScale * secondScale) : correlation * firstScale * secondScale;
			}
		} metric;
		const VectorN<Real, 2> point{0, 0}, coordinate{2, -3};
		const auto physical = VectorFieldOperations::ToPhysical(coordinate, point, metric);
		REQUIRE(IsWithinRel(physical[0], 2 * metric.firstScale, Precision::EqualityTolerance));
		REQUIRE(IsWithinRel(physical[1], -3 * metric.secondScale, Precision::EqualityTolerance));
		REQUIRE(VectorFieldOperations::ToContravariant(physical, point, metric).IsEqualTo(coordinate));
		metric.correlation = 4 * Precision::OrthogonalityTolerance;
		REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(coordinate, point, metric), DomainError);
		REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(physical, point, metric), DomainError);
		metric.correlation = 0;
		metric.firstScale = Precision::DivisionSafetyThreshold / 2;
		REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(coordinate, point, metric), DomainError);
		REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(physical, point, metric), DomainError);
	}

	TEST_CASE("Field components - invalid metrics and nonfinite conversions fail explicitly", "[fields][components][wewi.8.1]") {
		class SuppliedMetric : public MetricTensorField<2> {
		public:
			MatrixNM<Real, 2, 2> values;
			explicit SuppliedMetric(const MatrixNM<Real, 2, 2>& g)
					: values(g) {}
			Real component(int i, int j, const VectorN<Real, 2>&) const override { return values(i, j); }
		};
		const VectorN<Real, 2> point{1, 2}, vector{3, -4};
		const Real nan = std::numeric_limits<Real>::quiet_NaN();
		const Real infinity = std::numeric_limits<Real>::infinity();
		for (const auto& matrix :
				 {MatrixNM<Real, 2, 2>{1, Real(0.25), Real(0.25), 1}, MatrixNM<Real, 2, 2>{1, 0, Real(0.25), 1}, MatrixNM<Real, 2, 2>{0, 0, 0, 1},
					MatrixNM<Real, 2, 2>{-1, 0, 0, 1}, MatrixNM<Real, 2, 2>{infinity, 0, 0, 1}, MatrixNM<Real, 2, 2>{1, nan, 0, 1}}) {
			SuppliedMetric metric(matrix);
			REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(vector, point, metric), DomainError);
			REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(vector, point, metric), DomainError);
		}
		SuppliedMetric metric(MatrixNM<Real, 2, 2>{4, 0, 0, Real(0.25)});
		for (const auto& invalid : {VectorN<Real, 2>{nan, 1}, VectorN<Real, 2>{1, infinity}}) {
			REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(invalid, point, metric), DomainError);
			REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(invalid, point, metric), DomainError);
			REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(vector, invalid, metric), DomainError);
			REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(vector, invalid, metric), DomainError);
		}
		const Real maximum = std::numeric_limits<Real>::max();
		REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(VectorN<Real, 2>{maximum, 0}, point, metric), DomainError);
		REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(VectorN<Real, 2>{0, maximum}, point, metric), DomainError);
		MetricTensorSpherical spherical;
		for (const Vec3Sph& singular : {Vec3Sph{0, 1, 0}, {2, 0, 0}}) {
			REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(Vec3Sph{1, 2, 3}, singular, spherical), DomainError);
			REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(Vec3Sph{1, 2, 3}, singular, spherical), DomainError);
		}
		MetricTensorMinkowski minkowski;
		REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(VectorN<Real, 4>{1, 2, 3, 4}, VectorN<Real, 4>{0, 0, 0, 0}, minkowski),
											DomainError);
		REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(VectorN<Real, 4>{1, 2, 3, 4}, VectorN<Real, 4>{0, 0, 0, 0}, minkowski),
											DomainError);
		MetricTensorSphericalContravar inverseMetric;
		REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(Vec3Sph{1, 2, 3}, Vec3Sph{2, 1, 0}, inverseMetric), DomainError);
		REQUIRE_THROWS_AS(VectorFieldOperations::ToContravariant(Vec3Sph{1, 2, 3}, Vec3Sph{2, 1, 0}, inverseMetric), DomainError);
	}

	TEST_CASE("Mathematical identities - numerical helix curvature and signed torsion are parameter independent",
						"[identities][curves][wewi.2.6]") {
		class NumericalHelix : public Curves::ICurveCartesian3D {
		public:
			Real pitch;
			explicit NumericalHelix(Real b)
					: pitch(b) {}
			Real getMinT() const override { return -10; }
			Real getMaxT() const override { return 10; }
			VectorN<Real, 3> operator()(Real t) const override { return {2 * std::cos(t), 2 * std::sin(t), pitch * t}; }
		};
		for (Real pitch : {Real(-3), Real(0), Real(3)}) {
			NumericalHelix helix(pitch);
			for (Real t : {Real(-2), Real(0), Real(0.75), Real(4)}) {
				CAPTURE(pitch, t);
				REQUIRE(IsWithinAbs(helix.getCurvature(t), 2 / (4 + pitch * pitch), Precision::LooseTolerance));
				REQUIRE(IsWithinAbs(helix.getTorsion(t), pitch / (4 + pitch * pitch), Precision::CoarseTolerance));
			}
		}
	}


} // namespace MML::Tests::Core::MathematicalIdentityTests
