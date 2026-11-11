#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/Derivation/ADAdapters.h>
#include <mml/core/DifferentialGeometry/InducedMetric.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/Integration/PathIntegration.h>
#endif
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>

using namespace MML;

namespace {
	void Verify(Real actual, Real expected) {
		if (!IsNearlyEqual(actual, expected, Precision::DefaultTolerance))
			throw std::runtime_error("Automatic differentiation documentation example failed its reference check");
	}

	void CurveExample() {
		auto helix = AD::MakeADCurve(
			[](auto t) {
				using std::sin;
				using std::cos;
				using S = decltype(t);
				return std::array<S, 3>{2 * cos(t), 2 * sin(t), t};
			},
			0, 2 * Constants::PI);

		const Curves::ICurveCartesian3D& curve = helix;
		Real curvature = curve.getCurvature(0.7);
		Real torsion = curve.getTorsion(0.7);
		Real length = PathIntegration::ParametricCurveLength(curve, 0, 1);
		std::cout << "Helix: curvature=" << curvature << ", torsion=" << torsion << ", length=" << length << '\n';
		Verify(curvature, Real(0.4));
		Verify(torsion, Real(0.2));
		Verify(length, std::sqrt(Real(5)));
		Curves::CurveCartesian3D numerical([](Real t) { return VectorN<Real, 3>{2 * std::cos(t), 2 * std::sin(t), t}; });
		std::cout << std::scientific << "Helix curvature absolute error: AD=" << std::abs(curvature - Real(0.4))
							<< ", finite differences=" << std::abs(numerical.getCurvature(0.7) - Real(0.4)) << std::defaultfloat << '\n';
	}

	void SurfaceExample() {
		auto sphere = AD::MakeADSurface(
			[](auto u, auto w) {
				using std::sin;
				using std::cos;
				using S = decltype(u);
				return std::array<S, 3>{2 * sin(u) * cos(w), 2 * sin(u) * sin(w), 2 * cos(u)};
			},
			0, Constants::PI, 0, 2 * Constants::PI);

		Surfaces::ISurfaceCartesian& surface = sphere;
		Real extrinsicK = surface.GaussianCurvature(0.8, 0.4);
		Real intrinsicK = DifferentialGeometry::GaussianCurvatureIntrinsic(surface, 0.8, 0.4);
		DifferentialGeometry::InducedMetric2D induced(surface);
		Real scalarCurvature = induced.GetRicciScalar({0.8, 0.4});
		std::cout << "Sphere: extrinsic K=" << extrinsicK << ", intrinsic K=" << intrinsicK << ", scalar curvature=" << scalarCurvature << '\n';
		Verify(extrinsicK, Real(0.25));
		Verify(intrinsicK, Real(0.25));
		Verify(scalarCurvature, Real(0.5));
		Verify(std::abs(surface.MeanCurvature(0.8, 0.4)), Real(0.5));
		Surfaces::Sphere numerical(2);
		Real numericalIntrinsicK = DifferentialGeometry::GaussianCurvatureIntrinsic(numerical, 0.8, 0.4);
		std::cout << std::scientific << "Sphere intrinsic K absolute error: AD=" << std::abs(intrinsicK - Real(0.25))
							<< ", finite differences=" << std::abs(numericalIntrinsicK - Real(0.25)) << std::defaultfloat << '\n';
	}

	void MetricExample() {
		auto metric = AD::MakeADMetric<2>([](const auto& q) {
			using std::sin;
			using S = std::decay_t<decltype(q[0])>;
			return std::array<std::array<S, 2>, 2>{{{S(1), S(0)}, {S(0), sin(q[0]) * sin(q[0])}}};
		});

		const MetricTensorField<2>& geometry = metric;
		VectorN<Real, 2> point{0.8, 0.4};
		Real gamma = geometry.GetChristoffelSymbolSecondKind(0, 1, 1, point);
		Real riemann = geometry.GetRiemannCurvatureTensor(0, 1, 0, 1, point);
		Real scalar = geometry.GetRicciScalar(point);
		std::cout << "Unit sphere metric: Gamma=" << gamma << ", R^0_101=" << riemann << ", scalar curvature=" << scalar << '\n';
		Verify(gamma, -std::sin(point[0]) * std::cos(point[0]));
		Verify(riemann, std::sin(point[0]) * std::sin(point[0]));
		Verify(scalar, Real(2));
	}

	void FieldExample() {
		auto potential = AD::MakeADScalarField<3>([](const auto& p) { return p[0] * p[0] + p[1] * p[1] + p[2] * p[2]; });
		auto flow = AD::MakeADVectorField<3>([](const auto& p) {
			using S = std::decay_t<decltype(p[0])>;
			return std::array<S, 3>{-p[1], p[0], p[2] * p[2]};
		});

		VectorN<Real, 3> point{1, 2, 3};
		auto gradient = ScalarFieldOperations::GradientCart(potential, point);
		Real laplacian = ScalarFieldOperations::LaplacianCart(potential, point);
		Real divergence = VectorFieldOperations::DivCart(flow, point);
		auto curl = VectorFieldOperations::CurlCart(flow, point);
		auto hessian = potential.Hessian(point);
		auto jacobian = flow.Jacobian(point);
		auto detailed = ScalarFieldOperations::GradientCartDetailed(potential, point);
		std::cout << "Fields: gradient=" << gradient << ", Laplacian=" << laplacian << ", divergence=" << divergence << ", curl=" << curl
							<< '\n';
		Verify(gradient[0], Real(2));
		Verify(gradient[1], Real(4));
		Verify(gradient[2], Real(6));
		Verify(laplacian, Real(6));
		Verify(divergence, Real(6));
		Verify(curl[2], Real(2));
		Verify(hessian(0, 0), Real(2));
		Verify(jacobian(2, 2), Real(6));
		if (!detailed.IsSuccess() || !std::isnan(detailed.error[0]) || detailed.function_evaluations != 3)
			throw std::runtime_error("AD field result diagnostics do not match their documented contract");
	}

	void NumericalReferenceAndFailureExample() {
		auto potential = AD::MakeADScalarField<3>([](const auto& p) { return p[0] * p[0] + p[1] * p[1] + p[2] * p[2]; });
		VectorN<Real, 3> point{1, 2, 3};
		auto numericalGradient = ScalarFieldOperations::GradientCart(potential, point, 4);
		Verify(numericalGradient[0], Real(2));

		auto squareRoot = AD::MakeADScalarField<2>([](const auto& p) {
			using std::sqrt;
			return sqrt(p[0]) + p[1];
		});
		FieldOperationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		auto failed = ScalarFieldOperations::GradientCartDetailed(squareRoot, {0, 1}, config);
		if (failed.status != AlgorithmStatus::InvalidInput || failed.error_message.empty())
			throw std::runtime_error("Undefined AD derivative was not reported explicitly");
		std::cout << "Undefined derivative reported: " << failed.error_message << '\n';
	}
} // namespace

void Docs_Demo_AutomaticDifferentiation() {
	CurveExample();
	SurfaceExample();
	MetricExample();
	FieldExample();
	NumericalReferenceAndFailureExample();
}
