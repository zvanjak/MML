#include <catch2/catch_all.hpp>

#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/DifferentialGeometry/TypedFields.h>
#include <mml/core/DifferentialGeometry/FieldOperations.h>
#endif

#include <type_traits>
#include <utility>

using namespace MML;

namespace MML::Tests::Core::TypedFieldOperationsTests
{
	constexpr Real FieldDerivativeTolerance = Testing::Tol(REAL(1e-8), REAL(2e-6));
	constexpr Real FieldHigherOrderTolerance = Testing::Tol(REAL(1e-8), REAL(1e-5));
	constexpr Real FieldSecondDerivativeTolerance = Testing::Tol(REAL(1e-7), REAL(1e-3));

	template<class, class = void>
	struct HasExteriorDerivative : std::false_type { };

	template<class Field>
	struct HasExteriorDerivative<Field, std::void_t<decltype(exterior_derivative(std::declval<Field>()))>> : std::true_type { };

	template<class, class = void>
	struct HasMetricGradientWithoutMetric : std::false_type { };

	template<class Field>
	struct HasMetricGradientWithoutMetric<Field, std::void_t<decltype(gradient(std::declval<Field>()))>> : std::true_type { };

	class QuadraticScalarField : public IScalarField<2, Cartesian2>
	{
	public:
		Real operator()(const Point<Real, 2, Cartesian2>& point) const override
		{
			return point[0] * point[0] + REAL(3.0) * point[1];
		}
	};

	class PolynomialOneFormField : public IFormField<2, 1, Cartesian2>
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

	class LinearVectorField : public IVectorField<2, Cartesian2>
	{
	public:
		TangentVector<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>& point) const override
		{
			return TangentVector<2, Cartesian2>{ point[0] + point[1], point[0] - point[1] };
		}
	};

	class LinearVectorField3D : public IVectorField<3, Cartesian3>
	{
	public:
		TangentVector<3, Cartesian3> operator()(const Point<Real, 3, Cartesian3>& point) const override
		{
			return TangentVector<3, Cartesian3>{
				point[0] + point[1],
				REAL(2.0) * point[1] - point[2],
				REAL(3.0) * point[2] + point[0]
			};
		}
	};

	class SourceScalarFunction2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& x) const override
		{
			return x[0] + REAL(2.0) * x[1];
		}
	};

	class SourceVectorFunction2D : public IVectorFunction<2>
	{
	public:
		VectorN<Real, 2> operator()(const VectorN<Real, 2>& x) const override
		{
			return VectorN<Real, 2>{ REAL(2.0) * x[0], -x[1] };
		}
	};

	TEST_CASE("TypedFields - Interfaces evaluate typed points", "[TypedFields][DifferentialForms]")
	{
		Point<Real, 2, Cartesian2> point{ REAL(2.0), REAL(5.0) };
		QuadraticScalarField scalar;
		LinearVectorField vector;

		REQUIRE(scalar(point) == REAL(19.0));
		TangentVector<2, Cartesian2> value = vector(point);
		REQUIRE(value[0] == REAL(7.0));
		REQUIRE(value[1] == -REAL(3.0));
	}

	TEST_CASE("TypedFields - Function field adapters preserve values", "[TypedFields][DifferentialForms]")
	{
		Point<Real, 2, Cartesian2> point{ REAL(3.0), REAL(4.0) };
		SourceScalarFunction2D scalarFunc;
		SourceVectorFunction2D vectorFunc;
		ScalarFunctionFieldAdapter<2, Cartesian2> scalarField(scalarFunc);
		VectorFunctionFieldAdapter<2, Cartesian2> vectorField(vectorFunc);

		REQUIRE(scalarField(point) == REAL(11.0));
		TangentVector<2, Cartesian2> v = vectorField(point);
		REQUIRE(v[0] == REAL(6.0));
		REQUIRE(v[1] == -REAL(4.0));
	}

	TEST_CASE("TypedFields - Exterior derivative is a one-form directional derivative", "[TypedFields][DifferentialForms]")
	{
		QuadraticScalarField scalar;
		Point<Real, 2, Cartesian2> point{ REAL(2.0), REAL(5.0) };
		TangentVector<2, Cartesian2> direction{ REAL(0.5), -REAL(2.0) };

		Form1<2, Cartesian2> df = exterior_derivative(scalar)(point);

		REQUIRE(df.Component(0) == Catch::Approx(REAL(4.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(df.Component(1) == Catch::Approx(REAL(3.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(df(direction) == Catch::Approx(-REAL(4.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(DirectionalDerivative(scalar, point, direction) == Catch::Approx(df(direction)).epsilon(FieldDerivativeTolerance));
		REQUIRE(directional_derivative(scalar, point, direction) == Catch::Approx(df(direction)).epsilon(FieldDerivativeTolerance));
	}

	TEST_CASE("TypedFields - Metric gradient is sharp of exterior derivative", "[TypedFields][DifferentialForms]")
	{
		QuadraticScalarField scalar;
		Point<Real, 2, Cartesian2> point{ REAL(2.0), REAL(5.0) };

		TangentVector<2, Cartesian2> euclideanGradient = Gradient(scalar, Metric<2, Cartesian2>::Euclidean())(point);
		REQUIRE(euclideanGradient[0] == Catch::Approx(REAL(4.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(euclideanGradient[1] == Catch::Approx(REAL(3.0)).epsilon(FieldDerivativeTolerance));

		Metric<2, Cartesian2> diagonalMetric = Metric<2, Cartesian2>::Diagonal({ REAL(4.0), REAL(9.0) });
		TangentVector<2, Cartesian2> metricGradient = gradient(scalar, diagonalMetric)(point);
		REQUIRE(metricGradient[0] == Catch::Approx(REAL(1.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(metricGradient[1] == Catch::Approx(REAL(1.0) / REAL(3.0)).epsilon(FieldDerivativeTolerance));
	}

	TEST_CASE("TypedFields - Cartesian gradient and Laplacian wrappers preserve typed outputs", "[TypedFields][DifferentialForms]")
	{
		QuadraticScalarField scalar;
		Point<Real, 2, Cartesian2> point{ REAL(2.0), REAL(5.0) };

		TangentVector<2, Cartesian2> grad = GradientCart(scalar, point);
		TangentVector<2, Cartesian2> gradOrder6 = GradientCart(scalar, point, 6);

		REQUIRE(grad[0] == Catch::Approx(REAL(4.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(grad[1] == Catch::Approx(REAL(3.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(gradOrder6[0] == Catch::Approx(grad[0]).epsilon(FieldHigherOrderTolerance));
		REQUIRE(gradOrder6[1] == Catch::Approx(grad[1]).epsilon(FieldHigherOrderTolerance));
		REQUIRE(LaplacianCart(scalar, point) == Catch::Approx(REAL(2.0)).epsilon(FieldSecondDerivativeTolerance));
	}

	TEST_CASE("TypedFields - Cartesian divergence and curl wrappers preserve typed vector categories", "[TypedFields][DifferentialForms]")
	{
		LinearVectorField3D field;
		Point<Real, 3, Cartesian3> point{ REAL(1.0), REAL(2.0), REAL(3.0) };

		Real div = DivergenceCart(field, point);
		TangentVector<3, Cartesian3> curl = CurlCart(field, point);

		REQUIRE(div == Catch::Approx(REAL(6.0)).epsilon(FieldHigherOrderTolerance));
		REQUIRE(curl[0] == Catch::Approx(REAL(1.0)).epsilon(FieldHigherOrderTolerance));
		REQUIRE(curl[1] == Catch::Approx(-REAL(1.0)).epsilon(FieldHigherOrderTolerance));
		REQUIRE(curl[2] == Catch::Approx(-REAL(1.0)).epsilon(FieldHigherOrderTolerance));
		static_assert(std::is_same<decltype(curl), TangentVector<3, Cartesian3>>::value, "Typed curl preserves tangent-vector category");
	}

	TEST_CASE("TypedFields - Exterior derivative of one-form uses alternating component formula", "[TypedFields][DifferentialForms]")
	{
		PolynomialOneFormField alpha;
		Point<Real, 2, Cartesian2> point{ REAL(2.0), REAL(3.0) };

		Form2<2, Cartesian2> dAlpha = exterior_derivative(alpha)(point);

		REQUIRE(dAlpha.Component(0, 1) == Catch::Approx(REAL(2.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(dAlpha.Component(1, 0) == Catch::Approx(-REAL(2.0)).epsilon(FieldDerivativeTolerance));
		REQUIRE(dAlpha.Component(0, 0) == Catch::Approx(REAL(0.0)).margin(FieldDerivativeTolerance));
		REQUIRE(dAlpha.Component(1, 1) == Catch::Approx(REAL(0.0)).margin(FieldDerivativeTolerance));
	}

	TEST_CASE("TypedFields - Exterior derivative squares to zero for scalar fields", "[TypedFields][DifferentialForms]")
	{
		QuadraticScalarField scalar;
		Point<Real, 2, Cartesian2> point{ REAL(2.0), REAL(5.0) };

		auto df = exterior_derivative(scalar);
		Form2<2, Cartesian2> ddf = exterior_derivative(df)(point);

		REQUIRE(ddf.Component(0, 1) == Catch::Approx(REAL(0.0)).margin(FieldSecondDerivativeTolerance));
		REQUIRE(ddf.Component(1, 0) == Catch::Approx(REAL(0.0)).margin(FieldSecondDerivativeTolerance));
	}

	TEST_CASE("TypedFields - Exterior derivative rejects top-degree form fields", "[TypedFields][DifferentialForms]")
	{
		static_assert(HasExteriorDerivative<IFormField<2, 1, Cartesian2>>::value, "1-forms in 2D have exterior derivatives");
		static_assert(!HasExteriorDerivative<IFormField<2, 2, Cartesian2>>::value, "Top-degree forms do not have exterior derivatives");
		static_assert(HasExteriorDerivative<IScalarField<2, Cartesian2>>::value, "Scalar fields have a metric-free exterior derivative");
		static_assert(!HasMetricGradientWithoutMetric<IScalarField<2, Cartesian2>>::value, "Vector gradient requires an explicit metric");
		SUCCEED("Top-degree exterior derivative is rejected by the overload set.");
	}
}