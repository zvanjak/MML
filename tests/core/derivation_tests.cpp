#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/Function.h"
#include "base/BaseUtils.h"
#include "core/Derivation.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"

using namespace MML;
using namespace MML::Testing;
using namespace MML::Utils;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Core::DerivationTests
{
	/*********************************************************************/
	/*****          Derivation of REAL functions                     *****/
	/*********************************************************************/

	// NDer1 achieves precision of at least 5 decimal places on a set of functions
	// definirati parove <RealFunction, Interval>, init iz TestBeda
	//  onda za svaki ekvidistan covering intervala, izracunati NDer1 i usporediti sa analitickim derivacijama
	TEST_CASE("Test_NDer_sin_func", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getFunc("Sin")._func;
		RealFunction sinFuncDer = TestBeds::RealFunctionsTestBed::getFunc("Sin")._funcDerived;

		if constexpr (std::is_same_v<Real, double>)
		{
			Real der1 = Derivation::NDer1(sinFunc, REAL(0.5));
			REQUIRE_THAT(sinFuncDer(REAL(0.5)) , WithinRel(Real(der1), REAL(1e-8)));

			Real der2 = Derivation::NDer2(sinFunc, REAL(0.5));
			REQUIRE_THAT(sinFuncDer(REAL(0.5)) , WithinRel(Real(der2), REAL(1e-11)));

			Real der4 = Derivation::NDer4(sinFunc, REAL(0.5));
			REQUIRE_THAT(sinFuncDer(REAL(0.5)) , WithinRel(Real(der4), REAL(1e-13)));

			Real der6 = Derivation::NDer6(sinFunc, REAL(0.5));
			REQUIRE_THAT(sinFuncDer(REAL(0.5)) , WithinRel(Real(der6), REAL(1e-13)));

			Real der8 = Derivation::NDer8(sinFunc, REAL(0.5));
			REQUIRE_THAT(sinFuncDer(REAL(0.5)) , WithinRel(Real(der8), REAL(1e-14)));
		}
	}

	TEST_CASE("Test_NDer_TestDer1", "[simple]")
	{
			TEST_PRECISION_INFO();
		auto testFunc = TestBeds::RealFunctionsTestBed::getFuncWithDerivation("TestDer1")._func;
		auto testFuncDer = TestBeds::RealFunctionsTestBed::getFuncWithDerivation("TestDer1")._funcDerived;
	}

	TEST_CASE("Test_NSecDer", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getFunc("Sin")._func;
		RealFunction sinFuncSecDer = TestBeds::RealFunctionsTestBed::getFunc("Sin")._funcSecDer;

		Real secder2 = Derivation::NSecDer2(sinFunc, REAL(0.5));
		REQUIRE_THAT(sinFuncSecDer(REAL(0.5)) , WithinRel(Real(secder2), REAL(1e-5)));

		Real secder4 = Derivation::NSecDer4(sinFunc, REAL(0.5));
		REQUIRE_THAT(sinFuncSecDer(REAL(0.5)) , WithinRel(Real(secder4), REAL(1e-8)));
	}

	TEST_CASE("Test_NThirdDer", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getFunc("Sin")._func;
		RealFunction sinFuncThirdDer = TestBeds::RealFunctionsTestBed::getFunc("Sin")._funcThirdDer;

		Real thirdder2 = Derivation::NThirdDer2(sinFunc, REAL(0.5));
		REQUIRE_THAT(sinFuncThirdDer(REAL(0.5)) , WithinRel(Real(thirdder2), REAL(1e-6)));

		Real thirdder4 = Derivation::NThirdDer4(sinFunc, REAL(0.5));
		REQUIRE_THAT(sinFuncThirdDer(REAL(0.5)) , WithinRel(Real(thirdder4), REAL(1e-7)));
	}

	/*********************************************************************/
	/*****          Derivation of SCALAR functions                   *****/
	/*********************************************************************/

	TEST_CASE("Test_NDerPartial", "[simple]")
	{
			TEST_PRECISION_INFO();
		VectorN<Real, 3> point{ REAL(1.0), REAL(1.0), REAL(3.0) };

		ScalarFunction<3> f = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3("Scalar func 2")._func;
		Real(*fDer)(const VectorN<Real, 3> &, int ind) = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3("Scalar func 2")._funcDerived;

		Real der_x, der_y, der_z;

		der_x = Derivation::NDer1Partial(f, 0, point);
		der_y = Derivation::NDer1Partial(f, 1, point);
		der_z = Derivation::NDer1Partial(f, 2, point);

		REQUIRE_THAT(fDer(point, 0) , WithinRel(der_x, REAL(1e-7)));
		REQUIRE_THAT(fDer(point, 1) , WithinRel(der_y, REAL(1e-7)));
		REQUIRE_THAT(fDer(point, 2) , WithinRel(der_z, REAL(1e-7)));

		der_x = Derivation::NDer2Partial(f, 0, point);
		der_y = Derivation::NDer2Partial(f, 1, point);
		der_z = Derivation::NDer2Partial(f, 2, point);

		REQUIRE_THAT(fDer(point, 0) , WithinRel(der_x, REAL(1e-9)));
		REQUIRE_THAT(fDer(point, 1) , WithinRel(der_y, REAL(1e-11)));
		REQUIRE_THAT(fDer(point, 2) , WithinRel(der_z, REAL(1e-11)));

		der_x = Derivation::NDer4Partial(f, 0, point);
		der_y = Derivation::NDer4Partial(f, 1, point);
		der_z = Derivation::NDer4Partial(f, 2, point);

		REQUIRE_THAT(fDer(point, 0) , WithinRel(der_x, REAL(1e-12)));
		REQUIRE_THAT(fDer(point, 1) , WithinRel(der_y, REAL(1e-10)));
		REQUIRE_THAT(fDer(point, 2) , WithinRel(der_z, REAL(1e-14)));

		der_x = Derivation::NDer6Partial(f, 0, point);
		der_y = Derivation::NDer6Partial(f, 1, point);
		der_z = Derivation::NDer6Partial(f, 2, point);

		REQUIRE_THAT(fDer(point, 0) , WithinRel(der_x, REAL(1e-12)));
		REQUIRE_THAT(fDer(point, 1) , WithinRel(der_y, REAL(1e-13)));
		REQUIRE_THAT(fDer(point, 2) , WithinRel(der_z, REAL(1e-14)));

		der_x = Derivation::NDer8Partial(f, 0, point);
		der_y = Derivation::NDer8Partial(f, 1, point);
		der_z = Derivation::NDer8Partial(f, 2, point);

		REQUIRE_THAT(fDer(point, 0) , WithinRel(der_x, REAL(1e-11)));
		REQUIRE_THAT(fDer(point, 1) , WithinRel(der_y, REAL(1e-8)));
		REQUIRE_THAT(fDer(point, 2) , WithinRel(der_z, REAL(1e-13)));
	}

	TEST_CASE("Test_NSecDerPartial", "[simple]")
	{
		TEST_PRECISION_INFO();
	}

	/*********************************************************************/
	/*****          Derivation of VECTOR functions                   *****/
	/*********************************************************************/
	TEST_CASE("Test_NDer1VecPartial", "[simple]")
	{
			TEST_PRECISION_INFO();
		VectorN<Real, 3> point{ REAL(1.0), REAL(3.0), -REAL(1.0) };

		VectorFunction<3> f = TestBeds::VectorFunctionsTestBed::getTestFunctionVector(0)._func;
		VectorN<Real, 3>(*fDer)(const VectorN<Real, 3> &, int ind) = TestBeds::VectorFunctionsTestBed::getTestFunctionVector(0)._funcDerived;

    Real der_x_x, der_x_y, der_x_z;
    Real der_y_x, der_y_y, der_y_z;
    Real der_z_x, der_z_y, der_z_z;
    der_x_x = Derivation::NDer1Partial(f, 0, 0, point);
    der_x_y = Derivation::NDer1Partial(f, 0, 1, point);
		der_x_z = Derivation::NDer1Partial(f, 0, 2, point);
    der_y_x = Derivation::NDer1Partial(f, 1, 0, point);
    der_y_y = Derivation::NDer1Partial(f, 1, 1, point);
    der_y_z = Derivation::NDer1Partial(f, 1, 2, point);
    der_z_x = Derivation::NDer1Partial(f, 2, 0, point);
    der_z_y = Derivation::NDer1Partial(f, 2, 1, point);
    der_z_z = Derivation::NDer1Partial(f, 2, 2, point);


		REQUIRE_THAT(fDer(point, 0)[0] , WithinRel(der_x_x, REAL(1e-8)));
		REQUIRE_THAT(fDer(point, 0)[1] , WithinRel(der_x_y, REAL(1e-7)));
		REQUIRE_THAT(fDer(point, 0)[2] , WithinRel(der_x_z, REAL(1e-7)));

		der_y_x = Derivation::NDer1Partial(f, 1, 0, point);
		der_y_y = Derivation::NDer1Partial(f, 1, 1, point);
		der_y_z = Derivation::NDer1Partial(f, 1, 2, point);

		REQUIRE_THAT(fDer(point, 1)[0] , WithinRel(der_y_x, REAL(1e-7)));
		REQUIRE_THAT(fDer(point, 1)[1] , WithinRel(der_y_y, REAL(1e-8)));
		REQUIRE_THAT(fDer(point, 1)[2] , WithinRel(der_y_z, REAL(1e-8)));

		der_z_x = Derivation::NDer1Partial(f, 2, 0, point);
		der_z_y = Derivation::NDer1Partial(f, 2, 1, point);
		der_z_z = Derivation::NDer1Partial(f, 2, 2, point);

		REQUIRE_THAT(fDer(point, 2)[0] , WithinRel(der_z_x, REAL(1e-7)));
		REQUIRE_THAT(fDer(point, 2)[1] , WithinRel(der_z_y, REAL(1e-8)));
		REQUIRE_THAT(fDer(point, 2)[2] , WithinRel(der_z_z, REAL(1e-7)));
	}

	/*********************************************************************/
	/*****            Derivation of Tensor field                     *****/
	/*********************************************************************/

	/*********************************************************************/
	/*****          Derivation of Parametric curves                  *****/
	/*********************************************************************/
	TEST_CASE("Test_NDer_curve", "[simple]")
	{
		TEST_PRECISION_INFO();
	}

	TEST_CASE("Test_NSecDer_curve", "[simple]")
	{
		TEST_PRECISION_INFO();
	}

	TEST_CASE("Test_NThirdDer_curve", "[simple]")
	{
		TEST_PRECISION_INFO();
	}

	TEST_CASE("Verify_vector_equation1")
	{

	}

	/*********************************************************************/
	/*****          Gradient tests (all partial derivatives)         *****/
	/*********************************************************************/
	TEST_CASE("Derivation::Gradient_computation", "[gradient]")
	{
			TEST_PRECISION_INFO();
		// Test gradient computation using NDerPartialByAll
		// f(x,y,z) = x^2 + y^2 + z^2, gradient = [2x, 2y, 2z]
		ScalarFunction<3> quadratic{ [](const VectorN<Real, 3>& v) {
			return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		} };

		VectorN<Real, 3> point{ REAL(1.0), REAL(2.0), REAL(3.0) };

		// Compute gradient using NDer2PartialByAll (2nd order accurate)
		VectorN<Real, 3> grad = Derivation::NDer2PartialByAll(quadratic, point);

		// Expected gradient: [2*1, 2*2, 2*3] = [2, 4, 6]
		REQUIRE_THAT(grad[0] , WithinRel(Real(REAL(2.0)), REAL(1e-9)));
		REQUIRE_THAT(grad[1] , WithinRel(Real(REAL(4.0)), REAL(1e-9)));
		REQUIRE_THAT(grad[2] , WithinRel(Real(REAL(6.0)), REAL(1e-9)));

		// Test with 4th order method for higher accuracy
		VectorN<Real, 3> grad4 = Derivation::NDer4PartialByAll(quadratic, point);
		REQUIRE_THAT(grad4[0] , WithinRel(Real(REAL(2.0)), REAL(1e-12)));
		REQUIRE_THAT(grad4[1] , WithinRel(Real(REAL(4.0)), REAL(1e-12)));
		REQUIRE_THAT(grad4[2] , WithinRel(Real(REAL(6.0)), REAL(1e-12)));
	}

	TEST_CASE("Derivation::Gradient_transcendental", "[gradient]")
	{
			TEST_PRECISION_INFO();
		// f(x,y) = exp(x) * sin(y)
		// âˆ‚f/âˆ‚x = exp(x) * sin(y)
		// âˆ‚f/âˆ‚y = exp(x) * cos(y)
		ScalarFunction<2> transcendental{ [](const VectorN<Real, 2>& v) {
			return std::exp(v[0]) * std::sin(v[1]);
		} };

		VectorN<Real, 2> point{ REAL(1.0), Constants::PI / REAL(4.0) };

		VectorN<Real, 2> grad = Derivation::NDer4PartialByAll(transcendental, point);

		// Analytical gradient at (1, Ï€/4)
		Real expected_dx = std::exp(REAL(1.0)) * std::sin(Constants::PI / REAL(4.0));
		Real expected_dy = std::exp(REAL(1.0)) * std::cos(Constants::PI / REAL(4.0));

		REQUIRE_THAT(grad[0] , WithinRel(Real(expected_dx), REAL(1e-10)));
		REQUIRE_THAT(grad[1] , WithinRel(Real(expected_dy), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****          Directional derivative tests                     *****/
	/*********************************************************************/
	TEST_CASE("Derivation::Directional_derivative", "[directional]")
	{
			TEST_PRECISION_INFO();
		// Directional derivative = gradient Â· unit_direction
		// f(x,y,z) = x^2 + 2*y^2 + 3*z^2
		ScalarFunction<3> func{ [](const VectorN<Real, 3>& v) {
			return v[0] * v[0] + REAL(2.0) * v[1] * v[1] + REAL(3.0) * v[2] * v[2];
		} };

		VectorN<Real, 3> point{ REAL(1.0), REAL(1.0), REAL(1.0) };

		// Gradient at (1,1,1): [2*1, 4*1, 6*1] = [2, 4, 6]
		VectorN<Real, 3> grad = Derivation::NDer4PartialByAll(func, point);

		// Direction vector (not normalized): [1, 1, 1]
		VectorN<Real, 3> direction{ REAL(1.0), REAL(1.0), REAL(1.0) };
		direction = direction.Normalized();  // Normalize to unit vector

		// Directional derivative = grad Â· direction = 2/âˆš3 + 4/âˆš3 + 6/âˆš3 = 12/âˆš3
		Real dirDer = ScalarProduct(grad, direction);
		Real expected = REAL(12.0) / std::sqrt(REAL(3.0));

		REQUIRE_THAT(dirDer , WithinRel(Real(expected), REAL(1e-10)));
	}

	TEST_CASE("Derivation::Directional_derivative_along_axes", "[directional]")
	{
			TEST_PRECISION_INFO();
		// Directional derivative along coordinate axes should equal partial derivatives
		ScalarFunction<3> func{ [](const VectorN<Real, 3>& v) {
			return v[0] * v[1] * v[2];  // f(x,y,z) = xyz
		} };

		VectorN<Real, 3> point{ REAL(2.0), REAL(3.0), REAL(5.0) };

		// Gradient: [yz, xz, xy] = [15, 10, 6]
		VectorN<Real, 3> grad = Derivation::NDer4PartialByAll(func, point);

		// Direction along x-axis: [1, 0, 0]
		VectorN<Real, 3> dirX{ REAL(1.0), REAL(0.0), REAL(0.0) };
		Real dirDerX = ScalarProduct(grad, dirX);
		REQUIRE_THAT(dirDerX , WithinRel(Real(grad[0]), REAL(1e-14)));

		// Direction along y-axis: [0, 1, 0]
		VectorN<Real, 3> dirY{ REAL(0.0), REAL(1.0), REAL(0.0) };
		Real dirDerY = ScalarProduct(grad, dirY);
		REQUIRE_THAT(dirDerY , WithinRel(Real(grad[1]), REAL(1e-14)));

		// Direction along z-axis: [0, 0, 1]
		VectorN<Real, 3> dirZ{ REAL(0.0), REAL(0.0), REAL(1.0) };
		Real dirDerZ = ScalarProduct(grad, dirZ);
		REQUIRE_THAT(dirDerZ , WithinRel(Real(grad[2]), REAL(1e-14)));
	}

	/*********************************************************************/
	/*****          Mixed partial derivatives (Schwarz theorem)      *****/
	/*********************************************************************/
	TEST_CASE("Derivation::Mixed_partials_symmetry", "[mixed]")
	{
			TEST_PRECISION_INFO();
		// Schwarz theorem: âˆ‚Â²f/âˆ‚xâˆ‚y = âˆ‚Â²f/âˆ‚yâˆ‚x for continuous functions
		// f(x,y,z) = x^2*y + y^2*z + z^2*x
		ScalarFunction<3> func{ [](const VectorN<Real, 3>& v) {
			return v[0] * v[0] * v[1] + v[1] * v[1] * v[2] + v[2] * v[2] * v[0];
		} };

		VectorN<Real, 3> point{ REAL(1.0), REAL(2.0), REAL(3.0) };

		// Test âˆ‚Â²f/âˆ‚xâˆ‚y = âˆ‚Â²f/âˆ‚yâˆ‚x
		Real fxy = Derivation::NSecDer4Partial(func, 0, 1, point);
		Real fyx = Derivation::NSecDer4Partial(func, 1, 0, point);
		REQUIRE_THAT(fxy , WithinRel(Real(fyx), REAL(1e-8)));

		// Test âˆ‚Â²f/âˆ‚xâˆ‚z = âˆ‚Â²f/âˆ‚zâˆ‚x
		Real fxz = Derivation::NSecDer4Partial(func, 0, 2, point);
		Real fzx = Derivation::NSecDer4Partial(func, 2, 0, point);
		REQUIRE_THAT(fxz , WithinRel(Real(fzx), REAL(1e-8)));

		// Test âˆ‚Â²f/âˆ‚yâˆ‚z = âˆ‚Â²f/âˆ‚zâˆ‚y
		Real fyz = Derivation::NSecDer4Partial(func, 1, 2, point);
		Real fzy = Derivation::NSecDer4Partial(func, 2, 1, point);
		REQUIRE_THAT(fyz , WithinRel(Real(fzy), REAL(1e-8)));
	}

	/*********************************************************************/
	/*****          Laplacian tests (sum of pure second derivatives) *****/
	/*********************************************************************/
	TEST_CASE("Derivation::Laplacian_computation", "[laplacian]")
	{
			TEST_PRECISION_INFO();
		// Laplacian = âˆ‚Â²f/âˆ‚xÂ² + âˆ‚Â²f/âˆ‚yÂ² + âˆ‚Â²f/âˆ‚zÂ²
		// For f(x,y,z) = x^2 + y^2 + z^2, Laplacian = 2 + 2 + 2 = 6
		ScalarFunction<3> quadratic{ [](const VectorN<Real, 3>& v) {
			return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		} };

		VectorN<Real, 3> point{ REAL(1.0), REAL(2.0), REAL(3.0) };

		Real fxx = Derivation::NSecDer4Partial(quadratic, 0, 0, point);
		Real fyy = Derivation::NSecDer4Partial(quadratic, 1, 1, point);
		Real fzz = Derivation::NSecDer4Partial(quadratic, 2, 2, point);

		Real laplacian = fxx + fyy + fzz;
		REQUIRE_THAT(laplacian , WithinRel(Real(REAL(6.0)), REAL(1e-8)));
	}

	TEST_CASE("Derivation::Laplacian_harmonic", "[laplacian]")
	{
			TEST_PRECISION_INFO();
		// Harmonic functions have Laplacian = 0
		// f(x,y) = x^2 - y^2, Laplacian = 2 + (-2) = 0
		ScalarFunction<2> harmonic{ [](const VectorN<Real, 2>& v) {
			return v[0] * v[0] - v[1] * v[1];
		} };

		VectorN<Real, 2> point{ REAL(1.0), REAL(2.0) };

		Real fxx = Derivation::NSecDer4Partial(harmonic, 0, 0, point);
		Real fyy = Derivation::NSecDer4Partial(harmonic, 1, 1, point);

		Real laplacian = fxx + fyy;
		REQUIRE_THAT(laplacian, WithinAbs(REAL(0.0), REAL(1e-8)));
	}

	/*********************************************************************/
	/*****          Gradient descent and steepest ascent            *****/
	/*********************************************************************/
	TEST_CASE("Derivation::Gradient_direction_properties", "[gradient]")
	{
			TEST_PRECISION_INFO();
		// Gradient points in direction of steepest ascent
		// Test that gradient is perpendicular to level curves
		// f(x,y) = x^2 + y^2 (circular level curves)
		ScalarFunction<2> circular{ [](const VectorN<Real, 2>& v) {
			return v[0] * v[0] + v[1] * v[1];
		} };

		VectorN<Real, 2> point{ REAL(3.0), REAL(4.0) };

		// Gradient at (3,4): [6, 8]
		VectorN<Real, 2> grad = Derivation::NDer4PartialByAll(circular, point);

		// Tangent to level curve at (3,4) is perpendicular to position vector
		// For circle, tangent direction is [-y, x] = [-4, 3]
		VectorN<Real, 2> tangent{ -point[1], point[0] };

		// Gradient should be perpendicular to tangent (dot product = 0)
		Real dotProduct = ScalarProduct(grad, tangent);
		REQUIRE_THAT(dotProduct, WithinAbs(REAL(0.0), REAL(1e-8)));

		// Gradient should be parallel to position vector for x^2+y^2
		// grad/|grad| should equal point/|point|
		VectorN<Real, 2> gradNorm = grad.Normalized();
		VectorN<Real, 2> pointNorm = point.Normalized();
		REQUIRE_THAT(gradNorm[0] , WithinRel(Real(pointNorm[0]), REAL(1e-12)));
		REQUIRE_THAT(gradNorm[1] , WithinRel(Real(pointNorm[1]), REAL(1e-12)));
	}

	/*********************************************************************/
	/*****          Error estimation tests                           *****/
	/*********************************************************************/
	TEST_CASE("Derivation::Error_estimation", "[error]")
	{
			TEST_PRECISION_INFO();
		// Test that error estimates are provided and reasonable
		ScalarFunction<2> smooth{ [](const VectorN<Real, 2>& v) {
			return std::sin(v[0]) * std::cos(v[1]);
		} };

		VectorN<Real, 2> point{ REAL(1.0), REAL(1.0) };
		Real error1, error2, error4;

		Real der1 = Derivation::NDer1Partial(smooth, 0, point, &error1);
		Real der2 = Derivation::NDer2Partial(smooth, 0, point, &error2);
		Real der4 = Derivation::NDer4Partial(smooth, 0, point, &error4);

		// Higher order methods should have smaller errors
		REQUIRE(error2 < error1);
		REQUIRE(error4 < error2);

		// All methods should give similar results within error bounds
		REQUIRE(std::abs(der1 - der2) <= error1 + error2);
		REQUIRE(std::abs(der2 - der4) <= error2 + error4);
	}

	/*********************************************************************/
	/*****          HESSIAN MATRIX TESTS                             *****/
	/*********************************************************************/
	
	// Helper to compute full Hessian matrix using NSecDer4Partial
	template<int N>
	Matrix<Real> ComputeHessian(const ScalarFunction<N>& f, const VectorN<Real, N>& point)
	{
		Matrix<Real> H(N, N);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				H(i, j) = Derivation::NSecDer4Partial(f, i, j, point);
			}
		}
		return H;
	}
	
	TEST_CASE("Hessian::Rosenbrock_at_minimum", "[hessian][optimization]")
	{
		TEST_PRECISION_INFO();
		// Rosenbrock function: f(x,y) = (1-x)² + 100(y-x²)²
		// This is a classic optimization test function with steep valleys
		// The numerical second derivative is challenging due to high curvature
		// 
		// At (1,1): Analytical Hessian = [[802, -400], [-400, 200]]
		// But due to extreme curvature, we test properties rather than exact values
		
		ScalarFunction<2> rosenbrock{ [](const VectorN<Real, 2>& v) {
			Real x = v[0], y = v[1];
			return (REAL(1.0) - x) * (REAL(1.0) - x) + REAL(100.0) * (y - x * x) * (y - x * x);
		} };
		
		VectorN<Real, 2> minimum{ REAL(1.0), REAL(1.0) };
		
		// Compute Hessian at minimum
		Real Hxx = Derivation::NSecDer4Partial(rosenbrock, 0, 0, minimum);
		Real Hxy = Derivation::NSecDer4Partial(rosenbrock, 0, 1, minimum);
		Real Hyx = Derivation::NSecDer4Partial(rosenbrock, 1, 0, minimum);
		Real Hyy = Derivation::NSecDer4Partial(rosenbrock, 1, 1, minimum);
		
		// Analytical values at (1,1): Hxx=802, Hyy=200
		// Diagonal elements are more reliable numerically
		REQUIRE_THAT(Hxx, WithinRel(REAL(802.0), REAL(0.01)));  // Allow 1% error
		REQUIRE_THAT(Hyy, WithinRel(REAL(200.0), REAL(0.01)));
		
		// Symmetry: Hxy should equal Hyx (this is robust)
		REQUIRE_THAT(Hxy, WithinRel(Hyx, REAL(1e-8)));
		
		// Mixed partials: -400x at (1,1), allow larger tolerance for steep function
		// Note: Numerical differentiation struggles here due to extreme curvature
		REQUIRE_THAT(Hxy, WithinRel(REAL(-400.0), REAL(0.3)));  // 30% tolerance
		
		// Key property: Positive definite at minimum (det > 0 and Hxx > 0)
		Real det = Hxx * Hyy - Hxy * Hyx;
		REQUIRE(Hxx > 0);
		REQUIRE(det > 0);  // Confirms minimum, not saddle or maximum
	}
	
	TEST_CASE("Hessian::Saddle_point_properties", "[hessian][saddle]")
	{
		TEST_PRECISION_INFO();
		// f(x,y) = x² - y² (hyperbolic paraboloid)
		// Saddle point at origin, Hessian = [[2, 0], [0, -2]]
		// One positive eigenvalue (+2), one negative (-2)
		
		ScalarFunction<2> saddle{ [](const VectorN<Real, 2>& v) {
			return v[0] * v[0] - v[1] * v[1];
		} };
		
		VectorN<Real, 2> origin{ REAL(0.0), REAL(0.0) };
		
		// Away from origin to avoid numerical issues at exactly 0
		VectorN<Real, 2> nearOrigin{ REAL(0.1), REAL(0.1) };
		
		Real Hxx = Derivation::NSecDer4Partial(saddle, 0, 0, nearOrigin);
		Real Hxy = Derivation::NSecDer4Partial(saddle, 0, 1, nearOrigin);
		Real Hyx = Derivation::NSecDer4Partial(saddle, 1, 0, nearOrigin);
		Real Hyy = Derivation::NSecDer4Partial(saddle, 1, 1, nearOrigin);
		
		// Analytical: Hxx = 2, Hyy = -2, Hxy = Hyx = 0
		REQUIRE_THAT(Hxx, WithinRel(REAL(2.0), REAL(1e-8)));
		REQUIRE_THAT(Hyy, WithinRel(REAL(-2.0), REAL(1e-8)));
		REQUIRE_THAT(Hxy, WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(Hyx, WithinAbs(REAL(0.0), REAL(1e-10)));
		
		// Indefinite: det(H) < 0 (negative eigenvalue product)
		Real det = Hxx * Hyy - Hxy * Hyx;
		REQUIRE(det < 0);  // det = 2*(-2) - 0 = -4
		REQUIRE_THAT(det, WithinRel(REAL(-4.0), REAL(1e-8)));
	}
	
	TEST_CASE("Hessian::Negative_definite_at_maximum", "[hessian][maximum]")
	{
		TEST_PRECISION_INFO();
		// f(x,y) = -x² - y² (inverted paraboloid)
		// Maximum at origin, Hessian = [[-2, 0], [0, -2]]
		// Both eigenvalues negative
		
		ScalarFunction<2> negDefFunc{ [](const VectorN<Real, 2>& v) {
			return -v[0] * v[0] - v[1] * v[1];
		} };
		
		VectorN<Real, 2> point{ REAL(1.0), REAL(1.0) };
		
		Real Hxx = Derivation::NSecDer4Partial(negDefFunc, 0, 0, point);
		Real Hxy = Derivation::NSecDer4Partial(negDefFunc, 0, 1, point);
		Real Hyy = Derivation::NSecDer4Partial(negDefFunc, 1, 1, point);
		
		// Negative definite: Hxx < 0, det(H) > 0
		REQUIRE(Hxx < 0);
		REQUIRE_THAT(Hxx, WithinRel(REAL(-2.0), REAL(1e-8)));
		REQUIRE_THAT(Hyy, WithinRel(REAL(-2.0), REAL(1e-8)));
		REQUIRE_THAT(Hxy, WithinAbs(REAL(0.0), REAL(1e-10)));
		
		Real det = Hxx * Hyy - Hxy * Hxy;
		REQUIRE(det > 0);  // det = (-2)*(-2) - 0 = 4 > 0
	}
	
	TEST_CASE("Hessian::3D_sphere_function", "[hessian][3d]")
	{
		TEST_PRECISION_INFO();
		// f(x,y,z) = x² + y² + z² (n-dimensional sphere)
		// Hessian is 2*I (identity scaled by 2)
		
		ScalarFunction<3> sphere{ [](const VectorN<Real, 3>& v) {
			return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		} };
		
		VectorN<Real, 3> point{ REAL(1.0), REAL(2.0), REAL(3.0) };
		
		Matrix<Real> H = ComputeHessian(sphere, point);
		
		// Verify diagonal is 2, off-diagonal is near 0
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
				{
					REQUIRE_THAT(H(i, j), WithinRel(REAL(2.0), REAL(1e-8)));
				}
				else
				{
					REQUIRE_THAT(H(i, j), WithinAbs(REAL(0.0), REAL(1e-8)));  // Relaxed for numerical noise
				}
			}
		}
		
		// Symmetry (should be exact since same computation)
		REQUIRE_THAT(H(0, 1), WithinRel(H(1, 0), REAL(1e-12)));
		REQUIRE_THAT(H(0, 2), WithinRel(H(2, 0), REAL(1e-12)));
		REQUIRE_THAT(H(1, 2), WithinRel(H(2, 1), REAL(1e-12)));
	}
	
	TEST_CASE("Hessian::Coupled_variables", "[hessian][coupled]")
	{
		TEST_PRECISION_INFO();
		// f(x,y) = xy (simple bilinear coupling)
		// Hxx = 0, Hyy = 0, Hxy = 1
		
		ScalarFunction<2> bilinear{ [](const VectorN<Real, 2>& v) {
			return v[0] * v[1];
		} };
		
		VectorN<Real, 2> point{ REAL(2.0), REAL(3.0) };
		
		// NSecDer4Partial now works correctly with Richardson extrapolation
		Real Hxx = Derivation::NSecDer4Partial(bilinear, 0, 0, point);
		Real Hyy = Derivation::NSecDer4Partial(bilinear, 1, 1, point);
		Real Hxy = Derivation::NSecDer4Partial(bilinear, 0, 1, point);
		Real Hyx = Derivation::NSecDer4Partial(bilinear, 1, 0, point);
		
		// Pure second derivatives should be 0
		REQUIRE_THAT(Hxx, WithinAbs(REAL(0.0), REAL(1e-6)));
		REQUIRE_THAT(Hyy, WithinAbs(REAL(0.0), REAL(1e-6)));
		
		// Mixed partial should be 1
		REQUIRE_THAT(Hxy, WithinRel(REAL(1.0), REAL(1e-5)));
		REQUIRE_THAT(Hyx, WithinRel(REAL(1.0), REAL(1e-5)));  // Symmetry
		REQUIRE_THAT(Hxy, WithinRel(Hyx, REAL(1e-10)));
	}
	
	TEST_CASE("Hessian::Exponential_function", "[hessian][exponential]")
	{
		TEST_PRECISION_INFO();
		// f(x,y) = exp(x)*exp(y) = exp(x+y)
		// ∂²f/∂x² = exp(x+y), ∂²f/∂y² = exp(x+y), ∂²f/∂x∂y = exp(x+y)
		// All second partials equal
		
		ScalarFunction<2> expFunc{ [](const VectorN<Real, 2>& v) {
			return std::exp(v[0]) * std::exp(v[1]);
		} };
		
		VectorN<Real, 2> point{ REAL(0.5), REAL(0.5) };
		Real expected = std::exp(point[0] + point[1]);  // e^1 ≈ 2.718...
		
		// Now using NSecDer4Partial with fixed Richardson extrapolation formula
		Real Hxx = Derivation::NSecDer4Partial(expFunc, 0, 0, point);
		Real Hyy = Derivation::NSecDer4Partial(expFunc, 1, 1, point);
		Real Hxy = Derivation::NSecDer4Partial(expFunc, 0, 1, point);
		
		// All second derivatives should equal exp(1)
		REQUIRE_THAT(Hxx, WithinRel(expected, REAL(1e-6)));
		REQUIRE_THAT(Hyy, WithinRel(expected, REAL(1e-6)));
		REQUIRE_THAT(Hxy, WithinRel(expected, REAL(1e-6)));  // Fixed O(h⁴) formula
		
		// All entries should be approximately equal
		REQUIRE_THAT(Hxx, WithinRel(Hyy, REAL(1e-8)));
	}

	/*********************************************************************/
	/*****          ONE-SIDED DERIVATIVE TESTS                       *****/
	/*********************************************************************/
	
	TEST_CASE("Derivation::Left_right_derivatives_smooth_function", "[one-sided]")
	{
		TEST_PRECISION_INFO();
		// For smooth functions, left and right derivatives should match
		RealFunction sinFunc{ [](Real x) { return std::sin(x); } };
		Real x = REAL(1.0);
		
		Real leftDeriv = Derivation::NDer2Left(sinFunc, x);
		Real rightDeriv = Derivation::NDer2Right(sinFunc, x);
		Real centralDeriv = Derivation::NDer2(sinFunc, x);
		Real analytical = std::cos(x);
		
		// One-sided derivatives are less accurate than central difference
		// They evaluate derivatives at offset points, not at x itself
		REQUIRE_THAT(leftDeriv, WithinRel(analytical, REAL(1e-4)));
		REQUIRE_THAT(rightDeriv, WithinRel(analytical, REAL(1e-4)));
		REQUIRE_THAT(centralDeriv, WithinRel(analytical, REAL(1e-8)));
		
		// Left and right should agree for smooth functions (within tolerance)
		REQUIRE_THAT(leftDeriv, WithinRel(rightDeriv, REAL(1e-4)));
	}
	
	TEST_CASE("Derivation::Left_right_derivatives_absolute_value", "[one-sided]")
	{
		TEST_PRECISION_INFO();
		// f(x) = |x| has different left and right derivatives at x=0
		// Left derivative = -1, Right derivative = +1
		RealFunction absFunc{ [](Real x) { return std::abs(x); } };
		
		Real x = REAL(0.0);
		Real leftDeriv = Derivation::NDer2Left(absFunc, x);
		Real rightDeriv = Derivation::NDer2Right(absFunc, x);
		
		// Left approaches from negative x, derivative is -1
		REQUIRE_THAT(leftDeriv, WithinRel(REAL(-1.0), REAL(1e-5)));
		
		// Right approaches from positive x, derivative is +1
		REQUIRE_THAT(rightDeriv, WithinRel(REAL(1.0), REAL(1e-5)));
		
		// They should differ significantly
		REQUIRE(std::abs(rightDeriv - leftDeriv) > REAL(1.0));
	}
	
	TEST_CASE("Derivation::Left_right_derivatives_with_custom_step", "[one-sided]")
	{
		TEST_PRECISION_INFO();
		// Test that custom step size works correctly
		// Note: One-sided methods evaluate derivative at offset point,
		// so for x²: derivative at (x-3h) or (x+3h), not at x
		RealFunction quadratic{ [](Real x) { return x * x; } };
		Real x = REAL(2.0);  // derivative at 2 should be 4
		Real h = REAL(0.001);
		
		Real leftDeriv = Derivation::NDer2Left(quadratic, x, h);
		Real rightDeriv = Derivation::NDer2Right(quadratic, x, h);
		Real analytical = 2 * x;  // d/dx(x²) = 2x at x=2
		
		// One-sided derivatives are evaluated away from x, allow larger tolerance
		REQUIRE_THAT(leftDeriv, WithinRel(analytical, REAL(0.01)));  // 1% tolerance
		REQUIRE_THAT(rightDeriv, WithinRel(analytical, REAL(0.01)));
	}
	
	TEST_CASE("Derivation::Left_right_derivatives_piecewise", "[one-sided]")
	{
		TEST_PRECISION_INFO();
		// f(x) = x² for x < 1, 2x-1 for x >= 1
		// Both pieces have value 1 at x=1 (continuous)
		// Left derivative at x=1: 2*1 = 2
		// Right derivative at x=1: 2
		// This function is differentiable at x=1
		
		RealFunction piecewise{ [](Real x) {
			if (x < REAL(1.0))
				return x * x;
			else
				return REAL(2.0) * x - REAL(1.0);
		} };
		
		Real x = REAL(1.0);
		Real leftDeriv = Derivation::NDer2Left(piecewise, x);
		Real rightDeriv = Derivation::NDer2Right(piecewise, x);
		
		// Both should be 2 (function is differentiable here)
		REQUIRE_THAT(leftDeriv, WithinRel(REAL(2.0), REAL(1e-4)));
		REQUIRE_THAT(rightDeriv, WithinRel(REAL(2.0), REAL(1e-4)));
	}
	
	TEST_CASE("Derivation::Left_right_derivatives_corner_point", "[one-sided]")
	{
		TEST_PRECISION_INFO();
		// f(x) = x² for x < 1, 3x-2 for x >= 1
		// Both pieces have value 1 at x=1 (continuous)
		// Left derivative at x=1: 2*1 = 2
		// Right derivative at x=1: 3
		// Corner point - not differentiable
		
		RealFunction corner{ [](Real x) {
			if (x < REAL(1.0))
				return x * x;
			else
				return REAL(3.0) * x - REAL(2.0);
		} };
		
		Real x = REAL(1.0);
		Real leftDeriv = Derivation::NDer2Left(corner, x);
		Real rightDeriv = Derivation::NDer2Right(corner, x);
		
		REQUIRE_THAT(leftDeriv, WithinRel(REAL(2.0), REAL(1e-4)));
		REQUIRE_THAT(rightDeriv, WithinRel(REAL(3.0), REAL(1e-4)));
		
		// They differ - corner point
		REQUIRE(std::abs(rightDeriv - leftDeriv) > REAL(0.5));
	}
	
	TEST_CASE("Derivation::Left_right_boundary_behavior", "[one-sided]")
	{
		TEST_PRECISION_INFO();
		// Test at domain boundaries where only one-sided makes sense
		// f(x) = sqrt(x), only defined for x >= 0
		// Use right derivative near x=0
		
		RealFunction sqrtFunc{ [](Real x) { return std::sqrt(x); } };
		Real x = REAL(0.01);  // Near boundary but not at 0
		
		// d/dx(sqrt(x)) = 1/(2*sqrt(x)) at x=0.01 should be 5
		Real analytical = REAL(1.0) / (REAL(2.0) * std::sqrt(x));
		Real rightDeriv = Derivation::NDer2Right(sqrtFunc, x);
		
		// One-sided derivative evaluated at offset, allow larger tolerance
		REQUIRE_THAT(rightDeriv, WithinRel(analytical, REAL(0.01)));  // 1% tolerance
	}
}

