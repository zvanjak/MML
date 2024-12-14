#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Tensor.h"
#endif

using namespace MML;

namespace MML::Tests::Base::TensorsTests
{
	/*********************************************************************/
	/*****                      Tensors2 tests                       *****/
	/*********************************************************************/

	TEST_CASE("Tensors::Tensor2_init", "[simple]") {
		Tensor2<3> t1(1, 1);
		t1(0, 0) = 1;
		t1(0, 2) = -1;
		t1(1, 1) = 2;
		t1(2, 2) = 3;

		REQUIRE(t1(0, 0) == 1);
		REQUIRE(t1(0, 2) == -1);
		REQUIRE(t1(1, 1) == 2);
		REQUIRE(t1(2, 2) == 3);

		Tensor2<3> t2(1, 1, { 1.0, 2.0, 3.0,
													4.0, 5.0, 6.0,
													7.0, 8.0, 9.0 });

		REQUIRE(t2(0, 0) == 1);
		REQUIRE(t2(0, 1) == 2);
		REQUIRE(t2(0, 2) == 3);
		REQUIRE(t2(1, 0) == 4);
		REQUIRE(t2(1, 1) == 5);
		REQUIRE(t2(1, 2) == 6);
		REQUIRE(t2(2, 0) == 7);
		REQUIRE(t2(2, 1) == 8);
		REQUIRE(t2(2, 2) == 9);
	}

	TEST_CASE("Tensors::Tensor2_covar_contravar", "[simple]")
	{
		// check exceptions
		REQUIRE_THROWS_AS(Tensor2<3>(0, 3), TensorCovarContravarNumError);

		REQUIRE_THROWS_AS(Tensor2<3>(-1, 3), TensorCovarContravarNumError);

		// check NumContravar() & NumCovar
		Tensor2<3> t1(1, 1);

		REQUIRE(t1.NumContravar() == 1);
		REQUIRE(t1.NumCovar() == 1);

		// check _isContravar
	}

	TEST_CASE("Tensors::Tensor2_arithmetic_operations", "[simple]")
	{
		Tensor2<3> t1(1, 1, { 1.0, 2.0, 3.0,
													4.0, 5.0, 6.0,
													7.0, 8.0, 9.0 });
		Tensor2<3> t2(1, 1, { 1.0, 2.0, 3.0,
													4.0, 5.0, 6.0,
													7.0, 8.0, 9.0 });

		auto t3 = t1 + t2;

		REQUIRE(t3(0, 0) == 2);
		REQUIRE(t3(0, 1) == 4);
		REQUIRE(t3(0, 2) == 6);
		REQUIRE(t3(1, 0) == 8);
		REQUIRE(t3(1, 1) == 10);
		REQUIRE(t3(1, 2) == 12);
		REQUIRE(t3(2, 0) == 14);
		REQUIRE(t3(2, 1) == 16);
		REQUIRE(t3(2, 2) == 18);


		Vec3D v1{ 1, 2, 3 };
		Vec3D v2{ 1, 2, 3 };

		auto b = t2.Contract();
	}

	TEST_CASE("Tensors::Tensor2_Contraction", "[simple]")
	{
		Tensor2<3> t1(1, 1, { 1.0, 2.0, 3.0,
													4.0, 5.0, 6.0,
													7.0, 8.0, 9.0 });

		auto d = t1.Contract();

		REQUIRE(d == 15);
	}

	TEST_CASE("Tensors::Tensor2_operation()(Vector, Vector)", "[simple]")
	{
		Tensor2<3> t1(1, 1, { 1.0, 2.0, 3.0,
													4.0, 5.0, 6.0,
													7.0, 8.0, 9.0 });

		Vec3D v1{ 1, 2, 3 };
		Vec3D v2{ 1, 2, 3 };

		auto d = t1(v1, v2);

		REQUIRE(d == 228);

		// verifying directly
		MatrixNM<Real, 3, 3> mat({ 1.0, 2.0, 3.0,
															4.0, 5.0, 6.0,
															7.0, 8.0, 9.0 });
		VectorN<Real, 3> v3({ 1.0, 2.0, 3.0 });
		VectorN<Real, 3> v4({ 1.0, 2.0, 3.0 });

		auto a = v3 * mat;
		auto b = a.ScalarProductCartesian(v4);

		REQUIRE(b == 228);
	}

	/*********************************************************************/
	/*****                      Tensors3 tests                       *****/
	/*********************************************************************/
	TEST_CASE("Tensors::Tensor3_init", "[simple]") {
		Tensor3<3> t1(2, 1);
		t1(0, 0, 0) = 1;
		t1(0, 2, 1) = -1;
		t1(1, 1, 2) = 2;
		t1(2, 2, 0) = 3;

		REQUIRE(t1(0, 0, 0) == 1);
		REQUIRE(t1(0, 2, 1) == -1);
		REQUIRE(t1(1, 1, 2) == 2);
		REQUIRE(t1(2, 2, 0) == 3);
	}

	/*********************************************************************/
	/*****                      Tensors4 tests                       *****/
	/*********************************************************************/
	TEST_CASE("Tensors::Tensor4_init", "[simple]") {
		Tensor4<3> t1(2, 2);
		t1(0, 0, 0, 0) = 1;
		t1(0, 2, 1, 2) = -1;
		t1(1, 1, 2, 0) = 2;
		t1(2, 2, 0, 1) = 3;

		REQUIRE(t1(0, 0, 0, 0) == 1);
		REQUIRE(t1(0, 2, 1, 2) == -1);
		REQUIRE(t1(1, 1, 2, 0) == 2);
		REQUIRE(t1(2, 2, 0, 1) == 3);
	}

	/*********************************************************************/
	/*****                      Tensors5 tests                       *****/
	/*********************************************************************/
	TEST_CASE("Tensors::Tensor5_init", "[simple]") {
		Tensor5<3> t1(2, 3);
		t1(0, 0, 0, 0, 0) = 1;
		t1(0, 2, 1, 2, 1) = -1;
		t1(1, 1, 2, 0, 2) = 2;
		t1(2, 2, 0, 1, 0) = 3;

		REQUIRE(t1(0, 0, 0, 0, 0) == 1);
		REQUIRE(t1(0, 2, 1, 2, 1) == -1);
		REQUIRE(t1(1, 1, 2, 0, 2) == 2);
		REQUIRE(t1(2, 2, 0, 1, 0) == 3);
	}

}