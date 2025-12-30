#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Tensor.h"
#include "base/BaseUtils.h"
#endif

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
// WithinRel removed - use RealWithinRel instead

namespace MML::Tests::Base::TensorsTests
{
	/*********************************************************************/
	/*****                      Tensors2 tests                       *****/
	/*********************************************************************/

	TEST_CASE("Tensors::Tensor2_init", "[simple]") {
			TEST_PRECISION_INFO();
		Tensor2<3> t1(1, 1);
		t1(0, 0) = 1;
		t1(0, 2) = -1;
		t1(1, 1) = 2;
		t1(2, 2) = 3;

		REQUIRE(t1(0, 0) == 1);
		REQUIRE(t1(0, 2) == -1);
		REQUIRE(t1(1, 1) == 2);
		REQUIRE(t1(2, 2) == 3);

		Tensor2<3> t2(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
													REAL(4.0), REAL(5.0), REAL(6.0),
													REAL(7.0), REAL(8.0), REAL(9.0) });

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
			TEST_PRECISION_INFO();
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
			TEST_PRECISION_INFO();
		Tensor2<3> t1(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
													REAL(4.0), REAL(5.0), REAL(6.0),
													REAL(7.0), REAL(8.0), REAL(9.0) });
		Tensor2<3> t2(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
													REAL(4.0), REAL(5.0), REAL(6.0),
													REAL(7.0), REAL(8.0), REAL(9.0) });

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


	Vec3 v1{ REAL(1.0), REAL(2.0), REAL(3.0) };
	Vec3 v2{ REAL(1.0), REAL(2.0), REAL(3.0) };

		auto b = t2.Contract();
	}

	TEST_CASE("Tensors::Tensor2_Contraction", "[simple]")
	{
			TEST_PRECISION_INFO();
		Tensor2<3> t1(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
													REAL(4.0), REAL(5.0), REAL(6.0),
													REAL(7.0), REAL(8.0), REAL(9.0) });

		auto d = t1.Contract();

		REQUIRE(d == 15);
	}

	TEST_CASE("Tensors::Tensor2_operation()(Vector, Vector)", "[simple]")
	{
		Tensor2<3> t1(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
													REAL(4.0), REAL(5.0), REAL(6.0),
													REAL(7.0), REAL(8.0), REAL(9.0) });

Vec3 v1{ REAL(1.0), REAL(2.0), REAL(3.0) };
	Vec3 v2{ REAL(1.0), REAL(2.0), REAL(3.0) };

		auto d = t1(v1, v2);

		REQUIRE(d == 228);

		// verifying directly
		MatrixNM<Real, 3, 3> mat({ REAL(1.0), REAL(2.0), REAL(3.0),
															REAL(4.0), REAL(5.0), REAL(6.0),
															REAL(7.0), REAL(8.0), REAL(9.0) });
		VectorN<Real, 3> v3({ REAL(1.0), REAL(2.0), REAL(3.0) });
		VectorN<Real, 3> v4({ REAL(1.0), REAL(2.0), REAL(3.0) });

		auto a = v3 * mat;
		auto b = Utils::ScalarProduct<3>(a, v4);

		REQUIRE(b == 228);
	}

	TEST_CASE("Tensors::Tensor2_multiple_contractions", "[operations]")
	{
			TEST_PRECISION_INFO();
		// Test identity-like tensor (Kronecker delta) - requires 1 contravariant and 1 covariant
		Tensor2<3> identity(1, 1);
		identity(0, 0) = REAL(1.0);
		identity(1, 1) = REAL(1.0);
		identity(2, 2) = REAL(1.0);

		auto trace = identity.Contract();
		REQUIRE_THAT(trace , RealWithinRel(REAL(3.0), REAL(1e-5)));

		// Test symmetric mixed tensor
		Tensor2<3> symmetric(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
																 REAL(2.0), REAL(4.0), REAL(5.0),
																 REAL(3.0), REAL(5.0), REAL(6.0) });
		auto sym_trace = symmetric.Contract();
		REQUIRE_THAT(sym_trace , RealWithinRel(REAL(11.0), REAL(1e-5)));  // 1 + 4 + 6

		// Test antisymmetric mixed tensor (should have zero trace)
		Tensor2<3> antisym(1, 1, { REAL(0.0), REAL(1.0), -REAL(2.0),
															-REAL(1.0), REAL(0.0), REAL(3.0),
															REAL(2.0), -REAL(3.0), REAL(0.0) });
		auto antisym_trace = antisym.Contract();
		REQUIRE_THAT(antisym_trace , RealWithinRel(REAL(0.0), REAL(1e-5)));
	}

	TEST_CASE("Tensors::Tensor2_tensor_products", "[operations]")
	{
			TEST_PRECISION_INFO();
		// Test outer product behavior
Vec3 v1{ REAL(1.0), REAL(0.0), REAL(0.0) };
	Vec3 v2{ REAL(0.0), REAL(1.0), REAL(0.0) };

		Tensor2<3> t1(2, 0);
		// Manually construct outer product v1 ⊗ v2
		t1(0, 0) = v1[0] * v2[0];  // 0
		t1(0, 1) = v1[0] * v2[1];  // 1
		t1(0, 2) = v1[0] * v2[2];  // 0
		t1(1, 1) = v1[1] * v2[1];  // 0
		t1(2, 2) = v1[2] * v2[2];  // 0

		// Contract with v1 and v2
		auto result = t1(v1, v2);
		REQUIRE_THAT(result , RealWithinRel(REAL(1.0), REAL(1e-5)));
	}

	TEST_CASE("Tensors::Tensor2_tensor_symmetry", "[properties]")
	{
			TEST_PRECISION_INFO();
		// Test symmetric tensor properties (1 covariant, 1 contravariant)
		Tensor2<3> t_sym(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
														 REAL(2.0), REAL(4.0), REAL(5.0),
														 REAL(3.0), REAL(5.0), REAL(6.0) });

		// Verify symmetry: T_ij = T_ji
		REQUIRE(t_sym(0, 1) == t_sym(1, 0));
		REQUIRE(t_sym(0, 2) == t_sym(2, 0));
		REQUIRE(t_sym(1, 2) == t_sym(2, 1));

		// Test antisymmetric tensor properties (1 covariant, 1 contravariant)
		Tensor2<3> t_antisym(1, 1, { REAL(0.0), REAL(1.0), -REAL(2.0),
																-REAL(1.0), REAL(0.0), REAL(3.0),
																REAL(2.0), -REAL(3.0), REAL(0.0) });

		// Verify antisymmetry: T_ij = -T_ji
		REQUIRE(t_antisym(0, 1) == -t_antisym(1, 0));
		REQUIRE(t_antisym(0, 2) == -t_antisym(2, 0));
		REQUIRE(t_antisym(1, 2) == -t_antisym(2, 1));
		// Diagonal elements must be zero
		REQUIRE(t_antisym(0, 0) == REAL(0.0));
		REQUIRE(t_antisym(1, 1) == REAL(0.0));
		REQUIRE(t_antisym(2, 2) == REAL(0.0));
	}

	/*********************************************************************/
	/*****                      Tensors3 tests                       *****/
	/*********************************************************************/
	TEST_CASE("Tensors::Tensor3_init", "[simple]") {
			TEST_PRECISION_INFO();
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
			TEST_PRECISION_INFO();
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
			TEST_PRECISION_INFO();
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

	/*********************************************************************/
	/*****              Higher-order tensor operations               *****/
	/*********************************************************************/
	TEST_CASE("Tensors::Higher_order_tensor_operations", "[operations]")
	{
			TEST_PRECISION_INFO();
		// Test Tensor3 basic operations
		Tensor3<3> t3(2, 1);
		t3(0, 0, 0) = REAL(1.0);
		t3(1, 1, 1) = REAL(2.0);
		t3(2, 2, 2) = REAL(3.0);

		// Verify diagonal elements
		REQUIRE(t3(0, 0, 0) == REAL(1.0));
		REQUIRE(t3(1, 1, 1) == REAL(2.0));
		REQUIRE(t3(2, 2, 2) == REAL(3.0));

		// Test Tensor4 basic operations
		Tensor4<3> t4(2, 2);
		t4(0, 0, 0, 0) = REAL(1.0);
		t4(1, 1, 1, 1) = REAL(2.0);
		t4(2, 2, 2, 2) = REAL(3.0);

		// Verify diagonal elements
		REQUIRE(t4(0, 0, 0, 0) == REAL(1.0));
		REQUIRE(t4(1, 1, 1, 1) == REAL(2.0));
		REQUIRE(t4(2, 2, 2, 2) == REAL(3.0));

		// Test Tensor5 basic operations
		Tensor5<3> t5(3, 2);
		t5(0, 0, 0, 0, 0) = REAL(1.0);
		t5(1, 1, 1, 1, 1) = REAL(2.0);
		t5(2, 2, 2, 2, 2) = REAL(3.0);

		// Verify diagonal elements
		REQUIRE(t5(0, 0, 0, 0, 0) == REAL(1.0));
		REQUIRE(t5(1, 1, 1, 1, 1) == REAL(2.0));
		REQUIRE(t5(2, 2, 2, 2, 2) == REAL(3.0));
	}

	TEST_CASE("Tensors::Covariant_contravariant_consistency", "[properties]")
	{
			TEST_PRECISION_INFO();
		// Test that covariant and contravariant counts are consistent
		// Constructor signature: Tensor2(nCovar, nContraVar)
		Tensor2<3> t2_1(1, 1);  // 1 covariant, 1 contravariant
		REQUIRE(t2_1.NumContravar() == 1);
		REQUIRE(t2_1.NumCovar() == 1);
		REQUIRE(t2_1.NumContravar() + t2_1.NumCovar() == 2);

		Tensor2<3> t2_2(0, 2);  // 0 covariant, 2 contravariant
		REQUIRE(t2_2.NumContravar() == 2);
		REQUIRE(t2_2.NumCovar() == 0);

		Tensor2<3> t2_3(2, 0);  // 2 covariant, 0 contravariant
		REQUIRE(t2_3.NumContravar() == 0);
		REQUIRE(t2_3.NumCovar() == 2);

		// Test Tensor3: Tensor3(nCovar, nContraVar)
		Tensor3<3> t3(2, 1);  // 2 covariant, 1 contravariant
		REQUIRE(t3.NumContravar() == 1);
		REQUIRE(t3.NumCovar() == 2);
		REQUIRE(t3.NumContravar() + t3.NumCovar() == 3);

		// Test Tensor4: Tensor4(nCovar, nContraVar)
		Tensor4<3> t4(2, 2);  // 2 covariant, 2 contravariant
		REQUIRE(t4.NumContravar() == 2);
		REQUIRE(t4.NumCovar() == 2);
		REQUIRE(t4.NumContravar() + t4.NumCovar() == 4);
	}

}
