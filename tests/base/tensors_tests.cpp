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

	/*********************************************************************/
	/*****                  Contraction tests                        *****/
	/*********************************************************************/

	TEST_CASE("Tensors::Tensor3_Contract_all_index_pairs", "[contraction]")
	{
		TEST_PRECISION_INFO();

		// Create a simple Tensor3<3> with known values for easy verification
		// T[i][j][k] = i + 3*j + 9*k (unique values for each element)
		Tensor3<3> t3(1, 2);  // 1 covariant, 2 contravariant
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					t3(i, j, k) = Real(i + 3*j + 9*k);

		// Contract indices 0 and 1: result[k] = sum_i T[i][i][k]
		auto r01 = t3.Contract(0, 1);
		// sum_i (i + 3*i + 9*k) = sum_i (4*i + 9*k) = (0+4+8) + 3*9*k = 12 + 27*k
		REQUIRE_THAT(r01[0], RealWithinRel(REAL(12.0), REAL(1e-10)));   // k=0: 12
		REQUIRE_THAT(r01[1], RealWithinRel(REAL(39.0), REAL(1e-10)));   // k=1: 12+27=39
		REQUIRE_THAT(r01[2], RealWithinRel(REAL(66.0), REAL(1e-10)));   // k=2: 12+54=66

		// Contract indices 0 and 2: result[j] = sum_i T[i][j][i]
		auto r02 = t3.Contract(0, 2);
		// sum_i (i + 3*j + 9*i) = sum_i (10*i + 3*j) = (0+10+20) + 3*3*j = 30 + 9*j
		REQUIRE_THAT(r02[0], RealWithinRel(REAL(30.0), REAL(1e-10)));   // j=0: 30
		REQUIRE_THAT(r02[1], RealWithinRel(REAL(39.0), REAL(1e-10)));   // j=1: 30+9=39
		REQUIRE_THAT(r02[2], RealWithinRel(REAL(48.0), REAL(1e-10)));   // j=2: 30+18=48

		// Contract indices 1 and 2: result[i] = sum_j T[i][j][j]
		auto r12 = t3.Contract(1, 2);
		// sum_j (i + 3*j + 9*j) = sum_j (i + 12*j) = 3*i + (0+12+24) = 3*i + 36
		REQUIRE_THAT(r12[0], RealWithinRel(REAL(36.0), REAL(1e-10)));   // i=0: 36
		REQUIRE_THAT(r12[1], RealWithinRel(REAL(39.0), REAL(1e-10)));   // i=1: 36+3=39
		REQUIRE_THAT(r12[2], RealWithinRel(REAL(42.0), REAL(1e-10)));   // i=2: 36+6=42

		// Test reversed index order (should give same result)
		auto r10 = t3.Contract(1, 0);
		REQUIRE_THAT(r10[0], RealWithinRel(r01[0], REAL(1e-10)));
		REQUIRE_THAT(r10[1], RealWithinRel(r01[1], REAL(1e-10)));
		REQUIRE_THAT(r10[2], RealWithinRel(r01[2], REAL(1e-10)));
	}

	TEST_CASE("Tensors::Tensor3_Contract_error_handling", "[contraction]")
	{
		TEST_PRECISION_INFO();
		Tensor3<3> t3(1, 2);

		// Same index
		REQUIRE_THROWS_AS(t3.Contract(1, 1), TensorIndexError);

		// Out of range
		REQUIRE_THROWS_AS(t3.Contract(-1, 1), TensorIndexError);
		REQUIRE_THROWS_AS(t3.Contract(0, 3), TensorIndexError);
		REQUIRE_THROWS_AS(t3.Contract(5, 0), TensorIndexError);
	}

	TEST_CASE("Tensors::Tensor4_Contract_all_index_pairs", "[contraction]")
	{
		TEST_PRECISION_INFO();

		// Create Tensor4<3> with simple values: T[i][j][k][l] = i + 3*j + 9*k + 27*l
		Tensor4<3> t4(2, 2);  // 2 covariant, 2 contravariant
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++)
						t4(i, j, k, l) = Real(i + 3*j + 9*k + 27*l);

		// Contract indices 0 and 1: result[k][l] = sum_s T[s][s][k][l]
		auto r01 = t4.Contract(0, 1);
		// sum_s (s + 3*s + 9*k + 27*l) = 4*(0+1+2) + 3*(9*k + 27*l) = 12 + 27*k + 81*l
		REQUIRE_THAT(r01(0, 0), RealWithinRel(REAL(12.0), REAL(1e-10)));    // k=0,l=0
		REQUIRE_THAT(r01(1, 0), RealWithinRel(REAL(39.0), REAL(1e-10)));    // k=1,l=0
		REQUIRE_THAT(r01(0, 1), RealWithinRel(REAL(93.0), REAL(1e-10)));    // k=0,l=1

		// Contract indices 2 and 3: result[i][j] = sum_s T[i][j][s][s]
		auto r23 = t4.Contract(2, 3);
		// sum_s (i + 3*j + 9*s + 27*s) = 3*(i + 3*j) + 36*(0+1+2) = 3*i + 9*j + 108
		REQUIRE_THAT(r23(0, 0), RealWithinRel(REAL(108.0), REAL(1e-10)));   // i=0,j=0
		REQUIRE_THAT(r23(1, 0), RealWithinRel(REAL(111.0), REAL(1e-10)));   // i=1,j=0
		REQUIRE_THAT(r23(0, 1), RealWithinRel(REAL(117.0), REAL(1e-10)));   // i=0,j=1

		// Contract indices 0 and 2: result[j][l] = sum_s T[s][j][s][l]
		auto r02 = t4.Contract(0, 2);
		// sum_s (s + 3*j + 9*s + 27*l) = 10*(0+1+2) + 3*(3*j + 27*l) = 30 + 9*j + 81*l
		REQUIRE_THAT(r02(0, 0), RealWithinRel(REAL(30.0), REAL(1e-10)));    // j=0,l=0
		REQUIRE_THAT(r02(1, 0), RealWithinRel(REAL(39.0), REAL(1e-10)));    // j=1,l=0

		// Contract indices 1 and 3: result[i][k] = sum_s T[i][s][k][s]
		auto r13 = t4.Contract(1, 3);
		// sum_s (i + 3*s + 9*k + 27*s) = 3*(i + 9*k) + 30*(0+1+2) = 3*i + 27*k + 90
		REQUIRE_THAT(r13(0, 0), RealWithinRel(REAL(90.0), REAL(1e-10)));    // i=0,k=0
		REQUIRE_THAT(r13(1, 0), RealWithinRel(REAL(93.0), REAL(1e-10)));    // i=1,k=0
		REQUIRE_THAT(r13(0, 1), RealWithinRel(REAL(117.0), REAL(1e-10)));   // i=0,k=1
	}

	TEST_CASE("Tensors::Tensor4_Contract_error_handling", "[contraction]")
	{
		TEST_PRECISION_INFO();
		Tensor4<3> t4(2, 2);

		// Same index
		REQUIRE_THROWS_AS(t4.Contract(2, 2), TensorIndexError);

		// Out of range
		REQUIRE_THROWS_AS(t4.Contract(-1, 2), TensorIndexError);
		REQUIRE_THROWS_AS(t4.Contract(0, 4), TensorIndexError);
	}

	TEST_CASE("Tensors::Tensor5_Contract_selected_pairs", "[contraction]")
	{
		TEST_PRECISION_INFO();

		// Create Tensor5<2> (smaller dimension for faster test)
		// T[i][j][k][l][m] = i + 2*j + 4*k + 8*l + 16*m
		Tensor5<2> t5(2, 3);  // 2 covariant, 3 contravariant
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
					for (int l = 0; l < 2; l++)
						for (int m = 0; m < 2; m++)
							t5(i, j, k, l, m) = Real(i + 2*j + 4*k + 8*l + 16*m);

		// Contract indices 0 and 1: result[k][l][m] = sum_s T[s][s][k][l][m]
		auto r01 = t5.Contract(0, 1);
		// sum_s (s + 2*s + 4*k + 8*l + 16*m) = 3*(0+1) + 2*(4*k + 8*l + 16*m)
		// = 3 + 8*k + 16*l + 32*m
		REQUIRE_THAT(r01(0, 0, 0), RealWithinRel(REAL(3.0), REAL(1e-10)));   // k=0,l=0,m=0
		REQUIRE_THAT(r01(1, 0, 0), RealWithinRel(REAL(11.0), REAL(1e-10)));  // k=1,l=0,m=0
		REQUIRE_THAT(r01(0, 1, 0), RealWithinRel(REAL(19.0), REAL(1e-10)));  // k=0,l=1,m=0
		REQUIRE_THAT(r01(0, 0, 1), RealWithinRel(REAL(35.0), REAL(1e-10)));  // k=0,l=0,m=1

		// Contract indices 3 and 4: result[i][j][k] = sum_s T[i][j][k][s][s]
		auto r34 = t5.Contract(3, 4);
		// sum_s (i + 2*j + 4*k + 8*s + 16*s) = 2*(i + 2*j + 4*k) + 24*(0+1)
		// = 2*i + 4*j + 8*k + 24
		REQUIRE_THAT(r34(0, 0, 0), RealWithinRel(REAL(24.0), REAL(1e-10)));  // i=0,j=0,k=0
		REQUIRE_THAT(r34(1, 0, 0), RealWithinRel(REAL(26.0), REAL(1e-10)));  // i=1,j=0,k=0
		REQUIRE_THAT(r34(0, 1, 0), RealWithinRel(REAL(28.0), REAL(1e-10)));  // i=0,j=1,k=0
		REQUIRE_THAT(r34(0, 0, 1), RealWithinRel(REAL(32.0), REAL(1e-10)));  // i=0,j=0,k=1

		// Contract indices 1 and 3: result[i][k][m] = sum_s T[i][s][k][s][m]
		auto r13 = t5.Contract(1, 3);
		// sum_s (i + 2*s + 4*k + 8*s + 16*m) = 2*i + 10*(0+1) + 2*4*k + 2*16*m
		// = 2*i + 10 + 8*k + 32*m
		REQUIRE_THAT(r13(0, 0, 0), RealWithinRel(REAL(10.0), REAL(1e-10)));  // i=0,k=0,m=0
		REQUIRE_THAT(r13(1, 0, 0), RealWithinRel(REAL(12.0), REAL(1e-10)));  // i=1,k=0,m=0
		REQUIRE_THAT(r13(0, 1, 0), RealWithinRel(REAL(18.0), REAL(1e-10)));  // i=0,k=1,m=0
	}

	TEST_CASE("Tensors::Tensor5_Contract_error_handling", "[contraction]")
	{
		TEST_PRECISION_INFO();
		Tensor5<2> t5(2, 3);

		// Same index
		REQUIRE_THROWS_AS(t5.Contract(3, 3), TensorIndexError);

		// Out of range
		REQUIRE_THROWS_AS(t5.Contract(-1, 2), TensorIndexError);
		REQUIRE_THROWS_AS(t5.Contract(0, 5), TensorIndexError);
		REQUIRE_THROWS_AS(t5.Contract(6, 1), TensorIndexError);
	}

	TEST_CASE("Tensors::Tensor5_vector_application", "[operations]")
	{
		TEST_PRECISION_INFO();

		// Test operator() with 5 vectors
		Tensor5<2> t5(2, 3);
		// Set only diagonal: T[i][i][i][i][i] = i+1
		t5(0, 0, 0, 0, 0) = REAL(1.0);
		t5(1, 1, 1, 1, 1) = REAL(2.0);

		VectorN<Real, 2> v1({REAL(1.0), REAL(0.0)});
		VectorN<Real, 2> v2({REAL(1.0), REAL(0.0)});
		VectorN<Real, 2> v3({REAL(1.0), REAL(0.0)});
		VectorN<Real, 2> v4({REAL(1.0), REAL(0.0)});
		VectorN<Real, 2> v5({REAL(1.0), REAL(0.0)});

		// Only (0,0,0,0,0) contributes: 1 * 1 * 1 * 1 * 1 * 1 = 1
		auto result = t5(v1, v2, v3, v4, v5);
		REQUIRE_THAT(result, RealWithinRel(REAL(1.0), REAL(1e-10)));

		// All ones
		VectorN<Real, 2> ones({REAL(1.0), REAL(1.0)});
		// T[0,0,0,0,0]*1 + T[1,1,1,1,1]*1 = 1 + 2 = 3
		result = t5(ones, ones, ones, ones, ones);
		REQUIRE_THAT(result, RealWithinRel(REAL(3.0), REAL(1e-10)));
	}

	TEST_CASE("Tensors::Kronecker_delta_contraction", "[contraction]")
	{
		TEST_PRECISION_INFO();

		// Create Kronecker delta as Tensor2 (identity)
		Tensor2<3> delta(1, 1);  // mixed tensor
		delta(0, 0) = REAL(1.0);
		delta(1, 1) = REAL(1.0);
		delta(2, 2) = REAL(1.0);

		// Contraction of Kronecker delta = trace = dimension
		auto trace = delta.Contract();
		REQUIRE_THAT(trace, RealWithinRel(REAL(3.0), REAL(1e-10)));

		// Create a Tensor4 that is product of two Kronecker deltas
		// delta^i_j * delta^k_l
		Tensor4<3> delta4(2, 2);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++)
						delta4(i, j, k, l) = (i == j ? 1.0 : 0.0) * (k == l ? 1.0 : 0.0);

		// Contract first pair (0,1): should give delta^k_l * 3
		auto contracted = delta4.Contract(0, 1);
		for (int k = 0; k < 3; k++)
			for (int l = 0; l < 3; l++)
				REQUIRE_THAT(contracted(k, l), RealWithinRel(k == l ? REAL(3.0) : REAL(0.0), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****              Additional Coverage Tests                    *****/
	/*********************************************************************/

	TEST_CASE("Tensors::Tensor2_GetMatrix", "[Tensor2]")
	{
		TEST_PRECISION_INFO();
		
		Tensor2<3> t(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
		                     REAL(4.0), REAL(5.0), REAL(6.0),
		                     REAL(7.0), REAL(8.0), REAL(9.0) });
		
		auto mat = t.GetMatrix();
		
		REQUIRE(mat[0][0] == REAL(1.0));
		REQUIRE(mat[0][1] == REAL(2.0));
		REQUIRE(mat[1][0] == REAL(4.0));
		REQUIRE(mat[2][2] == REAL(9.0));
	}

	TEST_CASE("Tensors::Tensor2_IsContravar_IsCovar", "[Tensor2]")
	{
		TEST_PRECISION_INFO();
		
		// 1 covariant (index 0), 1 contravariant (index 1)
		Tensor2<3> t_mixed(1, 1);
		REQUIRE(t_mixed.IsCovar(0) == true);
		REQUIRE(t_mixed.IsContravar(0) == false);
		REQUIRE(t_mixed.IsCovar(1) == false);
		REQUIRE(t_mixed.IsContravar(1) == true);
		
		// 2 covariant, 0 contravariant
		Tensor2<3> t_covar(2, 0);
		REQUIRE(t_covar.IsCovar(0) == true);
		REQUIRE(t_covar.IsCovar(1) == true);
		REQUIRE(t_covar.IsContravar(0) == false);
		REQUIRE(t_covar.IsContravar(1) == false);
		
		// 0 covariant, 2 contravariant
		Tensor2<3> t_contra(0, 2);
		REQUIRE(t_contra.IsContravar(0) == true);
		REQUIRE(t_contra.IsContravar(1) == true);
	}

	TEST_CASE("Tensors::Tensor2_arithmetic_subtraction_scalar", "[Tensor2]")
	{
		TEST_PRECISION_INFO();
		
		Tensor2<3> t1(1, 1, { REAL(10.0), REAL(20.0), REAL(30.0),
		                      REAL(40.0), REAL(50.0), REAL(60.0),
		                      REAL(70.0), REAL(80.0), REAL(90.0) });
		Tensor2<3> t2(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0),
		                      REAL(4.0), REAL(5.0), REAL(6.0),
		                      REAL(7.0), REAL(8.0), REAL(9.0) });
		
		// Subtraction
		auto diff = t1 - t2;
		REQUIRE(diff(0, 0) == REAL(9.0));
		REQUIRE(diff(1, 1) == REAL(45.0));
		REQUIRE(diff(2, 2) == REAL(81.0));
		
		// Scalar multiplication
		auto scaled = t2 * REAL(3.0);
		REQUIRE(scaled(0, 0) == REAL(3.0));
		REQUIRE(scaled(1, 1) == REAL(15.0));
		REQUIRE(scaled(2, 2) == REAL(27.0));
		
		// Scalar multiplication (commutative)
		auto scaled2 = REAL(2.0) * t2;
		REQUIRE(scaled2(0, 0) == REAL(2.0));
		REQUIRE(scaled2(1, 1) == REAL(10.0));
		
		// Scalar division
		auto divided = t1 / REAL(10.0);
		REQUIRE(divided(0, 0) == REAL(1.0));
		REQUIRE(divided(1, 1) == REAL(5.0));
	}

	TEST_CASE("Tensors::Tensor2_type_mismatch_arithmetic", "[Tensor2][errors]")
	{
		TEST_PRECISION_INFO();
		
		Tensor2<3> t1(2, 0);  // 2 covariant
		Tensor2<3> t2(1, 1);  // 1 covariant, 1 contravariant
		
		REQUIRE_THROWS_AS(t1 + t2, TensorCovarContravarArithmeticError);
		REQUIRE_THROWS_AS(t1 - t2, TensorCovarContravarArithmeticError);
	}

	TEST_CASE("Tensors::Tensor2_Contract_wrong_type", "[Tensor2][errors]")
	{
		TEST_PRECISION_INFO();
		
		// Contract requires exactly 1 covariant and 1 contravariant
		Tensor2<3> t_covar(2, 0);
		REQUIRE_THROWS_AS(t_covar.Contract(), TensorCovarContravarNumError);
		
		Tensor2<3> t_contra(0, 2);
		REQUIRE_THROWS_AS(t_contra.Contract(), TensorCovarContravarNumError);
	}

	TEST_CASE("Tensors::Tensor2_stream_output", "[Tensor2][output]")
	{
		TEST_PRECISION_INFO();
		
		Tensor2<2> t(1, 1, { REAL(1.0), REAL(2.0),
		                     REAL(3.0), REAL(4.0) });
		
		std::ostringstream oss;
		oss << t;
		std::string output = oss.str();
		
		REQUIRE(!output.empty());
		REQUIRE(output.find("N = 2") != std::string::npos);
		REQUIRE(output.find("[") != std::string::npos);
	}

	TEST_CASE("Tensors::Tensor2_Print", "[Tensor2][output]")
	{
		TEST_PRECISION_INFO();
		
		Tensor2<2> t(1, 1, { REAL(1.0), REAL(2.0),
		                     REAL(3.0), REAL(4.0) });
		
		std::ostringstream oss;
		t.Print(oss, 8, 3);
		std::string output = oss.str();
		
		REQUIRE(!output.empty());
	}

	TEST_CASE("Tensors::Tensor2_checked_access_at", "[Tensor2][errors]")
	{
		TEST_PRECISION_INFO();
		
		Tensor2<3> t(1, 1);
		t(0, 0) = REAL(5.0);
		
		// Valid access
		REQUIRE(t.at(0, 0) == REAL(5.0));
		
		// Invalid access - should throw
		REQUIRE_THROWS_AS(t.at(-1, 0), TensorIndexError);
		REQUIRE_THROWS_AS(t.at(0, 3), TensorIndexError);
		REQUIRE_THROWS_AS(t.at(3, 0), TensorIndexError);
	}

	TEST_CASE("Tensors::Tensor3_initializer_list", "[Tensor3]")
	{
		TEST_PRECISION_INFO();
		
		// 2x2x2 tensor with initializer list
		Tensor3<2> t(1, 2, { REAL(1.0), REAL(2.0),
		                      REAL(3.0), REAL(4.0),
		                      REAL(5.0), REAL(6.0),
		                      REAL(7.0), REAL(8.0) });
		
		REQUIRE(t(0, 0, 0) == REAL(1.0));
		REQUIRE(t(0, 0, 1) == REAL(2.0));
		REQUIRE(t(0, 1, 0) == REAL(3.0));
		REQUIRE(t(0, 1, 1) == REAL(4.0));
		REQUIRE(t(1, 0, 0) == REAL(5.0));
		REQUIRE(t(1, 1, 1) == REAL(8.0));
	}

	TEST_CASE("Tensors::Tensor3_checked_access_at", "[Tensor3][errors]")
	{
		TEST_PRECISION_INFO();
		
		Tensor3<3> t(1, 2);
		t(0, 1, 2) = REAL(7.0);
		
		// Valid access
		REQUIRE(t.at(0, 1, 2) == REAL(7.0));
		
		// Invalid access
		REQUIRE_THROWS_AS(t.at(-1, 0, 0), TensorIndexError);
		REQUIRE_THROWS_AS(t.at(0, 3, 0), TensorIndexError);
		REQUIRE_THROWS_AS(t.at(0, 0, 3), TensorIndexError);
	}

	TEST_CASE("Tensors::Tensor3_vector_application", "[Tensor3]")
	{
		TEST_PRECISION_INFO();
		
		Tensor3<2> t(1, 2);
		// Set T[i][j][k] = 1 for all
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
					t(i, j, k) = REAL(1.0);
		
		VectorN<Real, 2> v1({REAL(1.0), REAL(1.0)});
		VectorN<Real, 2> v2({REAL(1.0), REAL(1.0)});
		VectorN<Real, 2> v3({REAL(1.0), REAL(1.0)});
		
		// Sum of all 8 elements times 1*1*1 = 8
		auto result = t(v1, v2, v3);
		REQUIRE_THAT(result, RealWithinRel(REAL(8.0), REAL(1e-10)));
	}

	TEST_CASE("Tensors::Tensor4_checked_access_at", "[Tensor4][errors]")
	{
		TEST_PRECISION_INFO();
		
		Tensor4<3> t(2, 2);
		t(0, 1, 2, 0) = REAL(11.0);
		
		// Valid access
		REQUIRE(t.at(0, 1, 2, 0) == REAL(11.0));
		
		// Invalid access
		REQUIRE_THROWS_AS(t.at(-1, 0, 0, 0), TensorIndexError);
		REQUIRE_THROWS_AS(t.at(0, 0, 0, 3), TensorIndexError);
	}

	TEST_CASE("Tensors::Tensor4_vector_application", "[Tensor4]")
	{
		TEST_PRECISION_INFO();
		
		Tensor4<2> t(2, 2);
		// Set only diagonal
		t(0, 0, 0, 0) = REAL(1.0);
		t(1, 1, 1, 1) = REAL(16.0);
		
		VectorN<Real, 2> e0({REAL(1.0), REAL(0.0)});
		VectorN<Real, 2> e1({REAL(0.0), REAL(1.0)});
		
		// Only (0,0,0,0) contributes
		auto r1 = t(e0, e0, e0, e0);
		REQUIRE_THAT(r1, RealWithinRel(REAL(1.0), REAL(1e-10)));
		
		// Only (1,1,1,1) contributes
		auto r2 = t(e1, e1, e1, e1);
		REQUIRE_THAT(r2, RealWithinRel(REAL(16.0), REAL(1e-10)));
	}

	TEST_CASE("Tensors::Tensor5_checked_access_at", "[Tensor5][errors]")
	{
		TEST_PRECISION_INFO();
		
		Tensor5<2> t(2, 3);
		t(0, 1, 0, 1, 0) = REAL(13.0);
		
		// Valid access
		REQUIRE(t.at(0, 1, 0, 1, 0) == REAL(13.0));
		
		// Invalid access
		REQUIRE_THROWS_AS(t.at(-1, 0, 0, 0, 0), TensorIndexError);
		REQUIRE_THROWS_AS(t.at(0, 0, 0, 0, 2), TensorIndexError);
	}

	TEST_CASE("Tensors::Tensor_covar_contravar_ctor_errors", "[errors]")
	{
		TEST_PRECISION_INFO();
		
		// Tensor3: must sum to 3
		REQUIRE_THROWS_AS(Tensor3<3>(1, 1), TensorCovarContravarNumError);  // 1+1=2
		REQUIRE_THROWS_AS(Tensor3<3>(0, 4), TensorCovarContravarNumError);  // 0+4=4
		
		// Tensor4: must sum to 4
		REQUIRE_THROWS_AS(Tensor4<3>(1, 1), TensorCovarContravarNumError);  // 1+1=2
		REQUIRE_THROWS_AS(Tensor4<3>(3, 3), TensorCovarContravarNumError);  // 3+3=6
		
		// Tensor5: must sum to 5
		REQUIRE_THROWS_AS(Tensor5<2>(2, 2), TensorCovarContravarNumError);  // 2+2=4
		REQUIRE_THROWS_AS(Tensor5<2>(0, 6), TensorCovarContravarNumError);  // 0+6=6
	}

	TEST_CASE("Tensors::Tensor2_initializer_list_overflow_throws", "[Tensor2]")
	{
		// 2x2 tensor = 4 elements, providing 5 should throw
		REQUIRE_THROWS_AS(
			Tensor2<2>(1, 1, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0) }),
			TensorCovarContravarNumError);
	}

	TEST_CASE("Tensors::Tensor3_initializer_list_overflow_throws", "[Tensor3]")
	{
		// 2x2x2 tensor = 8 elements, providing 9 should throw
		REQUIRE_THROWS_AS(
			Tensor3<2>(1, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0),
			                   REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0) }),
			TensorCovarContravarNumError);
	}

}

