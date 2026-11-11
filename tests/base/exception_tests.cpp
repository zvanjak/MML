///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        exception_tests.cpp                                                 ///
///  Description: Tests for MMLException base class and exception hierarchy           ///
///               Verifies catch(const MMLException&) catches all library exceptions  ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "../../mml/MMLExceptions.h"
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <type_traits>

using namespace MML;

TEST_CASE("Exception string constructors are explicit", "[exceptions][wewi.3.6]") {
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, IndexError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, SingularMatrixError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, IntegrationTooManySteps>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, RealFuncInterpInitError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, RealFuncInterpRuntimeError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, TensorIndexError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, RootFindingError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, ODESolverError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, StatisticsError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, NumericalMethodError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, MatrixNumericalError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, ConvergenceError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, GeometryError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, QuaternionError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, FileIOError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, FourierError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, VisualizerError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, CurveFittingError>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, DataError>);
	STATIC_REQUIRE(std::is_constructible_v<ConvergenceError, std::string>);
}

TEST_CASE("Vector checked access is an index error with correct metadata", "[exceptions][wewi.3.6]") {
	Vector<Real> vector(3);
	const Vector<Real>& constant = vector;
	VectorN<Real, 3> fixed;
	const VectorN<Real, 3>& constantFixed = fixed;
	for (int index : {-1, 3}) {
		REQUIRE_THROWS_AS(vector.at(index), VectorAccessBoundsError);
		REQUIRE_THROWS_AS(constant.at(index), IndexError);
		REQUIRE_THROWS_AS(fixed.at(index), IndexError);
		REQUIRE_THROWS_AS(constantFixed.at(index), std::out_of_range);
		try {
			(void)vector.at(index);
			FAIL("Expected index exception");
		} catch (const IndexError& error) {
			REQUIRE(error.index() == index);
			REQUIRE(error.size() == 3);
		}
	}
}

TEST_CASE("Cholesky distinguishes domain errors from dimension errors", "[exceptions][wewi.3.6]") {
	Matrix<Real> nonsymmetric(2, 2, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0)});
	REQUIRE_THROWS_AS(CholeskySolver<Real>(nonsymmetric), DomainError);
	Matrix<Real> nonsquare(2, 3);
	REQUIRE_THROWS_AS(CholeskySolver<Real>(nonsquare), MatrixDimensionError);
	Matrix<Complex> complexDiagonal(1, 1, {Complex(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0))});
	REQUIRE_THROWS_AS(CholeskySolver<Complex>(complexDiagonal), DomainError);
}

TEST_CASE("MMLException - catch all MML exceptions", "[exceptions]")
{
	SECTION("ArgumentError caught as MMLException") {
		bool caught = false;
		try {
			throw ArgumentError("test argument error");
		} catch (const MMLException& e) {
			caught = true;
			REQUIRE(std::string(e.message()) == "test argument error");
		}
		REQUIRE(caught);
	}

	SECTION("DomainError caught as MMLException") {
		bool caught = false;
		try {
			throw DomainError("domain violation");
		} catch (const MMLException& e) {
			caught = true;
			REQUIRE(std::string(e.message()) == "domain violation");
		}
		REQUIRE(caught);
	}

	SECTION("VectorDimensionError caught as MMLException") {
		bool caught = false;
		try {
			throw VectorDimensionError("dim mismatch", 3, 5);
		} catch (const MMLException& e) {
			caught = true;
			REQUIRE(std::string(e.message()) == "dim mismatch");
		}
		REQUIRE(caught);
	}

	SECTION("MatrixAccessBoundsError caught as MMLException") {
		bool caught = false;
		try {
			throw MatrixAccessBoundsError("out of bounds", 10, 5, 3, 3);
		} catch (const MMLException& e) {
			caught = true;
			REQUIRE(std::string(e.message()) == "out of bounds");
		}
		REQUIRE(caught);
	}

	SECTION("SingularMatrixError caught as MMLException") {
		bool caught = false;
		try {
			throw SingularMatrixError("singular matrix", 0.0, 2);
		} catch (const MMLException& e) {
			caught = true;
			REQUIRE(std::string(e.message()) == "singular matrix");
		}
		REQUIRE(caught);
	}

	SECTION("ConvergenceError caught as MMLException") {
		bool caught = false;
		try {
			throw ConvergenceError("did not converge", 1000, 1e-3);
		} catch (const MMLException& e) {
			caught = true;
			REQUIRE(std::string(e.message()) == "did not converge");
		}
		REQUIRE(caught);
	}

	SECTION("NotImplementedError caught as MMLException") {
		bool caught = false;
		try {
			throw NotImplementedError("not implemented");
		} catch (const MMLException& e) {
			caught = true;
			REQUIRE(std::string(e.message()) == "not implemented");
		}
		REQUIRE(caught);
	}
}

TEST_CASE("MMLException - interoperates with std exception types", "[exceptions]")
{
	SECTION("ArgumentError still caught as std::invalid_argument") {
		REQUIRE_THROWS_AS(throw ArgumentError("test"), std::invalid_argument);
	}

	SECTION("DomainError still caught as std::domain_error") {
		REQUIRE_THROWS_AS(throw DomainError("test"), std::domain_error);
	}

	SECTION("IndexError still caught as std::out_of_range") {
		REQUIRE_THROWS_AS(throw IndexError("test", 5, 3), std::out_of_range);
	}

	SECTION("NotImplementedError still caught as std::logic_error") {
		REQUIRE_THROWS_AS(throw NotImplementedError("test"), std::logic_error);
	}

	SECTION("RootFindingError still caught as std::runtime_error") {
		REQUIRE_THROWS_AS(throw RootFindingError("test"), std::runtime_error);
	}

	SECTION("All MML exceptions still caught as std::exception") {
		REQUIRE_THROWS_AS(throw ArgumentError("test"), std::exception);
		REQUIRE_THROWS_AS(throw DomainError("test"), std::exception);
		REQUIRE_THROWS_AS(throw ConvergenceError("test"), std::exception);
		REQUIRE_THROWS_AS(throw FileIOError("test"), std::exception);
	}
}

TEST_CASE("MMLException - enriched exception data preserved", "[exceptions]")
{
	SECTION("VectorDimensionError preserves dimensions") {
		try {
			throw VectorDimensionError("size mismatch", 3, 5);
		} catch (const VectorDimensionError& e) {
			REQUIRE(e.expected() == 3);
			REQUIRE(e.actual() == 5);
		}
	}

	SECTION("ConvergenceError preserves iteration info") {
		try {
			throw ConvergenceError("no convergence", 500, 1e-4);
		} catch (const ConvergenceError& e) {
			REQUIRE(e.iterations() == 500);
			REQUIRE(e.residual() == Catch::Approx(1e-4));
		}
	}

	SECTION("SingularMatrixError preserves pivot info") {
		try {
			throw SingularMatrixError("singular", 1e-15, 3);
		} catch (const SingularMatrixError& e) {
			REQUIRE(e.determinant() == Catch::Approx(1e-15));
			REQUIRE(e.pivot_row() == 3);
		}
	}
}
