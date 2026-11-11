#include <catch2/catch_all.hpp>
#include <mml/MMLExceptions.h>
#include <mml/algorithms/Fourier/Fourier.h>
#include <mml/algorithms/Statistics.h>
#include <mml/algorithms/CurveFitting.h>

#include <stdexcept>
#include <string>
#include <type_traits>

using namespace MML;

TEMPLATE_TEST_CASE("Centralized exceptions retain the shared exception contract",
                   "[exceptions][ownership]", FourierError, StatisticsError, CurveFittingError)
{
	STATIC_REQUIRE(std::is_base_of_v<MMLException, TestType>);
	STATIC_REQUIRE(std::is_base_of_v<std::exception, TestType>);
	STATIC_REQUIRE(std::is_constructible_v<TestType, std::string>);
	STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, TestType>);
	const std::string message = "shared domain exception";
	const TestType error(message);
	REQUIRE(std::string(error.what()) == message);
	REQUIRE(std::string(error.message()) == message);
	bool caughtAsMarker = false;
	try {
		throw error;
	} catch (const MMLException& shared) {
		caughtAsMarker = true;
		REQUIRE(dynamic_cast<const TestType*>(&shared) != nullptr);
		REQUIRE(std::string(shared.message()) == message);
	}
	REQUIRE(caughtAsMarker);
	bool caughtAsStandard = false;
	try {
		throw error;
	} catch (const std::exception& standard) {
		caughtAsStandard = true;
		REQUIRE(dynamic_cast<const TestType*>(&standard) != nullptr);
		REQUIRE(std::string(standard.what()) == message);
	}
	REQUIRE(caughtAsStandard);
}

TEST_CASE("Core algorithm error paths still throw their exact domain exceptions", "[exceptions][ownership]")
{
	STATIC_REQUIRE(std::is_base_of_v<std::invalid_argument, FourierError>);
	STATIC_REQUIRE(std::is_base_of_v<std::runtime_error, StatisticsError>);
	STATIC_REQUIRE(std::is_base_of_v<std::invalid_argument, CurveFittingError>);
	Vector<Complex> emptySpectrum;
	REQUIRE_THROWS_AS(Fourier::FFT::Transform(emptySpectrum), FourierError);
	REQUIRE_THROWS_AS(Statistics::Avg(Vector<Real>{}), StatisticsError);
	Real slope{}, intercept{};
	REQUIRE_THROWS_AS(LinearLeastSquares(Vector<Real>{}, Vector<Real>{}, slope, intercept), CurveFittingError);
}
