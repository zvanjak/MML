#include <catch2/catch_test_macros.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Vector/Vector.h>
#endif

#include <complex>
#include <memory>
#include <stdexcept>

namespace {
	template<typename Type>
	void CheckResizeContract() {
		for (bool preserve : {false, true}) {
			CAPTURE(preserve);
			MML::Vector<Type> values{Type{2}, Type{4}, Type{6}};
			const Type* original = &values[0];
			values.resize(2, preserve);
			REQUIRE(values.size() == 2);
			CHECK(&values[0] == original);
			CHECK(values[0] == Type{2});
			CHECK(values[1] == Type{4});
			values.resize(3, preserve);
			REQUIRE(values.size() == 3);
			CHECK(&values[0] == original);
			CHECK(values[2] == Type{0});
			values[2] = Type{8};
			values.resize(5, preserve);
			REQUIRE(values.size() == 5);
			CHECK(values[0] == Type{2});
			CHECK(values[1] == Type{4});
			CHECK(values[2] == Type{8});
			CHECK(values[3] == Type{0});
			CHECK(values[4] == Type{0});
			const Type* grown = &values[0];
			values.resize(std::size_t{5}, preserve);
			CHECK(&values[0] == grown);
			CHECK(values[2] == Type{8});
			CHECK_THROWS_AS(values.resize(-1, preserve), std::length_error);
			REQUIRE(values.size() == 5);
			CHECK(values[2] == Type{8});
			values.resize(std::size_t{1}, preserve);
			REQUIRE(values.size() == 1);
			CHECK(values[0] == Type{2});
			values.resize(0, preserve);
			CHECK(values.empty());
			values.resize(std::size_t{2}, preserve);
			REQUIRE(values.size() == 2);
			CHECK(values[0] == Type{0});
			CHECK(values[1] == Type{0});
		}
	}

	struct DefaultValue {
		int value = 7;
	};
}

TEST_CASE("Vector resize preserves prefixes and initializes every numeric tail", "[Vector][resize-contract]") {
	CheckResizeContract<MML::Real>();
	CheckResizeContract<int>();
	CheckResizeContract<unsigned>();
	CheckResizeContract<std::complex<MML::Real>>();
}

TEST_CASE("Vector resize retains default construction and move-only element support", "[Vector][resize-contract]") {
	MML::Vector<DefaultValue> defaults(2);
	defaults[0].value = 11;
	defaults.resize(5);
	REQUIRE(defaults.size() == 5);
	CHECK(defaults[0].value == 11);
	for (int i = 1; i < defaults.size(); ++i)
		CHECK(defaults[i].value == 7);

	MML::Vector<std::unique_ptr<int>> values(2);
	values[0] = std::make_unique<int>(7);
	values[1] = std::make_unique<int>(11);
	const int* first = values[0].get();
	values.resize(5);
	REQUIRE(values.size() == 5);
	CHECK(values[0].get() == first);
	CHECK(*values[1] == 11);
	for (int i = 2; i < values.size(); ++i)
		CHECK_FALSE(values[i]);
	values.resize(1, true);
	REQUIRE(values.size() == 1);
	CHECK(values[0].get() == first);
	values.resize(3, true);
	REQUIRE(values.size() == 3);
	CHECK(values[0].get() == first);
	CHECK_FALSE(values[1]);
	CHECK_FALSE(values[2]);
	values.resize(0);
	CHECK(values.empty());
}
