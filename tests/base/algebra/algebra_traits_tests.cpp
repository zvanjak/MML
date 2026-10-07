#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

#include <type_traits>

namespace MML::Tests::Base::AlgebraTests
{
	struct RuntimeGroup
	{
		using element_type = int;
	};

	struct ScalarAwareStructure
	{
		using element_type = int;
		using scalar_type = double;
	};

	TEST_CASE("AlgebraTraits expose conservative metadata defaults", "[Algebra][Base]")
	{
		using Traits = Algebra::AlgebraTraits<RuntimeGroup>;

		static_assert(std::is_same_v<Traits::element_type, int>);
		static_assert(std::is_same_v<Traits::scalar_type, int>);
		static_assert(Traits::static_order == Algebra::DynamicAlgebraExtent);
		static_assert(Traits::dimension == Algebra::DynamicAlgebraExtent);
		static_assert(Traits::equality == Algebra::AlgebraEquality::Exact);

		REQUIRE(Traits::static_order == Algebra::DynamicAlgebraExtent);
	}

	TEST_CASE("AlgebraTraits detect an explicit scalar type", "[Algebra][Base]")
	{
		using Traits = Algebra::AlgebraTraits<ScalarAwareStructure>;

		static_assert(std::is_same_v<Traits::element_type, int>);
		static_assert(std::is_same_v<Traits::scalar_type, double>);
		REQUIRE(Traits::equality == Algebra::AlgebraEquality::Exact);
	}
}