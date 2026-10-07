#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

#include <vector>

namespace MML::Tests::Core::AlgebraTests
{
	class CyclicGroup3
	{
	public:
		using element_type = int;

		const std::vector<int>& elements() const
		{
			static const std::vector<int> values{0, 1, 2};
			return values;
		}

		int identity() const { return 0; }
		int compose(int left, int right) const { return (left + right) % 3; }
		int inverse(int element) const { return (3 - element) % 3; }
	};

	class NonAssociativeStructure : public CyclicGroup3
	{
	public:
		int compose(int left, int right) const { return (left - right + 3) % 3; }
	};

	template<int Modulus>
	class ResidueRing
	{
	public:
		using element_type = int;

		const std::vector<int>& elements() const
		{
			static const std::vector<int> values = [] {
				std::vector<int> result;
				for (int value = 0; value < Modulus; value++)
					result.push_back(value);
				return result;
			}();
			return values;
		}

		int zero() const { return 0; }
		int one() const { return 1; }
		int add(int left, int right) const { return (left + right) % Modulus; }
		int negate(int element) const { return (Modulus - element) % Modulus; }
		int multiply(int left, int right) const { return (left * right) % Modulus; }

		int inverse(int element) const
		{
			for (int candidate = 1; candidate < Modulus; candidate++)
				if (multiply(element, candidate) == one())
					return candidate;
			return 0;
		}
	};

	TEST_CASE("Algebra law checks accept the cyclic group C3", "[Algebra][Core]")
	{
		CyclicGroup3 group;

		REQUIRE(Algebra::CheckClosure(group));
		REQUIRE(Algebra::CheckAssociativity(group));
		REQUIRE(Algebra::CheckIdentity(group));
		REQUIRE(Algebra::CheckInverses(group));
		REQUIRE(Algebra::CheckGroupLaws(group));
	}

	TEST_CASE("Algebra law checks reject a non-associative operation", "[Algebra][Core]")
	{
		NonAssociativeStructure structure;

		REQUIRE_FALSE(Algebra::CheckAssociativity(structure));
		REQUIRE_FALSE(Algebra::CheckGroupLaws(structure));
	}

	TEST_CASE("Algebra law checks distinguish rings from fields", "[Algebra][Core]")
	{
		ResidueRing<4> integersMod4;
		ResidueRing<5> integersMod5;

		REQUIRE(Algebra::CheckRingLaws(integersMod4));
		REQUIRE_FALSE(Algebra::CheckFieldLaws(integersMod4));
		REQUIRE(Algebra::CheckRingLaws(integersMod5));
		REQUIRE(Algebra::CheckFieldLaws(integersMod5));
	}

	TEST_CASE("Algebra law checks accept a caller-provided equality policy", "[Algebra][Core]")
	{
		CyclicGroup3 group;
		auto equal = [](int left, int right) { return left == right; };

		REQUIRE(Algebra::CheckGroupLaws(group, equal));
	}

	TEST_CASE("Algebra law checks validate the finite group primitive", "[Algebra][Core]")
	{
		Algebra::FiniteGroup<int> group(
			{0, 1, 2}, 0,
			[](const int& left, const int& right) { return (left + right) % 3; },
			[](const int& element) { return (3 - element) % 3; });

		REQUIRE(Algebra::CheckGroupLaws(group));
	}
}