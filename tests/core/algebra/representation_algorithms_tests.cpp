#include <catch2/catch_test_macros.hpp>

#include "../../TestPrecision.h"
#include <mml/base/Algebra_base.h>

#include <cmath>

namespace MML::Tests::Core::AlgebraTests
{
	constexpr Real RepresentationTolerance = Testing::Tol(REAL(1e-12), REAL(1e-6));

	bool MatrixNear(const MatrixNM<Real, 2, 2>& left, const MatrixNM<Real, 2, 2>& right)
	{
		return left.IsEqualTo(right, RepresentationTolerance);
	}

	bool VectorNear(const VectorN<Real, 2>& left, const VectorN<Real, 2>& right)
	{
		return left.IsEqualTo(right, RepresentationTolerance);
	}

	auto CyclicRotationRepresentation(const Algebra::CyclicGroup& group)
	{
		return Algebra::MakeRepresentation<Algebra::CyclicElement, Real, 2>(
			[degree = group.degree()](const Algebra::CyclicElement& element) {
				const Real angle = REAL(2.0) * Constants::PI * element.exponent / degree;
				const Real cosine = std::cos(angle);
				const Real sine = std::sin(angle);
				return MatrixNM<Real, 2, 2>{{cosine, -sine}, {sine, cosine}};
			});
	}

	auto DihedralPlaneRepresentation(const Algebra::DihedralGroup& group)
	{
		return Algebra::MakeRepresentation<Algebra::DihedralElement, Real, 2>(
			[degree = group.degree()](const Algebra::DihedralElement& element) {
				const Real angle = REAL(2.0) * Constants::PI * element.rotation / degree;
				const Real cosine = std::cos(angle);
				const Real sine = std::sin(angle);
				MatrixNM<Real, 2, 2> rotation{{cosine, -sine}, {sine, cosine}};
				MatrixNM<Real, 2, 2> reflection{{1, 0}, {0, -1}};
				return element.reflected ? rotation * reflection : rotation;
			});
	}

	TEST_CASE("C3 rotation representation satisfies homomorphism and characters", "[Algebra][Representation][Core]")
	{
		Algebra::CyclicGroup group(3);
		const auto representation = CyclicRotationRepresentation(group);

		REQUIRE(Algebra::VerifyRepresentation(group, representation, MatrixNear));
		REQUIRE(std::abs(Algebra::Character(representation, group.identity()) - REAL(2.0)) < RepresentationTolerance);
		REQUIRE(std::abs(Algebra::Character(representation, group.generator()) + REAL(1.0)) < RepresentationTolerance);
		REQUIRE(Algebra::IsCharacterConstantOnConjugacyClasses(group, representation,
			[](Real left, Real right) { return std::abs(left - right) < RepresentationTolerance; }));
	}

	TEST_CASE("D3 standard representation has class-constant character", "[Algebra][Representation][Core]")
	{
		Algebra::DihedralGroup group(3);
		const auto representation = DihedralPlaneRepresentation(group);

		REQUIRE(Algebra::VerifyRepresentation(group, representation, MatrixNear));
		REQUIRE(std::abs(Algebra::Character(representation, group.identity()) - REAL(2.0)) < RepresentationTolerance);
		REQUIRE(std::abs(Algebra::Character(representation, group.generator()) + REAL(1.0)) < RepresentationTolerance);
		REQUIRE(std::abs(Algebra::Character(representation, group.reflection())) < RepresentationTolerance);
		REQUIRE(Algebra::IsCharacterConstantOnConjugacyClasses(group, representation,
			[](Real left, Real right) { return std::abs(left - right) < RepresentationTolerance; }));
	}

	TEST_CASE("Invariant projection is idempotent and produces fixed vectors", "[Algebra][Representation][Core]")
	{
		Algebra::CyclicGroup group(2);
		const auto representation = Algebra::MakeRepresentation<Algebra::CyclicElement, Real, 2>(
			[](const Algebra::CyclicElement& element) {
				return element.exponent == 0
					? MatrixNM<Real, 2, 2>::Identity()
					: MatrixNM<Real, 2, 2>{{1, 0}, {0, -1}};
			});
		const auto projection = Algebra::InvariantProjection(group, representation);
		const VectorN<Real, 2> vector{2, -1};
		const auto symmetrized = Algebra::SymmetrizeVector(group, representation, vector);

		REQUIRE(MatrixNear(projection * projection, projection));
		REQUIRE(VectorNear(symmetrized, VectorN<Real, 2>{2, 0}));
		for (const auto& element : group.elements())
			REQUIRE(VectorNear(representation.apply(element, symmetrized), symmetrized));
	}

	TEST_CASE("GroupAverage symmetrizes additive sample values", "[Algebra][Representation][Core]")
	{
		Algebra::CyclicGroup group(3);
		const auto representation = CyclicRotationRepresentation(group);
		const VectorN<Real, 2> vector{3, 1};
		const auto averaged = Algebra::GroupAverage(group,
			[&representation](const Algebra::CyclicElement& element, const VectorN<Real, 2>& value) {
				return representation.apply(element, value);
			}, vector, VectorN<Real, 2>{0, 0},
			[](VectorN<Real, 2> sum, const VectorN<Real, 2>& value) { return sum + value; },
			[](VectorN<Real, 2> sum, std::size_t count) { return sum / static_cast<Real>(count); });

		REQUIRE(VectorNear(averaged, Algebra::SymmetrizeVector(group, representation, vector)));
	}
}