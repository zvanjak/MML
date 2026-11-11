#include <catch2/catch_all.hpp>
#include "../../TestMatchers.h"

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/base/Tensor.h>
#endif

#include <array>
#include <limits>
#include <type_traits>
#include <utility>

using namespace MML;

namespace {
	template<int Rank>
	using TensorOfRank = std::conditional_t<Rank == 2, Tensor2<2>,
		std::conditional_t<Rank == 3, Tensor3<2>, std::conditional_t<Rank == 4, Tensor4<2>, Tensor5<2>>>>;

	template<size_t Rank>
	std::array<int, Rank> Indices(int flat)
	{
		std::array<int, Rank> indices{};
		for (size_t slot = 0; slot < Rank; ++slot)
			indices[slot] = (flat >> slot) & 1;
		return indices;
	}

	template<class Tensor, size_t Rank, size_t... Slots>
	decltype(auto) Component(Tensor& tensor, const std::array<int, Rank>& indices, std::index_sequence<Slots...>)
	{
		return tensor(indices[Slots]...);
	}

	template<class Tensor, size_t Rank>
	decltype(auto) Component(Tensor& tensor, const std::array<int, Rank>& indices)
	{
		if constexpr (Rank == 1)
			return tensor[indices[0]];
		else
			return Component(tensor, indices, std::make_index_sequence<Rank>{});
	}

	template<size_t Rank>
	Real OracleComponent(const std::array<int, Rank>& indices)
	{
		Real value = 1, weight = 1;
		for (int index : indices) {
			value += weight * index;
			weight *= 3;
		}
		return value;
	}

	template<int Rank>
	void CheckMetadataBoundaries()
	{
		INFO("rank " << Rank);
		using Tensor = TensorOfRank<Rank>;
		const int maximum = std::numeric_limits<int>::max(), minimum = std::numeric_limits<int>::min();
		for (const auto& counts : {std::array<int, 2>{-1, Rank + 1}, {Rank + 1, -1}, {maximum, maximum},
			{minimum, minimum}, {0, Rank + 1}, {Rank, 1}}) {
			CHECK_THROWS_AS(Tensor(counts[0], counts[1]), TensorCovarContravarNumError);
			if constexpr (Rank == 2 || Rank == 3)
				CHECK_THROWS_AS(Tensor(counts[0], counts[1], {1}), TensorCovarContravarNumError);
		}
		for (int covariant = 0; covariant <= Rank; ++covariant) {
			Tensor tensor(covariant, Rank - covariant);
			CHECK(tensor.numCovar() == covariant);
			CHECK(tensor.numContravar() == Rank - covariant);
			if constexpr (Rank == 2 || Rank == 3) {
				const Tensor initialized(covariant, Rank - covariant, {1});
				CHECK(initialized.numCovar() == covariant);
				CHECK(initialized.numContravar() == Rank - covariant);
				for (int slot = 0; slot < Rank; ++slot)
					CHECK(initialized.isContravar(slot) == (slot >= covariant));
			}
			for (int slot = 0; slot < Rank; ++slot) {
				const bool original = slot >= covariant;
				CHECK(tensor.isContravar(slot) == original);
				tensor.setContravar(slot, !original);
				CHECK(tensor.numContravar() == Rank - covariant + (original ? -1 : 1));
				tensor.setContravar(slot, original);
			}
			const auto original = tensor;
			for (int invalid : {-1, Rank, maximum, minimum})
				CHECK_THROWS_AS(tensor.setContravar(invalid, true), TensorIndexError);
			CHECK(tensor == original);
		}
	}

	template<int Rank>
	void CheckContraction()
	{
		using Tensor = TensorOfRank<Rank>;
		using Interface = std::conditional_t<Rank == 3, ITensor3<2>,
			std::conditional_t<Rank == 4, ITensor4<2>, ITensor5<2>>>;
		for (int mask = 0; mask < (1 << Rank); ++mask) {
			INFO("rank " << Rank << ", variance mask " << mask);
			Tensor tensor(Rank, 0), transformed(Rank, 0);
			for (int slot = 0; slot < Rank; ++slot) {
				tensor.setContravar(slot, (mask & (1 << slot)) != 0);
				transformed.setContravar(slot, tensor.isContravar(slot));
			}
			for (int flat = 0; flat < (1 << Rank); ++flat) {
				const auto indices = Indices<Rank>(flat);
				Real factor = 1;
				for (int slot = 0; slot < Rank; ++slot) {
					const Real scale = indices[slot] == 0 ? Real(2) : Real(3);
					factor *= tensor.isContravar(slot) ? scale : Real(1) / scale;
				}
				Component(tensor, indices) = OracleComponent(indices);
				Component(transformed, indices) = OracleComponent(indices) * factor;
			}
			const Interface& interface = tensor;
			for (int first = 0; first < Rank; ++first)
				for (int second = first + 1; second < Rank; ++second) {
					const auto diagonal = tensor.EuclideanDiagonalSum(first, second);
					const auto reverse = tensor.EuclideanDiagonalSum(second, first);
					std::array<int, Rank - 2> remaining{};
					int outputSlot = 0;
					for (int slot = 0; slot < Rank; ++slot)
						if (slot != first && slot != second)
							remaining[outputSlot++] = slot;
					if constexpr (Rank > 3)
						for (int slot = 0; slot < Rank - 2; ++slot)
							CHECK(diagonal.isContravar(slot) == tensor.isContravar(remaining[slot]));
					for (int flat = 0; flat < (1 << (Rank - 2)); ++flat) {
						const auto outputIndices = Indices<Rank - 2>(flat);
						std::array<int, Rank> sourceIndices{};
						for (int slot = 0; slot < Rank - 2; ++slot)
							sourceIndices[remaining[slot]] = outputIndices[slot];
						Real expected = 0;
						for (int summed = 0; summed < 2; ++summed) {
							sourceIndices[first] = sourceIndices[second] = summed;
							expected += OracleComponent(sourceIndices);
						}
						CHECK(Component(diagonal, outputIndices) == expected);
						CHECK(Component(reverse, outputIndices) == expected);
					}
					if (tensor.isContravar(first) == tensor.isContravar(second)) {
						CHECK_THROWS_AS(interface.contract(first, second), TensorCovarContravarArithmeticError);
						CHECK_THROWS_AS(tensor.contract(second, first), TensorCovarContravarArithmeticError);
						continue;
					}
					const auto strict = interface.contract(first, second);
					const auto transformedStrict = transformed.contract(first, second);
					CHECK(strict == diagonal);
					CHECK(tensor.contract(second, first) == strict);
					for (int flat = 0; flat < (1 << (Rank - 2)); ++flat) {
						const auto indices = Indices<Rank - 2>(flat);
						Real factor = 1;
						for (int slot = 0; slot < Rank - 2; ++slot) {
							const Real scale = indices[slot] == 0 ? Real(2) : Real(3);
							factor *= tensor.isContravar(remaining[slot]) ? scale : Real(1) / scale;
						}
						CHECK_THAT(Component(transformedStrict, indices),
							RealWithinRelNoScale(Component(strict, indices) * factor, Precision::EqualityTolerance));
					}
				}
		}
		Tensor tensor(1, Rank - 1);
		for (int invalid : {-1, Rank, std::numeric_limits<int>::max(), std::numeric_limits<int>::min()}) {
			CHECK_THROWS_AS(tensor.contract(invalid, 0), TensorIndexError);
			CHECK_THROWS_AS(tensor.contract(0, invalid), TensorIndexError);
			CHECK_THROWS_AS(tensor.EuclideanDiagonalSum(invalid, 0), TensorIndexError);
			CHECK_THROWS_AS(tensor.EuclideanDiagonalSum(0, invalid), TensorIndexError);
		}
		for (int slot = 0; slot < Rank; ++slot) {
			CHECK_THROWS_AS(tensor.contract(slot, slot), TensorIndexError);
			CHECK_THROWS_AS(tensor.EuclideanDiagonalSum(slot, slot), TensorIndexError);
		}
	}
}

TEST_CASE("Tensor variance constructors and setters validate integer boundaries", "[TensorContract]")
{
	CheckMetadataBoundaries<2>();
	CheckMetadataBoundaries<3>();
	CheckMetadataBoundaries<4>();
	CheckMetadataBoundaries<5>();
}

TEST_CASE("Rank-three contraction obeys variance and basis-change laws", "[TensorContract]")
{
	CheckContraction<3>();
}

TEST_CASE("Rank-four contraction obeys variance and basis-change laws", "[TensorContract]")
{
	CheckContraction<4>();
}

TEST_CASE("Rank-five contraction obeys variance and basis-change laws", "[TensorContract]")
{
	CheckContraction<5>();
}

TEST_CASE("Tensor index raising and lowering preserve the exact transformed slot", "[TensorContract][TensorRaiseLower]")
{
	const Tensor2<2> metric(2, 0, {2, 1, 1, 3});
	const Tensor2<2> inverse(0, 2, {Real(3) / 5, -Real(1) / 5, -Real(1) / 5, Real(2) / 5});
	for (int mask = 0; mask < 4; ++mask) {
		Tensor2<2> tensor(2, 0, {1, 2, 3, 5});
		for (int slot = 0; slot < 2; ++slot)
			tensor.setContravar(slot, (mask & (1 << slot)) != 0);
		for (int slot = 0; slot < 2; ++slot) {
			const bool lowering = tensor.isContravar(slot);
			const auto changed = lowering ?
				(slot == 0 ? tensor.LowerFirstIndex(metric) : tensor.LowerSecondIndex(metric)) :
				(slot == 0 ? tensor.RaiseFirstIndex(inverse) : tensor.RaiseSecondIndex(inverse));
			for (int index = 0; index < 2; ++index)
				CHECK(changed.isContravar(index) == (index == slot ? !lowering : tensor.isContravar(index)));
			CHECK(changed.numContravar() == tensor.numContravar() + (lowering ? -1 : 1));
			const auto& transformation = lowering ? metric : inverse;
			for (int i = 0; i < 2; ++i)
				for (int j = 0; j < 2; ++j) {
					Real expected = 0, magnitude = 0;
					for (int k = 0; k < 2; ++k) {
						const Real term = slot == 0 ? transformation(i, k) * tensor(k, j) :
							tensor(i, k) * transformation(k, j);
						expected += term;
						magnitude += std::abs(term);
					}
					// Cancellation is bounded by the terms, not the nearly zero result.
					CHECK_THAT(changed(i, j), RealWithinAbsNoScale(expected, 4 * Precision::Eps * magnitude));
				}
			const auto restored = lowering ?
				(slot == 0 ? changed.RaiseFirstIndex(inverse) : changed.RaiseSecondIndex(inverse)) :
				(slot == 0 ? changed.LowerFirstIndex(metric) : changed.LowerSecondIndex(metric));
			for (int index = 0; index < 2; ++index)
				CHECK(restored.isContravar(index) == tensor.isContravar(index));
			for (int i = 0; i < 2; ++i)
				for (int j = 0; j < 2; ++j)
					CHECK_THAT(restored(i, j), RealWithinRelNoScale(tensor(i, j), Precision::EqualityTolerance));
		}
	}
}
