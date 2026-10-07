///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Permutation.h                                                       ///
///  Description: Fixed-size and runtime permutation value types                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERMUTATION_H
#define MML_PERMUTATION_H

#include <mml/MMLExceptions.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <numeric>
#include <utility>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		enum class PermutationParity
		{
			Even,
			Odd
		};

		namespace Detail
		{
			template<class Images>
			void ValidatePermutationImages(const Images& images)
			{
				std::vector<bool> seen(images.size(), false);
				for (std::size_t index = 0; index < images.size(); index++) {
					const int image = images[index];
					if (image < 0 || static_cast<std::size_t>(image) >= images.size())
						throw ArgumentError("Permutation image is out of range");
					if (seen[static_cast<std::size_t>(image)])
						throw ArgumentError("Permutation contains a duplicate image");
					seen[static_cast<std::size_t>(image)] = true;
				}
			}

			template<class Images>
			std::vector<std::vector<int>> PermutationCycles(const Images& images, bool includeFixedPoints)
			{
				std::vector<std::vector<int>> cycles;
				std::vector<bool> visited(images.size(), false);
				for (std::size_t start = 0; start < images.size(); start++) {
					if (visited[start])
						continue;

					std::vector<int> cycle;
					int current = static_cast<int>(start);
					do {
						cycle.push_back(current);
						visited[static_cast<std::size_t>(current)] = true;
						current = images[static_cast<std::size_t>(current)];
					} while (current != static_cast<int>(start));

					if (includeFixedPoints || cycle.size() > 1)
						cycles.push_back(std::move(cycle));
				}
				return cycles;
			}

			template<class Images>
			int PermutationTranspositionCount(const Images& images)
			{
				int count = 0;
				for (const auto& cycle : PermutationCycles(images, true))
					count += static_cast<int>(cycle.size()) - 1;
				return count;
			}

			template<class Images>
			std::size_t PermutationOrder(const Images& images)
			{
				std::size_t result = 1;
				for (const auto& cycle : PermutationCycles(images, false))
					result = std::lcm(result, cycle.size());
				return result;
			}

			template<class SetImage>
			void ApplyCycles(std::size_t degree,
				std::initializer_list<std::initializer_list<int>> cycles, SetImage setImage)
			{
				std::vector<bool> used(degree, false);
				for (const auto& cycle : cycles) {
					if (cycle.size() == 0)
						continue;

					for (int index : cycle) {
						if (index < 0 || static_cast<std::size_t>(index) >= degree)
							throw ArgumentError("Permutation cycle index is out of range");
						if (used[static_cast<std::size_t>(index)])
							throw ArgumentError("Permutation cycle notation repeats an index");
						used[static_cast<std::size_t>(index)] = true;
					}

					auto current = cycle.begin();
					auto next = current;
					++next;
					for (; next != cycle.end(); ++current, ++next)
						setImage(*current, *next);
					setImage(*current, *cycle.begin());
				}
			}
		}

		/// @brief A permutation of the indices [0, N).
		/// @details Images use the convention images[i] = p(i). Composition follows
		///          compose(left, right)(i) = left(right(i)).
		template<std::size_t N>
		class Permutation
		{
			static_assert(N > 0, "Permutation degree must be positive");
			std::array<int, N> _images{};

			explicit Permutation(const std::array<int, N>& images, bool validate)
				: _images(images)
			{
				if (validate)
					Detail::ValidatePermutationImages(_images);
			}

		public:
			using element_type = int;
			using storage_type = std::array<int, N>;
			static constexpr std::size_t Degree = N;

			Permutation()
			{
				for (std::size_t index = 0; index < N; index++)
					_images[index] = static_cast<int>(index);
			}

			explicit Permutation(const std::array<int, N>& images)
				: Permutation(images, true)
			{
			}

			Permutation(std::initializer_list<int> images)
			{
				if (images.size() != N)
					throw ArgumentError("Permutation initializer size must equal its degree");
				std::copy(images.begin(), images.end(), _images.begin());
				Detail::ValidatePermutationImages(_images);
			}

			static Permutation identity() { return Permutation(); }

			static Permutation from_cycles(std::initializer_list<std::initializer_list<int>> cycles)
			{
				Permutation result;
				Detail::ApplyCycles(N, cycles,
					[&result](int source, int destination) { result._images[static_cast<std::size_t>(source)] = destination; });
				return result;
			}

			int apply(int index) const
			{
				if (index < 0 || static_cast<std::size_t>(index) >= N)
					throw IndexError("Permutation index is out of range", index, static_cast<int>(N));
				return _images[static_cast<std::size_t>(index)];
			}

			template<class Value>
			std::array<Value, N> apply(const std::array<Value, N>& values) const
			{
				std::array<Value, N> result{};
				for (std::size_t source = 0; source < N; source++)
					result[static_cast<std::size_t>(_images[source])] = values[source];
				return result;
			}

			Permutation compose(const Permutation& right) const
			{
				std::array<int, N> result{};
				for (std::size_t index = 0; index < N; index++)
					result[index] = _images[static_cast<std::size_t>(right._images[index])];
				return Permutation(result, false);
			}

			Permutation inverse() const
			{
				std::array<int, N> result{};
				for (std::size_t index = 0; index < N; index++)
					result[static_cast<std::size_t>(_images[index])] = static_cast<int>(index);
				return Permutation(result, false);
			}

			std::vector<std::vector<int>> cycles(bool includeFixedPoints = false) const
			{
				return Detail::PermutationCycles(_images, includeFixedPoints);
			}

			int transposition_count() const { return Detail::PermutationTranspositionCount(_images); }
			PermutationParity parity() const
			{
				return transposition_count() % 2 == 0 ? PermutationParity::Even : PermutationParity::Odd;
			}
			int sign() const { return parity() == PermutationParity::Even ? 1 : -1; }
			std::size_t order() const { return Detail::PermutationOrder(_images); }

			const storage_type& images() const noexcept { return _images; }
			constexpr std::size_t size() const noexcept { return N; }

			bool operator==(const Permutation& other) const noexcept { return _images == other._images; }
			bool operator!=(const Permutation& other) const noexcept { return !(*this == other); }
		};

		template<std::size_t N>
		Permutation<N> compose(const Permutation<N>& left, const Permutation<N>& right)
		{
			return left.compose(right);
		}

		/// @brief A runtime-degree permutation with the same conventions as Permutation<N>.
		class DynamicPermutation
		{
			std::vector<int> _images;

			explicit DynamicPermutation(std::vector<int> images, bool validate)
				: _images(std::move(images))
			{
				if (validate)
					Detail::ValidatePermutationImages(_images);
			}

		public:
			using element_type = int;
			using storage_type = std::vector<int>;

			DynamicPermutation() = default;

			explicit DynamicPermutation(std::size_t degree)
				: _images(degree)
			{
				std::iota(_images.begin(), _images.end(), 0);
			}

			explicit DynamicPermutation(std::vector<int> images)
				: DynamicPermutation(std::move(images), true)
			{
			}

			DynamicPermutation(std::initializer_list<int> images)
				: DynamicPermutation(std::vector<int>(images), true)
			{
			}

			static DynamicPermutation identity(std::size_t degree) { return DynamicPermutation(degree); }

			static DynamicPermutation from_cycles(std::size_t degree,
				std::initializer_list<std::initializer_list<int>> cycles)
			{
				DynamicPermutation result(degree);
				Detail::ApplyCycles(degree, cycles,
					[&result](int source, int destination) { result._images[static_cast<std::size_t>(source)] = destination; });
				return result;
			}

			int apply(int index) const
			{
				if (index < 0 || static_cast<std::size_t>(index) >= _images.size())
					throw IndexError("DynamicPermutation index is out of range", index, static_cast<int>(_images.size()));
				return _images[static_cast<std::size_t>(index)];
			}

			template<class Value>
			std::vector<Value> apply(const std::vector<Value>& values) const
			{
				if (values.size() != _images.size())
					throw ArgumentError("DynamicPermutation value count must equal its degree");
				std::vector<Value> result(values.size());
				for (std::size_t source = 0; source < values.size(); source++)
					result[static_cast<std::size_t>(_images[source])] = values[source];
				return result;
			}

			DynamicPermutation compose(const DynamicPermutation& right) const
			{
				if (size() != right.size())
					throw ArgumentError("Cannot compose permutations of different degrees");
				std::vector<int> result(size());
				for (std::size_t index = 0; index < size(); index++)
					result[index] = _images[static_cast<std::size_t>(right._images[index])];
				return DynamicPermutation(std::move(result), false);
			}

			DynamicPermutation inverse() const
			{
				std::vector<int> result(size());
				for (std::size_t index = 0; index < size(); index++)
					result[static_cast<std::size_t>(_images[index])] = static_cast<int>(index);
				return DynamicPermutation(std::move(result), false);
			}

			std::vector<std::vector<int>> cycles(bool includeFixedPoints = false) const
			{
				return Detail::PermutationCycles(_images, includeFixedPoints);
			}

			int transposition_count() const { return Detail::PermutationTranspositionCount(_images); }
			PermutationParity parity() const
			{
				return transposition_count() % 2 == 0 ? PermutationParity::Even : PermutationParity::Odd;
			}
			int sign() const { return parity() == PermutationParity::Even ? 1 : -1; }
			std::size_t order() const { return Detail::PermutationOrder(_images); }

			const storage_type& images() const noexcept { return _images; }
			std::size_t size() const noexcept { return _images.size(); }

			bool operator==(const DynamicPermutation& other) const noexcept { return _images == other._images; }
			bool operator!=(const DynamicPermutation& other) const noexcept { return !(*this == other); }
		};

		inline DynamicPermutation compose(const DynamicPermutation& left, const DynamicPermutation& right)
		{
			return left.compose(right);
		}
	}
}

#endif // MML_PERMUTATION_H