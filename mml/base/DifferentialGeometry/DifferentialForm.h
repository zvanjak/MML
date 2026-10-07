///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/DifferentialForm.h                             ///
///  Description: Typed differential forms and exterior algebra                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DIFFERENTIAL_GEOMETRY_DIFFERENTIAL_FORM_H
#define MML_DIFFERENTIAL_GEOMETRY_DIFFERENTIAL_FORM_H

#include <mml/MMLBase.h>

#include <mml/base/Tensor/Rank1Tensor.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <type_traits>

namespace MML
{
	namespace Detail
	{
		constexpr int IntPow(int base, int exponent)
		{
			return exponent == 0 ? 1 : base * IntPow(base, exponent - 1);
		}

		constexpr int Factorial(int n)
		{
			return n <= 1 ? 1 : n * Factorial(n - 1);
		}

		template<int K>
		bool HasRepeatedIndex(const std::array<int, K == 0 ? 1 : K>& indices)
		{
			for (int i = 0; i < K; i++)
				for (int j = i + 1; j < K; j++)
					if (indices[i] == indices[j])
						return true;
			return false;
		}

		template<int K>
		int SortPermutationSign(std::array<int, K == 0 ? 1 : K>& indices)
		{
			int sign = 1;
			for (int i = 0; i < K; i++) {
				for (int j = i + 1; j < K; j++) {
					if (indices[j] < indices[i]) {
						std::swap(indices[i], indices[j]);
						sign = -sign;
					}
				}
			}
			return sign;
		}

		template<class Scalar>
		bool IsNearlyZeroForAlternating(Scalar value, Real tolerance)
		{
			return std::abs(value) <= tolerance;
		}
	}

	template<class Scalar, int N, int K, class FrameTag>
	class DifferentialForm
	{
		static_assert(N > 0, "DifferentialForm dimension must be positive");
		static_assert(K >= 0, "DifferentialForm degree cannot be negative");
		static_assert(K <= N, "DifferentialForm degree cannot exceed dimension");

		static constexpr int IndexStorageSize = K == 0 ? 1 : K;
		std::array<Scalar, Detail::IntPow(N, K)> _components{};

		static int FlattenIndex(const std::array<int, IndexStorageSize>& indices)
		{
			int flat = 0;
			for (int i = 0; i < K; i++)
				flat = flat * N + indices[i];
			return flat;
		}

		static std::array<int, IndexStorageSize> UnflattenIndex(int flat)
		{
			std::array<int, IndexStorageSize> indices{};
			for (int i = K - 1; i >= 0; i--) {
				indices[i] = flat % N;
				flat /= N;
			}
			return indices;
		}

	public:
		using value_type = Scalar;
		using frame_type = FrameTag;
		static constexpr int Dimension = N;
		static constexpr int Degree = K;
		static constexpr int ComponentCount = Detail::IntPow(N, K);

		DifferentialForm() = default;
		DifferentialForm(std::initializer_list<Scalar> values)
		{
			auto value = values.begin();
			for (int i = 0; i < ComponentCount; i++)
				_components[i] = value != values.end() ? *value++ : Scalar{};
		}

		Scalar& operator[](int flatIndex) noexcept { return _components[flatIndex]; }
		const Scalar& operator[](int flatIndex) const noexcept { return _components[flatIndex]; }

		template<class... Indices>
		Scalar& Component(Indices... indices)
		{
			static_assert(sizeof...(Indices) == K, "Component index count must match form degree");
			std::array<int, IndexStorageSize> multiIndex{ { static_cast<int>(indices)... } };
			return _components[FlattenIndex(multiIndex)];
		}

		template<class... Indices>
		const Scalar& Component(Indices... indices) const
		{
			static_assert(sizeof...(Indices) == K, "Component index count must match form degree");
			std::array<int, IndexStorageSize> multiIndex{ { static_cast<int>(indices)... } };
			return _components[FlattenIndex(multiIndex)];
		}

		Scalar& ComponentAt(const std::array<int, IndexStorageSize>& indices)
		{
			return _components[FlattenIndex(indices)];
		}

		const Scalar& ComponentAt(const std::array<int, IndexStorageSize>& indices) const
		{
			return _components[FlattenIndex(indices)];
		}

		bool IsAlternating(Real tolerance = PrecisionValues<Real>::DefaultTolerance) const
		{
			if constexpr (K <= 1) {
				return true;
			}
			else {
				for (int flat = 0; flat < ComponentCount; flat++) {
					auto indices = UnflattenIndex(flat);
					if (Detail::HasRepeatedIndex<K>(indices)) {
						if (!Detail::IsNearlyZeroForAlternating(_components[flat], tolerance))
							return false;
						continue;
					}

					auto sortedIndices = indices;
					int sign = Detail::SortPermutationSign<K>(sortedIndices);
					Scalar expected = Scalar(sign) * ComponentAt(sortedIndices);
					if (!Detail::IsNearlyZeroForAlternating(_components[flat] - expected, tolerance))
						return false;
				}
				return true;
			}
		}

		DifferentialForm AlternatingPart() const
		{
			if constexpr (K <= 1) {
				return *this;
			}
			else {
				DifferentialForm result;
				std::array<int, IndexStorageSize> permutation{};
				for (int flat = 0; flat < ComponentCount; flat++) {
					auto indices = UnflattenIndex(flat);
					if (Detail::HasRepeatedIndex<K>(indices)) {
						result.ComponentAt(indices) = Scalar{};
						continue;
					}

					for (int i = 0; i < K; i++)
						permutation[i] = i;

					Scalar sum{};
					do {
						std::array<int, IndexStorageSize> permuted{};
						for (int i = 0; i < K; i++)
							permuted[i] = indices[permutation[i]];

						int inversions = 0;
						for (int i = 0; i < K; i++)
							for (int j = i + 1; j < K; j++)
								if (permutation[i] > permutation[j])
									inversions++;

						Scalar sign = inversions % 2 == 0 ? Scalar{ 1 } : Scalar{ -1 };
						sum += sign * ComponentAt(permuted);
					} while (std::next_permutation(permutation.begin(), permutation.begin() + K));

					result.ComponentAt(indices) = sum / Scalar(Detail::Factorial(K));
				}
				return result;
			}
		}

		template<class... Indices>
		void SetAlternatingComponent(Scalar value, Indices... indices)
		{
			static_assert(sizeof...(Indices) == K, "Component index count must match form degree");
			std::array<int, IndexStorageSize> input{ { static_cast<int>(indices)... } };
			SetAlternatingComponentAt(value, input);
		}

		void SetAlternatingComponentAt(Scalar value, const std::array<int, IndexStorageSize>& indices)
		{
			if constexpr (K == 0) {
				_components[0] = value;
			}
			else if constexpr (K == 1) {
				ComponentAt(indices) = value;
			}
			else {
				if (Detail::HasRepeatedIndex<K>(indices)) {
					ComponentAt(indices) = Scalar{};
					return;
				}

				std::array<int, IndexStorageSize> permutation{};
				for (int i = 0; i < K; i++)
					permutation[i] = i;

				do {
					std::array<int, IndexStorageSize> permuted{};
					for (int i = 0; i < K; i++)
						permuted[i] = indices[permutation[i]];

					int inversions = 0;
					for (int i = 0; i < K; i++)
						for (int j = i + 1; j < K; j++)
							if (permutation[i] > permutation[j])
								inversions++;

					Scalar sign = inversions % 2 == 0 ? Scalar{ 1 } : Scalar{ -1 };
					ComponentAt(permuted) = sign * value;
				} while (std::next_permutation(permutation.begin(), permutation.begin() + K));
			}
		}

		template<class... Vectors>
		Scalar operator()(const Vectors&... vectors) const
			requires (sizeof...(Vectors) == K && (std::same_as<std::remove_cvref_t<Vectors>, TangentVector<N, FrameTag>> && ...))
		{
			std::array<const TangentVector<N, FrameTag>*, K> vectorRefs{ { &vectors... } };
			Scalar result{};

			for (int flat = 0; flat < ComponentCount; flat++) {
				auto indices = UnflattenIndex(flat);
				Scalar term = _components[flat];
				for (int i = 0; i < K; i++)
					term *= (*vectorRefs[i])[indices[i]];
				result += term;
			}

			return result;
		}
	};

	template<int N, class FrameTag>
	using ScalarForm = DifferentialForm<Real, N, 0, FrameTag>;

	template<int N, class FrameTag>
	using Form1 = DifferentialForm<Real, N, 1, FrameTag>;

	template<int N, class FrameTag>
	using OneForm = Form1<N, FrameTag>;

	template<int N, class FrameTag>
	using Form2 = DifferentialForm<Real, N, 2, FrameTag>;

	template<int N, class FrameTag>
	using Form3 = DifferentialForm<Real, N, 3, FrameTag>;

	template<int N, class FrameTag>
	using VolumeForm = DifferentialForm<Real, N, N, FrameTag>;

	template<int Index, int N, class FrameTag>
	Form1<N, FrameTag> BasisOneForm()
	{
		static_assert(Index >= 0 && Index < N, "Basis one-form index must be in range");
		Form1<N, FrameTag> result;
		result.Component(Index) = REAL(1.0);
		return result;
	}

	template<int N, class FrameTag>
	Form1<N, FrameTag> ToForm(const Covector<N, FrameTag>& alpha)
	{
		Form1<N, FrameTag> result;
		for (int i = 0; i < N; i++)
			result.Component(i) = alpha[i];
		return result;
	}

	template<int N, class FrameTag>
	Covector<N, FrameTag> ToCovector(const Form1<N, FrameTag>& form)
	{
		Covector<N, FrameTag> result;
		for (int i = 0; i < N; i++)
			result[i] = form.Component(i);
		return result;
	}

	template<class Scalar, int N, int P, int Q, class FrameTag>
	DifferentialForm<Scalar, N, P + Q, FrameTag> Wedge(const DifferentialForm<Scalar, N, P, FrameTag>& a,
		const DifferentialForm<Scalar, N, Q, FrameTag>& b)
		requires (P + Q <= N)
	{
		const int componentCount = Detail::IntPow(N, P + Q);
		constexpr int R = P + Q;
		constexpr int IndexStorageSize = R == 0 ? 1 : R;
		DifferentialForm<Scalar, N, R, FrameTag> result;
		std::array<int, IndexStorageSize> indices{};
		std::array<int, IndexStorageSize> permutation{};

		for (int flat = 0; flat < componentCount; flat++) {
			int remaining = flat;
			for (int i = R - 1; i >= 0; i--) {
				indices[i] = remaining % N;
				remaining /= N;
			}
			for (int i = 0; i < R; i++)
				permutation[i] = i;

			Scalar component{};
			do {
				std::array<int, P == 0 ? 1 : P> aIndices{};
				std::array<int, Q == 0 ? 1 : Q> bIndices{};
				for (int i = 0; i < P; i++)
					aIndices[i] = indices[permutation[i]];
				for (int i = 0; i < Q; i++)
					bIndices[i] = indices[permutation[P + i]];

				int inversions = 0;
				for (int i = 0; i < R; i++)
					for (int j = i + 1; j < R; j++)
						if (permutation[i] > permutation[j])
							inversions++;

				Scalar sign = inversions % 2 == 0 ? Scalar{ 1 } : Scalar{ -1 };
				component += sign * a.ComponentAt(aIndices) * b.ComponentAt(bIndices);
			} while (std::next_permutation(permutation.begin(), permutation.begin() + R));

			result[flat] = component / Scalar(Detail::Factorial(P) * Detail::Factorial(Q));
		}

		return result;
	}

	template<class Scalar, int N, int P, int Q, class FrameTag>
	DifferentialForm<Scalar, N, P + Q, FrameTag> wedge(const DifferentialForm<Scalar, N, P, FrameTag>& a,
		const DifferentialForm<Scalar, N, Q, FrameTag>& b)
		requires (P + Q <= N)
	{
		return Wedge(a, b);
	}
}

#endif // MML_DIFFERENTIAL_GEOMETRY_DIFFERENTIAL_FORM_H