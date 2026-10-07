///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLConcepts.h                                                       ///
///  Description: C++20 concepts and traits for MML scalar and shape constraints      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CONCEPTS_H
#define MML_CONCEPTS_H

#include <mml/MMLTypeDefs.h>

#include <complex>
#include <concepts>
#include <functional>
#include <type_traits>
#include <utility>

namespace MML
{
	template<typename T>
	struct is_std_complex : std::false_type { };

	template<typename T>
	struct is_std_complex<std::complex<T>> : std::true_type { };

	template<typename T>
	struct is_mml_complex_scalar : std::false_type { };

	template<std::floating_point T>
	struct is_mml_complex_scalar<std::complex<T>> : std::true_type { };

	template<typename T>
	struct is_simple_numeric : std::bool_constant<
		std::is_arithmetic_v<std::remove_cvref_t<T>> ||
		is_mml_complex_scalar<std::remove_cvref_t<T>>::value> { };

	template<typename T>
	inline constexpr bool is_MML_simple_numeric = is_simple_numeric<T>::value;

	template<typename T>
	inline constexpr bool is_MML_complex = is_mml_complex_scalar<std::remove_cvref_t<T>>::value;

	template<typename T>
	concept MMLArithmetic = std::is_arithmetic_v<std::remove_cvref_t<T>>;

	template<typename T>
	concept MMLReal = std::floating_point<std::remove_cvref_t<T>>;

	template<typename T>
	concept MMLComplex = is_MML_complex<T>;

	template<typename T>
	concept MMLNumeric = MMLArithmetic<T> || MMLComplex<T>;

	template<typename T>
	concept MMLScalar = MMLReal<T> || MMLComplex<T>;

	template<typename F>
	concept RealFunctionCallable = requires(F&& f, Real x) {
		{ std::invoke(std::forward<F>(f), x) } -> std::convertible_to<Real>;
	};

	template<typename T>
	concept Field = requires(T a, T b) {
		{ a + b } -> std::convertible_to<T>;
		{ a - b } -> std::convertible_to<T>;
		{ a * b } -> std::convertible_to<T>;
		{ a / b } -> std::convertible_to<T>;
		{ T{ 0 } } -> std::convertible_to<T>;
		{ T{ 1 } } -> std::convertible_to<T>;
	};

	template<typename V>
	concept VectorLike = requires(V v, int i) {
		typename std::remove_cvref_t<V>::value_type;
		{ v.size() } -> std::convertible_to<int>;
		{ v[i] } -> std::convertible_to<typename std::remove_cvref_t<V>::value_type>;
	};

	template<typename M>
	concept MatrixLike = requires(M matrix, int i, int j) {
		typename std::remove_cvref_t<M>::value_type;
		{ matrix.rows() } -> std::convertible_to<int>;
		{ matrix.cols() } -> std::convertible_to<int>;
		{ matrix(i, j) } -> std::convertible_to<typename std::remove_cvref_t<M>::value_type>;
	};
}

#endif // MML_CONCEPTS_H