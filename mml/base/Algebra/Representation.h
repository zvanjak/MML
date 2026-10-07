///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Representation.h                                                    ///
///  Description: Finite-dimensional matrix representation value adapter              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_REPRESENTATION_H
#define MML_REPRESENTATION_H

#include <mml/MMLExceptions.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Vector/VectorN.h>

#include <functional>
#include <utility>

namespace MML
{
	namespace Algebra
	{
		/// @brief Matrix-valued map rho: G -> GL(N, Scalar).
		template<class GroupElement, class Scalar, int N>
		class Representation
		{
			static_assert(N > 0, "Representation dimension must be positive");

		public:
			using group_element_type = GroupElement;
			using scalar_type = Scalar;
			using matrix_type = MatrixNM<Scalar, N, N>;
			using vector_type = VectorN<Scalar, N>;
			using evaluate_type = std::function<matrix_type(const GroupElement&)>;
			static constexpr int dimension = N;

		private:
			evaluate_type _evaluate;

		public:
			explicit Representation(evaluate_type evaluate)
				: _evaluate(std::move(evaluate))
			{
				if (!_evaluate)
					throw ArgumentError("Representation evaluator must be callable");
			}

			matrix_type matrix(const GroupElement& element) const { return _evaluate(element); }
			matrix_type operator()(const GroupElement& element) const { return matrix(element); }
			vector_type apply(const GroupElement& element, const vector_type& vector) const
			{
				return matrix(element) * vector;
			}
		};

		template<class GroupElement, class Scalar, int N, class Evaluate>
		Representation<GroupElement, Scalar, N> MakeRepresentation(Evaluate evaluate)
		{
			return Representation<GroupElement, Scalar, N>(std::move(evaluate));
		}
	}
}

#endif // MML_REPRESENTATION_H