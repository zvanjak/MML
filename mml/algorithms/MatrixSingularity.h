#if !defined MML_MATRIX_SINGULARITY_H
#define MML_MATRIX_SINGULARITY_H

#include <mml/MMLSingularityHandling.h>
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>

namespace MML::Singularity {
		/// Apply policy only to rank failure; invalid data and overflow still throw.
		template<class TMatrix>
		TMatrix InverseMatrix(const TMatrix& matrix, SingularityPolicy policy = DEFAULT_POLICY) {
			RequireMatrixPolicy(policy);
			RequireFiniteMatrix(matrix, "Matrix inversion");
			try {
				auto inverse = MatrixAlg::Inverse(matrix);
				RequireFiniteMatrix(inverse, "Matrix inversion result");
				return inverse;
			} catch (const SingularMatrixError&) {
				if (policy == SingularityPolicy::Throw)
					throw;
				return MatrixAtSingularity<TMatrix>(policy, "Matrix inversion: singular matrix");
			}
		}
}

#endif // MML_MATRIX_SINGULARITY_H
