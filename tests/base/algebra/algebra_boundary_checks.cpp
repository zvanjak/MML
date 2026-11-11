#include <mml/base/Algebra_base.h>
#include <mml/core/Fields/VectorFieldOperations.h>
#include <mml/core/VectorSpaces/Basis.h>

namespace MML::Algebra {
	class CyclicGroup;
	class DihedralGroup;
}
template<class Type>
concept CompleteOptionalAlgebraType = requires { sizeof(Type); };
static_assert(!CompleteOptionalAlgebraType<MML::Algebra::CyclicGroup>);
static_assert(!CompleteOptionalAlgebraType<MML::Algebra::DihedralGroup>);
static_assert(sizeof(MML::Algebra::SO3) > 0);
static_assert(sizeof(MML::Algebra::SE3) > 0);
