#include <MML.h>

namespace MML::Algebra {
	class CyclicGroup;
	class DihedralGroup;
}
template<class Type>
concept CompleteOptionalSingleHeaderAlgebraType = requires { sizeof(Type); };
static_assert(!CompleteOptionalSingleHeaderAlgebraType<MML::Algebra::CyclicGroup>);
static_assert(!CompleteOptionalSingleHeaderAlgebraType<MML::Algebra::DihedralGroup>);
static_assert(sizeof(MML::Algebra::SO3) > 0);
static_assert(sizeof(MML::Algebra::SE3) > 0);
