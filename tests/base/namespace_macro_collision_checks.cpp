constexpr int ApplicationRealLiteral(int value) { return value + 1; }
#define REAL(value) ApplicationRealLiteral(value)

#include <mml/MMLBase.h>
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>

static_assert(REAL(3) == 4);
static_assert(std::is_same_v<MML::Complex, std::complex<MML::Real>>);

#undef REAL
