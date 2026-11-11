#include <mml/MMLTypeDefs.h>
#include <mml/MMLPrecision.h>
#include <type_traits>

#ifdef REAL
#error MML scalar configuration must not inject REAL.
#endif

#if defined(MML_USE_FLOAT)
static_assert(std::is_same_v<MML::Real, float>);
#elif defined(MML_USE_LONG_DOUBLE)
static_assert(std::is_same_v<MML::Real, long double>);
#else
static_assert(std::is_same_v<MML::Real, double>);
#endif
static_assert(std::is_same_v<MML::Complex, std::complex<MML::Real>>);
static_assert(MML::Precision::Eps == std::numeric_limits<MML::Real>::epsilon());
