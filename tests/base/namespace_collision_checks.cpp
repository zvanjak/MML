struct Real { int applicationValue; };
struct Complex { int applicationValue; };
inline int GetEnv(const char*) { return 7; }

#ifdef MML_NAMESPACE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLTypeDefs.h>
#include <mml/MMLBase.h>
#include <mml/tools/visualizers/VisualizerConfig.h>
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>
#endif

#include <type_traits>

#ifdef REAL
#error MML must not define a global REAL macro.
#endif

static_assert(std::is_same_v<MML::Complex, std::complex<MML::Real>>);
static_assert(!std::is_same_v<::Real, MML::Real>);
static_assert(!std::is_same_v<::Complex, MML::Complex>);

void MMLNamespaceCollisionCompileCheck()
{
	::Real application{1};
	::Complex applicationComplex{2};
	MML::Vector<MML::Real> values{MML::Real(1), MML::Real(2)};
	MML::Matrix<MML::Real> matrix(1, 1, {MML::Real(2)});
	const auto inverse = MML::MatrixAlg::Inverse(matrix);
	(void)application;
	(void)applicationComplex;
	(void)values;
	(void)inverse;
	(void)::GetEnv("application-owned");
}
