#include <mml/core/Curves.h>
#include <mml/core/Surfaces.h>

#if defined(MML_PARAMETRIC_CURVES_DEFS_H) || defined(MML_PARAMETRIC_SURFACES_DEFS_H)
#error Core geometry must not depend on development fixture definitions.
#endif

namespace MML::Curves {
	class LemniscateCurve;
	class DeltoidCurve;
	class AstroidCurve;
	class EpitrochoidCurve;
	class ButterflyCurve;
	class TwistedCubicCurve;
	class ToroidalSpiralCurve;
}
namespace MML::Surfaces {
	class MonkeySaddle;
	class MobiusStrip;
	class Helicoid;
	class EnneperSurface;
	class KleinBottle;
	class DiniSurface;
}

template<class Type>
concept CompleteCoreShowcaseType = requires { sizeof(Type); };

static_assert(!CompleteCoreShowcaseType<MML::Curves::LemniscateCurve>);
static_assert(!CompleteCoreShowcaseType<MML::Curves::DeltoidCurve>);
static_assert(!CompleteCoreShowcaseType<MML::Curves::AstroidCurve>);
static_assert(!CompleteCoreShowcaseType<MML::Curves::EpitrochoidCurve>);
static_assert(!CompleteCoreShowcaseType<MML::Curves::ButterflyCurve>);
static_assert(!CompleteCoreShowcaseType<MML::Curves::TwistedCubicCurve>);
static_assert(!CompleteCoreShowcaseType<MML::Curves::ToroidalSpiralCurve>);
static_assert(!CompleteCoreShowcaseType<MML::Surfaces::MonkeySaddle>);
static_assert(!CompleteCoreShowcaseType<MML::Surfaces::MobiusStrip>);
static_assert(!CompleteCoreShowcaseType<MML::Surfaces::Helicoid>);
static_assert(!CompleteCoreShowcaseType<MML::Surfaces::EnneperSurface>);
static_assert(!CompleteCoreShowcaseType<MML::Surfaces::KleinBottle>);
static_assert(!CompleteCoreShowcaseType<MML::Surfaces::DiniSurface>);
