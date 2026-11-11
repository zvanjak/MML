#include <mml/base/Geometry/Geometry2D.h>
#include <mml/base/Geometry/Geometry3DBodies.h>
#include <mml/algorithms/ComputationalGeometry.h>

static_assert(MML::Precision::GeometryEpsilon == MML::Precision::EqualityTolerance);
static_assert(MML::PrecisionValues<float>::GeometryEpsilon == MML::PrecisionValues<float>::EqualityTolerance);
static_assert(MML::PrecisionValues<double>::GeometryEpsilon == MML::PrecisionValues<double>::EqualityTolerance);
static_assert(MML::PrecisionValues<long double>::GeometryEpsilon == MML::PrecisionValues<long double>::EqualityTolerance);
