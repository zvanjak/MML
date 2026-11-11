#ifdef MML_RESERVED_MEMBER_SINGLE_HEADER
#define MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/DiracDeltaFunction.h>
#include <mml/base/Geometry/Geometry3D/Geometry3DPlane.h>
#include <mml/base/Geometry/GeometryBase/Geometry3DShapes.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DSpherePyramid.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DTorusCylinder.h>
#include <mml/core/Surfaces.h>
#include <mml/core/LinAlgEqSolvers/LinAlgQR.h>
#include <mml/systems/DiscreteMaps.h>
#include <mml/systems/LinearSystem.h>
#endif
#include <test_data/stiff_ode_defs.h>
#include <test_data/diff_eq_systems_defs.h>

// Clang's parser checks declarations rather than matching equations, strings, or platform macros.
static_assert(sizeof(MML::Plane3D) > 0);
