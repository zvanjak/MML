///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Integration.h                                                       ///
///  Description: Numerical integration umbrella header                               ///
///               Includes 1D, 2D, 3D, improper, and Gaussian quadrature              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_INTEGRATION_H
#define MML_INTEGRATION_H

/// @file Integration1D.h - 1D numerical integration (Simpson, Romberg, adaptive methods)
#include "Integration/Integration1D.h"
/// @file Integration2D.h - 2D numerical integration over rectangular and arbitrary regions
#include "Integration/Integration2D.h"
/// @file Integration3D.h - 3D numerical integration over volumetric domains
#include "Integration/Integration3D.h"
/// @file IntegrationImproper.h - Improper integrals with infinite limits or singularities
#include "Integration/IntegrationImproper.h"
/// @file GaussianQuadrature.h - Gaussian quadrature with Legendre, Chebyshev, Hermite, Laguerre polynomials
#include "Integration/GaussianQuadrature.h"
/// @file GaussKronrod.h - Gauss-Kronrod quadrature with adaptive error estimation
#include "Integration/GaussKronrod.h"
/// @file Integration2DAdaptive.h - Adaptive 2D integration with quadtree subdivision
#include "Integration/Integration2DAdaptive.h"
/// @file Integration3DAdaptive.h - Adaptive 3D integration with octree subdivision
#include "Integration/Integration3DAdaptive.h"


#endif