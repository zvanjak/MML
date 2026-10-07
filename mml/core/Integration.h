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

/// @file Integration.h
/// @brief Aggregate header for numerical integration.
/// This header includes one-, two-, and three-dimensional integration,
/// improper integrals, Gaussian quadrature, Gauss-Kronrod quadrature, and
/// adaptive multidimensional integration.
///
/// Main types and functions included here:
/// - IntegrationResult, IntegrationDetailedResult, IntegrationConfig, and IntegrationMethod
/// - IntegrateTrap, IntegrateSimpson, IntegrateRomberg, IntegrateGauss10, IntegrateGK21, and Integrate
/// - Integrate2D and Integrate3D - nested quadrature over rectangular and bounded domains
/// - IntegrateOpen, IntegrateUpperInf, IntegrateLowerInf, IntegrateInf, and singular-integral helpers
/// - GaussQuadratureRule, GaussLegendre, GaussLaguerre, GaussHermite, GaussJacobi, and Chebyshev rules
/// - GaussKronrod::GKResult, GaussKronrod::GKRule, and adaptive Gauss-Kronrod helpers
/// - AdaptiveConfig2D, AdaptiveResult2D, AdaptiveConfig3D, AdaptiveResult3D, IntegrateAdaptive2D, and IntegrateAdaptive3D
///
/// For more focused includes, use the individual headers under Integration/.

#if !defined MML_INTEGRATION_H
#define MML_INTEGRATION_H

#include "Integration/Integration1D.h"
#include "Integration/Integration2D.h"
#include "Integration/Integration3D.h"
#include "Integration/IntegrationImproper.h"
#include "Integration/GaussianQuadrature.h"
#include "Integration/GaussKronrod.h"
#include "Integration/Integration2DAdaptive.h"
#include "Integration/Integration3DAdaptive.h"

#endif