///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/FrameTags.h                                    ///
///  Description: Compile-time coordinate frame tags for typed geometry APIs          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DIFFERENTIAL_GEOMETRY_FRAME_TAGS_H
#define MML_DIFFERENTIAL_GEOMETRY_FRAME_TAGS_H

namespace MML
{
	struct Cartesian2 { };
	struct Cartesian3 { };
	struct Polar2 { };
	struct Cylindrical3 { };
	struct Spherical3 { };
	struct Parameter2 { };
}

#endif // MML_DIFFERENTIAL_GEOMETRY_FRAME_TAGS_H