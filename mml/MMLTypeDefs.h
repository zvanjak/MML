///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLTypeDefs.h                                                       ///
///  Description: Core type definition - the Real floating-point type                 ///
///               Minimal header with zero dependencies, included by MMLPrecision.h   ///
///               and MMLExceptions.h before MMLBase.h orchestrates everything        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_TYPEDEFS_H
#define MML_TYPEDEFS_H

///////////////////////////////////////////////////////////////////////////////////////////
// Real Type Configuration
///////////////////////////////////////////////////////////////////////////////////////////
//
// CONFIGURATION MODEL: Per-Library Build (NOT Per-Translation-Unit)
//
// The Real typedef defines the floating-point precision for the entire library.
// This is a **build-time configuration** that must be consistent across ALL translation
// units that use MML together in the same program.
//
// USAGE:
//   1. Choose ONE Real type for your build (double, float, long double, or __float128)
//   2. Ensure ALL source files in your project see the SAME definition
//   3. Do NOT mix different Real types in the same executable
//
// CURRENT CONFIGURATION:
typedef double Real; // default real type
//
// OTHER SUPPORTED OPTIONS (uncomment ONE, comment others):
// typedef float						 Real;    // Lower precision, faster, smaller memory
// typedef long double       Real;    // Extended precision (80-bit on x86, 128-bit on some platforms)
// typedef __float128				 Real;    // Quad precision (GCC only, 128-bit IEEE 754)
//
// ABI AND COMPATIBILITY GUARANTEES:
//
// 1. BINARY COMPATIBILITY:
//    - Changing Real breaks ABI compatibility
//    - All libraries/object files must be recompiled with the same Real type
//    - Linking objects built with different Real types causes undefined behavior
//
// 2. SERIALIZATION:
//    - Data serialized with one Real type may NOT be readable with another
//    - Binary format depends on Real's size and representation
//    - Always document which Real type was used when saving data
//    - Consider text-based formats for cross-precision compatibility
//
// 3. PRECISION CONSTANTS:
//    - Constants in MML::Constants are stored as 'double' literals
//    - For float builds: implicit conversion (acceptable precision loss)
//    - For long double/__float128: constants may have less precision than Real
//    - Use REAL() macro for literals that should match Real precision
//
// 4. API CONTRACTS:
//    - All MML APIs use Real consistently
//    - Template parameters like Vector<Real>, Matrix<Real> must match the global Real
//    - Mixing Real with explicit float/double types is supported but requires care
//
// 5. THREAD SAFETY:
//    - Real type is compile-time constant (no runtime changes)
//    - Safe to use across multiple threads
//    - Thread-local contexts (Defaults) are independent of Real choice
//
// RECOMMENDATIONS:
//   - Use 'double' (default) for most applications (good balance of speed/precision)
//   - Use 'float' when memory or performance is critical and precision allows
//   - Use 'long double' for extended precision needs (note: compiler/platform dependent)
//   - Use '__float128' for maximum precision (GCC only, slower, requires libquadmath)
//
// PORTABILITY NOTES:
//   - 'double' and 'float' are fully portable (C++ standard)
//   - 'long double' precision varies by platform (80-bit x86, 128-bit ARM/PowerPC, 64-bit MSVC)
//   - '__float128' requires GCC and libquadmath linkage
//
///////////////////////////////////////////////////////////////////////////////////////////

#endif

