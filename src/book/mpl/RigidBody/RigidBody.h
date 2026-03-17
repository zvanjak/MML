///////////////////////////////////////////////////////////////////////////////////////////
///                    Minimal Physics Library (MPL)                                    ///
///                                                                                     ///
///  File:        RigidBody.h                                                           ///
///  Description: LEGACY COMPATIBILITY HEADER - forwards to new inheritance hierarchy  ///
///                                                                                     ///
///  DEPRECATED: This header exists for backward compatibility only.                    ///
///              New code should use RigidBodies.h and the derived classes:             ///
///              - RigidBodyBox for boxes/parallelepipeds                               ///
///              - RigidBodySphere for spheres                                          ///
///                                                                                     ///
///  Part of:     Rigid Body Collision Simulator (Epic MinimalMathLibrary-28bb)         ///
///  Refactored:  Epic 28bb.14 - Inheritance-based hierarchy                            ///
///                                                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MPL_RIGID_BODY_H
#define MPL_RIGID_BODY_H

// Forward to new inheritance-based headers
#include "RigidBodies.h"

// NOTE: The old monolithic RigidBody class with std::variant has been replaced
// by an inheritance hierarchy:
//
//   RigidBody (abstract base in RigidBodyBase.h)
//   ├── RigidBodyBox (RigidBodyBox.h)
//   └── RigidBodySphere (RigidBodySphere.h)
//
// The old constructors are NO LONGER AVAILABLE:
//   OLD: RigidBody(mass, halfA, halfB, halfC)  -- REMOVED
//   OLD: RigidBody(mass, radius)               -- REMOVED (was ambiguous!)
//
// Use the new unambiguous constructors:
//   NEW: RigidBodyBox(mass, halfA, halfB, halfC)
//   NEW: RigidBodySphere(mass, radius)
//
// Or use factory functions:
//   MPL::CreateBox(mass, halfA, halfB, halfC)
//   MPL::CreateSphere(mass, radius)

#endif // MPL_RIGID_BODY_H
