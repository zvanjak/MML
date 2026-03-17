///////////////////////////////////////////////////////////////////////////////////////////
///                    Minimal Physics Library (MPL)                                    ///
///                                                                                     ///
///  File:        RigidBodies.h                                                         ///
///  Description: Convenience header that includes all rigid body types                 ///
///                                                                                     ///
///  Part of:     Rigid Body Collision Simulator (Epic MinimalMathLibrary-28bb)         ///
///  Refactored:  Epic 28bb.14 - Inheritance-based hierarchy                            ///
///                                                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MPL_RIGID_BODIES_H
#define MPL_RIGID_BODIES_H

// Abstract base class
#include "RigidBodyBase.h"

// Concrete implementations
#include "RigidBodyBox.h"
#include "RigidBodySphere.h"

// Support types
#include "RigidBodyState.h"
#include "RigidBodyShape.h"  // For ShapeType enum (deprecated variant types)

namespace MPL
{
    // ========================= Factory Functions =========================
    // 
    // These provide a convenient way to create rigid bodies without
    // remembering which derived class to use.

    /// @brief Create a box rigid body
    /// @param mass Mass in kg
    /// @param halfA Half-width along X axis [m]
    /// @param halfB Half-width along Y axis [m]
    /// @param halfC Half-width along Z axis [m]
    /// @return Unique pointer to new RigidBodyBox
    inline std::unique_ptr<RigidBody> CreateBox(Real mass, Real halfA, Real halfB, Real halfC)
    {
        return std::make_unique<RigidBodyBox>(mass, halfA, halfB, halfC);
    }

    /// @brief Create a box rigid body with initial state
    /// @param mass Mass in kg
    /// @param halfExtents Half-widths (a, b, c)
    /// @param state Initial dynamic state
    /// @return Unique pointer to new RigidBodyBox
    inline std::unique_ptr<RigidBody> CreateBox(Real mass, const Vec3Cart& halfExtents, 
                                                 const RigidBodyState& state)
    {
        return std::make_unique<RigidBodyBox>(mass, halfExtents, state);
    }

    /// @brief Create a sphere rigid body
    /// @param mass Mass in kg
    /// @param radius Sphere radius [m]
    /// @return Unique pointer to new RigidBodySphere
    inline std::unique_ptr<RigidBody> CreateSphere(Real mass, Real radius)
    {
        return std::make_unique<RigidBodySphere>(mass, radius);
    }

    /// @brief Create a sphere rigid body with initial state
    /// @param mass Mass in kg
    /// @param radius Sphere radius [m]
    /// @param state Initial dynamic state
    /// @return Unique pointer to new RigidBodySphere
    inline std::unique_ptr<RigidBody> CreateSphere(Real mass, Real radius, 
                                                    const RigidBodyState& state)
    {
        return std::make_unique<RigidBodySphere>(mass, radius, state);
    }

    // ========================= Type Checking Utilities =========================

    /// @brief Check if a rigid body is a box
    /// @param body Pointer to rigid body
    /// @return true if body is a RigidBodyBox
    inline bool IsBox(const RigidBody* body)
    {
        return body && body->GetShapeType() == ShapeType::Box;
    }

    /// @brief Check if a rigid body is a box
    /// @param body Reference to rigid body
    /// @return true if body is a RigidBodyBox
    inline bool IsBox(const RigidBody& body)
    {
        return body.GetShapeType() == ShapeType::Box;
    }

    /// @brief Check if a rigid body is a sphere
    /// @param body Pointer to rigid body
    /// @return true if body is a RigidBodySphere
    inline bool IsSphere(const RigidBody* body)
    {
        return body && body->GetShapeType() == ShapeType::Sphere;
    }

    /// @brief Check if a rigid body is a sphere
    /// @param body Reference to rigid body
    /// @return true if body is a RigidBodySphere
    inline bool IsSphere(const RigidBody& body)
    {
        return body.GetShapeType() == ShapeType::Sphere;
    }

    /// @brief Cast to RigidBodyBox (returns nullptr if not a box)
    /// @param body Pointer to rigid body
    /// @return Pointer to RigidBodyBox, or nullptr
    inline const RigidBodyBox* AsBox(const RigidBody* body)
    {
        if (body && body->GetShapeType() == ShapeType::Box)
            return static_cast<const RigidBodyBox*>(body);
        return nullptr;
    }

    /// @brief Cast to RigidBodyBox (returns nullptr if not a box)
    /// @param body Pointer to rigid body
    /// @return Pointer to RigidBodyBox, or nullptr
    inline RigidBodyBox* AsBox(RigidBody* body)
    {
        if (body && body->GetShapeType() == ShapeType::Box)
            return static_cast<RigidBodyBox*>(body);
        return nullptr;
    }

    /// @brief Cast to RigidBodySphere (returns nullptr if not a sphere)
    /// @param body Pointer to rigid body
    /// @return Pointer to RigidBodySphere, or nullptr
    inline const RigidBodySphere* AsSphere(const RigidBody* body)
    {
        if (body && body->GetShapeType() == ShapeType::Sphere)
            return static_cast<const RigidBodySphere*>(body);
        return nullptr;
    }

    /// @brief Cast to RigidBodySphere (returns nullptr if not a sphere)
    /// @param body Pointer to rigid body
    /// @return Pointer to RigidBodySphere, or nullptr
    inline RigidBodySphere* AsSphere(RigidBody* body)
    {
        if (body && body->GetShapeType() == ShapeType::Sphere)
            return static_cast<RigidBodySphere*>(body);
        return nullptr;
    }

} // namespace MPL

#endif // MPL_RIGID_BODIES_H
