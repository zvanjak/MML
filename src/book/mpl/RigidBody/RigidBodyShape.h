///////////////////////////////////////////////////////////////////////////////////////////
/// @file RigidBodyShape.h
/// @brief Shape abstraction for rigid bodies - DEPRECATED in favor of inheritance
///
/// @deprecated This file contains the old std::variant-based shape abstraction.
///             Use the new inheritance-based types instead:
///             - RigidBody (abstract base)  <- RigidBodyBase.h
///             - RigidBodyBox               <- RigidBodyBox.h
///             - RigidBodySphere            <- RigidBodySphere.h
///             - Convenience header         <- RigidBodies.h
///
/// RETAINED TYPES (still valid):
/// - ShapeType enum (Box, Sphere)
/// - ShapeTypeName() function
///
/// DEPRECATED TYPES (will be removed in v2.0):
/// - BoxShape struct
/// - SphereShape struct
/// - RigidBodyShape variant
/// - GetShapeType(), IsBox(), IsSphere(), GetBox(), GetSphere() for variants
///
/// @author Generated for MinimalMathLibrary
/// @date February 2026
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MPL_RIGID_BODY_SHAPE_H
#define MPL_RIGID_BODY_SHAPE_H

#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"

#include <variant>
#include <stdexcept>

namespace MPL
{
    using namespace MML;

    // ========================= Shape Type Enum =========================

    /// @brief Enumeration of supported rigid body shape types
    enum class ShapeType
    {
        Box,        ///< Rectangular parallelepiped (6 faces, 8 vertices, 12 edges)
        Sphere      ///< Perfect sphere (infinite symmetry)
    };

    // ========================= Shape Structs (DEPRECATED) =========================

    /// @brief Box shape defined by half-extents along each axis
    /// @deprecated Use RigidBodyBox instead (from RigidBodies.h)
    /// @details A box centered at origin with corners at (±a, ±b, ±c)
    ///          Full dimensions are 2a × 2b × 2c
    struct [[deprecated("Use RigidBodyBox from RigidBodies.h instead")]] BoxShape
    {
        Vec3Cart halfExtents;   ///< Half-widths (a, b, c) along X, Y, Z axes [m]

        /// @brief Default constructor - unit cube (half-extents = 0.5)
        BoxShape() : halfExtents(0.5, 0.5, 0.5) {}

        /// @brief Construct with specified half-extents
        explicit BoxShape(const Vec3Cart& extents) : halfExtents(extents) {}

        /// @brief Construct with individual half-extents
        BoxShape(Real halfA, Real halfB, Real halfC) 
            : halfExtents(halfA, halfB, halfC) {}

        /// @brief Get half-extent for a specific axis (0=X, 1=Y, 2=Z)
        Real HalfExtent(int axis) const { return halfExtents[axis]; }

        /// @brief Compute bounding sphere radius (for broad-phase collision)
        /// @return Distance from center to farthest corner
        Real BoundingRadius() const
        {
            return halfExtents.NormL2();
        }
    };

    /// @brief Sphere shape defined by radius
    /// @deprecated Use RigidBodySphere instead (from RigidBodies.h)
    /// @details A sphere centered at origin with given radius
    struct [[deprecated("Use RigidBodySphere from RigidBodies.h instead")]] SphereShape
    {
        Real radius;    ///< Radius of the sphere [m]

        /// @brief Default constructor - unit sphere (radius = 1.0)
        SphereShape() : radius(1.0) {}

        /// @brief Construct with specified radius
        explicit SphereShape(Real r) : radius(r) {}

        /// @brief Compute bounding sphere radius (trivial for sphere)
        /// @return The radius itself
        Real BoundingRadius() const
        {
            return radius;
        }
    };

    // ========================= Variant Type (DEPRECATED) =========================

    /// @brief Type-safe union of all supported shapes
    /// @deprecated Use polymorphic RigidBody hierarchy instead (RigidBodyBox, RigidBodySphere)
    using RigidBodyShape = std::variant<BoxShape, SphereShape>;

    // ========================= Utility Functions (DEPRECATED) =========================

    /// @brief Get the shape type from a RigidBodyShape variant
    /// @deprecated Use body.GetShapeType() on polymorphic RigidBody instead
    /// @param shape The shape variant to query
    /// @return ShapeType::Box or ShapeType::Sphere
    [[deprecated("Use body.GetShapeType() on polymorphic RigidBody instead")]]
    inline ShapeType GetShapeType(const RigidBodyShape& shape)
    {
        return std::visit([](const auto& s) -> ShapeType {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, BoxShape>)
                return ShapeType::Box;
            else if constexpr (std::is_same_v<T, SphereShape>)
                return ShapeType::Sphere;
        }, shape);
    }

    /// @brief Check if shape is a box
    /// @deprecated Use IsBox(rigidBody) from RigidBodies.h instead
    [[deprecated("Use IsBox(rigidBody) from RigidBodies.h instead")]]
    inline bool IsBox(const RigidBodyShape& shape)
    {
        return std::holds_alternative<BoxShape>(shape);
    }

    /// @brief Check if shape is a sphere
    /// @deprecated Use IsSphere(rigidBody) from RigidBodies.h instead
    [[deprecated("Use IsSphere(rigidBody) from RigidBodies.h instead")]]
    inline bool IsSphere(const RigidBodyShape& shape)
    {
        return std::holds_alternative<SphereShape>(shape);
    }

    /// @brief Get box shape (throws if not a box)
    /// @deprecated Use AsBox(rigidBody) from RigidBodies.h instead
    [[deprecated("Use AsBox(rigidBody) from RigidBodies.h instead")]]
    inline const BoxShape& GetBox(const RigidBodyShape& shape)
    {
        if (!IsBox(shape))
            throw std::runtime_error("Shape is not a box");
        return std::get<BoxShape>(shape);
    }

    /// @brief Get box shape mutable reference (throws if not a box)
    /// @deprecated Use AsBox(rigidBody) from RigidBodies.h instead
    [[deprecated("Use AsBox(rigidBody) from RigidBodies.h instead")]]
    inline BoxShape& GetBox(RigidBodyShape& shape)
    {
        if (!IsBox(shape))
            throw std::runtime_error("Shape is not a box");
        return std::get<BoxShape>(shape);
    }

    /// @brief Get sphere shape (throws if not a sphere)
    /// @deprecated Use AsSphere(rigidBody) from RigidBodies.h instead
    [[deprecated("Use AsSphere(rigidBody) from RigidBodies.h instead")]]
    inline const SphereShape& GetSphere(const RigidBodyShape& shape)
    {
        if (!IsSphere(shape))
            throw std::runtime_error("Shape is not a sphere");
        return std::get<SphereShape>(shape);
    }

    /// @brief Get sphere shape mutable reference (throws if not a sphere)
    /// @deprecated Use AsSphere(rigidBody) from RigidBodies.h instead
    [[deprecated("Use AsSphere(rigidBody) from RigidBodies.h instead")]]
    inline SphereShape& GetSphere(RigidBodyShape& shape)
    {
        if (!IsSphere(shape))
            throw std::runtime_error("Shape is not a sphere");
        return std::get<SphereShape>(shape);
    }

    /// @brief Get bounding radius for any shape (for broad-phase collision detection)
    /// @deprecated Use body.BoundingRadius() on polymorphic RigidBody instead
    /// @param shape The shape to query
    /// @return Radius of smallest enclosing sphere centered at shape origin
    [[deprecated("Use body.BoundingRadius() on polymorphic RigidBody instead")]]
    inline Real GetBoundingRadius(const RigidBodyShape& shape)
    {
        return std::visit([](const auto& s) -> Real {
            return s.BoundingRadius();
        }, shape);
    }

    // ========================= Utility Functions (NOT DEPRECATED) =========================

    /// @brief Get shape type as string (for debugging/output)
    /// @note This function is NOT deprecated - still useful for the ShapeType enum
    inline const char* ShapeTypeName(ShapeType type)
    {
        switch (type)
        {
            case ShapeType::Box:    return "Box";
            case ShapeType::Sphere: return "Sphere";
            default:                return "Unknown";
        }
    }

    /// @brief Get shape type name from shape variant
    /// @deprecated Use ShapeTypeName(body.GetShapeType()) instead
    [[deprecated("Use ShapeTypeName(body.GetShapeType()) instead")]]
    inline const char* GetShapeTypeName(const RigidBodyShape& shape)
    {
        return ShapeTypeName(GetShapeType(shape));
    }

} // namespace MPL

#endif // MPL_RIGID_BODY_SHAPE_H
