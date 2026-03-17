///////////////////////////////////////////////////////////////////////////////////////////
///                    Minimal Physics Library (MPL)                                    ///
///                                                                                     ///
///  File:        RigidBodyBox.h                                                        ///
///  Description: Box (parallelepiped) rigid body implementation                        ///
///                                                                                     ///
///  Part of:     Rigid Body Collision Simulator (Epic MinimalMathLibrary-28bb)         ///
///  Refactored:  Epic 28bb.14 - Inheritance-based hierarchy                            ///
///                                                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MPL_RIGID_BODY_BOX_H
#define MPL_RIGID_BODY_BOX_H

#include "RigidBodyBase.h"
#include "base/Geometry/Geometry3DBodies.h"  // For MML::Box3D integration

namespace MPL
{
    using namespace MML;

    /// @brief Box (parallelepiped) rigid body
    /// @details Axis-aligned box in body frame with half-extents (a, b, c).
    ///          Inertia tensor: I = m/3 * diag(b²+c², a²+c², a²+b²)
    /// 
    /// Construction is now unambiguous:
    /// @code
    /// RigidBodyBox box(5.0, 0.5, 0.5, 0.5);  // 5kg box with half-extents 0.5m
    /// @endcode
    class RigidBodyBox final : public RigidBody
    {
    private:
        Vec3Cart _halfExtents;              ///< Half-widths (a, b, c) along body axes

    public:
        // ========================= Construction =========================

        /// @brief Construct box with half-extents
        /// @param mass Mass in kg
        /// @param halfA Half-width along X axis [m]
        /// @param halfB Half-width along Y axis [m]
        /// @param halfC Half-width along Z axis [m]
        RigidBodyBox(Real mass, Real halfA, Real halfB, Real halfC)
            : RigidBody(mass)
            , _halfExtents(halfA, halfB, halfC)
        {
            ComputeInertiaTensor();
        }

        /// @brief Construct box with half-extents vector
        /// @param mass Mass in kg
        /// @param halfExtents Half-widths (a, b, c) as Vec3Cart
        RigidBodyBox(Real mass, const Vec3Cart& halfExtents)
            : RigidBody(mass)
            , _halfExtents(halfExtents)
        {
            ComputeInertiaTensor();
        }

        /// @brief Construct box with half-extents and initial state
        /// @param mass Mass in kg
        /// @param halfExtents Half-widths (a, b, c)
        /// @param initialState Initial dynamic state (position, velocity, orientation, angular vel)
        RigidBodyBox(Real mass, const Vec3Cart& halfExtents, const RigidBodyState& initialState)
            : RigidBody(mass, initialState)
            , _halfExtents(halfExtents)
        {
            ComputeInertiaTensor();
        }

        /// @brief Construct box with separate half-extents and initial state
        /// @param mass Mass in kg
        /// @param halfA Half-width along X axis [m]
        /// @param halfB Half-width along Y axis [m]
        /// @param halfC Half-width along Z axis [m]
        /// @param initialState Initial dynamic state
        RigidBodyBox(Real mass, Real halfA, Real halfB, Real halfC, const RigidBodyState& initialState)
            : RigidBody(mass, initialState)
            , _halfExtents(halfA, halfB, halfC)
        {
            ComputeInertiaTensor();
        }

        // ========================= RigidBody Interface =========================

        /// @brief Get shape type
        /// @return ShapeType::Box
        [[nodiscard]] ShapeType GetShapeType() const override 
        { 
            return ShapeType::Box; 
        }

        /// @brief Get bounding radius (half-diagonal of the box)
        /// @return √(a² + b² + c²)
        [[nodiscard]] Real BoundingRadius() const override
        {
            return _halfExtents.NormL2();  // Diagonal half-length
        }

        /// @brief Create a polymorphic copy
        /// @return Unique pointer to cloned RigidBodyBox
        [[nodiscard]] std::unique_ptr<RigidBody> Clone() const override
        {
            return std::make_unique<RigidBodyBox>(*this);
        }

        // ========================= Box-Specific Accessors =========================

        /// @brief Get half-extents (a, b, c)
        /// @return Vec3Cart containing half-widths along each axis
        [[nodiscard]] const Vec3Cart& HalfExtents() const { return _halfExtents; }

        /// @brief Get half-extent for a specific axis (0=X, 1=Y, 2=Z)
        /// @param axis Axis index (0, 1, or 2)
        /// @return Half-width along specified axis
        [[nodiscard]] Real HalfExtent(int axis) const { return _halfExtents[axis]; }

        // ========================= Box Geometry =========================

        /// @brief Get the 8 corner vertices of the box in world coordinates
        /// @return Vector of 8 corner positions
        [[nodiscard]] std::vector<Vec3Cart> GetWorldVertices() const
        {
            std::vector<Vec3Cart> vertices;
            vertices.reserve(8);
            
            Real a = _halfExtents[0];
            Real b = _halfExtents[1];
            Real c = _halfExtents[2];
            
            // 8 corners in local coordinates: (±a, ±b, ±c)
            std::array<Vec3Cart, 8> localCorners = {{
                Vec3Cart(-a, -b, -c),
                Vec3Cart(+a, -b, -c),
                Vec3Cart(-a, +b, -c),
                Vec3Cart(+a, +b, -c),
                Vec3Cart(-a, -b, +c),
                Vec3Cart(+a, -b, +c),
                Vec3Cart(-a, +b, +c),
                Vec3Cart(+a, +b, +c)
            }};
            
            // Transform each corner to world frame
            for (const auto& local : localCorners)
            {
                Vec3Cart world = _state.orientation.Rotate(local);
                world = world + _state.position;
                vertices.push_back(world);
            }
            
            return vertices;
        }

        /// @brief Get the 3 edge direction vectors in world frame (scaled by half-extents)
        /// @return Vector of 3 edge vectors
        [[nodiscard]] std::vector<Vec3Cart> GetWorldEdgeVectors() const
        {
            std::vector<Vec3Cart> edges;
            edges.reserve(3);
            
            edges.push_back(_state.orientation.Rotate(Vec3Cart(_halfExtents[0], 0, 0)));
            edges.push_back(_state.orientation.Rotate(Vec3Cart(0, _halfExtents[1], 0)));
            edges.push_back(_state.orientation.Rotate(Vec3Cart(0, 0, _halfExtents[2])));
            
            return edges;
        }

        /// @brief Get the 6 face centers in world coordinates
        /// @return Vector of 6 face center positions
        [[nodiscard]] std::vector<Vec3Cart> GetWorldFaceCenters() const
        {
            std::vector<Vec3Cart> centers;
            centers.reserve(6);
            
            Real a = _halfExtents[0];
            Real b = _halfExtents[1];
            Real c = _halfExtents[2];
            
            // 6 face centers in local coordinates
            std::array<Vec3Cart, 6> localCenters = {{
                Vec3Cart(-a, 0, 0),  // -X face
                Vec3Cart(+a, 0, 0),  // +X face
                Vec3Cart(0, -b, 0),  // -Y face
                Vec3Cart(0, +b, 0),  // +Y face
                Vec3Cart(0, 0, -c),  // -Z face
                Vec3Cart(0, 0, +c)   // +Z face
            }};
            
            for (const auto& local : localCenters)
            {
                Vec3Cart world = _state.orientation.Rotate(local);
                world = world + _state.position;
                centers.push_back(world);
            }
            
            return centers;
        }

        // ========================= MML Geometry Integration =========================

        /// @brief Convert to MML AABB geometry for visualization
        /// @details Creates an axis-aligned Box3D at current position.
        ///          Note: Rotation is NOT applied (Box3D is axis-aligned).
        ///          For full OBB visualization, use GetWorldVertices() to create a mesh.
        /// @return Box3D representing the bounding box at current position
        [[nodiscard]] MML::Box3D ToMMLGeometry() const
        {
            MML::Pnt3Cart center(Position()[0], Position()[1], Position()[2]);
            return MML::Box3D::FromCenterAndHalfExtents(center, 
                                                        _halfExtents[0], 
                                                        _halfExtents[1], 
                                                        _halfExtents[2]);
        }

        /// @brief Get MML points for the 8 vertices (for mesh export)
        /// @return Vector of 8 MML::Pnt3Cart in world coordinates
        [[nodiscard]] std::vector<MML::Pnt3Cart> ToMMLVertices() const
        {
            auto worldVerts = GetWorldVertices();
            std::vector<MML::Pnt3Cart> points;
            points.reserve(8);
            for (const auto& v : worldVerts)
            {
                points.emplace_back(v[0], v[1], v[2]);
            }
            return points;
        }

        // ========================= Printing =========================

        /// @brief Print box information
        void Print(std::ostream& os) const override
        {
            os << "RigidBodyBox:\n"
               << "  Mass: " << _mass << " kg\n"
               << "  Half-extents: (" << _halfExtents[0] << ", " 
                                      << _halfExtents[1] << ", " 
                                      << _halfExtents[2] << ") m\n"
               << "  Bounding radius: " << BoundingRadius() << " m\n"
               << "  " << _state;
        }

    protected:
        /// @brief Compute box inertia tensor
        /// @details I = m/3 * diag(b²+c², a²+c², a²+b²) for solid box
        void ComputeInertiaTensor() override
        {
            Real a = _halfExtents[0];
            Real b = _halfExtents[1];
            Real c = _halfExtents[2];
            
            Real a2 = a * a;
            Real b2 = b * b;
            Real c2 = c * c;
            Real coeff = _mass / 3.0;
            
            // Diagonal elements
            _inertiaTensorBody(0, 0) = coeff * (b2 + c2);
            _inertiaTensorBody(1, 1) = coeff * (a2 + c2);
            _inertiaTensorBody(2, 2) = coeff * (a2 + b2);
            
            // Off-diagonal elements (zero for symmetric box about principal axes)
            _inertiaTensorBody(0, 1) = 0;
            _inertiaTensorBody(0, 2) = 0;
            _inertiaTensorBody(1, 0) = 0;
            _inertiaTensorBody(1, 2) = 0;
            _inertiaTensorBody(2, 0) = 0;
            _inertiaTensorBody(2, 1) = 0;
            
            // Inverse (trivial for diagonal matrix)
            _inertiaTensorBodyInv(0, 0) = 1.0 / _inertiaTensorBody(0, 0);
            _inertiaTensorBodyInv(1, 1) = 1.0 / _inertiaTensorBody(1, 1);
            _inertiaTensorBodyInv(2, 2) = 1.0 / _inertiaTensorBody(2, 2);
            _inertiaTensorBodyInv(0, 1) = 0;
            _inertiaTensorBodyInv(0, 2) = 0;
            _inertiaTensorBodyInv(1, 0) = 0;
            _inertiaTensorBodyInv(1, 2) = 0;
            _inertiaTensorBodyInv(2, 0) = 0;
            _inertiaTensorBodyInv(2, 1) = 0;
        }
    };

} // namespace MPL

#endif // MPL_RIGID_BODY_BOX_H
