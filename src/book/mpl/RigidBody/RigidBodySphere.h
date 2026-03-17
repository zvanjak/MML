///////////////////////////////////////////////////////////////////////////////////////////
///                    Minimal Physics Library (MPL)                                    ///
///                                                                                     ///
///  File:        RigidBodySphere.h                                                     ///
///  Description: Spherical rigid body implementation                                   ///
///                                                                                     ///
///  Part of:     Rigid Body Collision Simulator (Epic MinimalMathLibrary-28bb)         ///
///  Refactored:  Epic 28bb.14 - Inheritance-based hierarchy                            ///
///                                                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MPL_RIGID_BODY_SPHERE_H
#define MPL_RIGID_BODY_SPHERE_H

#include "RigidBodyBase.h"
#include "base/Geometry/Geometry3DBodies.h"  // For MML::Sphere3D integration

namespace MPL
{
    using namespace MML;

    /// @brief Spherical rigid body
    /// @details Uniform solid sphere with radius r.
    ///          Inertia tensor: I = 2/5 * m * r² * Identity
    /// 
    /// Construction is now unambiguous:
    /// @code
    /// RigidBodySphere sphere(5.0, 1.0);  // 5kg sphere with radius 1.0m
    /// @endcode
    /// 
    /// No more confusion with RigidBodyBox(5.0, 1.0, 1.0, 1.0)!
    class RigidBodySphere final : public RigidBody
    {
    private:
        Real _radius;                       ///< Sphere radius [m]

    public:
        // ========================= Construction =========================

        /// @brief Construct sphere with radius
        /// @param mass Mass in kg
        /// @param radius Sphere radius [m]
        RigidBodySphere(Real mass, Real radius)
            : RigidBody(mass)
            , _radius(radius)
        {
            ComputeInertiaTensor();
        }

        /// @brief Construct sphere with initial state
        /// @param mass Mass in kg
        /// @param radius Sphere radius [m]
        /// @param initialState Initial dynamic state (position, velocity, orientation, angular vel)
        RigidBodySphere(Real mass, Real radius, const RigidBodyState& initialState)
            : RigidBody(mass, initialState)
            , _radius(radius)
        {
            ComputeInertiaTensor();
        }

        // ========================= RigidBody Interface =========================

        /// @brief Get shape type
        /// @return ShapeType::Sphere
        [[nodiscard]] ShapeType GetShapeType() const override 
        { 
            return ShapeType::Sphere; 
        }

        /// @brief Get bounding radius (same as radius for sphere)
        /// @return Sphere radius
        [[nodiscard]] Real BoundingRadius() const override
        {
            return _radius;
        }

        /// @brief Create a polymorphic copy
        /// @return Unique pointer to cloned RigidBodySphere
        [[nodiscard]] std::unique_ptr<RigidBody> Clone() const override
        {
            return std::make_unique<RigidBodySphere>(*this);
        }

        // ========================= Sphere-Specific Accessors =========================

        /// @brief Get radius (type-safe, no throwing!)
        /// @return Sphere radius in meters
        /// @note Unlike the old RigidBody::Radius(), this never throws!
        [[nodiscard]] Real Radius() const { return _radius; }

        // ========================= Sphere Geometry =========================

        /// @brief Get surface area of the sphere
        /// @return 4πr² [m²]
        [[nodiscard]] Real SurfaceArea() const
        {
            return 4.0 * Constants::PI * _radius * _radius;
        }

        /// @brief Get volume of the sphere
        /// @return 4/3 πr³ [m³]
        [[nodiscard]] Real Volume() const
        {
            return (4.0 / 3.0) * Constants::PI * _radius * _radius * _radius;
        }

        /// @brief Get density (assuming uniform density)
        /// @return Mass / Volume [kg/m³]
        [[nodiscard]] Real Density() const
        {
            return _mass / Volume();
        }

        /// @brief Get a point on the sphere surface in world coordinates
        /// @param theta Polar angle from +Z axis [0, π]
        /// @param phi Azimuthal angle from +X axis [0, 2π]
        /// @return Point on sphere surface in world frame
        [[nodiscard]] Vec3Cart GetSurfacePoint(Real theta, Real phi) const
        {
            // Spherical to Cartesian (local frame)
            Real sinTheta = std::sin(theta);
            Vec3Cart local(
                _radius * sinTheta * std::cos(phi),
                _radius * sinTheta * std::sin(phi),
                _radius * std::cos(theta)
            );
            
            // Transform to world frame
            Vec3Cart world = _state.orientation.Rotate(local);
            return world + _state.position;
        }

        // ========================= MML Geometry Integration =========================

        /// @brief Convert to MML Sphere3D geometry for visualization
        /// @details Creates a Sphere3D at current position.
        /// @param numLatitude Latitude divisions for mesh generation (default 16)
        /// @param numLongitude Longitude divisions for mesh generation (default 20)
        /// @return Sphere3D representing the sphere at current position
        [[nodiscard]] MML::Sphere3D ToMMLGeometry(int numLatitude = 16, int numLongitude = 20) const
        {
            MML::Pnt3Cart center(Position()[0], Position()[1], Position()[2]);
            return MML::Sphere3D(_radius, center, numLatitude, numLongitude);
        }

        // ========================= Printing =========================

        /// @brief Print sphere information
        void Print(std::ostream& os) const override
        {
            os << "RigidBodySphere:\n"
               << "  Mass: " << _mass << " kg\n"
               << "  Radius: " << _radius << " m\n"
               << "  Volume: " << Volume() << " m³\n"
               << "  Density: " << Density() << " kg/m³\n"
               << "  " << _state;
        }

    protected:
        /// @brief Compute sphere inertia tensor
        /// @details I = 2/5 * m * r² * Identity for solid sphere
        void ComputeInertiaTensor() override
        {
            Real r2 = _radius * _radius;
            Real I = (2.0 / 5.0) * _mass * r2;
            
            // Diagonal elements (all equal for sphere)
            _inertiaTensorBody(0, 0) = I;
            _inertiaTensorBody(1, 1) = I;
            _inertiaTensorBody(2, 2) = I;
            
            // Off-diagonal elements (zero for sphere)
            _inertiaTensorBody(0, 1) = 0;
            _inertiaTensorBody(0, 2) = 0;
            _inertiaTensorBody(1, 0) = 0;
            _inertiaTensorBody(1, 2) = 0;
            _inertiaTensorBody(2, 0) = 0;
            _inertiaTensorBody(2, 1) = 0;
            
            // Inverse (trivial for diagonal matrix)
            Real invI = 1.0 / I;
            _inertiaTensorBodyInv(0, 0) = invI;
            _inertiaTensorBodyInv(1, 1) = invI;
            _inertiaTensorBodyInv(2, 2) = invI;
            _inertiaTensorBodyInv(0, 1) = 0;
            _inertiaTensorBodyInv(0, 2) = 0;
            _inertiaTensorBodyInv(1, 0) = 0;
            _inertiaTensorBodyInv(1, 2) = 0;
            _inertiaTensorBodyInv(2, 0) = 0;
            _inertiaTensorBodyInv(2, 1) = 0;
        }
    };

} // namespace MPL

#endif // MPL_RIGID_BODY_SPHERE_H
