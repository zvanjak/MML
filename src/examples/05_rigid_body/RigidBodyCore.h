///////////////////////////////////////////////////////////////////////////////////////////
/// @file RigidBodyCore.h
/// @brief Core rigid body types for Example 05: Rigid Body Collision Simulator
/// @details Self-contained header with all core rigid body types:
///          - ShapeType enum (Box, Sphere)
///          - RigidBodyState struct
///          - SafeNormalize/SafeNormalized utilities
///          - RigidBody abstract base class
///          - RigidBodyBox concrete class
///          - RigidBodySphere concrete class
///          - Factory functions and type checking utilities
///
/// This file is extracted from the MPL library for use in Example 05.
/// Namespace changed from MPL to RigidBodySim to avoid conflicts.
///
/// @author MinimalMathLibrary - Example 05
/// @date February 2026
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef RIGID_BODY_CORE_H
#define RIGID_BODY_CORE_H

#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#include "base/Quaternions.h"
#include "base/Matrix/MatrixNM.h"
#include "base/Vector/VectorN.h"
#include "base/Geometry/Geometry3DBodies.h"

#include <memory>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>

namespace RigidBodySim
{
    using namespace MML;

    // ========================= Shape Type Enum =========================

    /// @brief Enumeration of supported rigid body shape types
    enum class ShapeType
    {
        Box,        ///< Rectangular parallelepiped (6 faces, 8 vertices, 12 edges)
        Sphere      ///< Perfect sphere (infinite symmetry)
    };

    // ========================= Safe Quaternion Normalization =========================

    /// @brief Safely normalize quaternion, returning identity if near-zero
    /// @details Prevents crashes during ODE integration when quaternion becomes
    ///          very small due to numerical precision issues.
    inline Quaternion SafeNormalized(const Quaternion& q)
    {
        Real normSq = q.w()*q.w() + q.x()*q.x() + q.y()*q.y() + q.z()*q.z();
        constexpr Real threshold = 1e-12;
        
        if (normSq < threshold)
            return Quaternion::Identity();
        
        Real invNorm = 1.0 / std::sqrt(normSq);
        return Quaternion(q.w() * invNorm, q.x() * invNorm, 
                         q.y() * invNorm, q.z() * invNorm);
    }

    /// @brief Safely normalize quaternion in-place
    inline void SafeNormalize(Quaternion& q)
    {
        q = SafeNormalized(q);
    }

    // ========================= Rigid Body State =========================

    /// @brief Complete state of a rigid body in 3D space
    /// @details Represents 13 degrees of freedom:
    ///          - 3 DOF: position (x, y, z)
    ///          - 3 DOF: velocity (vx, vy, vz)
    ///          - 4 DOF: orientation quaternion (qw, qx, qy, qz) - normalized, so effectively 3 DOF
    ///          - 3 DOF: angular velocity (ωx, ωy, ωz) in body frame
    struct RigidBodyState
    {
        Vec3Cart position;       ///< Center of mass position in world frame [m]
        Vec3Cart velocity;       ///< Linear velocity of COM in world frame [m/s]
        Quaternion orientation;  ///< Orientation as unit quaternion (body-to-world rotation)
        Vec3Cart angularVel;     ///< Angular velocity in body frame [rad/s]

        /// @brief Default constructor - body at rest at origin
        RigidBodyState()
            : position(0, 0, 0)
            , velocity(0, 0, 0)
            , orientation(Quaternion::Identity())
            , angularVel(0, 0, 0)
        {}

        /// @brief Full state constructor
        /// @param pos Position of center of mass
        /// @param vel Linear velocity
        /// @param orient Orientation quaternion (will be normalized)
        /// @param angVel Angular velocity in body frame
        RigidBodyState(const Vec3Cart& pos, const Vec3Cart& vel, 
                       const Quaternion& orient, const Vec3Cart& angVel)
            : position(pos)
            , velocity(vel)
            , orientation(SafeNormalized(orient))
            , angularVel(angVel)
        {}

        /// @brief Get angular velocity in world frame
        /// @return ω_world = q * ω_body * q^(-1)
        Vec3Cart GetAngularVelWorld() const
        {
            return orientation.Rotate(angularVel);
        }

        /// @brief Set angular velocity from world frame value
        /// @param angVelWorld Angular velocity in world frame
        void SetAngularVelFromWorld(const Vec3Cart& angVelWorld)
        {
            // Transform from world to body: ω_body = q^(-1) * ω_world * q
            angularVel = orientation.Conjugate().Rotate(angVelWorld);
        }

        /// @brief Normalize the orientation quaternion (fixes numerical drift)
        void NormalizeOrientation()
        {
            SafeNormalize(orientation);
        }

        /// @brief Get rotation matrix from body to world frame
        /// @return 3x3 rotation matrix R such that v_world = R * v_body
        MatrixNM<Real, 3, 3> GetRotationMatrix() const
        {
            return orientation.ToRotationMatrix();
        }

        /// @brief Pack state into a 13-component vector for ODE integration
        /// @return [x, y, z, vx, vy, vz, qw, qx, qy, qz, ωx, ωy, ωz]
        VectorN<Real, 13> ToVector() const
        {
            VectorN<Real, 13> v;
            // Position
            v[0] = position[0];
            v[1] = position[1];
            v[2] = position[2];
            // Velocity
            v[3] = velocity[0];
            v[4] = velocity[1];
            v[5] = velocity[2];
            // Orientation quaternion
            v[6] = orientation.w();
            v[7] = orientation.x();
            v[8] = orientation.y();
            v[9] = orientation.z();
            // Angular velocity (body frame)
            v[10] = angularVel[0];
            v[11] = angularVel[1];
            v[12] = angularVel[2];
            return v;
        }

        /// @brief Unpack state from a 13-component vector
        /// @param v State vector [x, y, z, vx, vy, vz, qw, qx, qy, qz, ωx, ωy, ωz]
        static RigidBodyState FromVector(const VectorN<Real, 13>& v)
        {
            RigidBodyState state;
            // Position
            state.position = Vec3Cart(v[0], v[1], v[2]);
            // Velocity
            state.velocity = Vec3Cart(v[3], v[4], v[5]);
            // Orientation quaternion
            state.orientation = Quaternion(v[6], v[7], v[8], v[9]);
            SafeNormalize(state.orientation);  // Ensure unit quaternion (safe)
            // Angular velocity (body frame)
            state.angularVel = Vec3Cart(v[10], v[11], v[12]);
            return state;
        }

        /// @brief Print state to output stream
        friend std::ostream& operator<<(std::ostream& os, const RigidBodyState& state)
        {
            os << "RigidBodyState:\n"
               << "  Position: (" << state.position[0] << ", " 
                                  << state.position[1] << ", " 
                                  << state.position[2] << ") m\n"
               << "  Velocity: (" << state.velocity[0] << ", " 
                                  << state.velocity[1] << ", " 
                                  << state.velocity[2] << ") m/s\n"
               << "  Orientation: [" << state.orientation.w() << ", " 
                                     << state.orientation.x() << ", "
                                     << state.orientation.y() << ", "
                                     << state.orientation.z() << "]\n"
               << "  Angular vel (body): (" << state.angularVel[0] << ", " 
                                            << state.angularVel[1] << ", " 
                                            << state.angularVel[2] << ") rad/s";
            return os;
        }
    };

    // ========================= Forward Declarations =========================
    
    class RigidBodyBox;
    class RigidBodySphere;

    // ========================= Rigid Body Base Class =========================

    /// @brief Abstract base class for all rigid bodies
    /// @details Contains common state and physics shared by all shape types.
    ///          Derived classes must implement shape-specific behavior.
    class RigidBody
    {
    public:
        RigidBodyState _state;              ///< Dynamic state (position, velocity, orientation, angular vel)
        
    protected:
        Real _mass;                         ///< Mass [kg]
        MatrixNM<Real, 3, 3> _inertiaTensorBody;     ///< Inertia tensor in body frame
        MatrixNM<Real, 3, 3> _inertiaTensorBodyInv;  ///< Inverse inertia tensor in body frame

        // ========================= Protected Construction =========================
        
        /// @brief Construct with mass (called by derived classes)
        /// @param mass Mass in kg
        explicit RigidBody(Real mass)
            : _mass(mass)
        {
            // Zero-initialize inertia tensors (derived class fills in)
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    _inertiaTensorBody(i, j) = 0.0;
                    _inertiaTensorBodyInv(i, j) = 0.0;
                }
            }
        }

        /// @brief Construct with mass and initial state (called by derived classes)
        /// @param mass Mass in kg
        /// @param initialState Initial dynamic state
        RigidBody(Real mass, const RigidBodyState& initialState)
            : _state(initialState)
            , _mass(mass)
        {
            // Zero-initialize inertia tensors (derived class fills in)
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    _inertiaTensorBody(i, j) = 0.0;
                    _inertiaTensorBodyInv(i, j) = 0.0;
                }
            }
        }

        /// @brief Compute shape-specific inertia tensor
        /// @details Called by derived class constructors after setting shape parameters
        virtual void ComputeInertiaTensor() = 0;

    public:
        // ========================= Virtual Destructor =========================
        
        /// @brief Virtual destructor for proper cleanup of derived classes
        virtual ~RigidBody() = default;

        // ========================= Pure Virtual Interface =========================
        
        /// @brief Get the shape type (Box, Sphere, etc.)
        [[nodiscard]] virtual ShapeType GetShapeType() const = 0;
        
        /// @brief Get bounding radius (smallest sphere containing the body)
        [[nodiscard]] virtual Real BoundingRadius() const = 0;
        
        /// @brief Create a polymorphic copy
        [[nodiscard]] virtual std::unique_ptr<RigidBody> Clone() const = 0;

        // ========================= Common Accessors =========================

        [[nodiscard]] Real Mass() const { return _mass; }
        
        [[nodiscard]] Vec3Cart& Position() { return _state.position; }
        [[nodiscard]] const Vec3Cart& Position() const { return _state.position; }
        
        [[nodiscard]] Vec3Cart& Velocity() { return _state.velocity; }
        [[nodiscard]] const Vec3Cart& Velocity() const { return _state.velocity; }
        
        [[nodiscard]] Quaternion& Orientation() { return _state.orientation; }
        [[nodiscard]] const Quaternion& Orientation() const { return _state.orientation; }
        
        [[nodiscard]] Vec3Cart& AngularVel() { return _state.angularVel; }
        [[nodiscard]] const Vec3Cart& AngularVel() const { return _state.angularVel; }
        
        /// @brief Get angular velocity in world frame
        [[nodiscard]] Vec3Cart GetAngularVelWorld() const { return _state.GetAngularVelWorld(); }

        [[nodiscard]] const MatrixNM<Real, 3, 3>& InertiaTensorBody() const { return _inertiaTensorBody; }
        [[nodiscard]] const MatrixNM<Real, 3, 3>& InertiaTensorBodyInv() const { return _inertiaTensorBodyInv; }

        // ========================= Inertia Tensor =========================

        /// @brief Get inertia tensor transformed to world frame
        /// @return I_world = R * I_body * R^T
        [[nodiscard]] MatrixNM<Real, 3, 3> GetWorldInertiaTensor() const
        {
            MatrixNM<Real, 3, 3> R = _state.orientation.ToRotationMatrix();
            MatrixNM<Real, 3, 3> Rt = R.transpose();
            return R * _inertiaTensorBody * Rt;
        }

        /// @brief Get inverse inertia tensor in world frame
        /// @return I_world^(-1) = R * I_body^(-1) * R^T
        [[nodiscard]] MatrixNM<Real, 3, 3> GetWorldInertiaTensorInv() const
        {
            MatrixNM<Real, 3, 3> R = _state.orientation.ToRotationMatrix();
            MatrixNM<Real, 3, 3> Rt = R.transpose();
            return R * _inertiaTensorBodyInv * Rt;
        }

        // ========================= Geometry =========================

        /// @brief Get center of the body (same as position for uniform density)
        [[nodiscard]] Vec3Cart GetCenter() const
        {
            return _state.position;
        }

        /// @brief Get the 3 face normal directions in world frame
        /// @return Vector of 3 unit vectors (principal axes in world frame)
        [[nodiscard]] std::vector<Vec3Cart> GetWorldAxes() const
        {
            std::vector<Vec3Cart> axes;
            axes.reserve(3);
            
            axes.push_back(_state.orientation.Rotate(Vec3Cart(1, 0, 0)));
            axes.push_back(_state.orientation.Rotate(Vec3Cart(0, 1, 0)));
            axes.push_back(_state.orientation.Rotate(Vec3Cart(0, 0, 1)));
            
            return axes;
        }

        // ========================= Kinematics =========================

        /// @brief Get velocity of a point on the rigid body (in world frame)
        /// @param worldPoint A point attached to the body, in world coordinates
        /// @return Velocity of that point = v_COM + ω × r
        [[nodiscard]] Vec3Cart GetPointVelocity(const Vec3Cart& worldPoint) const
        {
            Vec3Cart r = worldPoint - _state.position;
            Vec3Cart omega_world = _state.GetAngularVelWorld();
            return _state.velocity + VectorProduct(omega_world, r);
        }

        // ========================= Energy =========================

        /// @brief Compute translational kinetic energy (center of mass motion)
        /// @return KE_trans = ½ m v² [J]
        [[nodiscard]] Real TranslationalKineticEnergy() const
        {
            return 0.5 * _mass * (_state.velocity * _state.velocity);
        }

        /// @brief Compute rotational kinetic energy (spinning motion)
        /// @return KE_rot = ½ ω · I · ω [J]
        [[nodiscard]] Real RotationalKineticEnergy() const
        {
            Vec3Cart Iw = MatVecMult(_inertiaTensorBody, _state.angularVel);
            return 0.5 * (_state.angularVel * Iw);
        }

        /// @brief Compute kinetic energy (translational + rotational)
        /// @return Total kinetic energy [J]
        [[nodiscard]] Real KineticEnergy() const
        {
            return TranslationalKineticEnergy() + RotationalKineticEnergy();
        }

        /// @brief Compute angular momentum in world frame
        /// @return L = I_world · ω_world
        [[nodiscard]] Vec3Cart AngularMomentum() const
        {
            Vec3Cart omega_world = _state.GetAngularVelWorld();
            MatrixNM<Real, 3, 3> I_world = GetWorldInertiaTensor();
            return MatVecMult(I_world, omega_world);
        }

        /// @brief Compute linear momentum
        /// @return p = m · v
        [[nodiscard]] Vec3Cart LinearMomentum() const
        {
            return _state.velocity * _mass;
        }

        // ========================= Helpers =========================

        /// @brief Normalize orientation quaternion (call periodically to fix numerical drift)
        void NormalizeOrientation()
        {
            _state.orientation.Normalize();
        }

        /// @brief Print body information (virtual for derived customization)
        virtual void Print(std::ostream& os) const
        {
            os << "RigidBody (abstract):\n"
               << "  Mass: " << _mass << " kg\n"
               << "  " << _state;
        }

        /// @brief Stream output operator
        friend std::ostream& operator<<(std::ostream& os, const RigidBody& body)
        {
            body.Print(os);
            return os;
        }

    protected:
        /// @brief Matrix-vector multiplication helper (3x3 matrix * 3-vector)
        static Vec3Cart MatVecMult(const MatrixNM<Real, 3, 3>& M, const Vec3Cart& v)
        {
            return Vec3Cart(
                M(0, 0) * v[0] + M(0, 1) * v[1] + M(0, 2) * v[2],
                M(1, 0) * v[0] + M(1, 1) * v[1] + M(1, 2) * v[2],
                M(2, 0) * v[0] + M(2, 1) * v[1] + M(2, 2) * v[2]
            );
        }
    };

    // ========================= Rigid Body Box =========================

    /// @brief Box (parallelepiped) rigid body
    /// @details Axis-aligned box in body frame with half-extents (a, b, c).
    ///          Inertia tensor: I = m/3 * diag(b²+c², a²+c², a²+b²)
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

    // ========================= Rigid Body Sphere =========================

    /// @brief Spherical rigid body
    /// @details Uniform solid sphere with radius r.
    ///          Inertia tensor: I = 2/5 * m * r² * Identity
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

    // ========================= Factory Functions =========================

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

} // namespace RigidBodySim

#endif // RIGID_BODY_CORE_H
