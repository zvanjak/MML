///////////////////////////////////////////////////////////////////////////////////////////
///                    Minimal Physics Library (MPL)                                    ///
///                                                                                     ///
///  File:        RigidBodyBase.h                                                       ///
///  Description: Abstract base class for rigid body dynamics                           ///
///               Derived classes: RigidBodyBox, RigidBodySphere                        ///
///                                                                                     ///
///  Part of:     Rigid Body Collision Simulator (Epic MinimalMathLibrary-28bb)         ///
///  Refactored:  Epic 28bb.14 - Inheritance-based hierarchy                            ///
///                                                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MPL_RIGID_BODY_BASE_H
#define MPL_RIGID_BODY_BASE_H

#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#include "base/Quaternions.h"
#include "base/Matrix/MatrixNM.h"

#include "RigidBodyState.h"
#include "RigidBodyShape.h"  // For ShapeType enum

#include <memory>
#include <vector>
#include <array>

namespace MPL
{
    using namespace MML;

    /// @brief Abstract base class for all rigid bodies
    /// @details Contains common state and physics shared by all shape types.
    ///          Derived classes must implement shape-specific behavior.
    /// 
    /// @note Replaces the old std::variant-based RigidBody class.
    ///       Use RigidBodyBox or RigidBodySphere for construction.
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

} // namespace MPL

#endif // MPL_RIGID_BODY_BASE_H
