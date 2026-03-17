///////////////////////////////////////////////////////////////////////////////////////////
///                    Minimal Physics Library (MPL)                                ///
///                                                                                 ///
///  File:        RigidBodyState.h                                                  ///
///  Description: State representation for rigid body dynamics                       ///
///               Position, velocity, orientation (quaternion), angular velocity    ///
///                                                                                 ///
///  Part of:     Rigid Body Collision Simulator (Epic MinimalMathLibrary-28bb)     ///
///                                                                                 ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MPL_RIGID_BODY_STATE_H
#define MPL_RIGID_BODY_STATE_H

#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#include "base/Quaternions.h"

#include <cmath>

namespace MPL
{
    using namespace MML;

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

} // namespace MPL

#endif // MPL_RIGID_BODY_STATE_H
