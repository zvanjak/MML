///////////////////////////////////////////////////////////////////////////////////////////
/// @file RigidBodyDynamics.h
/// @brief ODE systems for rigid body dynamics - Example 05: Rigid Body Collision Simulator
/// @details Self-contained header with Euler equations ODE systems:
///          - SingleBodyODESystem (13 DOF)
///          - TwoBodyODESystem (26 DOF)
///          - ExtractBodyState helper function
///          - PackBodyState helper function
///
/// Physics:
///   Euler equations: I·dω/dt = -ω × (I·ω)   (torque-free case)
///   Quaternion kinematics: dq/dt = ½ [0, ω] * q
///
/// This file is extracted from the MPL library for use in Example 05.
/// Namespace changed from MPL to RigidBodySim to avoid conflicts.
///
/// @author MinimalMathLibrary - Example 05
/// @date February 2026
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef RIGID_BODY_DYNAMICS_H
#define RIGID_BODY_DYNAMICS_H

#include "RigidBodyCore.h"
#include "interfaces/IODESystem.h"
#include "base/Vector/Vector.h"

#include <array>

namespace RigidBodySim
{
    using namespace MML;

    // ========================= Single Body ODE System =========================

    /// @brief ODE system for a single rigid body in torque-free motion
    /// @details State vector (13 components):
    ///          [0-2]  : position (x, y, z)
    ///          [3-5]  : velocity (vx, vy, vz)
    ///          [6-9]  : orientation quaternion (qw, qx, qy, qz)
    ///          [10-12]: angular velocity in body frame (ωx, ωy, ωz)
    ///
    /// Equations of motion (torque-free, zero gravity):
    ///   dx/dt = v                           (position derivative)
    ///   dv/dt = 0                           (no forces)
    ///   dq/dt = ½ [0, ω_world] * q          (quaternion kinematics)
    ///   dω/dt = I⁻¹(-ω × (I·ω))             (Euler equations)
    class SingleBodyODESystem : public IODESystem
    {
    private:
        MatrixNM<Real, 3, 3> _inertiaTensorBody;     ///< Inertia tensor in body frame
        MatrixNM<Real, 3, 3> _inertiaTensorBodyInv;  ///< Inverse inertia tensor in body frame

    public:
        /// @brief Construct ODE system for a specific rigid body
        /// @param body The rigid body to simulate
        SingleBodyODESystem(const RigidBody& body)
            : _inertiaTensorBody(body.InertiaTensorBody())
            , _inertiaTensorBodyInv(body.InertiaTensorBodyInv())
        {}

        /// @brief Get system dimension (13 DOF for single body)
        int getDim() const override { return 13; }

        /// @brief Compute state derivatives
        /// @param t Current time (unused for autonomous system)
        /// @param x Current state vector
        /// @param dxdt Output derivative vector
        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
        {
            // Extract state components
            // Position: x[0], x[1], x[2]
            // Velocity: x[3], x[4], x[5]
            Vec3Cart velocity(x[3], x[4], x[5]);
            
            // Orientation quaternion: x[6], x[7], x[8], x[9]
            Quaternion q(x[6], x[7], x[8], x[9]);
            SafeNormalize(q);  // Ensure unit quaternion (safe for near-zero)
            
            // Angular velocity in body frame: x[10], x[11], x[12]
            Vec3Cart omega_body(x[10], x[11], x[12]);

            // ==================== Position derivative ====================
            // dx/dt = v
            dxdt[0] = velocity[0];
            dxdt[1] = velocity[1];
            dxdt[2] = velocity[2];

            // ==================== Velocity derivative ====================
            // dv/dt = F/m = 0 (zero gravity, no external forces)
            dxdt[3] = 0.0;
            dxdt[4] = 0.0;
            dxdt[5] = 0.0;

            // ==================== Quaternion derivative ====================
            // dq/dt = ½ ω ⊗ q  where ω is pure quaternion [0, ω_body]
            // Using the formula: dq/dt = 0.5 * [0, ω_body] * q
            //
            // For body-frame angular velocity:
            // dq/dt = 0.5 * q * [0, ω_body]
            Quaternion omega_quat(0, omega_body[0], omega_body[1], omega_body[2]);
            Quaternion dqdt = q * omega_quat * 0.5;
            
            dxdt[6] = dqdt.w();
            dxdt[7] = dqdt.x();
            dxdt[8] = dqdt.y();
            dxdt[9] = dqdt.z();

            // ==================== Angular velocity derivative ====================
            // Euler's equations for torque-free rotation:
            // I · dω/dt = -ω × (I · ω)
            // dω/dt = I⁻¹ · (-ω × (I · ω))
            //
            // For diagonal inertia tensor (principal axes):
            // I₁ dω₁/dt = (I₂ - I₃) ω₂ ω₃
            // I₂ dω₂/dt = (I₃ - I₁) ω₃ ω₁
            // I₃ dω₃/dt = (I₁ - I₂) ω₁ ω₂
            
            Vec3Cart Iw = MatVecMult(_inertiaTensorBody, omega_body);
            Vec3Cart omega_cross_Iw = VectorProduct(omega_body, Iw);
            Vec3Cart domega = MatVecMult(_inertiaTensorBodyInv, omega_cross_Iw * (-1.0));
            
            dxdt[10] = domega[0];
            dxdt[11] = domega[1];
            dxdt[12] = domega[2];
        }

        /// @brief Get variable name for output
        std::string getVarName(int ind) const override
        {
            static const std::array<std::string, 13> names = {
                "x", "y", "z", "vx", "vy", "vz",
                "qw", "qx", "qy", "qz",
                "wx", "wy", "wz"
            };
            if (ind >= 0 && ind < 13)
                return names[ind];
            return "var" + std::to_string(ind);
        }

    private:
        /// @brief Matrix-vector multiplication helper
        static Vec3Cart MatVecMult(const MatrixNM<Real, 3, 3>& M, const Vec3Cart& v)
        {
            return Vec3Cart(
                M(0, 0) * v[0] + M(0, 1) * v[1] + M(0, 2) * v[2],
                M(1, 0) * v[0] + M(1, 1) * v[1] + M(1, 2) * v[2],
                M(2, 0) * v[0] + M(2, 1) * v[1] + M(2, 2) * v[2]
            );
        }
    };

    // ========================= Two Body ODE System =========================

    /// @brief ODE system for two rigid bodies in torque-free motion
    /// @details State vector (26 components):
    ///          [0-12] : Body 1 state (position, velocity, quaternion, angular vel)
    ///          [13-25]: Body 2 state (position, velocity, quaternion, angular vel)
    class TwoBodyODESystem : public IODESystem
    {
    private:
        std::array<MatrixNM<Real, 3, 3>, 2> _inertiaTensorBody;
        std::array<MatrixNM<Real, 3, 3>, 2> _inertiaTensorBodyInv;

    public:
        /// @brief Construct ODE system for two rigid bodies
        TwoBodyODESystem(const RigidBody& body1, const RigidBody& body2)
        {
            _inertiaTensorBody[0] = body1.InertiaTensorBody();
            _inertiaTensorBody[1] = body2.InertiaTensorBody();
            _inertiaTensorBodyInv[0] = body1.InertiaTensorBodyInv();
            _inertiaTensorBodyInv[1] = body2.InertiaTensorBodyInv();
        }

        int getDim() const override { return 26; }

        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
        {
            // Process each body
            for (int bodyIdx = 0; bodyIdx < 2; bodyIdx++)
            {
                int offset = bodyIdx * 13;
                
                // Velocity
                Vec3Cart velocity(x[offset + 3], x[offset + 4], x[offset + 5]);
                
                // Quaternion
                Quaternion q(x[offset + 6], x[offset + 7], x[offset + 8], x[offset + 9]);
                SafeNormalize(q);  // Safe normalization for ODE integration
                
                // Angular velocity (body frame)
                Vec3Cart omega_body(x[offset + 10], x[offset + 11], x[offset + 12]);

                // Position derivative
                dxdt[offset + 0] = velocity[0];
                dxdt[offset + 1] = velocity[1];
                dxdt[offset + 2] = velocity[2];

                // Velocity derivative (zero for torque-free, gravity-free)
                dxdt[offset + 3] = 0.0;
                dxdt[offset + 4] = 0.0;
                dxdt[offset + 5] = 0.0;

                // Quaternion derivative
                Quaternion omega_quat(0, omega_body[0], omega_body[1], omega_body[2]);
                Quaternion dqdt = q * omega_quat * 0.5;
                
                dxdt[offset + 6] = dqdt.w();
                dxdt[offset + 7] = dqdt.x();
                dxdt[offset + 8] = dqdt.y();
                dxdt[offset + 9] = dqdt.z();

                // Angular velocity derivative (Euler equations)
                Vec3Cart Iw = MatVecMult(_inertiaTensorBody[bodyIdx], omega_body);
                Vec3Cart omega_cross_Iw = VectorProduct(omega_body, Iw);
                Vec3Cart domega = MatVecMult(_inertiaTensorBodyInv[bodyIdx], omega_cross_Iw * (-1.0));
                
                dxdt[offset + 10] = domega[0];
                dxdt[offset + 11] = domega[1];
                dxdt[offset + 12] = domega[2];
            }
        }

        std::string getVarName(int ind) const override
        {
            static const std::array<std::string, 13> baseNames = {
                "x", "y", "z", "vx", "vy", "vz",
                "qw", "qx", "qy", "qz",
                "wx", "wy", "wz"
            };
            
            int bodyIdx = ind / 13;
            int varIdx = ind % 13;
            
            if (varIdx >= 0 && varIdx < 13)
                return "body" + std::to_string(bodyIdx + 1) + "_" + baseNames[varIdx];
            
            return "var" + std::to_string(ind);
        }

    private:
        static Vec3Cart MatVecMult(const MatrixNM<Real, 3, 3>& M, const Vec3Cart& v)
        {
            return Vec3Cart(
                M(0, 0) * v[0] + M(0, 1) * v[1] + M(0, 2) * v[2],
                M(1, 0) * v[0] + M(1, 1) * v[1] + M(1, 2) * v[2],
                M(2, 0) * v[0] + M(2, 1) * v[1] + M(2, 2) * v[2]
            );
        }
    };

    // ========================= State Extraction/Packing Helpers =========================

    /// @brief Helper to extract rigid body state from ODE solution vector
    inline RigidBodyState ExtractBodyState(const Vector<Real>& x, int bodyIndex = 0)
    {
        int offset = bodyIndex * 13;
        
        RigidBodyState state;
        state.position = Vec3Cart(x[offset + 0], x[offset + 1], x[offset + 2]);
        state.velocity = Vec3Cart(x[offset + 3], x[offset + 4], x[offset + 5]);
        state.orientation = Quaternion(x[offset + 6], x[offset + 7], x[offset + 8], x[offset + 9]);
        SafeNormalize(state.orientation);  // Safe normalization
        state.angularVel = Vec3Cart(x[offset + 10], x[offset + 11], x[offset + 12]);
        
        return state;
    }

    /// @brief Helper to pack rigid body state into ODE solution vector
    inline void PackBodyState(Vector<Real>& x, const RigidBodyState& state, int bodyIndex = 0)
    {
        int offset = bodyIndex * 13;
        
        x[offset + 0] = state.position[0];
        x[offset + 1] = state.position[1];
        x[offset + 2] = state.position[2];
        x[offset + 3] = state.velocity[0];
        x[offset + 4] = state.velocity[1];
        x[offset + 5] = state.velocity[2];
        x[offset + 6] = state.orientation.w();
        x[offset + 7] = state.orientation.x();
        x[offset + 8] = state.orientation.y();
        x[offset + 9] = state.orientation.z();
        x[offset + 10] = state.angularVel[0];
        x[offset + 11] = state.angularVel[1];
        x[offset + 12] = state.angularVel[2];
    }

} // namespace RigidBodySim

#endif // RIGID_BODY_DYNAMICS_H
