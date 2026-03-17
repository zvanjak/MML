///////////////////////////////////////////////////////////////////////////////////////////
/// @file CollisionResponse.h
/// @brief Impulse-based collision response for rigid body dynamics
/// @details Implements physically correct elastic collision response between rigid bodies
///          and between rigid bodies and walls.
///
/// PHYSICS:
/// For collision at point P with normal n (pointing from A to B):
/// 
/// Relative velocity at contact point:
///   v_rel = (v_B + ω_B × r_BP) - (v_A + ω_A × r_AP)
/// 
/// Normal component of relative velocity:
///   v_n = v_rel · n
/// 
/// If v_n > 0, objects are separating → no impulse needed
/// 
/// Impulse magnitude (coefficient of restitution e):
///   j = -(1 + e) * v_n / (1/m_A + 1/m_B + n·(I_A^(-1)(r_AP × n) × r_AP) + n·(I_B^(-1)(r_BP × n) × r_BP))
/// 
/// Apply impulse:
///   v_A -= j * n / m_A
///   v_B += j * n / m_B
///   ω_A -= I_A^(-1) * (r_AP × j*n)
///   ω_B += I_B^(-1) * (r_BP × j*n)
///
/// REFERENCE: Physics for Game Developers by David M. Bourg, Chapter 5
///            Real-Time Collision Detection by Christer Ericson, Chapter 5
///
/// @author Generated for MinimalMathLibrary
/// @date January 2026
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MPL_COLLISION_RESPONSE_H
#define MPL_COLLISION_RESPONSE_H

#include "mpl/RigidBody/RigidBodies.h"
#include "mpl/RigidBody/CollisionDetection.h"

namespace MPL {

    using namespace MML;

    /// @brief Collision response parameters
    struct CollisionParams
    {
        Real coeffRestitution = 1.0;    ///< 1.0 = perfectly elastic, 0.0 = inelastic
        Real minSeparationSpeed = 1e-6; ///< Below this, objects are considered resting
    };

    /// @brief Impulse-based collision response for rigid bodies
    class CollisionResponse
    {
    public:
        /// @brief Apply collision response between two rigid bodies
        /// @param A First rigid body (modified)
        /// @param B Second rigid body (modified)
        /// @param collision Collision info from detection
        /// @param params Collision parameters (restitution)
        /// @return True if impulse was applied
        static bool ApplyBoxBoxImpulse(
            RigidBody& A, 
            RigidBody& B, 
            const CollisionInfo& collision,
            const CollisionParams& params = CollisionParams())
        {
            if (!collision.hasCollision)
                return false;
            
            Vec3Cart n = collision.contactNormal;  // Points from A to B
            Vec3Cart P = collision.contactPoint;
            
            // Vectors from centers to contact point
            Vec3Cart rAP = P - A.Position();
            Vec3Cart rBP = P - B.Position();
            
            // Velocities at contact point
            // v_contact = v_center + ω × r
            Vec3Cart vAP = A.Velocity() + VectorProduct(A.GetAngularVelWorld(), rAP);
            Vec3Cart vBP = B.Velocity() + VectorProduct(B.GetAngularVelWorld(), rBP);
            
            // Relative velocity (B relative to A) at contact point
            Vec3Cart vRel = vBP - vAP;
            
            // Normal component of relative velocity
            Real vn = ScalarProduct(vRel, n);
            
            // If objects are separating, no impulse needed
            if (vn > params.minSeparationSpeed)
                return false;
            
            // Compute impulse denominator
            Real invMassA = 1.0 / A.Mass();
            Real invMassB = 1.0 / B.Mass();
            
            // Get inverse inertia tensors in world frame
            MatrixNM<Real, 3, 3> invIA = A.GetWorldInertiaTensorInv();
            MatrixNM<Real, 3, 3> invIB = B.GetWorldInertiaTensorInv();
            
            // Angular contribution: n · (I^(-1) * (r × n)) × r
            Vec3Cart rAxN = VectorProduct(rAP, n);
            Vec3Cart rBxN = VectorProduct(rBP, n);
            
            // I^(-1) * (r × n)
            Vec3Cart invIA_rAxN = MatrixVectorProduct(invIA, rAxN);
            Vec3Cart invIB_rBxN = MatrixVectorProduct(invIB, rBxN);
            
            // (I^(-1) * (r × n)) × r
            Vec3Cart angularTermA = VectorProduct(invIA_rAxN, rAP);
            Vec3Cart angularTermB = VectorProduct(invIB_rBxN, rBP);
            
            // n · angular_term
            Real angularA = ScalarProduct(n, angularTermA);
            Real angularB = ScalarProduct(n, angularTermB);
            
            // Full denominator
            Real denom = invMassA + invMassB + angularA + angularB;
            
            // Compute impulse magnitude
            Real e = params.coeffRestitution;
            Real j = -(1.0 + e) * vn / denom;
            
            // Apply impulse
            Vec3Cart impulse = n * j;
            
            // Linear velocity changes
            A.Velocity() = A.Velocity() - impulse * invMassA;
            B.Velocity() = B.Velocity() + impulse * invMassB;
            
            // Angular velocity changes (in world frame, then convert to body frame)
            Vec3Cart deltaOmegaA_world = MatrixVectorProduct(invIA, VectorProduct(rAP, impulse));
            Vec3Cart deltaOmegaB_world = MatrixVectorProduct(invIB, VectorProduct(rBP, impulse));
            
            // Convert to body frame (ω_body = R^T * ω_world)
            MatrixNM<Real, 3, 3> RA = A.Orientation().ToRotationMatrix();
            MatrixNM<Real, 3, 3> RB = B.Orientation().ToRotationMatrix();
            MatrixNM<Real, 3, 3> RAt = RA.transpose();
            MatrixNM<Real, 3, 3> RBt = RB.transpose();
            
            Vec3Cart deltaOmegaA_body = MatrixVectorProduct(RAt, deltaOmegaA_world);
            Vec3Cart deltaOmegaB_body = MatrixVectorProduct(RBt, deltaOmegaB_world);
            
            A.AngularVel() = A.AngularVel() - deltaOmegaA_body;
            B.AngularVel() = B.AngularVel() + deltaOmegaB_body;
            
            // Position correction (separate overlapping objects)
            Real correctionRatio = 0.5;  // Split correction between both bodies
            A.Position() = A.Position() - n * (collision.penetrationDepth * correctionRatio);
            B.Position() = B.Position() + n * (collision.penetrationDepth * correctionRatio);
            
            return true;
        }

        /// @brief Apply collision response between any two rigid bodies (unified)
        /// @details The impulse physics are shape-independent - this method works for
        ///          any combination of shapes (box-box, sphere-sphere, sphere-box).
        ///          It is functionally identical to ApplyBoxBoxImpulse but named 
        ///          generically to reflect its universal applicability.
        ///
        /// This is the RECOMMENDED method for body-body collision response.
        ///
        /// @param A First rigid body (modified)
        /// @param B Second rigid body (modified)
        /// @param collision Collision info from detection
        /// @param params Collision parameters (restitution)
        /// @return True if impulse was applied
        static bool ApplyBodyBodyImpulse(
            RigidBody& A, 
            RigidBody& B, 
            const CollisionInfo& collision,
            const CollisionParams& params = CollisionParams())
        {
            // Physics are identical regardless of shape - delegate to the implementation
            return ApplyBoxBoxImpulse(A, B, collision, params);
        }

        /// @brief Apply collision response between a rigid body and a wall
        /// @param box Rigid body (modified)
        /// @param collision Wall collision info
        /// @param params Collision parameters
        /// @return True if impulse was applied
        static bool ApplyWallImpulse(
            RigidBody& box,
            const WallCollisionInfo& collision,
            const CollisionParams& params = CollisionParams())
        {
            if (!collision.hasCollision)
                return false;
            
            Vec3Cart n = collision.contactNormal;  // Points into container
            Vec3Cart P = collision.contactPoint;
            
            // Vector from center to contact point
            Vec3Cart rP = P - box.Position();
            
            // Velocity at contact point
            Vec3Cart vP = box.Velocity() + VectorProduct(box.GetAngularVelWorld(), rP);
            
            // Normal component of velocity
            // Since n points INTO container, vP·n < 0 means moving toward wall
            Real vn = ScalarProduct(vP, n);
            
            // If moving away from wall (into container), no impulse needed
            if (vn > -params.minSeparationSpeed)
                return false;
            
            // Compute impulse (wall has infinite mass, zero velocity)
            Real invMass = 1.0 / box.Mass();
            
            MatrixNM<Real, 3, 3> invI = box.GetWorldInertiaTensorInv();
            
            Vec3Cart rPxN = VectorProduct(rP, n);
            Vec3Cart invI_rPxN = MatrixVectorProduct(invI, rPxN);
            Vec3Cart angularTerm = VectorProduct(invI_rPxN, rP);
            Real angularContrib = ScalarProduct(n, angularTerm);
            
            Real denom = invMass + angularContrib;
            
            Real e = params.coeffRestitution;
            Real j = -(1.0 + e) * vn / denom;
            
            // Apply impulse (wall normal points inward, so positive j means bouncing)
            Vec3Cart impulse = n * j;
            
            box.Velocity() = box.Velocity() + impulse * invMass;
            
            Vec3Cart deltaOmega_world = MatrixVectorProduct(invI, VectorProduct(rP, impulse));
            
            MatrixNM<Real, 3, 3> R = box.Orientation().ToRotationMatrix();
            MatrixNM<Real, 3, 3> Rt = R.transpose();
            Vec3Cart deltaOmega_body = MatrixVectorProduct(Rt, deltaOmega_world);
            
            box.AngularVel() = box.AngularVel() + deltaOmega_body;
            
            // Position correction (push box back inside container)
            box.Position() = box.Position() + n * collision.penetrationDepth;
            
            return true;
        }

    private:
        /// @brief Multiply 3x3 matrix by 3-vector
        static Vec3Cart MatrixVectorProduct(const MatrixNM<Real, 3, 3>& M, const Vec3Cart& v)
        {
            return Vec3Cart(
                M(0, 0) * v.X() + M(0, 1) * v.Y() + M(0, 2) * v.Z(),
                M(1, 0) * v.X() + M(1, 1) * v.Y() + M(1, 2) * v.Z(),
                M(2, 0) * v.X() + M(2, 1) * v.Y() + M(2, 2) * v.Z()
            );
        }
    };

} // namespace MPL

#endif // MPL_COLLISION_RESPONSE_H
