///////////////////////////////////////////////////////////////////////////////////////////
/// @file CollisionDetection.h
/// @brief Collision detection for oriented bounding boxes, spheres, and walls
/// @details Implements collision detection algorithms for multiple primitive types:
///          - OBB-OBB: Separating Axis Theorem (SAT) with 15 axes
///          - Sphere-Sphere: Distance-based detection (O(1))
///          - Sphere-Wall: Center distance to wall vs radius (O(1))
///          - Box-Wall: Vertex-based detection
///          - Sphere-Box: Closest point on OBB algorithm (TODO)
///
/// ALGORITHM OVERVIEW (OBB-OBB):
/// Test 15 potential separating axes:
/// - 3 face normals of box A
/// - 3 face normals of box B
/// - 9 cross-products of edges (3 × 3)
///
/// For each axis:
/// - Project both boxes onto the axis
/// - Check if projections overlap
/// - If ANY axis has no overlap → no collision
/// - If ALL axes overlap → collision detected (use minimum penetration axis)
///
/// ALGORITHM OVERVIEW (Sphere-Sphere):
/// - Distance between centers vs sum of radii
/// - Contact point at midpoint along connecting line
/// - Contact normal = (center_B - center_A).Normalized()
///
/// REFERENCE: Real-Time Collision Detection by Christer Ericson, Chapters 4 & 5
///
/// @author Generated for MinimalMathLibrary
/// @date January 2026
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MPL_COLLISION_DETECTION_H
#define MPL_COLLISION_DETECTION_H

#include "mpl/RigidBody/RigidBodies.h"

namespace MPL {

    using namespace MML;

    /// @brief Result of collision detection between two bodies
    struct CollisionInfo
    {
        bool hasCollision = false;          ///< Whether collision was detected
        Vec3Cart contactPoint;              ///< World-space contact point
        Vec3Cart contactNormal;             ///< Normal pointing from A to B
        Real penetrationDepth = 0.0;        ///< How deep the objects overlap
        
        // Additional info for debugging/visualization
        int separatingAxisType = -1;        ///< Which axis type detected collision
                                            ///< 0-2: face of A, 3-5: face of B, 6-14: edge-edge
    };

    /// @brief Result of collision detection with a wall
    struct WallCollisionInfo
    {
        bool hasCollision = false;          ///< Whether collision was detected
        Vec3Cart contactPoint;              ///< Contact point on box surface
        Vec3Cart contactNormal;             ///< Wall normal (points into container)
        Real penetrationDepth = 0.0;        ///< How far past wall boundary
        int wallIndex = -1;                 ///< Which wall: 0=-X, 1=+X, 2=-Y, 3=+Y, 4=-Z, 5=+Z
    };

    /// @brief Collision detection for rigid body dynamics
    /// @details Implements SAT algorithm for OBB-OBB and box-wall detection
    class CollisionDetector
    {
    public:
        /// Small epsilon for numerical stability
        static constexpr Real EPSILON = 1e-9;

        // ========================= OBB-OBB Collision =========================

        /// @brief Detect collision between two oriented bounding boxes using SAT
        /// @param A First rigid body (must be a box)
        /// @param B Second rigid body (must be a box)
        /// @return CollisionInfo with contact details if collision exists
        static CollisionInfo DetectBoxBoxCollision(const RigidBodyBox& A, const RigidBodyBox& B)
        {
            CollisionInfo result;
            result.hasCollision = false;
            
            // Get world-space properties
            Vec3Cart posA = A.Position();
            Vec3Cart posB = B.Position();
            Vec3Cart centerVec = posB - posA;  // Vector from A center to B center
            
            // Get rotation matrices and axes
            MatrixNM<Real, 3, 3> RA = A.Orientation().ToRotationMatrix();
            MatrixNM<Real, 3, 3> RB = B.Orientation().ToRotationMatrix();
            
            // Extract axis vectors (columns of rotation matrices)
            Vec3Cart axisA[3], axisB[3];
            for (int i = 0; i < 3; i++)
            {
                axisA[i] = Vec3Cart(RA(0, i), RA(1, i), RA(2, i));
                axisB[i] = Vec3Cart(RB(0, i), RB(1, i), RB(2, i));
            }
            
            // Half-extents
            Real aExtent[3] = { A.HalfExtent(0), A.HalfExtent(1), A.HalfExtent(2) };
            Real bExtent[3] = { B.HalfExtent(0), B.HalfExtent(1), B.HalfExtent(2) };
            
            // Compute rotation matrix expressing B in A's coordinate frame: R = RA^T * RB
            // And absolute values for stability
            MatrixNM<Real, 3, 3> R, AbsR;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    R(i, j) = ScalarProduct(axisA[i], axisB[j]);
                    AbsR(i, j) = std::abs(R(i, j)) + EPSILON;
                }
            }
            
            // Translation in A's coordinate frame
            Real t[3];
            for (int i = 0; i < 3; i++)
                t[i] = ScalarProduct(centerVec, axisA[i]);
            
            // Track minimum penetration for contact info
            Real minPenetration = std::numeric_limits<Real>::max();
            int minPenetrationAxis = -1;
            Vec3Cart minPenetrationNormal;
            
            // ===== Test 15 separating axes =====
            
            // Test axes L = A0, A1, A2 (face normals of A)
            for (int i = 0; i < 3; i++)
            {
                Real ra = aExtent[i];
                Real rb = bExtent[0] * AbsR(i, 0) + bExtent[1] * AbsR(i, 1) + bExtent[2] * AbsR(i, 2);
                Real separation = std::abs(t[i]) - (ra + rb);
                
                if (separation > 0)
                    return result;  // Separating axis found, no collision
                
                if (-separation < minPenetration)
                {
                    minPenetration = -separation;
                    minPenetrationAxis = i;
                    minPenetrationNormal = axisA[i];
                    if (t[i] < 0) minPenetrationNormal = minPenetrationNormal * (-1.0);
                }
            }
            
            // Test axes L = B0, B1, B2 (face normals of B)
            for (int i = 0; i < 3; i++)
            {
                Real s = ScalarProduct(centerVec, axisB[i]);
                Real ra = aExtent[0] * AbsR(0, i) + aExtent[1] * AbsR(1, i) + aExtent[2] * AbsR(2, i);
                Real rb = bExtent[i];
                Real separation = std::abs(s) - (ra + rb);
                
                if (separation > 0)
                    return result;  // Separating axis found, no collision
                
                if (-separation < minPenetration)
                {
                    minPenetration = -separation;
                    minPenetrationAxis = 3 + i;
                    minPenetrationNormal = axisB[i];
                    if (s < 0) minPenetrationNormal = minPenetrationNormal * (-1.0);
                }
            }
            
            // Test 9 edge-edge axes: L = Ai × Bj
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    Vec3Cart axis = VectorProduct(axisA[i], axisB[j]);
                    Real len = axis.NormL2();
                    
                    // Skip parallel edges (degenerate axis)
                    if (len < EPSILON)
                        continue;
                    
                    axis = axis / len;  // Normalize
                    
                    // Compute projected radii
                    Real s = ScalarProduct(centerVec, axis);
                    
                    // Projected half-extents of A onto axis
                    Real ra = 0;
                    for (int k = 0; k < 3; k++)
                        ra += aExtent[k] * std::abs(ScalarProduct(axisA[k], axis));
                    
                    // Projected half-extents of B onto axis
                    Real rb = 0;
                    for (int k = 0; k < 3; k++)
                        rb += bExtent[k] * std::abs(ScalarProduct(axisB[k], axis));
                    
                    Real separation = std::abs(s) - (ra + rb);
                    
                    if (separation > 0)
                        return result;  // Separating axis found, no collision
                    
                    if (-separation < minPenetration)
                    {
                        minPenetration = -separation;
                        minPenetrationAxis = 6 + i * 3 + j;
                        minPenetrationNormal = axis;
                        if (s < 0) minPenetrationNormal = minPenetrationNormal * (-1.0);
                    }
                }
            }
            
            // All 15 axes overlap - collision detected!
            result.hasCollision = true;
            result.penetrationDepth = minPenetration;
            result.contactNormal = minPenetrationNormal;
            result.separatingAxisType = minPenetrationAxis;
            
            // Compute contact point (approximation: midpoint along line of minimum separation)
            // For face-face: use center - penetration/2 along normal
            // This is a simplified contact point - could be refined for edge-edge cases
            result.contactPoint = posA + (posB - posA) * 0.5;
            
            // Refine contact point based on collision type
            if (minPenetrationAxis < 3)
            {
                // Face of A - contact point on A's face
                result.contactPoint = posB - result.contactNormal * (minPenetration / 2.0);
            }
            else if (minPenetrationAxis < 6)
            {
                // Face of B - contact point on B's face
                result.contactPoint = posA + result.contactNormal * (minPenetration / 2.0);
            }
            // For edge-edge (axis >= 6), the midpoint approximation is reasonable
            
            return result;
        }

        // ========================= Sphere-Sphere Collision =========================

        /// @brief Detect collision between two spheres
        /// @details O(1) algorithm: checks if distance between centers < sum of radii
        /// @param A First sphere (type-safe)
        /// @param B Second sphere (type-safe)
        /// @return CollisionInfo with contact details if collision exists
        static CollisionInfo DetectSphereSphereCollision(const RigidBodySphere& A, const RigidBodySphere& B)
        {
            CollisionInfo result;
            result.hasCollision = false;
            
            // Get world-space centers
            Vec3Cart centerA = A.Position();
            Vec3Cart centerB = B.Position();
            
            // Vector from A to B
            Vec3Cart centerVec = centerB - centerA;
            Real distanceSquared = ScalarProduct(centerVec, centerVec);
            
            // Sum of radii
            Real radiusSum = A.Radius() + B.Radius();
            Real radiusSumSquared = radiusSum * radiusSum;
            
            // No collision if distance > sum of radii
            if (distanceSquared > radiusSumSquared)
            {
                return result;
            }
            
            // Collision detected!
            Real distance = std::sqrt(distanceSquared);
            
            // Handle degenerate case: centers coincide (should be very rare)
            if (distance < EPSILON)
            {
                // Use arbitrary normal (positive X)
                result.hasCollision = true;
                result.contactNormal = Vec3Cart(1.0, 0.0, 0.0);
                result.penetrationDepth = radiusSum;
                result.contactPoint = centerA;
                result.separatingAxisType = -1;  // Special case marker
                return result;
            }
            
            // Contact normal: points from A to B (normalized direction)
            result.hasCollision = true;
            result.contactNormal = centerVec * (1.0 / distance);
            result.penetrationDepth = radiusSum - distance;
            
            // Contact point: on the surface of A, toward B
            // Or equivalently: midpoint between the two surface points
            // Using midpoint approach for symmetry:
            //   Point on A's surface toward B: centerA + normalAB * radiusA
            //   Point on B's surface toward A: centerB - normalAB * radiusB
            //   Midpoint: (surfA + surfB) / 2
            Vec3Cart surfA = centerA + result.contactNormal * A.Radius();
            Vec3Cart surfB = centerB - result.contactNormal * B.Radius();
            result.contactPoint = (surfA + surfB) * 0.5;
            
            result.separatingAxisType = -2;  // Marker for sphere-sphere collision
            
            return result;
        }

        // ========================= Sphere-Wall Collision =========================

        /// @brief Detect collision between a sphere and axis-aligned container walls
        /// @details O(1) per wall: checks center distance to each wall plane vs radius
        /// @param sphere The sphere to test
        /// @param containerHalfSize Half-size of the cubic container (walls at ±containerHalfSize)
        /// @return Vector of all wall collisions (sphere can touch multiple walls simultaneously)
        static std::vector<WallCollisionInfo> DetectSphereWallCollisions(
            const RigidBodySphere& sphere, 
            Real containerHalfSize)
        {
            std::vector<WallCollisionInfo> collisions;
            
            Vec3Cart center = sphere.Position();
            Real radius = sphere.Radius();
            
            // Wall normals (pointing INTO the container)
            static const Vec3Cart wallNormals[6] = {
                Vec3Cart( 1,  0,  0),  // -X wall (at -containerHalfSize)
                Vec3Cart(-1,  0,  0),  // +X wall (at +containerHalfSize)
                Vec3Cart( 0,  1,  0),  // -Y wall
                Vec3Cart( 0, -1,  0),  // +Y wall
                Vec3Cart( 0,  0,  1),  // -Z wall
                Vec3Cart( 0,  0, -1)   // +Z wall
            };
            
            // Sphere center coordinates
            Real coords[3] = { center.X(), center.Y(), center.Z() };
            
            // Check each of the 6 walls
            for (int wall = 0; wall < 6; wall++)
            {
                int axis = wall / 2;  // 0=X, 1=Y, 2=Z
                bool isNegativeWall = (wall % 2 == 0);
                
                Real penetration = 0.0;
                
                if (isNegativeWall)
                {
                    // Negative wall: at -containerHalfSize
                    // Sphere penetrates if center - radius < -containerHalfSize
                    // Penetration = -containerHalfSize - (center - radius) = -containerHalfSize - center + radius
                    penetration = -containerHalfSize - coords[axis] + radius;
                }
                else
                {
                    // Positive wall: at +containerHalfSize
                    // Sphere penetrates if center + radius > containerHalfSize
                    // Penetration = (center + radius) - containerHalfSize
                    penetration = coords[axis] + radius - containerHalfSize;
                }
                
                if (penetration > 0)
                {
                    WallCollisionInfo info;
                    info.hasCollision = true;
                    info.wallIndex = wall;
                    info.contactNormal = wallNormals[wall];
                    info.penetrationDepth = penetration;
                    
                    // Contact point: on sphere surface, at the wall
                    // This is center projected onto wall + offset along normal
                    Vec3Cart contactPoint = center;
                    if (isNegativeWall)
                    {
                        // Project center's axis coordinate to wall, add radius toward center
                        // Contact is on sphere surface: center - normal * radius
                        contactPoint = center - wallNormals[wall] * radius;
                        // Actually for negative wall with normal pointing IN (+),
                        // the contact on sphere is center - normal * radius
                    }
                    else
                    {
                        // For positive wall with normal pointing IN (-),
                        // contact on sphere is center - normal * radius (but normal is negative)
                        contactPoint = center - wallNormals[wall] * radius;
                    }
                    info.contactPoint = contactPoint;
                    
                    collisions.push_back(info);
                }
            }
            
            return collisions;
        }

        // ========================= Sphere-Box Collision =========================

        /// @brief Detect collision between a sphere and an oriented bounding box
        /// @details Algorithm from Ericson "Real-Time Collision Detection" Ch. 5.2.4:
        ///          1. Transform sphere center to box's local coordinate frame
        ///          2. Clamp to box extents to find closest point on box
        ///          3. Transform closest point back to world space
        ///          4. Check distance from closest point to sphere center vs radius
        /// 
        /// This is the most geometrically complex of the collision primitives.
        /// Time complexity: O(1)
        ///
        /// @param sphere The sphere rigid body
        /// @param box The box (OBB) rigid body
        /// @return CollisionInfo with contact details if collision exists
        static CollisionInfo DetectSphereBoxCollision(const RigidBodySphere& sphere, const RigidBodyBox& box)
        {
            CollisionInfo result;
            result.hasCollision = false;
            
            // Get sphere properties
            Vec3Cart sphereCenter = sphere.Position();
            Real radius = sphere.Radius();
            
            // Get box properties
            Vec3Cart boxCenter = box.Position();
            MatrixNM<Real, 3, 3> boxRotation = box.Orientation().ToRotationMatrix();
            
            // Extract box axes (columns of rotation matrix)
            Vec3Cart boxAxisX(boxRotation(0, 0), boxRotation(1, 0), boxRotation(2, 0));
            Vec3Cart boxAxisY(boxRotation(0, 1), boxRotation(1, 1), boxRotation(2, 1));
            Vec3Cart boxAxisZ(boxRotation(0, 2), boxRotation(1, 2), boxRotation(2, 2));
            
            // Box half-extents
            Real halfExtentX = box.HalfExtent(0);
            Real halfExtentY = box.HalfExtent(1);
            Real halfExtentZ = box.HalfExtent(2);
            
            // Vector from box center to sphere center
            Vec3Cart d = sphereCenter - boxCenter;
            
            // Transform sphere center to box's local coordinate frame
            // by projecting onto each axis
            Real localX = ScalarProduct(d, boxAxisX);
            Real localY = ScalarProduct(d, boxAxisY);
            Real localZ = ScalarProduct(d, boxAxisZ);
            
            // Find closest point on box by clamping to extents
            Real closestX = std::clamp(localX, -halfExtentX, halfExtentX);
            Real closestY = std::clamp(localY, -halfExtentY, halfExtentY);
            Real closestZ = std::clamp(localZ, -halfExtentZ, halfExtentZ);
            
            // Transform closest point back to world space
            Vec3Cart closestPoint = boxCenter 
                + boxAxisX * closestX 
                + boxAxisY * closestY 
                + boxAxisZ * closestZ;
            
            // Vector from closest point on box to sphere center
            Vec3Cart toSphere = sphereCenter - closestPoint;
            Real distanceSquared = ScalarProduct(toSphere, toSphere);
            Real radiusSquared = radius * radius;
            
            // No collision if distance > radius
            if (distanceSquared > radiusSquared)
            {
                return result;
            }
            
            // Collision detected!
            Real distance = std::sqrt(distanceSquared);
            
            // Handle special case: sphere center is inside the box
            // (closest point equals sphere center)
            if (distance < EPSILON)
            {
                // Sphere center is inside the box
                // Find the closest face of the box and use its normal
                
                // Calculate signed distances to each face in local coordinates
                Real distToFaces[6] = {
                    halfExtentX - localX,   // +X face
                    halfExtentX + localX,   // -X face
                    halfExtentY - localY,   // +Y face
                    halfExtentY + localY,   // -Y face
                    halfExtentZ - localZ,   // +Z face
                    halfExtentZ + localZ    // -Z face
                };
                
                // Find minimum distance (closest face)
                int minFace = 0;
                Real minDist = distToFaces[0];
                for (int i = 1; i < 6; i++)
                {
                    if (distToFaces[i] < minDist)
                    {
                        minDist = distToFaces[i];
                        minFace = i;
                    }
                }
                
                // Set contact normal based on closest face
                // Face normals point OUTWARD from box, but we need normal from A (sphere) to B (box)
                // So we negate them to point INTO the box (from sphere toward box center)
                Vec3Cart faceNormals[6] = {
                    boxAxisX * (-1.0),      // +X face -> point toward -X (into box)
                    boxAxisX,               // -X face -> point toward +X (into box)
                    boxAxisY * (-1.0),      // +Y face -> point toward -Y (into box)
                    boxAxisY,               // -Y face -> point toward +Y (into box)
                    boxAxisZ * (-1.0),      // +Z face -> point toward -Z (into box)
                    boxAxisZ                // -Z face -> point toward +Z (into box)
                };
                
                result.hasCollision = true;
                result.contactNormal = faceNormals[minFace];
                result.penetrationDepth = radius + minDist;  // radius + distance to face
                result.contactPoint = sphereCenter + result.contactNormal * (result.penetrationDepth / 2.0);
                result.separatingAxisType = -3;  // Marker for sphere-inside-box
                
                return result;
            }
            
            // Normal case: sphere center is outside box
            result.hasCollision = true;
            // Contact normal must point from A (sphere) to B (box) for collision response
            // toSphere points from closestPoint on box TO sphere, so negate it
            result.contactNormal = toSphere * (-1.0 / distance);  // Normalized, points from sphere to box
            result.penetrationDepth = radius - distance;
            
            // Contact point: midway between sphere surface and box surface
            result.contactPoint = closestPoint + toSphere * (0.5 * result.penetrationDepth / distance);
            
            result.separatingAxisType = -4;  // Marker for sphere-box collision
            
            return result;
        }

        // ========================= Box-Wall Collision =========================

        /// @brief Detect collision between a box and axis-aligned container walls
        /// @param box The box rigid body to test
        /// @param containerHalfSize Half-size of the cubic container (walls at ±containerHalfSize)
        /// @return Vector of all wall collisions (box can touch multiple walls)
        static std::vector<WallCollisionInfo> DetectBoxWallCollisions(
            const RigidBodyBox& box, 
            Real containerHalfSize)
        {
            std::vector<WallCollisionInfo> collisions;
            
            // Get all 8 vertices of the box in world coordinates
            std::vector<Vec3Cart> vertices = box.GetWorldVertices();
            
            // Wall normals (pointing INTO the container)
            static const Vec3Cart wallNormals[6] = {
                Vec3Cart( 1,  0,  0),  // -X wall (at -containerHalfSize)
                Vec3Cart(-1,  0,  0),  // +X wall (at +containerHalfSize)
                Vec3Cart( 0,  1,  0),  // -Y wall
                Vec3Cart( 0, -1,  0),  // +Y wall
                Vec3Cart( 0,  0,  1),  // -Z wall
                Vec3Cart( 0,  0, -1)   // +Z wall
            };
            
            // For each wall, check all vertices
            for (int wall = 0; wall < 6; wall++)
            {
                WallCollisionInfo info;
                info.hasCollision = false;
                info.wallIndex = wall;
                info.contactNormal = wallNormals[wall];
                info.penetrationDepth = 0.0;
                
                Real wallPos = (wall % 2 == 0) ? -containerHalfSize : containerHalfSize;
                int axis = wall / 2;  // 0=X, 1=Y, 2=Z
                
                // Find deepest penetrating vertex for this wall
                for (const auto& vertex : vertices)
                {
                    Real coord = (axis == 0) ? vertex.X() : ((axis == 1) ? vertex.Y() : vertex.Z());
                    Real penetration = 0.0;
                    
                    if (wall % 2 == 0)
                    {
                        // Negative wall (at -containerHalfSize)
                        penetration = -containerHalfSize - coord;
                    }
                    else
                    {
                        // Positive wall (at +containerHalfSize)
                        penetration = coord - containerHalfSize;
                    }
                    
                    if (penetration > info.penetrationDepth)
                    {
                        info.hasCollision = true;
                        info.penetrationDepth = penetration;
                        info.contactPoint = vertex;
                    }
                }
                
                if (info.hasCollision)
                    collisions.push_back(info);
            }
            
            return collisions;
        }

        /// @brief Check if a single vertex is outside the container
        /// @param vertex The vertex position
        /// @param containerHalfSize Half-size of the cubic container
        /// @return Wall collision info (hasCollision will be true if outside)
        static WallCollisionInfo CheckVertexWallCollision(
            const Vec3Cart& vertex,
            Real containerHalfSize)
        {
            WallCollisionInfo result;
            result.hasCollision = false;
            result.penetrationDepth = 0.0;
            
            // Check each axis
            const Real coords[3] = { vertex.X(), vertex.Y(), vertex.Z() };
            const Vec3Cart normals[6] = {
                Vec3Cart( 1,  0,  0), Vec3Cart(-1,  0,  0),
                Vec3Cart( 0,  1,  0), Vec3Cart( 0, -1,  0),
                Vec3Cart( 0,  0,  1), Vec3Cart( 0,  0, -1)
            };
            
            for (int axis = 0; axis < 3; axis++)
            {
                // Check negative wall
                Real penNeg = -containerHalfSize - coords[axis];
                if (penNeg > result.penetrationDepth)
                {
                    result.hasCollision = true;
                    result.penetrationDepth = penNeg;
                    result.contactPoint = vertex;
                    result.contactNormal = normals[axis * 2];
                    result.wallIndex = axis * 2;
                }
                
                // Check positive wall
                Real penPos = coords[axis] - containerHalfSize;
                if (penPos > result.penetrationDepth)
                {
                    result.hasCollision = true;
                    result.penetrationDepth = penPos;
                    result.contactPoint = vertex;
                    result.contactNormal = normals[axis * 2 + 1];
                    result.wallIndex = axis * 2 + 1;
                }
            }
            
            return result;
        }

        // ========================= UNIFIED COLLISION DISPATCHER =========================
        
        /// @brief Detect collision between any two rigid bodies (unified dispatcher)
        /// @details Automatically selects the correct algorithm based on shape types:
        ///          - Box-Box: SAT algorithm (15 axes)
        ///          - Sphere-Sphere: Distance-based (O(1))
        ///          - Sphere-Box: Closest point on OBB (Ericson)
        ///
        /// This is the RECOMMENDED entry point for body-body collision detection
        /// as it handles all shape type combinations transparently.
        ///
        /// @param A First rigid body
        /// @param B Second rigid body
        /// @return CollisionInfo with contact details if collision exists
        static CollisionInfo DetectBodyBodyCollision(const RigidBody& A, const RigidBody& B)
        {
            // Dispatch based on shape types using type-safe casting
            const auto* boxA = AsBox(&A);
            const auto* boxB = AsBox(&B);
            const auto* sphereA = AsSphere(&A);
            const auto* sphereB = AsSphere(&B);
            
            if (boxA && boxB)
            {
                return DetectBoxBoxCollision(*boxA, *boxB);
            }
            else if (sphereA && sphereB)
            {
                return DetectSphereSphereCollision(*sphereA, *sphereB);
            }
            else if (sphereA && boxB)
            {
                return DetectSphereBoxCollision(*sphereA, *boxB);
            }
            else if (boxA && sphereB)
            {
                // Flip the order: sphere first, box second
                // Then flip the contact normal
                CollisionInfo result = DetectSphereBoxCollision(*sphereB, *boxA);
                if (result.hasCollision)
                {
                    result.contactNormal = result.contactNormal * (-1.0);  // Reverse direction
                }
                return result;
            }
            
            // Unknown shape combination (shouldn't happen with current types)
            CollisionInfo result;
            result.hasCollision = false;
            return result;
        }

        /// @brief Detect collision between a body and container walls (unified dispatcher)
        /// @details Automatically selects the correct algorithm based on shape type:
        ///          - Box: Vertex-based detection
        ///          - Sphere: Center distance detection (O(1))
        ///
        /// This is the RECOMMENDED entry point for body-wall collision detection.
        ///
        /// @param body The rigid body to test
        /// @param containerHalfSize Half-size of the cubic container
        /// @return Vector of all wall collisions
        static std::vector<WallCollisionInfo> DetectBodyWallCollisions(
            const RigidBody& body, 
            Real containerHalfSize)
        {
            // Type-safe dispatch
            if (const auto* box = AsBox(&body))
            {
                return DetectBoxWallCollisions(*box, containerHalfSize);
            }
            else if (const auto* sphere = AsSphere(&body))
            {
                return DetectSphereWallCollisions(*sphere, containerHalfSize);
            }
            
            // Unknown shape type (shouldn't happen)
            return {};
        }
    };

} // namespace MPL

#endif // MPL_COLLISION_DETECTION_H
