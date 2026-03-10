#ifndef MML_CONVEX_HULL_3D_H
#define MML_CONVEX_HULL_3D_H

/////////////////////////////////////////////////////////////////////////////////////////
///  @file      ConvexHull3D.h
///  @brief     3D Convex Hull computation using incremental algorithm
///  @details   Part of the MML Computational Geometry module
/////////////////////////////////////////////////////////////////////////////////////////

#include "mml/MMLBase.h"
#include "mml/base/Geometry/Geometry3D.h"
#include "mml/base/Vector/VectorN.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>

namespace MML::CompGeometry {

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief 3D Convex Hull structure
/// 
/// Represents the convex hull of a set of 3D points as a collection of
/// triangular faces. Each face is defined by three vertex indices with
/// outward-pointing normal (CCW winding when viewed from outside).
//////////////////////////////////////////////////////////////////////////////////////////
struct ConvexHull3D {
	std::vector<Point3Cartesian> vertices;           ///< Hull vertices
	std::vector<std::array<int, 3>> faces;           ///< Triangular faces (vertex indices, CCW from outside)
	std::vector<std::array<int, 3>> neighbors;       ///< Adjacent face indices for each face edge

	int NumVertices() const { return static_cast<int>(vertices.size()); }
	int NumFaces() const { return static_cast<int>(faces.size()); }

	/// @brief Check if a point is inside or on the convex hull
	bool Contains(const Point3Cartesian& p, Real tolerance = 1e-10) const {
		// Point is inside if it's on the negative side of all face planes
		for (int f = 0; f < NumFaces(); f++) {
			const auto& face = faces[f];
			const Point3Cartesian& a = vertices[face[0]];
			const Point3Cartesian& b = vertices[face[1]];
			const Point3Cartesian& c = vertices[face[2]];
			
			// Compute face normal (outward pointing)
			Vector3Cartesian ab(b.X() - a.X(), b.Y() - a.Y(), b.Z() - a.Z());
			Vector3Cartesian ac(c.X() - a.X(), c.Y() - a.Y(), c.Z() - a.Z());
			Vector3Cartesian n = VectorProduct(ab, ac);
			
			// Signed distance from plane
			Real dist = n.X() * (p.X() - a.X()) + 
			            n.Y() * (p.Y() - a.Y()) + 
			            n.Z() * (p.Z() - a.Z());
			
			if (dist > tolerance) return false;  // Outside this face
		}
		return true;
	}

	/// @brief Compute the volume of the convex hull
	Real Volume() const {
		if (NumFaces() == 0) return 0.0;
		
		// Use signed tetrahedra with origin
		Real volume = 0.0;
		for (const auto& face : faces) {
			const Point3Cartesian& a = vertices[face[0]];
			const Point3Cartesian& b = vertices[face[1]];
			const Point3Cartesian& c = vertices[face[2]];
			
			// Signed volume of tetrahedron with origin
			volume += (a.X() * (b.Y() * c.Z() - b.Z() * c.Y()) +
			           a.Y() * (b.Z() * c.X() - b.X() * c.Z()) +
			           a.Z() * (b.X() * c.Y() - b.Y() * c.X())) / 6.0;
		}
		return std::abs(volume);
	}

	/// @brief Compute the surface area of the convex hull
	Real SurfaceArea() const {
		Real area = 0.0;
		for (const auto& face : faces) {
			const Point3Cartesian& a = vertices[face[0]];
			const Point3Cartesian& b = vertices[face[1]];
			const Point3Cartesian& c = vertices[face[2]];
			
			Vector3Cartesian ab(b.X() - a.X(), b.Y() - a.Y(), b.Z() - a.Z());
			Vector3Cartesian ac(c.X() - a.X(), c.Y() - a.Y(), c.Z() - a.Z());
			Vector3Cartesian cross = VectorProduct(ab, ac);
			
			area += 0.5 * std::sqrt(cross.X() * cross.X() + 
			                        cross.Y() * cross.Y() + 
			                        cross.Z() * cross.Z());
		}
		return area;
	}

	/// @brief Get a face as Triangle3D
	Triangle3D GetFace(int index) const {
		if (index < 0 || index >= NumFaces()) {
			return Triangle3D(Point3Cartesian(0,0,0), Point3Cartesian(1,0,0), Point3Cartesian(0,1,0));
		}
		return Triangle3D(vertices[faces[index][0]], 
		                  vertices[faces[index][1]], 
		                  vertices[faces[index][2]]);
	}

	/// @brief Get the centroid of the convex hull
	Point3Cartesian Centroid() const {
		if (vertices.empty()) return Point3Cartesian(0, 0, 0);
		
		Real cx = 0, cy = 0, cz = 0;
		for (const auto& v : vertices) {
			cx += v.X();
			cy += v.Y();
			cz += v.Z();
		}
		int n = NumVertices();
		return Point3Cartesian(cx / n, cy / n, cz / n);
	}
};

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief 3D Convex Hull computation
/// @details Uses an incremental algorithm: O(n²) worst case, typically better
//////////////////////////////////////////////////////////////////////////////////////////
class ConvexHull3DComputer {
	// Use centralized geometry epsilon from Constants
	static constexpr Real EPSILON = Constants::GEOMETRY_EPSILON;

public:
	/// @brief Compute 3D convex hull using incremental algorithm
	/// 
	/// This implementation uses an incremental approach:
	/// 1. Start with an initial tetrahedron from 4 non-coplanar points
	/// 2. For each remaining point:
	///    - If inside hull, skip
	///    - If outside, find visible faces and replace with new faces
	/// 
	/// Time complexity: O(n²) worst case, typically better
	/// 
	/// @param points Input points
	/// @return ConvexHull3D structure with vertices and faces
	static ConvexHull3D Compute(std::vector<Point3Cartesian> points) {
		ConvexHull3D hull;
		
		int n = static_cast<int>(points.size());
		if (n < 4) {
			// Degenerate cases
			hull.vertices = points;
			if (n == 3) {
				hull.faces.push_back({0, 1, 2});
				hull.faces.push_back({0, 2, 1});  // Both orientations
			}
			return hull;
		}
		
		// Find extreme points to form initial tetrahedron
		int minX = 0, maxX = 0, minY = 0, maxY = 0, minZ = 0, maxZ = 0;
		for (int i = 1; i < n; i++) {
			if (points[i].X() < points[minX].X()) minX = i;
			if (points[i].X() > points[maxX].X()) maxX = i;
			if (points[i].Y() < points[minY].Y()) minY = i;
			if (points[i].Y() > points[maxY].Y()) maxY = i;
			if (points[i].Z() < points[minZ].Z()) minZ = i;
			if (points[i].Z() > points[maxZ].Z()) maxZ = i;
		}
		
		// Find the two most distant extreme points
		std::vector<int> extremes = {minX, maxX, minY, maxY, minZ, maxZ};
		int p0 = -1, p1 = -1;
		Real maxDist2 = -1;
		for (size_t i = 0; i < extremes.size(); i++) {
			for (size_t j = i + 1; j < extremes.size(); j++) {
				Real d2 = Distance2(points[extremes[i]], points[extremes[j]]);
				if (d2 > maxDist2) {
					maxDist2 = d2;
					p0 = extremes[i];
					p1 = extremes[j];
				}
			}
		}
		
		if (maxDist2 < 1e-20) {
			// All points coincident
			hull.vertices.push_back(points[0]);
			return hull;
		}
		
		// Find point furthest from line p0-p1
		int p2 = -1;
		Real maxLineDist2 = -1;
		Vector3Cartesian lineDir(points[p1].X() - points[p0].X(),
		                         points[p1].Y() - points[p0].Y(),
		                         points[p1].Z() - points[p0].Z());
		Real lineLen2 = lineDir.X() * lineDir.X() + lineDir.Y() * lineDir.Y() + lineDir.Z() * lineDir.Z();
		
		for (int i = 0; i < n; i++) {
			if (i == p0 || i == p1) continue;
			Real d2 = PointToLineDistance2(points[i], points[p0], lineDir, lineLen2);
			if (d2 > maxLineDist2) {
				maxLineDist2 = d2;
				p2 = i;
			}
		}
		
		if (p2 < 0 || maxLineDist2 < 1e-20) {
			// Points are collinear
			hull.vertices.push_back(points[p0]);
			hull.vertices.push_back(points[p1]);
			return hull;
		}
		
		// Find point furthest from plane p0-p1-p2
		Vector3Cartesian v01(points[p1].X() - points[p0].X(),
		                     points[p1].Y() - points[p0].Y(),
		                     points[p1].Z() - points[p0].Z());
		Vector3Cartesian v02(points[p2].X() - points[p0].X(),
		                     points[p2].Y() - points[p0].Y(),
		                     points[p2].Z() - points[p0].Z());
		Vector3Cartesian planeNormal = VectorProduct(v01, v02);
		
		int p3 = -1;
		Real maxPlaneDist = 0;
		for (int i = 0; i < n; i++) {
			if (i == p0 || i == p1 || i == p2) continue;
			Real d = std::abs(planeNormal.X() * (points[i].X() - points[p0].X()) +
			                  planeNormal.Y() * (points[i].Y() - points[p0].Y()) +
			                  planeNormal.Z() * (points[i].Z() - points[p0].Z()));
			if (d > maxPlaneDist) {
				maxPlaneDist = d;
				p3 = i;
			}
		}
		
		if (p3 < 0 || maxPlaneDist < 1e-14) {
			// Points are coplanar - return 2D hull as degenerate 3D hull
			hull.vertices.push_back(points[p0]);
			hull.vertices.push_back(points[p1]);
			hull.vertices.push_back(points[p2]);
			hull.faces.push_back({0, 1, 2});
			return hull;
		}
		
		// Build initial tetrahedron
		// Map original indices to hull vertex indices
		std::map<int, int> pointToVertex;
		pointToVertex[p0] = 0; hull.vertices.push_back(points[p0]);
		pointToVertex[p1] = 1; hull.vertices.push_back(points[p1]);
		pointToVertex[p2] = 2; hull.vertices.push_back(points[p2]);
		pointToVertex[p3] = 3; hull.vertices.push_back(points[p3]);
		
		// Determine orientation: is p3 above or below plane p0-p1-p2?
		Real signedDist = planeNormal.X() * (points[p3].X() - points[p0].X()) +
		                  planeNormal.Y() * (points[p3].Y() - points[p0].Y()) +
		                  planeNormal.Z() * (points[p3].Z() - points[p0].Z());
		
		// Create faces with outward normals (CCW from outside)
		if (signedDist > 0) {
			// p3 is above plane p0-p1-p2
			hull.faces.push_back({0, 2, 1});  // Bottom face (looking from below)
			hull.faces.push_back({0, 1, 3});  // Side faces
			hull.faces.push_back({1, 2, 3});
			hull.faces.push_back({2, 0, 3});
		} else {
			// p3 is below plane p0-p1-p2
			hull.faces.push_back({0, 1, 2});  // Top face
			hull.faces.push_back({0, 3, 1});  // Side faces
			hull.faces.push_back({1, 3, 2});
			hull.faces.push_back({2, 3, 0});
		}
		
		// Assign remaining points to faces they're outside of
		struct FaceData {
			std::array<int, 3> vertices;
			Vector3Cartesian normal;
			Point3Cartesian centroid;
			std::vector<int> outsidePoints;  // Original point indices
			bool valid = true;
		};
		
		auto computeFaceData = [&hull](const std::array<int, 3>& face) -> FaceData {
			FaceData fd;
			fd.vertices = face;
			const Point3Cartesian& a = hull.vertices[face[0]];
			const Point3Cartesian& b = hull.vertices[face[1]];
			const Point3Cartesian& c = hull.vertices[face[2]];
			
			Vector3Cartesian ab(b.X() - a.X(), b.Y() - a.Y(), b.Z() - a.Z());
			Vector3Cartesian ac(c.X() - a.X(), c.Y() - a.Y(), c.Z() - a.Z());
			fd.normal = VectorProduct(ab, ac);
			
			// Normalize
			Real len = std::sqrt(fd.normal.X() * fd.normal.X() + 
			                     fd.normal.Y() * fd.normal.Y() + 
			                     fd.normal.Z() * fd.normal.Z());
			if (len > 1e-14) {
				fd.normal = Vector3Cartesian(fd.normal.X() / len, fd.normal.Y() / len, fd.normal.Z() / len);
			}
			
			fd.centroid = Point3Cartesian((a.X() + b.X() + c.X()) / 3.0,
			                              (a.Y() + b.Y() + c.Y()) / 3.0,
			                              (a.Z() + b.Z() + c.Z()) / 3.0);
			return fd;
		};
		
		std::vector<FaceData> faceDataList;
		for (const auto& face : hull.faces) {
			faceDataList.push_back(computeFaceData(face));
		}
		
		// Assign points to faces
		std::set<int> usedPoints = {p0, p1, p2, p3};
		for (int i = 0; i < n; i++) {
			if (usedPoints.count(i)) continue;
			
			for (auto& fd : faceDataList) {
				const Point3Cartesian& a = hull.vertices[fd.vertices[0]];
				Real dist = fd.normal.X() * (points[i].X() - a.X()) +
				            fd.normal.Y() * (points[i].Y() - a.Y()) +
				            fd.normal.Z() * (points[i].Z() - a.Z());
				if (dist > 1e-10) {
					fd.outsidePoints.push_back(i);
					break;  // Each point assigned to at most one face
				}
			}
		}
		
		// Process faces with outside points
		int maxIterations = n * n;  // Safety limit
		int iterations = 0;
		while (iterations++ < maxIterations) {
			// Find a face with outside points
			int faceIdx = -1;
			for (size_t i = 0; i < faceDataList.size(); i++) {
				if (faceDataList[i].valid && !faceDataList[i].outsidePoints.empty()) {
					faceIdx = static_cast<int>(i);
					break;
				}
			}
			
			if (faceIdx < 0) break;  // No more faces with outside points
			
			FaceData& currentFace = faceDataList[faceIdx];
			
			// Find furthest outside point
			int furthestPt = -1;
			Real maxDist = 0;
			const Point3Cartesian& a = hull.vertices[currentFace.vertices[0]];
			for (int pt : currentFace.outsidePoints) {
				Real dist = currentFace.normal.X() * (points[pt].X() - a.X()) +
				            currentFace.normal.Y() * (points[pt].Y() - a.Y()) +
				            currentFace.normal.Z() * (points[pt].Z() - a.Z());
				if (dist > maxDist) {
					maxDist = dist;
					furthestPt = pt;
				}
			}
			
			if (furthestPt < 0) {
				currentFace.outsidePoints.clear();
				continue;
			}
			
			// Find all faces visible from furthestPt
			std::vector<int> visibleFaces;
			for (size_t i = 0; i < faceDataList.size(); i++) {
				if (!faceDataList[i].valid) continue;
				
				const Point3Cartesian& fa = hull.vertices[faceDataList[i].vertices[0]];
				Real dist = faceDataList[i].normal.X() * (points[furthestPt].X() - fa.X()) +
				            faceDataList[i].normal.Y() * (points[furthestPt].Y() - fa.Y()) +
				            faceDataList[i].normal.Z() * (points[furthestPt].Z() - fa.Z());
				if (dist > 1e-10) {
					visibleFaces.push_back(static_cast<int>(i));
				}
			}
			
			if (visibleFaces.empty()) {
				currentFace.outsidePoints.erase(
					std::remove(currentFace.outsidePoints.begin(), currentFace.outsidePoints.end(), furthestPt),
					currentFace.outsidePoints.end());
				continue;
			}
			
			// Collect orphaned points from visible faces
			std::vector<int> orphanedPoints;
			for (int vf : visibleFaces) {
				for (int pt : faceDataList[vf].outsidePoints) {
					if (pt != furthestPt) {
						orphanedPoints.push_back(pt);
					}
				}
			}
			
			// Find horizon edges
			std::map<std::pair<int,int>, int> edgeCount;
			std::map<std::pair<int,int>, std::pair<int,int>> edgeOrientation;
			
			for (int vf : visibleFaces) {
				const auto& face = faceDataList[vf].vertices;
				for (int e = 0; e < 3; e++) {
					int v1 = face[e];
					int v2 = face[(e + 1) % 3];
					auto canonicalEdge = (v1 < v2) ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
					edgeCount[canonicalEdge]++;
					if (edgeCount[canonicalEdge] == 1) {
						edgeOrientation[canonicalEdge] = {v1, v2};
					}
				}
			}
			
			std::vector<std::pair<int, int>> horizonEdges;
			for (const auto& [edge, count] : edgeCount) {
				if (count == 1) {
					auto [v1, v2] = edgeOrientation[edge];
					horizonEdges.push_back({v2, v1});
				}
			}
			
			if (horizonEdges.size() < 3) {
				currentFace.outsidePoints.erase(
					std::remove(currentFace.outsidePoints.begin(), currentFace.outsidePoints.end(), furthestPt),
					currentFace.outsidePoints.end());
				continue;
			}
			
			// Add new vertex
			int newVertexIdx = static_cast<int>(hull.vertices.size());
			hull.vertices.push_back(points[furthestPt]);
			usedPoints.insert(furthestPt);
			
			// Mark visible faces as invalid
			for (int vf : visibleFaces) {
				faceDataList[vf].valid = false;
			}
			
			// Compute hull centroid for orientation
			Point3Cartesian hullCentroid(0, 0, 0);
			int validCount = 0;
			for (const auto& fd : faceDataList) {
				if (fd.valid) {
					hullCentroid = Point3Cartesian(
						hullCentroid.X() + fd.centroid.X(),
						hullCentroid.Y() + fd.centroid.Y(),
						hullCentroid.Z() + fd.centroid.Z());
					validCount++;
				}
			}
			if (validCount > 0) {
				hullCentroid = Point3Cartesian(
					hullCentroid.X() / validCount,
					hullCentroid.Y() / validCount,
					hullCentroid.Z() / validCount);
			}
			
			// Create new faces from horizon edges
			std::vector<int> newFaceIndices;
			for (const auto& edge : horizonEdges) {
				std::array<int, 3> newFace = {edge.first, edge.second, newVertexIdx};
				FaceData fd = computeFaceData(newFace);
				
				// Verify orientation
				const Point3Cartesian& faceVertex = hull.vertices[newFace[0]];
				Real towardCentroid = fd.normal.X() * (hullCentroid.X() - faceVertex.X()) +
				                      fd.normal.Y() * (hullCentroid.Y() - faceVertex.Y()) +
				                      fd.normal.Z() * (hullCentroid.Z() - faceVertex.Z());
				
				if (towardCentroid > 0) {
					newFace = {edge.second, edge.first, newVertexIdx};
					fd = computeFaceData(newFace);
				}
				
				fd.valid = true;
				
				int newFaceIdx = static_cast<int>(faceDataList.size());
				faceDataList.push_back(fd);
				newFaceIndices.push_back(newFaceIdx);
			}
			
			// Reassign orphaned points
			for (int pt : orphanedPoints) {
				for (int nfi : newFaceIndices) {
					FaceData& fd = faceDataList[nfi];
					const Point3Cartesian& fa = hull.vertices[fd.vertices[0]];
					Real dist = fd.normal.X() * (points[pt].X() - fa.X()) +
					            fd.normal.Y() * (points[pt].Y() - fa.Y()) +
					            fd.normal.Z() * (points[pt].Z() - fa.Z());
					if (dist > 1e-10) {
						fd.outsidePoints.push_back(pt);
						break;
					}
				}
			}
		}
		
		// Collect valid faces
		hull.faces.clear();
		for (const auto& fd : faceDataList) {
			if (fd.valid) {
				hull.faces.push_back(fd.vertices);
			}
		}
		
		return hull;
	}

private:
	/// @brief Squared distance between two 3D points
	static Real Distance2(const Point3Cartesian& a, const Point3Cartesian& b) {
		Real dx = b.X() - a.X();
		Real dy = b.Y() - a.Y();
		Real dz = b.Z() - a.Z();
		return dx * dx + dy * dy + dz * dz;
	}

	/// @brief Squared distance from point to line
	static Real PointToLineDistance2(const Point3Cartesian& p, 
	                                  const Point3Cartesian& linePoint,
	                                  const Vector3Cartesian& lineDir,
	                                  Real lineDirLen2) {
		Vector3Cartesian v(p.X() - linePoint.X(), p.Y() - linePoint.Y(), p.Z() - linePoint.Z());
		Vector3Cartesian cross = VectorProduct(v, lineDir);
		Real crossLen2 = cross.X() * cross.X() + cross.Y() * cross.Y() + cross.Z() * cross.Z();
		return crossLen2 / lineDirLen2;
	}
};

} // namespace MML::CompGeometry

#endif // MML_CONVEX_HULL_3D_H
