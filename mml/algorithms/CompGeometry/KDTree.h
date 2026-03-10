///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        KDTree.h                                                            ///
///  Description: KD-Tree spatial data structures for efficient neighbor queries     ///
///               Supports 2D, 3D, and N-dimensional point sets                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_KDTREE_H
#define MML_KDTREE_H

#include "mml/MMLBase.h"
#include "mml/core/AlgorithmTypes.h"
#include "mml/base/Vector/Vector.h"
#include "mml/base/Geometry/Geometry2D.h"
#include "mml/base/Geometry/Geometry3D.h"

#include <vector>
#include <memory>
#include <algorithm>
#include <queue>
#include <cmath>

namespace MML {

	/******************************************************************************/
	/*****                    Configuration                                   *****/
	/******************************************************************************/

	/// @brief Configuration for KD-tree construction
	struct KDTreeConfig {
		int max_leaf_size = 10;           ///< Maximum points per leaf node before splitting
		bool balanced = true;             ///< Use median split (balanced) vs midpoint split
		
		/// Default configuration
		static KDTreeConfig Default() { return {}; }
		
		/// High performance configuration (larger leaves, faster build)
		static KDTreeConfig HighPerformance() {
			KDTreeConfig cfg;
			cfg.max_leaf_size = 32;
			return cfg;
		}
		
		/// Memory efficient configuration (smaller leaves, slower build)
		static KDTreeConfig LowMemory() {
			KDTreeConfig cfg;
			cfg.max_leaf_size = 4;
			return cfg;
		}
	};

	/******************************************************************************/
	/*****                    Statistics                                      *****/
	/******************************************************************************/

	/// @brief Statistics about KD-tree construction
	struct KDTreeStats {
		size_t num_points = 0;            ///< Total points in tree
		size_t num_nodes = 0;             ///< Total nodes (internal + leaf)
		size_t num_leaf_nodes = 0;        ///< Number of leaf nodes
		size_t tree_depth = 0;            ///< Maximum depth of tree
		double build_time_ms = 0.0;       ///< Time to build tree
	};

	/******************************************************************************/
	/*****                    Result Structures                               *****/
	/******************************************************************************/

	/// @brief Result of nearest neighbor query
	struct NearestNeighborResult {
		size_t index = SIZE_MAX;          ///< Index of nearest point (SIZE_MAX if not found)
		Real distance = 0.0;              ///< Distance to nearest point
		bool found = false;               ///< True if a neighbor was found
		AlgorithmStatus status = AlgorithmStatus::Success;
		int nodes_visited = 0;            ///< Number of nodes visited (for diagnostics)
	};

	/// @brief Result of k-nearest neighbors query
	struct KNearestResult {
		std::vector<size_t> indices;      ///< Indices of k nearest points (sorted by distance)
		std::vector<Real> distances;      ///< Corresponding distances
		AlgorithmStatus status = AlgorithmStatus::Success;
		int nodes_visited = 0;            ///< Number of nodes visited
	};

	/// @brief Result of radius search query
	struct RadiusSearchResult {
		std::vector<size_t> indices;      ///< Indices of all points within radius
		std::vector<Real> distances;      ///< Corresponding distances
		AlgorithmStatus status = AlgorithmStatus::Success;
	};

	/******************************************************************************/
	/*****                    Generic KD-Tree                                 *****/
	/******************************************************************************/

	/// @brief N-dimensional KD-Tree for efficient spatial queries
	///
	/// A KD-tree is a space-partitioning data structure that organizes points
	/// in k-dimensional space. It enables efficient nearest neighbor queries,
	/// range queries, and k-nearest neighbor searches.
	///
	/// @tparam Dim Number of dimensions (typically 2 or 3)
	///
	/// @example
	/// // Build a 2D tree
	/// std::vector<VectorN<Real, 2>> points = { {0,0}, {1,1}, {2,2}, {3,3} };
	/// KDTree<2> tree;
	/// tree.build(points);
	///
	/// // Find nearest neighbor
	/// auto result = tree.findNearest({1.5, 1.5});
	/// if (result.found) {
	///     std::cout << "Nearest point index: " << result.index << "\n";
	/// }
	template<int Dim>
	class KDTree {
	public:
		using PointType = VectorN<Real, Dim>;

	private:
		/// Internal node structure
		struct Node {
			int split_axis = -1;                    ///< Axis used for splitting (0 to Dim-1)
			Real split_value = 0.0;                 ///< Split value on the axis
			std::unique_ptr<Node> left;             ///< Left child (points < split_value)
			std::unique_ptr<Node> right;            ///< Right child (points >= split_value)
			std::vector<size_t> point_indices;      ///< Point indices (only for leaf nodes)
			
			bool isLeaf() const { return left == nullptr && right == nullptr; }
		};

		std::unique_ptr<Node> _root;
		std::vector<PointType> _points;             ///< Copy of all points
		KDTreeConfig _config;
		KDTreeStats _stats;

	public:
		/// Default constructor (empty tree)
		KDTree() = default;
		
		/// Move constructor
		KDTree(KDTree&&) = default;
		
		/// Move assignment
		KDTree& operator=(KDTree&&) = default;
		
		// Disable copy (tree can be large)
		KDTree(const KDTree&) = delete;
		KDTree& operator=(const KDTree&) = delete;

		/**************************************************************************/
		/*****                    Construction                                *****/
		/**************************************************************************/

		/// @brief Build KD-tree from a set of points
		///
		/// @param points Vector of N-dimensional points
		/// @param config Tree configuration
		///
		/// Complexity: O(n log n) average case
		void build(const std::vector<PointType>& points, 
		           const KDTreeConfig& config = KDTreeConfig::Default()) {
			AlgorithmTimer timer;
			
			_config = config;
			_points = points;  // Store copy
			_stats = KDTreeStats{};
			_stats.num_points = points.size();
			
			if (points.empty()) {
				_root = nullptr;
				_stats.build_time_ms = timer.elapsed_ms();
				return;
			}
			
			// Create index array
			std::vector<size_t> indices(points.size());
			for (size_t i = 0; i < points.size(); ++i) {
				indices[i] = i;
			}
			
			// Build tree recursively
			size_t depth = 0;
			_root = buildRecursive(indices, 0, depth);
			_stats.tree_depth = depth;
			_stats.build_time_ms = timer.elapsed_ms();
		}

		/**************************************************************************/
		/*****                    Query Operations                            *****/
		/**************************************************************************/

		/// @brief Find the nearest neighbor to a query point
		///
		/// @param query The query point
		/// @return NearestNeighborResult with index and distance of nearest point
		///
		/// Complexity: O(log n) average case, O(n) worst case
		NearestNeighborResult findNearest(const PointType& query) const {
			NearestNeighborResult result;
			
			if (!_root) {
				result.status = AlgorithmStatus::InvalidInput;
				return result;
			}
			
			Real bestDistSq = std::numeric_limits<Real>::max();
			size_t bestIndex = SIZE_MAX;
			
			nearestRecursive(_root.get(), query, bestIndex, bestDistSq, result.nodes_visited);
			
			if (bestIndex != SIZE_MAX) {
				result.index = bestIndex;
				result.distance = std::sqrt(bestDistSq);
				result.found = true;
			}
			
			return result;
		}

		/// @brief Find k nearest neighbors to a query point
		///
		/// @param query The query point
		/// @param k Number of neighbors to find
		/// @return KNearestResult with indices and distances (sorted by distance)
		///
		/// Complexity: O(k log n) average case
		KNearestResult findKNearest(const PointType& query, int k) const {
			KNearestResult result;
			
			if (!_root || k <= 0) {
				if (k <= 0) result.status = AlgorithmStatus::InvalidInput;
				return result;
			}
			
			// Max-heap: (distance_sq, index)
			using HeapEntry = std::pair<Real, size_t>;
			std::priority_queue<HeapEntry> maxHeap;
			
			kNearestRecursive(_root.get(), query, k, maxHeap, result.nodes_visited);
			
			// Extract results (in reverse order since it's a max-heap)
			result.indices.resize(maxHeap.size());
			result.distances.resize(maxHeap.size());
			
			int i = static_cast<int>(maxHeap.size()) - 1;
			while (!maxHeap.empty()) {
				result.indices[i] = maxHeap.top().second;
				result.distances[i] = std::sqrt(maxHeap.top().first);
				maxHeap.pop();
				--i;
			}
			
			return result;
		}

		/// @brief Find all points within a given radius
		///
		/// @param query The query point (center of search sphere)
		/// @param radius Search radius
		/// @return RadiusSearchResult with all points within radius
		///
		/// Complexity: O(n^(1-1/d) + m) where m is the number of reported points
		RadiusSearchResult findInRadius(const PointType& query, Real radius) const {
			RadiusSearchResult result;
			
			if (!_root) {
				return result;
			}
			
			if (radius < 0) {
				result.status = AlgorithmStatus::InvalidInput;
				return result;
			}
			
			Real radiusSq = radius * radius;
			radiusRecursive(_root.get(), query, radiusSq, result.indices, result.distances);
			
			// Convert squared distances to actual distances
			for (auto& d : result.distances) {
				d = std::sqrt(d);
			}
			
			return result;
		}

		/// @brief Find all points within an axis-aligned bounding box
		///
		/// @param minCorner Minimum corner of the box
		/// @param maxCorner Maximum corner of the box
		/// @return Vector of point indices within the box
		std::vector<size_t> findInBox(const PointType& minCorner, 
		                               const PointType& maxCorner) const {
			std::vector<size_t> result;
			
			if (!_root) {
				return result;
			}
			
			boxRecursive(_root.get(), minCorner, maxCorner, result);
			return result;
		}

		/**************************************************************************/
		/*****                    Accessors                                   *****/
		/**************************************************************************/

		/// Check if tree is empty
		bool empty() const { return _root == nullptr; }
		
		/// Number of points in tree
		size_t size() const { return _stats.num_points; }
		
		/// Get build statistics
		const KDTreeStats& stats() const { return _stats; }
		
		/// Access a point by index
		const PointType& point(size_t index) const { return _points[index]; }
		
		/// Clear tree and release memory
		void clear() {
			_root.reset();
			_points.clear();
			_stats = KDTreeStats{};
		}

	private:
		/**************************************************************************/
		/*****                    Build Helpers                               *****/
		/**************************************************************************/

		std::unique_ptr<Node> buildRecursive(std::vector<size_t>& indices, 
		                                      int depth, size_t& maxDepth) {
			if (indices.empty()) {
				return nullptr;
			}
			
			maxDepth = std::max(maxDepth, static_cast<size_t>(depth));
			
			auto node = std::make_unique<Node>();
			_stats.num_nodes++;
			
			// Create leaf if few points
			if (static_cast<int>(indices.size()) <= _config.max_leaf_size) {
				node->point_indices = std::move(indices);
				_stats.num_leaf_nodes++;
				return node;
			}
			
			// Choose split axis (cycle through dimensions)
			int axis = depth % Dim;
			node->split_axis = axis;
			
			// Find median along this axis
			size_t mid = indices.size() / 2;
			std::nth_element(indices.begin(), indices.begin() + mid, indices.end(),
				[this, axis](size_t a, size_t b) {
					return _points[a][axis] < _points[b][axis];
				});
			
			node->split_value = _points[indices[mid]][axis];
			
			// Split indices
			std::vector<size_t> leftIndices(indices.begin(), indices.begin() + mid);
			std::vector<size_t> rightIndices(indices.begin() + mid, indices.end());
			
			// Recurse
			node->left = buildRecursive(leftIndices, depth + 1, maxDepth);
			node->right = buildRecursive(rightIndices, depth + 1, maxDepth);
			
			return node;
		}

		/**************************************************************************/
		/*****                    Query Helpers                               *****/
		/**************************************************************************/

		/// Squared Euclidean distance
		Real distanceSquared(const PointType& a, const PointType& b) const {
			Real sum = 0.0;
			for (int i = 0; i < Dim; ++i) {
				Real diff = a[i] - b[i];
				sum += diff * diff;
			}
			return sum;
		}

		void nearestRecursive(const Node* node, const PointType& query,
		                      size_t& bestIndex, Real& bestDistSq,
		                      int& nodesVisited) const {
			if (!node) return;
			
			nodesVisited++;
			
			if (node->isLeaf()) {
				// Check all points in leaf
				for (size_t idx : node->point_indices) {
					Real distSq = distanceSquared(query, _points[idx]);
					if (distSq < bestDistSq) {
						bestDistSq = distSq;
						bestIndex = idx;
					}
				}
				return;
			}
			
			// Determine which child to search first
			int axis = node->split_axis;
			Real diff = query[axis] - node->split_value;
			
			Node* firstChild = (diff < 0) ? node->left.get() : node->right.get();
			Node* secondChild = (diff < 0) ? node->right.get() : node->left.get();
			
			// Search the closer child first
			nearestRecursive(firstChild, query, bestIndex, bestDistSq, nodesVisited);
			
			// Check if we need to search the other child
			// (only if the splitting plane is closer than current best)
			if (diff * diff < bestDistSq) {
				nearestRecursive(secondChild, query, bestIndex, bestDistSq, nodesVisited);
			}
		}

		void kNearestRecursive(const Node* node, const PointType& query, int k,
		                       std::priority_queue<std::pair<Real, size_t>>& maxHeap,
		                       int& nodesVisited) const {
			if (!node) return;
			
			nodesVisited++;
			
			if (node->isLeaf()) {
				for (size_t idx : node->point_indices) {
					Real distSq = distanceSquared(query, _points[idx]);
					
					if (static_cast<int>(maxHeap.size()) < k) {
						maxHeap.push({distSq, idx});
					} else if (distSq < maxHeap.top().first) {
						maxHeap.pop();
						maxHeap.push({distSq, idx});
					}
				}
				return;
			}
			
			int axis = node->split_axis;
			Real diff = query[axis] - node->split_value;
			
			Node* firstChild = (diff < 0) ? node->left.get() : node->right.get();
			Node* secondChild = (diff < 0) ? node->right.get() : node->left.get();
			
			kNearestRecursive(firstChild, query, k, maxHeap, nodesVisited);
			
			// Check other side if heap not full or plane is closer than k-th best
			Real worstDistSq = maxHeap.empty() ? std::numeric_limits<Real>::max() 
			                                    : maxHeap.top().first;
			if (static_cast<int>(maxHeap.size()) < k || diff * diff < worstDistSq) {
				kNearestRecursive(secondChild, query, k, maxHeap, nodesVisited);
			}
		}

		void radiusRecursive(const Node* node, const PointType& query, Real radiusSq,
		                     std::vector<size_t>& indices, 
		                     std::vector<Real>& distancesSq) const {
			if (!node) return;
			
			if (node->isLeaf()) {
				for (size_t idx : node->point_indices) {
					Real distSq = distanceSquared(query, _points[idx]);
					if (distSq <= radiusSq) {
						indices.push_back(idx);
						distancesSq.push_back(distSq);
					}
				}
				return;
			}
			
			int axis = node->split_axis;
			Real diff = query[axis] - node->split_value;
			
			// Search child containing query point
			if (diff < 0) {
				radiusRecursive(node->left.get(), query, radiusSq, indices, distancesSq);
				// Search other child if sphere crosses splitting plane
				if (diff * diff <= radiusSq) {
					radiusRecursive(node->right.get(), query, radiusSq, indices, distancesSq);
				}
			} else {
				radiusRecursive(node->right.get(), query, radiusSq, indices, distancesSq);
				if (diff * diff <= radiusSq) {
					radiusRecursive(node->left.get(), query, radiusSq, indices, distancesSq);
				}
			}
		}

		void boxRecursive(const Node* node, const PointType& minCorner,
		                  const PointType& maxCorner, std::vector<size_t>& result) const {
			if (!node) return;
			
			if (node->isLeaf()) {
				for (size_t idx : node->point_indices) {
					const PointType& p = _points[idx];
					bool inside = true;
					for (int i = 0; i < Dim && inside; ++i) {
						if (p[i] < minCorner[i] || p[i] > maxCorner[i]) {
							inside = false;
						}
					}
					if (inside) {
						result.push_back(idx);
					}
				}
				return;
			}
			
			int axis = node->split_axis;
			
			// Check if we need to search left child
			if (minCorner[axis] <= node->split_value) {
				boxRecursive(node->left.get(), minCorner, maxCorner, result);
			}
			
			// Check if we need to search right child
			if (maxCorner[axis] >= node->split_value) {
				boxRecursive(node->right.get(), minCorner, maxCorner, result);
			}
		}
	};

	/******************************************************************************/
	/*****                    2D Specialization                               *****/
	/******************************************************************************/

	/// @brief 2D KD-Tree with Point2Cartesian support
	///
	/// Convenient wrapper for 2D point queries using MML's Point2Cartesian type.
	class KDTree2D {
	public:
		using PointType = Point2Cartesian;
		
	private:
		KDTree<2> _tree;
		std::vector<Point2Cartesian> _originalPoints;
		
		static VectorN<Real, 2> toVectorN(const Point2Cartesian& p) {
			return VectorN<Real, 2>({p.X(), p.Y()});
		}
		
	public:
		/// Build tree from Point2Cartesian points
		void build(const std::vector<Point2Cartesian>& points,
		           const KDTreeConfig& config = KDTreeConfig::Default()) {
			_originalPoints = points;
			
			std::vector<VectorN<Real, 2>> vecPoints;
			vecPoints.reserve(points.size());
			for (const auto& p : points) {
				vecPoints.push_back(toVectorN(p));
			}
			
			_tree.build(vecPoints, config);
		}
		
		/// Find nearest neighbor
		NearestNeighborResult findNearest(const Point2Cartesian& query) const {
			return _tree.findNearest(toVectorN(query));
		}
		
		/// Find k nearest neighbors
		KNearestResult findKNearest(const Point2Cartesian& query, int k) const {
			return _tree.findKNearest(toVectorN(query), k);
		}
		
		/// Find all points within radius
		RadiusSearchResult findInRadius(const Point2Cartesian& query, Real radius) const {
			return _tree.findInRadius(toVectorN(query), radius);
		}
		
		/// Find all points within bounding box
		std::vector<size_t> findInBox(const Point2Cartesian& minCorner,
		                               const Point2Cartesian& maxCorner) const {
			return _tree.findInBox(toVectorN(minCorner), toVectorN(maxCorner));
		}
		
		bool empty() const { return _tree.empty(); }
		size_t size() const { return _tree.size(); }
		const KDTreeStats& stats() const { return _tree.stats(); }
		const Point2Cartesian& point(size_t index) const { return _originalPoints[index]; }
		void clear() { _tree.clear(); _originalPoints.clear(); }
	};

	/******************************************************************************/
	/*****                    3D Specialization                               *****/
	/******************************************************************************/

	/// @brief 3D KD-Tree with Point3Cartesian support
	///
	/// Convenient wrapper for 3D point queries using MML's Point3Cartesian type.
	class KDTree3D {
	public:
		using PointType = Point3Cartesian;
		
	private:
		KDTree<3> _tree;
		std::vector<Point3Cartesian> _originalPoints;
		
		static VectorN<Real, 3> toVectorN(const Point3Cartesian& p) {
			return VectorN<Real, 3>({p.X(), p.Y(), p.Z()});
		}
		
	public:
		/// Build tree from Point3Cartesian points
		void build(const std::vector<Point3Cartesian>& points,
		           const KDTreeConfig& config = KDTreeConfig::Default()) {
			_originalPoints = points;
			
			std::vector<VectorN<Real, 3>> vecPoints;
			vecPoints.reserve(points.size());
			for (const auto& p : points) {
				vecPoints.push_back(toVectorN(p));
			}
			
			_tree.build(vecPoints, config);
		}
		
		/// Find nearest neighbor
		NearestNeighborResult findNearest(const Point3Cartesian& query) const {
			return _tree.findNearest(toVectorN(query));
		}
		
		/// Find k nearest neighbors
		KNearestResult findKNearest(const Point3Cartesian& query, int k) const {
			return _tree.findKNearest(toVectorN(query), k);
		}
		
		/// Find all points within radius
		RadiusSearchResult findInRadius(const Point3Cartesian& query, Real radius) const {
			return _tree.findInRadius(toVectorN(query), radius);
		}
		
		/// Find all points within bounding box
		std::vector<size_t> findInBox(const Point3Cartesian& minCorner,
		                               const Point3Cartesian& maxCorner) const {
			return _tree.findInBox(toVectorN(minCorner), toVectorN(maxCorner));
		}
		
		bool empty() const { return _tree.empty(); }
		size_t size() const { return _tree.size(); }
		const KDTreeStats& stats() const { return _tree.stats(); }
		const Point3Cartesian& point(size_t index) const { return _originalPoints[index]; }
		void clear() { _tree.clear(); _originalPoints.clear(); }
	};

	/******************************************************************************/
	/*****                    Utility Functions                               *****/
	/******************************************************************************/

	/// @brief Brute-force nearest neighbor search (for testing/comparison)
	inline NearestNeighborResult FindNearestBruteForce2D(
		const std::vector<Point2Cartesian>& points,
		const Point2Cartesian& query) {
		
		NearestNeighborResult result;
		
		if (points.empty()) {
			return result;
		}
		
		Real bestDistSq = std::numeric_limits<Real>::max();
		
		for (size_t i = 0; i < points.size(); ++i) {
			Real dx = points[i].X() - query.X();
			Real dy = points[i].Y() - query.Y();
			Real distSq = dx * dx + dy * dy;
			
			if (distSq < bestDistSq) {
				bestDistSq = distSq;
				result.index = i;
			}
		}
		
		result.distance = std::sqrt(bestDistSq);
		result.found = true;
		result.nodes_visited = static_cast<int>(points.size());  // All points checked
		
		return result;
	}

	/// @brief Brute-force nearest neighbor search for 3D (for testing/comparison)
	inline NearestNeighborResult FindNearestBruteForce3D(
		const std::vector<Point3Cartesian>& points,
		const Point3Cartesian& query) {
		
		NearestNeighborResult result;
		
		if (points.empty()) {
			return result;
		}
		
		Real bestDistSq = std::numeric_limits<Real>::max();
		
		for (size_t i = 0; i < points.size(); ++i) {
			Real dx = points[i].X() - query.X();
			Real dy = points[i].Y() - query.Y();
			Real dz = points[i].Z() - query.Z();
			Real distSq = dx * dx + dy * dy + dz * dz;
			
			if (distSq < bestDistSq) {
				bestDistSq = distSq;
				result.index = i;
			}
		}
		
		result.distance = std::sqrt(bestDistSq);
		result.found = true;
		result.nodes_visited = static_cast<int>(points.size());
		
		return result;
	}

} // namespace MML

#endif // MML_KDTREE_H
