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

#include <mml/MMLBase.h>
#include <mml/MMLNumericValidation.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Geometry/Geometry2D.h>
#include <mml/base/Geometry/Geometry3D.h>

#include <vector>
#include <memory>
#include <algorithm>
#include <queue>
#include <cmath>
#include <limits>
#include <numeric>

namespace MML {

	namespace KDTreeDetail {
		template<int Dim>
		bool IsFinitePoint(const VectorN<Real, Dim>& point) {
			for (int d = 0; d < Dim; ++d)
				if (!std::isfinite(point[d])) return false;
			return true;
		}

		template<int Dim>
		Real Distance(const VectorN<Real, Dim>& a, const VectorN<Real, Dim>& b) {
			Real distance = 0;
			for (int d = 0; d < Dim; ++d)
				distance = std::hypot(distance, a[d] - b[d]);
			// An overflowing coordinate difference already implies an unrepresentable distance.
			// Keep it infinite internally; only a selected output distance makes the query fail.
			return distance;
		}

		inline VectorN<Real, 2> Coordinates(const Point2Cartesian& point) { return {point.X(), point.Y()}; }
		inline VectorN<Real, 3> Coordinates(const Point3Cartesian& point) { return {point.X(), point.Y(), point.Z()}; }
	}

	/******************************************************************************/
	/*****                    Configuration                                   *****/
	/******************************************************************************/

	/// @brief Configuration for KD-tree construction
	struct KDTreeConfig {
		int max_leaf_size = 10;           ///< Positive maximum points per leaf before splitting
		bool balanced = true;             ///< Median vs midpoint split, with median fallback for empty partitions
		
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
		static_assert(Dim > 0, "KDTree requires a positive dimension");
	public:
		using PointType = VectorN<Real, Dim>;

	private:
		/// Internal node structure
		struct Node {
			int split_axis = -1;                    ///< Axis used for splitting (0 to Dim-1)
			Real split_value = 0.0;                 ///< Split value on the axis
			std::unique_ptr<Node> left;             ///< Left child (points <= split_value; ties may straddle a median)
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
		/// Rejects nonpositive leaf limits and nonfinite coordinates before replacing the tree.
		/// Complexity: O(n log n) average case for balanced construction
		void build(const std::vector<PointType>& points, 
		           const KDTreeConfig& config = KDTreeConfig::Default()) {
			AlgorithmTimer timer;
			if (config.max_leaf_size <= 0)
				throw ArgumentError("KDTree: max_leaf_size must be positive");
			for (const auto& point : points)
				for (int d = 0; d < Dim; ++d) ValidateFinite(point[d], "KDTree point coordinate");

			KDTree next;
			next._config = config;
			next._points = points;
			next._stats.num_points = points.size();
			
			// Create index array
			std::vector<size_t> indices(points.size());
			for (size_t i = 0; i < points.size(); ++i) {
				indices[i] = i;
			}
			
			// Build tree recursively
			size_t depth = 0;
			next._root = next.buildRecursive(indices, 0, depth);
			next._stats.tree_depth = depth;
			next._stats.build_time_ms = timer.elapsed_ms();
			*this = std::move(next);
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
			
			if (!_root || !KDTreeDetail::IsFinitePoint(query)) {
				result.status = AlgorithmStatus::InvalidInput;
				return result;
			}
			
			Real bestDistance = std::numeric_limits<Real>::infinity();
			size_t bestIndex = SIZE_MAX;
			
			nearestRecursive(_root.get(), query, bestIndex, bestDistance, result.nodes_visited);
			if (!std::isfinite(bestDistance)) {
				result.status = AlgorithmStatus::NumericalInstability;
				return result;
			}
			
			if (bestIndex != SIZE_MAX) {
				result.index = bestIndex;
				result.distance = bestDistance;
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
			
			if (k <= 0 || !KDTreeDetail::IsFinitePoint(query)) {
				result.status = AlgorithmStatus::InvalidInput;
				return result;
			}
			if (!_root) return result;
			
			// Max-heap: (distance, index)
			using HeapEntry = std::pair<Real, size_t>;
			std::priority_queue<HeapEntry> maxHeap;
			
			kNearestRecursive(_root.get(), query, k, maxHeap, result.nodes_visited);
			if (!maxHeap.empty() && !std::isfinite(maxHeap.top().first)) {
				result.status = AlgorithmStatus::NumericalInstability;
				return result;
			}
			
			// Extract results (in reverse order since it's a max-heap)
			result.indices.resize(maxHeap.size());
			result.distances.resize(maxHeap.size());
			
			int i = static_cast<int>(maxHeap.size()) - 1;
			while (!maxHeap.empty()) {
				result.indices[i] = maxHeap.top().second;
				result.distances[i] = maxHeap.top().first;
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
			
			if (!KDTreeDetail::IsFinitePoint(query) || !std::isfinite(radius) || radius < 0) {
				result.status = AlgorithmStatus::InvalidInput;
				return result;
			}
			
			radiusRecursive(_root.get(), query, radius, result.indices, result.distances);
			
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
			for (int d = 0; d < Dim; ++d) {
				ValidateFinite(minCorner[d], "KDTree box lower coordinate");
				ValidateFinite(maxCorner[d], "KDTree box upper coordinate");
				if (minCorner[d] > maxCorner[d])
					throw ArgumentError("KDTree: box lower corner must not exceed upper corner");
			}
			
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
		                                      size_t depth, size_t& maxDepth) {
			if (indices.empty()) {
				return nullptr;
			}
			
			maxDepth = std::max(maxDepth, depth);
			
			auto node = std::make_unique<Node>();
			_stats.num_nodes++;
			
			// Create leaf if few points
			if (indices.size() <= static_cast<size_t>(_config.max_leaf_size)) {
				node->point_indices = std::move(indices);
				_stats.num_leaf_nodes++;
				return node;
			}
			
			// Choose split axis (cycle through dimensions)
			const int axis = static_cast<int>(depth % static_cast<size_t>(Dim));
			node->split_axis = axis;
			
			size_t mid = indices.size() / 2;
			const auto compare = [this, axis](size_t a, size_t b) {
				return _points[a][axis] < _points[b][axis];
			};
			bool useMedian = true;
			if (!_config.balanced) {
				const auto bounds = std::minmax_element(indices.begin(), indices.end(), compare);
				const Real split = std::midpoint(_points[*bounds.first][axis], _points[*bounds.second][axis]);
				const auto boundary = std::partition(indices.begin(), indices.end(), [this, axis, split](size_t index) {
					return _points[index][axis] < split;
				});
				const size_t count = static_cast<size_t>(boundary - indices.begin());
				if (count != 0 && count != indices.size()) {
					mid = count;
					node->split_value = split;
					useMedian = false;
				}
			}
			if (useMedian) {
				std::nth_element(indices.begin(), indices.begin() + mid, indices.end(), compare);
				node->split_value = _points[indices[mid]][axis];
			}
			
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

		void nearestRecursive(const Node* node, const PointType& query,
		                      size_t& bestIndex, Real& bestDistance,
		                      int& nodesVisited) const {
			if (!node) return;
			
			nodesVisited++;
			
			if (node->isLeaf()) {
				// Check all points in leaf
				for (size_t idx : node->point_indices) {
					const Real distance = KDTreeDetail::Distance(query, _points[idx]);
					if (bestIndex == SIZE_MAX || distance < bestDistance) {
						bestDistance = distance;
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
			nearestRecursive(firstChild, query, bestIndex, bestDistance, nodesVisited);
			
			// Check if we need to search the other child
			// (only if the splitting plane is closer than current best)
			if (std::abs(diff) < bestDistance) {
				nearestRecursive(secondChild, query, bestIndex, bestDistance, nodesVisited);
			}
		}

		void kNearestRecursive(const Node* node, const PointType& query, int k,
		                       std::priority_queue<std::pair<Real, size_t>>& maxHeap,
		                       int& nodesVisited) const {
			if (!node) return;
			
			nodesVisited++;
			
			if (node->isLeaf()) {
				for (size_t idx : node->point_indices) {
					const Real distance = KDTreeDetail::Distance(query, _points[idx]);
					
					if (static_cast<int>(maxHeap.size()) < k) {
						maxHeap.push({distance, idx});
					} else if (distance < maxHeap.top().first) {
						maxHeap.pop();
						maxHeap.push({distance, idx});
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
			Real worstDistance = maxHeap.empty() ? std::numeric_limits<Real>::infinity()
			                                    : maxHeap.top().first;
			if (static_cast<int>(maxHeap.size()) < k || std::abs(diff) < worstDistance) {
				kNearestRecursive(secondChild, query, k, maxHeap, nodesVisited);
			}
		}

		void radiusRecursive(const Node* node, const PointType& query, Real radius,
		                     std::vector<size_t>& indices, 
		                     std::vector<Real>& distances) const {
			if (!node) return;
			
			if (node->isLeaf()) {
				for (size_t idx : node->point_indices) {
					const Real distance = KDTreeDetail::Distance(query, _points[idx]);
					if (distance <= radius) {
						indices.push_back(idx);
						distances.push_back(distance);
					}
				}
				return;
			}
			
			int axis = node->split_axis;
			Real diff = query[axis] - node->split_value;
			
			// Search child containing query point
			if (diff < 0) {
				radiusRecursive(node->left.get(), query, radius, indices, distances);
				// Search other child if sphere crosses splitting plane
				if (std::abs(diff) <= radius) {
					radiusRecursive(node->right.get(), query, radius, indices, distances);
				}
			} else {
				radiusRecursive(node->right.get(), query, radius, indices, distances);
				if (std::abs(diff) <= radius) {
					radiusRecursive(node->left.get(), query, radius, indices, distances);
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
			auto originalPoints = points;
			
			std::vector<VectorN<Real, 2>> vecPoints;
			vecPoints.reserve(points.size());
			for (const auto& p : points) {
				vecPoints.push_back(toVectorN(p));
			}
			
			_tree.build(vecPoints, config);
			_originalPoints = std::move(originalPoints);
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
			auto originalPoints = points;
			
			std::vector<VectorN<Real, 3>> vecPoints;
			vecPoints.reserve(points.size());
			for (const auto& p : points) {
				vecPoints.push_back(toVectorN(p));
			}
			
			_tree.build(vecPoints, config);
			_originalPoints = std::move(originalPoints);
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

	namespace KDTreeDetail {
		template<class CartesianPoint>
		NearestNeighborResult FindNearestBruteForce(const std::vector<CartesianPoint>& points, const CartesianPoint& query) {
			NearestNeighborResult result;
			const auto queryCoordinates = Coordinates(query);
			if (!IsFinitePoint(queryCoordinates)) {
				result.status = AlgorithmStatus::InvalidInput;
				return result;
			}
			Real bestDistance = std::numeric_limits<Real>::infinity();
			size_t bestIndex = SIZE_MAX;
			for (size_t i = 0; i < points.size(); ++i) {
				const auto coordinates = Coordinates(points[i]);
				if (!IsFinitePoint(coordinates)) {
					result.status = AlgorithmStatus::InvalidInput;
					return result;
				}
				const Real distance = Distance(coordinates, queryCoordinates);
				if (bestIndex == SIZE_MAX || distance < bestDistance) {
					bestDistance = distance;
					bestIndex = i;
				}
			}
			result.nodes_visited = static_cast<int>(points.size());
			if (bestIndex == SIZE_MAX) return result;
			if (!std::isfinite(bestDistance)) {
				result.status = AlgorithmStatus::NumericalInstability;
				return result;
			}
			result.index = bestIndex;
			result.distance = bestDistance;
			result.found = true;
			return result;
		}
	}

	/// @brief Brute-force nearest neighbor search (for testing/comparison)
	inline NearestNeighborResult FindNearestBruteForce2D(
		const std::vector<Point2Cartesian>& points,
		const Point2Cartesian& query) {
		return KDTreeDetail::FindNearestBruteForce(points, query);
	}

	/// @brief Brute-force nearest neighbor search for 3D (for testing/comparison)
	inline NearestNeighborResult FindNearestBruteForce3D(
		const std::vector<Point3Cartesian>& points,
		const Point3Cartesian& query) {
		return KDTreeDetail::FindNearestBruteForce(points, query);
	}

} // namespace MML

#endif // MML_KDTREE_H
