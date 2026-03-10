///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        comp_geometry/Circles.h                                             ///
///  Description: Circle-related algorithms                                           ///
///               Minimum Enclosing Circle (Welzl's Algorithm)                        ///
///               Closest Pair of Points (Divide and Conquer)                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMP_GEOMETRY_CIRCLES_H
#define MML_COMP_GEOMETRY_CIRCLES_H

#include "mml/algorithms/CompGeometry/CompGeometryBase.h"
#include "mml/MMLBase.h"
#include "mml/base/Geometry/Geometry2D.h"

#include <algorithm>
#include <random>
#include <vector>

namespace MML {
namespace CompGeometry {

	/// @brief Circle-related computational geometry algorithms
	class Circles {
	private:
		// Use centralized geometry epsilon from Constants
		static constexpr Real EPSILON = Constants::GEOMETRY_EPSILON;

		// ========================================================================
		// Helper functions for Minimum Enclosing Circle
		// ========================================================================

		static void CircleFrom2Points(const Point2Cartesian& p1, const Point2Cartesian& p2, 
		                              Point2Cartesian& center, Real& radius) {
			center = Point2Cartesian((p1.X() + p2.X()) / 2.0, (p1.Y() + p2.Y()) / 2.0);
			radius = p1.Dist(p2) / 2.0;
		}

		static bool CircleFrom3Points(const Point2Cartesian& p1, const Point2Cartesian& p2, 
		                              const Point2Cartesian& p3, Point2Cartesian& center, Real& radius) {
			Real ax = p2.X() - p1.X();
			Real ay = p2.Y() - p1.Y();
			Real bx = p3.X() - p1.X();
			Real by = p3.Y() - p1.Y();

			Real d = 2.0 * (ax * by - ay * bx);
			if (std::abs(d) < EPSILON)
				return false; // Collinear points

			Real aSq = ax * ax + ay * ay;
			Real bSq = bx * bx + by * by;

			Real cx = (by * aSq - ay * bSq) / d;
			Real cy = (ax * bSq - bx * aSq) / d;

			center = Point2Cartesian(p1.X() + cx, p1.Y() + cy);
			radius = std::sqrt(cx * cx + cy * cy);
			return true;
		}

		static bool IsInsideCircle(const Point2Cartesian& p, const Point2Cartesian& center, Real radius) {
			return p.Dist(center) <= radius + EPSILON;
		}

		static void WelzlHelper(std::vector<Point2Cartesian>& points, std::vector<Point2Cartesian>& boundary, 
		                        int n, Point2Cartesian& center, Real& radius, std::mt19937& rng) {
			if (n == 0 || boundary.size() == 3) {
				if (boundary.empty()) {
					center = Point2Cartesian(0, 0);
					radius = 0;
				} else if (boundary.size() == 1) {
					center = boundary[0];
					radius = 0;
				} else if (boundary.size() == 2) {
					CircleFrom2Points(boundary[0], boundary[1], center, radius);
				} else {
					if (!CircleFrom3Points(boundary[0], boundary[1], boundary[2], center, radius)) {
						Real d01 = boundary[0].Dist(boundary[1]);
						Real d02 = boundary[0].Dist(boundary[2]);
						Real d12 = boundary[1].Dist(boundary[2]);

						if (d01 >= d02 && d01 >= d12)
							CircleFrom2Points(boundary[0], boundary[1], center, radius);
						else if (d02 >= d12)
							CircleFrom2Points(boundary[0], boundary[2], center, radius);
						else
							CircleFrom2Points(boundary[1], boundary[2], center, radius);
					}
				}
				return;
			}

			std::uniform_int_distribution<int> dist(0, n - 1);
			int idx = dist(rng);
			std::swap(points[idx], points[n - 1]);
			Point2Cartesian p = points[n - 1];

			WelzlHelper(points, boundary, n - 1, center, radius, rng);

			if (IsInsideCircle(p, center, radius))
				return;

			boundary.push_back(p);
			WelzlHelper(points, boundary, n - 1, center, radius, rng);
			boundary.pop_back();
		}

		// ========================================================================
		// Helper functions for Closest Pair
		// ========================================================================

		static ClosestPairResult ClosestPairBruteForce(const std::vector<Point2Cartesian>& points, 
		                                               int left, int right) {
			ClosestPairResult best;
			for (int i = left; i < right; i++) {
				for (int j = i + 1; j <= right; j++) {
					Real d = points[i].Dist(points[j]);
					if (d < best.distance)
						best = ClosestPairResult(points[i], points[j], d);
				}
			}
			return best;
		}

		static ClosestPairResult ClosestPairStrip(std::vector<Point2Cartesian>& strip, Real delta) {
			ClosestPairResult best;
			best.distance = delta;

			std::sort(strip.begin(), strip.end(), 
				[](const Point2Cartesian& a, const Point2Cartesian& b) { return a.Y() < b.Y(); });

			int n = static_cast<int>(strip.size());
			for (int i = 0; i < n; i++) {
				for (int j = i + 1; j < n && (strip[j].Y() - strip[i].Y()) < best.distance; j++) {
					Real d = strip[i].Dist(strip[j]);
					if (d < best.distance)
						best = ClosestPairResult(strip[i], strip[j], d);
				}
			}
			return best;
		}

		static ClosestPairResult ClosestPairHelper(std::vector<Point2Cartesian>& pointsX, int left, int right) {
			if (right - left <= 3)
				return ClosestPairBruteForce(pointsX, left, right);

			int mid = (left + right) / 2;
			Point2Cartesian midPoint = pointsX[mid];

			ClosestPairResult leftResult = ClosestPairHelper(pointsX, left, mid);
			ClosestPairResult rightResult = ClosestPairHelper(pointsX, mid + 1, right);

			ClosestPairResult best = (leftResult.distance < rightResult.distance) ? leftResult : rightResult;

			std::vector<Point2Cartesian> strip;
			for (int i = left; i <= right; i++) {
				if (std::abs(pointsX[i].X() - midPoint.X()) < best.distance)
					strip.push_back(pointsX[i]);
			}

			ClosestPairResult stripResult = ClosestPairStrip(strip, best.distance);
			if (stripResult.distance < best.distance)
				best = stripResult;

			return best;
		}

	public:
		// ========================================================================
		// Minimum Enclosing Circle (Welzl's Algorithm)
		// ========================================================================
		// Computes the smallest circle containing all given points.
		// Complexity: Expected O(n) time
		// ========================================================================

		static Circle2D MinimumEnclosingCircle(std::vector<Point2Cartesian> points) {
			int n = static_cast<int>(points.size());

			if (n == 0)
				return Circle2D(Point2Cartesian(0, 0), 0);
			if (n == 1)
				return Circle2D(points[0], 0);

			std::random_device rd;
			std::mt19937 rng(rd());
			std::shuffle(points.begin(), points.end(), rng);

			Point2Cartesian center;
			Real radius;
			std::vector<Point2Cartesian> boundary;

			WelzlHelper(points, boundary, n, center, radius, rng);

			return Circle2D(center, radius);
		}

		static Circle2D MinimumEnclosingCircle(const Polygon2D& polygon) { 
			return MinimumEnclosingCircle(polygon.Vertices()); 
		}

		// ========================================================================
		// Closest Pair of Points - Divide and Conquer
		// ========================================================================
		// Finds the two closest points in a set using divide-and-conquer.
		// Complexity: O(n log n)
		// ========================================================================

		static ClosestPairResult ClosestPair(std::vector<Point2Cartesian> points) {
			int n = static_cast<int>(points.size());
			if (n < 2)
				return ClosestPairResult();

			if (n == 2)
				return ClosestPairResult(points[0], points[1], points[0].Dist(points[1]));

			std::sort(points.begin(), points.end(), 
				[](const Point2Cartesian& a, const Point2Cartesian& b) { return a.X() < b.X(); });

			return ClosestPairHelper(points, 0, n - 1);
		}
	};

} // namespace CompGeometry
} // namespace MML

#endif // MML_COMP_GEOMETRY_CIRCLES_H
