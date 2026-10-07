///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        detail/GraphAlgorithmUtils.h                                        ///
///  Description: Shared graph algorithm utility helpers                             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_ALGORITHMS_DETAIL_GRAPHALGORITHMUTILS_H
#define MML_GRAPH_ALGORITHMS_DETAIL_GRAPHALGORITHMUTILS_H

#include <mml/base/Graph.h>
#include <chrono>
#include <string>

namespace MML
{
	namespace detail
	{
		inline double ElapsedMilliseconds(const std::chrono::steady_clock::time_point& start)
		{
			return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
		}

		template<typename V, typename E>
		bool HasNegativeWeight(const Graph<V, E>& graph)
		{
			for (const auto& [from, to, weight] : graph.directedEdges())
			{
				(void)from;
				(void)to;
				if (static_cast<Real>(weight) < Real{0})
					return true;
			}
			return false;
		}

		inline void SetGraphMessage(std::string& message, std::string& diagnostics, const std::string& text)
		{
			message = text;
			diagnostics = text;
		}

		inline ShortestPathTreeResult MakeShortestPathTree(const TraversalResult& traversal, size_t source, const std::string& algorithmName)
		{
			ShortestPathTreeResult result;
			result.status = traversal.status;
			result.graphStatus = traversal.graphStatus;
			result.source = source;
			result.parent = traversal.parent;
			result.distance = traversal.distance;
			result.nodesVisited = traversal.nodesVisited;
			result.algorithm_name = algorithmName;
			result.message = traversal.message;
			result.diagnostics = traversal.diagnostics;
			result.elapsed_time_ms = traversal.elapsed_time_ms;
			return result;
		}
	}  // namespace detail
}

#endif
