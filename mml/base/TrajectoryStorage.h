#ifndef MML_TRAJECTORY_STORAGE_H
#define MML_TRAJECTORY_STORAGE_H

#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>
#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MML::detail {
	inline int TrajectoryMaxPoints(std::initializer_list<int> dimensions)
	{
		int limit = static_cast<int>(std::min(static_cast<size_t>(Matrix<Real>::MaxDimension()),
			Matrix<Real>::MaxElementCount()));
		for (int rows : dimensions) {
			if (rows < 0)
				throw MatrixDimensionError("Trajectory - invalid state dimension", rows, 0,
					Matrix<Real>::MaxDimension(), Matrix<Real>::MaxDimension());
			if (rows > 0)
				limit = static_cast<int>(std::min(static_cast<size_t>(limit),
					Matrix<Real>::MaxElementCount() / static_cast<size_t>(rows)));
		}
		return limit;
	}

	inline int TrajectoryCapacity(int current, std::int64_t required, int limit, bool geometric)
	{
		if (required < 0)
			throw ArgumentError("Trajectory capacity must be non-negative");
		if (required > limit)
			throw ResourceLimitError("Trajectory storage exceeds its element or index limit");
		int capacity = static_cast<int>(required);
		if (geometric) {
			capacity = current;
			while (capacity < required)
				capacity += std::min(limit - capacity, capacity / 2 + 1);
		}
		return capacity;
	}

	inline void ResizeTrajectoryStorage(int capacity, Vector<Real>& times,
		std::initializer_list<Matrix<Real>*> matrices)
	{
		if (capacity == times.size()) return;
		try {
			Vector<Real> newTimes(capacity);
			std::copy_n(times.begin(), std::min(times.size(), capacity), newTimes.begin());
			std::vector<Matrix<Real>> replacements;
			replacements.reserve(matrices.size());
			for (const auto* matrix : matrices) {
				auto& replacement = replacements.emplace_back(matrix->rows(), capacity);
				for (int row = 0; row < matrix->rows(); ++row)
					std::copy_n((*matrix)[row], std::min(matrix->cols(), capacity), replacement[row]);
			}
			// Every allocation and copy succeeds before the noexcept commit.
			times = std::move(newTimes);
			size_t index = 0;
			for (auto* matrix : matrices)
				*matrix = std::move(replacements[index++]);
		}
		catch (const std::bad_alloc&) {
			throw ResourceLimitError("Unable to allocate trajectory storage");
		}
		catch (const std::length_error&) {
			throw ResourceLimitError("Trajectory allocation exceeds the allocator's size limit");
		}
	}

	template<class Solution>
	bool PrepareTrajectoryPoint(Solution& solution, int index, AlgorithmStatus& status,
		std::string& message)
	{
		try {
			solution.reserveForPoint(index);
			return true;
		}
		catch (const ResourceLimitError& error) {
			status = AlgorithmStatus::ResourceLimitExceeded;
			message = error.what();
			return false;
		}
	}
}
#endif
