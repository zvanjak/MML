#ifndef MML_PARALLEL_EXECUTION_H
#define MML_PARALLEL_EXECUTION_H

#include <mml/MMLBase.h>
#include <mml/base/ThreadPool.h>
#include <algorithm>
#include <exception>
#include <vector>

namespace MML {
	struct ParallelExecutionConfig {
		size_t thread_count = 0; ///< Zero selects hardware concurrency, or one when unavailable.
		size_t min_parallel_tasks = ParallelDefaults::MinTasks; ///< Smaller batches execute inline.
	};

	namespace ParallelDetail {
		inline size_t WorkerCount(size_t requested, size_t tasks) {
			const size_t available = requested ? requested : std::max(1u, std::thread::hardware_concurrency());
			return std::min(available, tasks);
		}

		class ContextGuard {
			AlgorithmContext previous = AlgorithmContext::Get();
		public:
			explicit ContextGuard(const AlgorithmContext& context) { AlgorithmContext::Get() = context; }
			~ContextGuard() { AlgorithmContext::Get() = previous; }
		};

		/// A private, fully drained indexed batch. Failures are rethrown by index, not completion order.
		template<class Function>
		size_t Execute(size_t count, const ParallelExecutionConfig& config, Function&& function) {
			if (count == 0) return 0;
			const size_t workers = WorkerCount(config.thread_count, count);
			const auto context = AlgorithmContext::Get();
			std::vector<std::exception_ptr> failures(count);
			auto run = [&](size_t index) {
				ContextGuard guard(context);
				try { function(index); }
				catch (...) { failures[index] = std::current_exception(); }
			};
			const bool parallel = workers > 1 && count >= config.min_parallel_tasks;
			if (parallel) {
				// Pool lifetime ends before the referenced task state, including on enqueue failure.
				ThreadPool pool(workers);
				for (size_t i = 0; i < count; ++i)
					pool.enqueue([&, i] { run(i); });
				pool.wait_for_tasks();
			} else {
				for (size_t i = 0; i < count; ++i) run(i);
			}
			for (const auto& failure : failures)
				if (failure) std::rethrow_exception(failure);
			return parallel ? workers : 1;
		}
	}
}
#endif
