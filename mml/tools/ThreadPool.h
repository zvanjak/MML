///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ThreadPool.h                                                        ///
///  Description: Thread pool for parallel computation                                ///
///               Task queue with worker threads for concurrent execution             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_THREAD_POOL_H
#define MML_THREAD_POOL_H

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>

namespace MML
{
	class ThreadPool {
	private:
		std::vector<std::thread> workers;
		std::queue<std::function<void()>> tasks;
		std::mutex queue_mutex;
		std::condition_variable condition;
		std::condition_variable tasks_done_condition;  // P0-2: For wait_for_tasks()
		std::atomic<bool> stop{false};                 // P0-3: Atomic to prevent data race
		std::atomic<size_t> in_flight_tasks{0};        // P0-2: Track executing tasks

	public:
		ThreadPool(size_t numThreads)
		{
			for (size_t i = 0; i < numThreads; ++i)
				workers.emplace_back([this]
					{
						while (true) {
							std::function<void()> task;
							{
								std::unique_lock<std::mutex> lock(this->queue_mutex);

								this->condition.wait(lock, [this] { return this->stop.load() || !this->tasks.empty(); });

								if (this->stop.load() && this->tasks.empty())
									return;

								task = std::move(this->tasks.front());
								this->tasks.pop();
								++this->in_flight_tasks;  // P0-2: Track task start
							}
							task();
							
							// P0-2: Track task completion and notify waiters
							if (--this->in_flight_tasks == 0) {
								std::lock_guard<std::mutex> lock(this->queue_mutex);
								if (this->tasks.empty()) {
									this->tasks_done_condition.notify_all();
								}
							}
						}
					});
		}
		~ThreadPool()
		{
			// P0-1: Set stop flag and notify under lock, then release BEFORE joining
			{
				std::lock_guard<std::mutex> lock(queue_mutex);
				stop.store(true);
			}
			condition.notify_all();
			
			// P0-1: Now safe to join - workers can acquire mutex to wake up
			for (std::thread& worker : workers)
				worker.join();
		}

		// Delete copy constructor and assignment operator
		ThreadPool(const ThreadPool&) = delete;
		ThreadPool& operator=(const ThreadPool&) = delete;

		ThreadPool(ThreadPool&&) = delete;
		ThreadPool& operator=(ThreadPool&&) = delete;

		void enqueue(std::function<void()> f) 
		{
			{
				std::lock_guard<std::mutex> lock(queue_mutex);
				tasks.push(std::move(f));
			}
			condition.notify_one();
		}
		
		// P0-2: Fixed - waits for queue empty AND all tasks completed
		void wait_for_tasks() {
			std::unique_lock<std::mutex> lock(queue_mutex);
			tasks_done_condition.wait(lock, [this] { 
				return tasks.empty() && in_flight_tasks.load() == 0; 
			});
		}
		
		bool has_tasks() {
			std::lock_guard<std::mutex> lock(queue_mutex);
			return !tasks.empty() || in_flight_tasks.load() > 0;
		}
		
		size_t pending_tasks() const {
			return in_flight_tasks.load();
		}

	};
} // namespace MML

#endif // MML_THREAD_POOL_H