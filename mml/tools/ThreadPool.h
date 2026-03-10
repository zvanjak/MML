///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ThreadPool.h                                                        ///
///  Description: Thread pool for parallel computation                                ///
///               Task queue with worker threads for concurrent execution             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
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
#include <exception>
#include <stdexcept>

namespace MML
{
	/// @brief A thread pool for parallel task execution
	/// @details ThreadPool manages a fixed number of worker threads that process tasks from a queue.
	/// Tasks are submitted via enqueue() and executed asynchronously. The pool handles thread lifecycle,
	/// task scheduling, and provides synchronization primitives for waiting on task completion.
	/// 
	/// **Thread Safety**: All public methods are thread-safe and can be called from multiple threads.
	/// 
	/// **Exception Handling**: Exceptions thrown in tasks are captured and can be retrieved via
	/// get_last_exception(). An optional error callback can be set for immediate notification.
	/// 
	/// **Usage Patterns**:
	/// 
	/// 1. **Fire-and-forget tasks** (simplest):
	/// @code
	/// ThreadPool pool(4);  // 4 worker threads
	/// pool.enqueue([] { compute_something(); });
	/// pool.enqueue([] { compute_something_else(); });
	/// pool.wait_for_tasks();  // Block until all done
	/// @endcode
	/// 
	/// 2. **Parallel map over data**:
	/// @code
	/// ThreadPool pool(std::thread::hardware_concurrency());
	/// std::vector<double> data(1000);
	/// std::atomic<size_t> index{0};
	/// for (size_t i = 0; i < pool.thread_count(); ++i) {
	///     pool.enqueue([&] {
	///         size_t i;
	///         while ((i = index.fetch_add(1)) < data.size()) {
	///             data[i] = expensive_computation(i);
	///         }
	///     });
	/// }
	/// pool.wait_for_tasks();
	/// @endcode
	/// 
	/// 3. **Task batching with synchronization**:
	/// @code
	/// ThreadPool pool(4);
	/// // Phase 1
	/// for (auto& chunk : phase1_chunks) pool.enqueue([&] { process(chunk); });
	/// pool.wait_for_tasks();  // Barrier
	/// // Phase 2 (depends on phase 1)
	/// for (auto& chunk : phase2_chunks) pool.enqueue([&] { process(chunk); });
	/// pool.wait_for_tasks();
	/// @endcode
	/// 
	/// 4. **With exception handling**:
	/// @code
	/// ThreadPool pool(4);
	/// pool.set_error_callback([](std::exception_ptr e) {
	///     try { std::rethrow_exception(e); }
	///     catch (const std::exception& ex) { std::cerr << "Task failed: " << ex.what() << "\n"; }
	/// });
	/// pool.enqueue([] { throw std::runtime_error("oops"); });
	/// pool.wait_for_tasks();
	/// if (pool.has_exceptions()) {
	///     // Handle accumulated exceptions
	/// }
	/// @endcode
	/// 
	/// @note The pool is non-copyable and non-movable. Destructor waits for all tasks to complete.
	/// @see wait_for_tasks(), has_tasks(), enqueue()
	class ThreadPool {
	private:
		/// Worker threads that process tasks from the queue
		std::vector<std::thread> workers;
		/// Queue of pending tasks awaiting execution
		std::queue<std::function<void()>> tasks;
		/// Mutex protecting the task queue and related state
		mutable std::mutex queue_mutex;
		/// Condition variable for worker thread wake-up
		std::condition_variable condition;
		/// Condition variable for wait_for_tasks() synchronization
		std::condition_variable tasks_done_condition;
		/// Flag indicating pool is shutting down (atomic to prevent data race)
		std::atomic<bool> stop{false};
		/// Count of tasks currently being executed by worker threads
		std::atomic<size_t> in_flight_tasks{0};
		/// Queue of exceptions caught from failed tasks
		std::queue<std::exception_ptr> exceptions;
		/// Mutex protecting the exceptions queue
		mutable std::mutex exceptions_mutex;
		/// Optional callback invoked when a task throws an exception
		std::function<void(std::exception_ptr)> error_callback;
		/// Number of worker threads in the pool
		size_t num_threads;

	public:
		/// @brief Construct a thread pool with the specified number of worker threads
		/// @param numThreads Number of worker threads to create (typically std::thread::hardware_concurrency())
		/// @throws std::invalid_argument if numThreads is 0
		explicit ThreadPool(size_t numThreads) : num_threads(numThreads)
		{
			if (numThreads == 0)
				throw std::invalid_argument("ThreadPool requires at least 1 thread");
				
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
								++this->in_flight_tasks;
							}
							
							// Execute task with exception capture
							try {
								task();
							}
							catch (...) {
								// Capture exception
								std::exception_ptr eptr = std::current_exception();
								{
									std::lock_guard<std::mutex> lock(this->exceptions_mutex);
									this->exceptions.push(eptr);
								}
								// Invoke error callback if set
								if (this->error_callback) {
									try {
										this->error_callback(eptr);
									}
									catch (...) {
										// Ignore exceptions from callback
									}
								}
							}
							
							// Track task completion and notify waiters
							if (--this->in_flight_tasks == 0) {
								std::lock_guard<std::mutex> lock(this->queue_mutex);
								if (this->tasks.empty()) {
									this->tasks_done_condition.notify_all();
								}
							}
						}
					});
		}
		
		/// @brief Destructor - waits for all pending tasks to complete before joining threads
		/// @note Blocks until all queued and in-flight tasks finish execution
		~ThreadPool()
		{
			// Set stop flag and notify under lock, then release BEFORE joining
			{
				std::lock_guard<std::mutex> lock(queue_mutex);
				stop.store(true);
			}
			condition.notify_all();
			
			// Now safe to join - workers can acquire mutex to wake up
			for (std::thread& worker : workers)
				worker.join();
		}

		// Non-copyable and non-movable
		ThreadPool(const ThreadPool&) = delete;
		ThreadPool& operator=(const ThreadPool&) = delete;
		ThreadPool(ThreadPool&&) = delete;
		ThreadPool& operator=(ThreadPool&&) = delete;

		/// @brief Submit a task for asynchronous execution
		/// @param f Callable to execute (typically a lambda or std::function<void()>)
		/// @note Thread-safe. Task will be executed by an available worker thread.
		/// @note Exceptions thrown by the task are captured and can be retrieved via get_exception()
		void enqueue(std::function<void()> f) 
		{
			{
				std::lock_guard<std::mutex> lock(queue_mutex);
				tasks.push(std::move(f));
			}
			condition.notify_one();
		}
		
		/// @brief Block until all queued tasks have completed execution
		/// @details Waits for both the task queue to be empty AND all in-flight tasks to finish.
		/// This provides a synchronization barrier - when wait_for_tasks() returns, all tasks
		/// that were queued before the call have completed.
		/// @note Thread-safe. Can be called from multiple threads simultaneously.
		/// @note Does not prevent new tasks from being enqueued while waiting.
		void wait_for_tasks() {
			std::unique_lock<std::mutex> lock(queue_mutex);
			tasks_done_condition.wait(lock, [this] { 
				return tasks.empty() && in_flight_tasks.load() == 0; 
			});
		}
		
		/// @brief Check if there are pending or executing tasks
		/// @return true if queue is non-empty or tasks are currently executing
		bool has_tasks() const {
			std::lock_guard<std::mutex> lock(queue_mutex);
			return !tasks.empty() || in_flight_tasks.load() > 0;
		}
		
		/// @brief Get the number of tasks currently being executed
		/// @return Count of in-flight tasks (not including queued tasks)
		size_t pending_tasks() const noexcept {
			return in_flight_tasks.load();
		}
		
		/// @brief Get the number of tasks waiting in the queue
		/// @return Count of queued tasks not yet started
		size_t queued_tasks() const {
			std::lock_guard<std::mutex> lock(queue_mutex);
			return tasks.size();
		}
		
		/// @brief Get the number of worker threads in the pool
		/// @return Thread count specified at construction
		size_t thread_count() const noexcept {
			return num_threads;
		}
		
		//=====================================================================================
		// Exception Handling
		//=====================================================================================
		
		/// @brief Set a callback to be invoked when a task throws an exception
		/// @param callback Function taking std::exception_ptr, called immediately when exception occurs
		/// @note The callback is invoked from the worker thread. Keep it short to avoid blocking.
		/// @note Exceptions thrown by the callback itself are silently ignored.
		void set_error_callback(std::function<void(std::exception_ptr)> callback) {
			std::lock_guard<std::mutex> lock(exceptions_mutex);
			error_callback = std::move(callback);
		}
		
		/// @brief Check if any task has thrown an exception
		/// @return true if there are captured exceptions
		bool has_exceptions() const {
			std::lock_guard<std::mutex> lock(exceptions_mutex);
			return !exceptions.empty();
		}
		
		/// @brief Get the number of captured exceptions
		/// @return Count of exceptions from failed tasks
		size_t exception_count() const {
			std::lock_guard<std::mutex> lock(exceptions_mutex);
			return exceptions.size();
		}
		
		/// @brief Retrieve and remove the oldest captured exception
		/// @return The exception_ptr, or nullptr if no exceptions are pending
		/// @note Use std::rethrow_exception() to handle the exception
		std::exception_ptr get_exception() {
			std::lock_guard<std::mutex> lock(exceptions_mutex);
			if (exceptions.empty())
				return nullptr;
			std::exception_ptr eptr = exceptions.front();
			exceptions.pop();
			return eptr;
		}
		
		/// @brief Clear all captured exceptions without processing them
		void clear_exceptions() {
			std::lock_guard<std::mutex> lock(exceptions_mutex);
			while (!exceptions.empty())
				exceptions.pop();
		}
		
		/// @brief Rethrow the first captured exception, clearing it from the queue
		/// @throws The captured exception if one exists
		/// @note Does nothing if no exceptions are pending
		void rethrow_if_exception() {
			std::exception_ptr eptr = get_exception();
			if (eptr)
				std::rethrow_exception(eptr);
		}

	};
} // namespace MML

#endif // MML_THREAD_POOL_H