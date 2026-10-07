///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        threadpool_tests.cpp                                                ///
///  Description: Unit tests for ThreadPool                                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include <atomic>
#include <vector>
#include <thread>
#include <chrono>
#include <stdexcept>
#include <stop_token>

#include <mml/MMLBase.h>
#include <mml/tools/ThreadPool.h>

using namespace MML;

#define TEST_PRECISION_INFO() \
	INFO("Test precision: " << Constants::Eps << " (" << (sizeof(Real) == 8 ? "double" : "float") << ")")

namespace MML::Tests::Tools::ThreadPoolTests {

	/////////////////////////////////////////////////////////////////////////////////////
	///                         BASIC THREADPOOL OPERATIONS                           ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ThreadPool - Basic construction", "[threadpool][basic]") {
		TEST_PRECISION_INFO();

		SECTION("ThreadPool can be constructed with various thread counts") {
			{
				ThreadPool pool1(1);
				// Should not throw
				REQUIRE(true);
			}
			{
				ThreadPool pool2(2);
				REQUIRE(true);
			}
			{
				ThreadPool pool4(4);
				REQUIRE(true);
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         TASK EXECUTION                                        ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ThreadPool - Task execution", "[threadpool][execution]") {
		TEST_PRECISION_INFO();

		SECTION("Single task executes correctly") {
			ThreadPool pool(2);
			
			std::atomic<int> result{0};
			
			pool.enqueue([&result]() {
				result.store(42);
			});
			
			pool.wait_for_tasks();
			
			REQUIRE(result.load() == 42);
		}

		SECTION("Multiple tasks execute correctly") {
			ThreadPool pool(4);
			
			std::atomic<int> counter{0};
			constexpr int NUM_TASKS = 100;
			
			for (int i = 0; i < NUM_TASKS; ++i) {
				pool.enqueue([&counter]() {
					counter.fetch_add(1);
				});
			}
			
			pool.wait_for_tasks();
			
			REQUIRE(counter.load() == NUM_TASKS);
		}

		SECTION("Tasks with computation") {
			ThreadPool pool(2);
			
			std::atomic<Real> sum{0.0};
			constexpr int NUM_TASKS = 50;
			
			for (int i = 0; i < NUM_TASKS; ++i) {
				pool.enqueue([&sum, i]() {
					// Each task adds its index
					Real current = sum.load();
					while (!sum.compare_exchange_weak(current, current + static_cast<Real>(i))) {
						// Retry
					}
				});
			}
			
			pool.wait_for_tasks();
			
			// Sum of 0..49 = 49*50/2 = 1225
			REQUIRE(std::abs(sum.load() - 1225.0) < 0.1);
		}

		SECTION("Batch tasks execute correctly") {
			ThreadPool pool(4);
			std::atomic<int> counter{0};
			std::vector<std::function<void()>> batch;
			constexpr int NUM_TASKS = 100;

			for (int i = 0; i < NUM_TASKS; ++i) {
				batch.emplace_back([&counter]() {
					counter.fetch_add(1);
				});
			}

			pool.enqueue_batch(std::move(batch));
			pool.wait_for_tasks();

			REQUIRE(counter.load() == NUM_TASKS);
		}

		SECTION("Empty batch is accepted") {
			ThreadPool pool(2);
			std::vector<std::function<void()>> batch;

			REQUIRE_NOTHROW(pool.enqueue_batch(std::move(batch)));
			REQUIRE_NOTHROW(pool.wait_for_tasks());
			REQUIRE_FALSE(pool.has_tasks());
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         WAIT AND STATUS                                       ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ThreadPool - Wait and status", "[threadpool][wait]") {
		TEST_PRECISION_INFO();

		SECTION("wait_for_tasks blocks until complete") {
			ThreadPool pool(2);
			
			std::atomic<bool> taskCompleted{false};
			
			pool.enqueue([&taskCompleted]() {
				std::this_thread::sleep_for(std::chrono::milliseconds(50));
				taskCompleted.store(true);
			});
			
			pool.wait_for_tasks();
			
			// After wait returns, task must be complete
			REQUIRE(taskCompleted.load() == true);
		}

		SECTION("has_tasks reflects pending work") {
			ThreadPool pool(1);  // Single thread to control execution order
			
			std::atomic<bool> taskStarted{false};
			std::atomic<bool> taskCanFinish{false};
			
			pool.enqueue([&taskStarted, &taskCanFinish]() {
				taskStarted.store(true);
				while (!taskCanFinish.load()) {
					std::this_thread::sleep_for(std::chrono::milliseconds(1));
				}
			});
			
			// Wait for task to start
			while (!taskStarted.load()) {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
			}
			
			// While task is running, has_tasks should be true
			REQUIRE(pool.has_tasks() == true);
			
			// Let task finish
			taskCanFinish.store(true);
			pool.wait_for_tasks();
			
			// After completion, no pending tasks
			REQUIRE(pool.has_tasks() == false);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         CONCURRENT EXECUTION                                  ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ThreadPool - Concurrent execution", "[threadpool][concurrent]") {
		TEST_PRECISION_INFO();

		SECTION("Tasks execute concurrently") {
			ThreadPool pool(4);
			
			std::atomic<int> maxConcurrent{0};
			std::atomic<int> currentActive{0};
			constexpr int NUM_TASKS = 20;
			
			for (int i = 0; i < NUM_TASKS; ++i) {
				pool.enqueue([&maxConcurrent, &currentActive]() {
					int active = currentActive.fetch_add(1) + 1;
					
					// Track max concurrent
					int prevMax = maxConcurrent.load();
					while (active > prevMax && !maxConcurrent.compare_exchange_weak(prevMax, active)) {
						// Retry
					}
					
					// Simulate work
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
					
					currentActive.fetch_sub(1);
				});
			}
			
			pool.wait_for_tasks();
			
			// With 4 threads and enough tasks, we should have seen some concurrency
			REQUIRE(maxConcurrent.load() >= 2);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         EDGE CASES                                            ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ThreadPool - Edge cases", "[threadpool][edge]") {
		TEST_PRECISION_INFO();

		SECTION("wait_for_tasks with no tasks is immediate") {
			ThreadPool pool(2);
			
			auto start = std::chrono::steady_clock::now();
			pool.wait_for_tasks();
			auto elapsed = std::chrono::steady_clock::now() - start;
			
			// Should be nearly immediate (< 10ms)
			REQUIRE(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() < 10);
		}

		SECTION("Pool handles many small tasks") {
			ThreadPool pool(4);
			
			std::atomic<int> counter{0};
			constexpr int NUM_TASKS = 1000;
			
			for (int i = 0; i < NUM_TASKS; ++i) {
				pool.enqueue([&counter]() {
					counter.fetch_add(1);
				});
			}
			
			pool.wait_for_tasks();
			
			REQUIRE(counter.load() == NUM_TASKS);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         EXCEPTION CALLBACKS                                    ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ThreadPool - Error callback can re-enter pool without deadlock", "[threadpool][exceptions]") {
		TEST_PRECISION_INFO();

		ThreadPool pool(2);
		std::atomic<bool> callbackObservedException{false};
		std::atomic<bool> recoveryTaskRan{false};

		pool.set_error_callback([&pool, &callbackObservedException, &recoveryTaskRan](std::exception_ptr) {
			callbackObservedException.store(pool.exception_count() > 0);
			pool.enqueue([&recoveryTaskRan]() {
				recoveryTaskRan.store(true);
			});
		});

		pool.enqueue([]() {
			throw std::runtime_error("expected ThreadPool test failure");
		});

		pool.wait_for_tasks();

		REQUIRE(callbackObservedException.load());
		REQUIRE(recoveryTaskRan.load());
		REQUIRE(pool.exception_count() == 1);
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         CANCELLATION                                            ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ThreadPool - Cancellable queued task is skipped when stop requested", "[threadpool][cancellation]") {
		TEST_PRECISION_INFO();

		ThreadPool pool(1);
		std::stop_source stopSource;
		std::atomic<bool> blockerCanFinish{false};
		std::atomic<bool> cancelledTaskRan{false};

		pool.enqueue([&blockerCanFinish]() {
			while (!blockerCanFinish.load()) {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
			}
		});

		pool.enqueue(stopSource.get_token(), [&cancelledTaskRan](std::stop_token) {
			cancelledTaskRan.store(true);
		});

		stopSource.request_stop();
		blockerCanFinish.store(true);
		pool.wait_for_tasks();

		REQUIRE_FALSE(cancelledTaskRan.load());
	}

	TEST_CASE("ThreadPool - Cancellable batch skips tasks when stop requested", "[threadpool][cancellation]") {
		TEST_PRECISION_INFO();

		ThreadPool pool(1);
		std::stop_source stopSource;
		std::atomic<bool> blockerCanFinish{false};
		std::atomic<int> cancelledTasksRan{0};
		std::vector<std::function<void(std::stop_token)>> batch;

		pool.enqueue([&blockerCanFinish]() {
			while (!blockerCanFinish.load()) {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
			}
		});

		for (int i = 0; i < 3; ++i) {
			batch.emplace_back([&cancelledTasksRan](std::stop_token) {
				cancelledTasksRan.fetch_add(1);
			});
		}

		pool.enqueue_batch(stopSource.get_token(), std::move(batch));
		stopSource.request_stop();
		blockerCanFinish.store(true);
		pool.wait_for_tasks();

		REQUIRE(cancelledTasksRan.load() == 0);
	}

	TEST_CASE("ThreadPool - Running cancellable task observes stop token", "[threadpool][cancellation]") {
		TEST_PRECISION_INFO();

		ThreadPool pool(1);
		std::stop_source stopSource;
		std::atomic<bool> taskStarted{false};
		std::atomic<bool> taskObservedStop{false};

		pool.enqueue(stopSource.get_token(), [&taskStarted, &taskObservedStop](std::stop_token token) {
			taskStarted.store(true);
			while (!token.stop_requested()) {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
			}
			taskObservedStop.store(true);
		});

		while (!taskStarted.load()) {
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}

		stopSource.request_stop();
		pool.wait_for_tasks();

		REQUIRE(taskObservedStop.load());
	}

} // namespace MML::Tests::Tools::ThreadPoolTests

