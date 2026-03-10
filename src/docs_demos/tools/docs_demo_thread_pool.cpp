///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: ThreadPool - Parallel Task Execution
///////////////////////////////////////////////////////////////////////////////////////////

#include "MMLBase.h"
#include "mml/tools/ThreadPool.h"

#include <iostream>
#include <vector>
#include <atomic>
#include <cmath>
#include <chrono>
#include <thread>
#include <numeric>

using namespace MML;

namespace MML::docs_demos::thread_pool
{
    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              BASIC USAGE                                         ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_basic_usage()
    {
        std::cout << "\n=== Basic ThreadPool Usage ===\n";

        // Create pool with hardware threads
        size_t numThreads = std::min(size_t(4), size_t(std::thread::hardware_concurrency()));
        ThreadPool pool(numThreads);

        std::cout << "Created pool with " << pool.thread_count() << " threads\n";

        // Enqueue simple tasks
        std::atomic<int> counter{0};

        for (int i = 0; i < 10; ++i) {
            pool.enqueue([&counter, i]() {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
                counter.fetch_add(1);
            });
        }

        // Wait for all tasks to complete
        pool.wait_for_tasks();

        std::cout << "All tasks completed. Counter = " << counter.load() << "\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              PARALLEL COMPUTATION                                ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_parallel_computation()
    {
        std::cout << "\n=== Parallel Computation ===\n";

        const size_t dataSize = 10000;
        std::vector<double> data(dataSize);
        std::vector<double> results(dataSize);

        // Initialize data
        for (size_t i = 0; i < dataSize; ++i) {
            data[i] = static_cast<double>(i);
        }

        size_t numThreads = std::min(size_t(4), size_t(std::thread::hardware_concurrency()));
        ThreadPool pool(numThreads);

        // Parallel map: compute sqrt for each element
        std::atomic<size_t> index{0};

        auto start = std::chrono::steady_clock::now();

        for (size_t t = 0; t < numThreads; ++t) {
            pool.enqueue([&]() {
                size_t i;
                while ((i = index.fetch_add(1)) < dataSize) {
                    results[i] = std::sin(data[i] * 0.001) * std::cos(data[i] * 0.001);
                }
            });
        }

        pool.wait_for_tasks();

        auto end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double, std::milli>(end - start).count();

        // Verify
        double sum = 0.0;
        for (const auto& r : results) sum += r;

        std::cout << "Processed " << dataSize << " elements in " << elapsed << " ms\n";
        std::cout << "Sum of results: " << sum << "\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              TASK BATCHING                                       ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_task_batching()
    {
        std::cout << "\n=== Task Batching with Synchronization ===\n";

        size_t numThreads = std::min(size_t(4), size_t(std::thread::hardware_concurrency()));
        ThreadPool pool(numThreads);

        std::vector<double> phase1_results(4);
        std::vector<double> phase2_results(4);

        // Phase 1: Compute initial values
        std::cout << "Starting Phase 1...\n";
        for (size_t i = 0; i < 4; ++i) {
            pool.enqueue([&phase1_results, i]() {
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
                phase1_results[i] = std::sqrt(static_cast<double>(i + 1));
            });
        }
        pool.wait_for_tasks();  // Barrier
        std::cout << "Phase 1 complete. Results: ";
        for (const auto& r : phase1_results) std::cout << r << " ";
        std::cout << "\n";

        // Phase 2: Use Phase 1 results
        std::cout << "Starting Phase 2 (depends on Phase 1)...\n";
        for (size_t i = 0; i < 4; ++i) {
            pool.enqueue([&phase1_results, &phase2_results, i]() {
                std::this_thread::sleep_for(std::chrono::milliseconds(15));
                phase2_results[i] = phase1_results[i] * phase1_results[i];
            });
        }
        pool.wait_for_tasks();  // Barrier
        std::cout << "Phase 2 complete. Results: ";
        for (const auto& r : phase2_results) std::cout << r << " ";
        std::cout << "\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              EXCEPTION HANDLING                                  ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_exception_handling()
    {
        std::cout << "\n=== Exception Handling ===\n";

        size_t numThreads = 2;
        ThreadPool pool(numThreads);

        // Set error callback
        pool.set_error_callback([](std::exception_ptr eptr) {
            try {
                std::rethrow_exception(eptr);
            } catch (const std::exception& e) {
                std::cout << "  [Error callback] Task failed: " << e.what() << "\n";
            }
        });

        // Enqueue tasks, some will throw
        pool.enqueue([]() { 
            std::cout << "  Task 1: Success\n"; 
        });

        pool.enqueue([]() {
            throw std::runtime_error("Simulated error in task 2");
        });

        pool.enqueue([]() { 
            std::cout << "  Task 3: Success\n"; 
        });

        pool.wait_for_tasks();

        // Check for accumulated exceptions
        if (pool.has_exceptions()) {
            std::cout << "Pool has unhandled exceptions\n";
            pool.clear_exceptions();
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              NUMERICAL INTEGRATION                               ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_parallel_integration()
    {
        std::cout << "\n=== Parallel Numerical Integration ===\n";

        // Integrate sin(x) from 0 to π using parallel trapezoidal rule
        const int totalPoints = 100000;
        const double a = 0.0, b = Constants::PI;
        const double h = (b - a) / totalPoints;

        size_t numThreads = std::min(size_t(4), size_t(std::thread::hardware_concurrency()));
        ThreadPool pool(numThreads);

        std::vector<double> partialSums(numThreads, 0.0);
        int pointsPerThread = totalPoints / static_cast<int>(numThreads);

        auto start = std::chrono::steady_clock::now();

        for (size_t t = 0; t < numThreads; ++t) {
            pool.enqueue([t, numThreads, pointsPerThread, totalPoints, a, h, &partialSums]() {
                int startIdx = static_cast<int>(t) * pointsPerThread;
                int endIdx = (t == numThreads - 1) ? totalPoints : startIdx + pointsPerThread;
                
                double sum = 0.0;
                for (int i = startIdx; i < endIdx; ++i) {
                    double x = a + i * h;
                    sum += std::sin(x);
                }
                partialSums[t] = sum;
            });
        }

        pool.wait_for_tasks();

        auto end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double, std::milli>(end - start).count();

        // Combine results
        double totalSum = 0.0;
        for (const auto& ps : partialSums) totalSum += ps;
        double integral = (totalSum - 0.5 * (std::sin(a) + std::sin(b))) * h;

        std::cout << "Integral of sin(x) from 0 to π:\n";
        std::cout << "  Computed: " << integral << "\n";
        std::cout << "  Exact:    2.0\n";
        std::cout << "  Error:    " << std::abs(integral - 2.0) << "\n";
        std::cout << "  Time:     " << elapsed << " ms\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              POOL STATUS                                         ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_pool_status()
    {
        std::cout << "\n=== Pool Status Queries ===\n";

        size_t numThreads = 2;
        ThreadPool pool(numThreads);

        std::cout << "Thread count: " << pool.thread_count() << "\n";
        std::cout << "Has tasks (before): " << (pool.has_tasks() ? "yes" : "no") << "\n";

        // Add some tasks
        for (int i = 0; i < 5; ++i) {
            pool.enqueue([i]() {
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            });
        }

        std::cout << "Has tasks (during): " << (pool.has_tasks() ? "yes" : "no") << "\n";

        pool.wait_for_tasks();

        std::cout << "Has tasks (after): " << (pool.has_tasks() ? "yes" : "no") << "\n";
    }

    void Run()
    {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "       MML Documentation Demo: ThreadPool\n";
        std::cout << std::string(70, '=') << "\n";

        demo_basic_usage();
        demo_parallel_computation();
        demo_task_batching();
        demo_exception_handling();
        demo_parallel_integration();
        demo_pool_status();

        std::cout << "\n" << std::string(70, '=') << "\n";
    }
} // namespace

void Docs_Demo_ThreadPool() {
    MML::docs_demos::thread_pool::Run();
}
