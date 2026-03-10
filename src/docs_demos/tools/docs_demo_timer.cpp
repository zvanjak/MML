///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: Timer - Performance Measurement
///////////////////////////////////////////////////////////////////////////////////////////

#include "MMLBase.h"
#include "tools/Timer.h"

#include <iostream>
#include <cmath>
#include <thread>
#include <chrono>
#include <vector>

using namespace MML;

namespace MML::docs_demos::timer
{
    // Simulated work functions
    void simulate_work(int ms) {
        std::this_thread::sleep_for(std::chrono::milliseconds(ms));
    }

    double compute_intensive() {
        double sum = 0.0;
        for (int i = 0; i < 100000; ++i) {
            sum += std::sin(static_cast<double>(i) * 0.001);
        }
        return sum;
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              TIMER CLASS                                         ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_basic_timer()
    {
        std::cout << "\n=== Basic Timer Usage ===\n";

        Timer timer;
        timer.Start();

        // Phase 1: Initialization
        simulate_work(50);
        timer.MarkTime("Initialization");

        // Phase 2: Computation
        double result = compute_intensive();
        timer.MarkTime("Computation");

        // Phase 3: Cleanup
        simulate_work(30);
        timer.MarkTime("Cleanup");

        // Print results
        std::cout << "\nResult: " << result << "\n";
        timer.Print();
    }

    void demo_interval_access()
    {
        std::cout << "\n=== Accessing Individual Intervals ===\n";

        Timer timer;
        timer.Start();

        simulate_work(25);
        timer.MarkTime("Step 1");

        simulate_work(75);
        timer.MarkTime("Step 2");

        simulate_work(50);
        timer.MarkTime("Step 3");

        // Access individual intervals
        std::cout << "Interval 0 (Start->Step1): " << timer.GetIntervalTime(0) * 1000 << " ms\n";
        std::cout << "Interval 1 (Step1->Step2): " << timer.GetIntervalTime(1) * 1000 << " ms\n";
        std::cout << "Interval 2 (Step2->Step3): " << timer.GetIntervalTime(2) * 1000 << " ms\n";
        std::cout << "Total elapsed: " << timer.GetTotalTime() * 1000 << " ms\n";
        std::cout << "Mark count: " << timer.GetMarkCount() << "\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              SCOPED TIMER                                        ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_scoped_timer()
    {
        std::cout << "\n=== Scoped Timer (RAII) ===\n";

        // ScopedTimer automatically prints on destruction
        {
            ScopedTimer st("Matrix multiplication");
            compute_intensive();
        } // Timer prints here

        {
            ScopedTimer st("I/O simulation");
            simulate_work(100);
        } // Timer prints here
    }

    void demo_nested_timing()
    {
        std::cout << "\n=== Nested Timing ===\n";

        ScopedTimer outer("Total operation");

        {
            ScopedTimer inner1("Phase 1");
            simulate_work(40);
        }

        {
            ScopedTimer inner2("Phase 2");
            simulate_work(60);
        }

        {
            ScopedTimer inner3("Phase 3");
            simulate_work(50);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              PERFORMANCE COMPARISON                              ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_algorithm_comparison()
    {
        std::cout << "\n=== Algorithm Performance Comparison ===\n";

        const int N = 1000;
        std::vector<double> data(N);

        Timer timer;
        timer.Start();

        // Algorithm 1: Direct summation
        double sum1 = 0.0;
        for (int i = 0; i < N; ++i) {
            sum1 += std::sin(static_cast<double>(i) * 0.01);
        }
        timer.MarkTime("Direct sum");

        // Algorithm 2: Kahan summation
        double sum2 = 0.0, c = 0.0;
        for (int i = 0; i < N; ++i) {
            double y = std::sin(static_cast<double>(i) * 0.01) - c;
            double t = sum2 + y;
            c = (t - sum2) - y;
            sum2 = t;
        }
        timer.MarkTime("Kahan sum");

        // Algorithm 3: Pairwise summation (simplified)
        std::vector<double> vals(N);
        for (int i = 0; i < N; ++i) {
            vals[i] = std::sin(static_cast<double>(i) * 0.01);
        }
        double sum3 = 0.0;
        for (const auto& v : vals) sum3 += v;
        timer.MarkTime("Pairwise sum");

        std::cout << "\nResults:\n";
        std::cout << "  Direct:   " << sum1 << "\n";
        std::cout << "  Kahan:    " << sum2 << "\n";
        std::cout << "  Pairwise: " << sum3 << "\n";

        timer.Print();
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              TIMING UTILITIES                                    ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_timing_wrapper()
    {
        std::cout << "\n=== Timing a Function ===\n";

        // Time a lambda function
        auto task = []() -> double {
            double sum = 0.0;
            for (int i = 0; i < 50000; ++i) {
                sum += std::sqrt(static_cast<double>(i));
            }
            return sum;
        };

        Timer timer;
        timer.Start();
        double result = task();
        timer.MarkTime("Task execution");

        std::cout << "Result: " << result << "\n";
        std::cout << "Time: " << timer.GetIntervalTime(0) * 1000 << " ms\n";
    }

    void demo_repeated_timing()
    {
        std::cout << "\n=== Repeated Timing (Averaging) ===\n";

        const int iterations = 5;
        std::vector<double> times;
        times.reserve(iterations);

        for (int i = 0; i < iterations; ++i) {
            Timer timer;
            timer.Start();
            
            compute_intensive();
            
            timer.MarkTime("Iteration");
            times.push_back(timer.GetIntervalTime(0) * 1000);
        }

        // Calculate statistics
        double total = 0.0;
        for (double t : times) total += t;
        double avg = total / iterations;

        double variance = 0.0;
        for (double t : times) variance += (t - avg) * (t - avg);
        double stddev = std::sqrt(variance / iterations);

        std::cout << "Iterations: " << iterations << "\n";
        std::cout << "Average time: " << avg << " ms\n";
        std::cout << "Std deviation: " << stddev << " ms\n";
        std::cout << "Min: " << *std::min_element(times.begin(), times.end()) << " ms\n";
        std::cout << "Max: " << *std::max_element(times.begin(), times.end()) << " ms\n";
    }

    void Run()
    {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "       MML Documentation Demo: Timer\n";
        std::cout << std::string(70, '=') << "\n";

        demo_basic_timer();
        demo_interval_access();
        demo_scoped_timer();
        demo_nested_timing();
        demo_algorithm_comparison();
        demo_timing_wrapper();
        demo_repeated_timing();

        std::cout << "\n" << std::string(70, '=') << "\n";
    }
} // namespace

void Docs_Demo_Timer() {
    MML::docs_demos::timer::Run();
}
