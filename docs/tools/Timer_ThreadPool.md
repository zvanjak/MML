# Timer & ThreadPool - Performance and Parallelization

## Overview

Two essential utility classes for performance measurement and parallel computation:

- **Timer**: High-precision performance measurement with multiple time markers
- **ThreadPool**: Thread pool for parallel task execution

---

## Highlights
- `Timer`: High-resolution timing, labeled checkpoints, interval and total times, formatted output.
- `ThreadPool`: Fixed-size pool, task queue with futures, safe shutdown and synchronization.
- Designed for profiling numerical algorithms and parallelizing workloads.
- Clean APIs integrate with ConsolePrinter/Serializer for reporting.

## Timer Class

### Purpose

Measure execution time with multiple checkpoints for profiling algorithms and comparing performance.

### Core Features

- High-resolution timing using `std::chrono::steady_clock`
- Multiple time markers with labels
- Interval and cumulative time measurement
- Automatic formatted output

### API Reference

#### Construction & Initialization

```cpp
Timer timer;              // Create timer
timer.Start();            // Start timing
```

#### Marking Time Points

```cpp
void MarkTime();                        // Unnamed mark
void MarkTime(const std::string& name); // Named mark
```

#### Retrieving Time Data

```cpp
double GetIntervalTime(int index = 0) const;      // Time since previous mark
double GetMarkTimeFromStart(int index = 0) const; // Cumulative time to mark
double GetTotalTime() const;                      // Total elapsed time
```

**Throws:** `std::out_of_range` if index invalid, `std::runtime_error` if no marks exist.

#### Display Results

```cpp
void Print();  // Formatted console output
```

### Usage Examples

#### Basic Timing

```cpp
#include "tools/Timer.h"

Timer timer;
timer.Start();

// Your code here
performComputation();

timer.MarkTime("Computation");
timer.Print();

// Output:
// Timer Results:
// Mark 0 (            Computation): 1.23456 seconds
// Total Time: 1.23456 seconds
```

#### Multi-Stage Profiling

```cpp
Timer timer;
timer.Start();

loadData();
timer.MarkTime("Data Loading");

preprocessData();
timer.MarkTime("Preprocessing");

runAlgorithm();
timer.MarkTime("Algorithm");

saveResults();
timer.MarkTime("Save Results");

timer.Print();

// Output:
// Timer Results:
// Mark 0 (          Data Loading): 0.523 seconds
// Mark 1 (         Preprocessing): 0.145 seconds
// Mark 2 (              Algorithm): 2.891 seconds
// Mark 3 (          Save Results): 0.312 seconds
// Total Time: 3.871 seconds
```

#### Interval vs Cumulative Time

```cpp
timer.Start();

step1();
timer.MarkTime("Step 1");

step2();
timer.MarkTime("Step 2");

step3();
timer.MarkTime("Step 3");

// Interval times (each step duration)
double t1 = timer.GetIntervalTime(0);  // Step 1 duration
double t2 = timer.GetIntervalTime(1);  // Step 2 duration  
double t3 = timer.GetIntervalTime(2);  // Step 3 duration

// Cumulative times (from start)
double cum1 = timer.GetMarkTimeFromStart(0);  // After step 1
double cum2 = timer.GetMarkTimeFromStart(1);  // After step 1+2
double cum3 = timer.GetMarkTimeFromStart(2);  // After step 1+2+3
```

#### Algorithm Comparison

```cpp
#include "algorithms/RootFinding.h"

Timer timer;

// Test bisection
timer.Start();
for (int i = 0; i < 1000; ++i) {
    FindRootBisection(f, a, b, tol);
}
timer.MarkTime("Bisection (1000x)");

// Test Newton-Raphson
for (int i = 0; i < 1000; ++i) {
    FindRootNewton(f, df, x0, tol);
}
timer.MarkTime("Newton (1000x)");

timer.Print();

// Compare performance
double bisectionTime = timer.GetIntervalTime(0);
double newtonTime = timer.GetIntervalTime(1);

std::cout << "Newton is " << (bisectionTime / newtonTime) 
          << "x faster\n";
```

#### Matrix Operation Profiling

```cpp
Timer timer;
timer.Start();

Matrix<Real> A = generateRandomMatrix(1000, 1000);
timer.MarkTime("Matrix Generation");

Matrix<Real> B = A.Transpose();
timer.MarkTime("Transpose");

Matrix<Real> C = A * B;
timer.MarkTime("Multiplication");

auto [L, U] = A.LUDecomposition();
timer.MarkTime("LU Decomposition");

timer.Print();
```

### Best Practices

**1. Always call Start() before MarkTime()**
```cpp
Timer timer;
timer.Start();              // Required!
timer.MarkTime("Work");
```

**2. Use descriptive names**
```cpp
// Poor
timer.MarkTime("1");
timer.MarkTime("2");

// Good
timer.MarkTime("Load Data");
timer.MarkTime("Solve System");
```

**3. Create new Timer for separate measurements**
```cpp
// Don't reuse same timer for unrelated measurements
Timer parseTimer;
parseTimer.Start();
parseFile();
parseTimer.MarkTime("Parse");

Timer computeTimer;  // New timer for different task
computeTimer.Start();
compute();
computeTimer.MarkTime("Compute");
```

**4. Exception safety**
```cpp
try {
    double interval = timer.GetIntervalTime(5);
} catch (const std::out_of_range& e) {
    std::cerr << "Invalid mark index: " << e.what() << '\n';
}
```

### Performance Characteristics

- **Overhead**: ~100-500 nanoseconds per MarkTime() call
- **Resolution**: Nanosecond precision (platform-dependent)
- **Accuracy**: Suitable for operations > 1 microsecond
- **Memory**: ~40 bytes per mark (name + timestamp)

---

## ThreadPool Class

### Purpose

Manage a pool of worker threads for parallel task execution, avoiding thread creation overhead.

### Core Features

- Fixed-size thread pool
- Task queue with FIFO execution
- Thread-safe enqueueing
- Wait for completion
- RAII design (automatic cleanup)

### API Reference

#### Construction

```cpp
ThreadPool(size_t numThreads);  // Create pool with N threads
```

**Note:** Non-copyable, non-movable (deleted copy/move constructors).

#### Task Management

```cpp
void enqueue(std::function<void()> task);  // Add task to queue
void wait_for_tasks();                      // Block until queue empty
bool has_tasks();                           // Check if tasks pending
```

#### Destruction

```cpp
~ThreadPool();  // Joins all threads, waits for tasks to complete
```

### Usage Examples

#### Basic Parallel Execution

```cpp
#include "tools/ThreadPool.h"

ThreadPool pool(4);  // 4 worker threads

// Enqueue tasks
for (int i = 0; i < 100; ++i) {
    pool.enqueue([i]() {
        processData(i);
    });
}

// ThreadPool destructor waits for all tasks
```

#### Parallel Matrix Operations

```cpp
void parallelMatrixMultiply(const Matrix<Real>& A, 
                           const Matrix<Real>& B,
                           Matrix<Real>& result) {
    int rows = A.RowNum();
    int cols = B.ColNum();
    
    ThreadPool pool(std::thread::hardware_concurrency());
    
    // Each thread computes several rows
    for (int i = 0; i < rows; ++i) {
        pool.enqueue([&A, &B, &result, i, cols]() {
            for (int j = 0; j < cols; ++j) {
                Real sum = 0.0;
                for (int k = 0; k < A.ColNum(); ++k) {
                    sum += A(i, k) * B(k, j);
                }
                result(i, j) = sum;
            }
        });
    }
    // Destructor waits for completion
}
```

#### Parallel Function Evaluation

```cpp
// Evaluate function at many points in parallel
std::vector<Real> parallelEvaluate(
    const IRealFunction& f,
    const std::vector<Real>& points) 
{
    std::vector<Real> results(points.size());
    ThreadPool pool(8);
    
    for (size_t i = 0; i < points.size(); ++i) {
        pool.enqueue([&f, &points, &results, i]() {
            results[i] = f(points[i]);
        });
    }
    
    return results;  // Destructor ensures completion
}
```

#### Parallel Integration (Monte Carlo)

```cpp
Real parallelMonteCarloIntegration(
    const IScalarFunction<2>& f,
    Real xMin, Real xMax,
    Real yMin, Real yMax,
    size_t numSamples) 
{
    const int numThreads = 8;
    ThreadPool pool(numThreads);
    
    std::vector<Real> partialSums(numThreads, 0.0);
    size_t samplesPerThread = numSamples / numThreads;
    
    for (int t = 0; t < numThreads; ++t) {
        pool.enqueue([&, t]() {
            std::mt19937 rng(t);  // Separate RNG per thread
            std::uniform_real_distribution<Real> distX(xMin, xMax);
            std::uniform_real_distribution<Real> distY(yMin, yMax);
            
            Real sum = 0.0;
            for (size_t i = 0; i < samplesPerThread; ++i) {
                VectorN<Real, 2> point{distX(rng), distY(rng)};
                sum += f(point);
            }
            partialSums[t] = sum;
        });
    }
    
    pool.wait_for_tasks();  // Explicit wait
    
    Real totalSum = std::accumulate(partialSums.begin(), 
                                     partialSums.end(), 0.0);
    Real area = (xMax - xMin) * (yMax - yMin);
    return area * totalSum / numSamples;
}
```

#### Task Dependencies (Advanced)

```cpp
// Use wait_for_tasks() for dependencies
ThreadPool pool(4);

// Phase 1: Independent tasks
pool.enqueue([]() { taskA(); });
pool.enqueue([]() { taskB(); });
pool.enqueue([]() { taskC(); });

pool.wait_for_tasks();  // Wait for phase 1

// Phase 2: Depends on phase 1 results
pool.enqueue([]() { taskD(); });
pool.enqueue([]() { taskE(); });

// Destructor waits for phase 2
```

#### Progress Monitoring

```cpp
ThreadPool pool(4);
std::atomic<int> completed{0};
const int totalTasks = 100;

for (int i = 0; i < totalTasks; ++i) {
    pool.enqueue([&completed, i]() {
        processTask(i);
        completed++;
    });
}

// Monitor progress
while (completed < totalTasks) {
    std::cout << "Progress: " << completed << "/" << totalTasks 
              << "\r" << std::flush;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
}
std::cout << "\nDone!\n";
```

### Best Practices

**1. Choose thread count wisely**
```cpp
// CPU-bound tasks: Use hardware concurrency
ThreadPool pool(std::thread::hardware_concurrency());

// I/O-bound tasks: Can use more threads
ThreadPool ioPool(2 * std::thread::hardware_concurrency());

// Limited resources: Fixed small count
ThreadPool limitedPool(4);
```

**2. Avoid data races**
```cpp
// BAD: Race condition
int counter = 0;
pool.enqueue([&counter]() { counter++; });  // Unsafe!

// GOOD: Use atomic or mutex
std::atomic<int> counter{0};
pool.enqueue([&counter]() { counter++; });  // Safe

// Or use separate storage per thread
std::vector<int> results(numThreads);
pool.enqueue([&results, threadId]() {
    results[threadId] = compute();  // No race
});
```

**3. Capture by value when possible**
```cpp
// Risky: Capture reference to temporary
for (int i = 0; i < 10; ++i) {
    pool.enqueue([&i]() { work(i); });  // BAD: i changes!
}

// Safe: Capture by value
for (int i = 0; i < 10; ++i) {
    pool.enqueue([i]() { work(i); });   // GOOD: i copied
}
```

**4. Exception handling in tasks**
```cpp
pool.enqueue([]() {
    try {
        riskyOperation();
    } catch (const std::exception& e) {
        std::cerr << "Task failed: " << e.what() << '\n';
    }
});
```

**5. Avoid pool in pool**
```cpp
// DON'T create nested thread pools
pool.enqueue([]() {
    ThreadPool innerPool(4);  // BAD: Oversubscription!
    // ...
});

// Instead, reuse existing pool or use sequential processing
```

### Performance Characteristics

- **Thread creation cost**: Amortized over many tasks
- **Task overhead**: ~1-10 microseconds per enqueue
- **Scalability**: Linear speedup for CPU-bound, independent tasks
- **Optimal pool size**: 
  - CPU-bound: `std::thread::hardware_concurrency()`
  - I/O-bound: 2-4× hardware concurrency
  - Mixed workload: Benchmark and tune

### Common Patterns

#### Map-Reduce

```cpp
// Map phase
std::vector<int> data(1000);
std::vector<int> mapped(1000);

ThreadPool pool(8);
for (size_t i = 0; i < data.size(); ++i) {
    pool.enqueue([&, i]() {
        mapped[i] = transform(data[i]);
    });
}
pool.wait_for_tasks();

// Reduce phase (sequential or parallel chunks)
int result = std::accumulate(mapped.begin(), mapped.end(), 0);
```

#### Pipeline Processing

```cpp
std::queue<Data> stage1Queue;
std::queue<Data> stage2Queue;
std::mutex queueMutex;

ThreadPool pool(6);

// Stage 1 workers
for (int i = 0; i < 2; ++i) {
    pool.enqueue([&]() {
        while (hasWork()) {
            auto data = fetchInput();
            auto processed = stage1Process(data);
            
            std::lock_guard<std::mutex> lock(queueMutex);
            stage1Queue.push(processed);
        }
    });
}

// Stage 2 workers
for (int i = 0; i < 4; ++i) {
    pool.enqueue([&]() {
        while (hasWork()) {
            Data data;
            {
                std::lock_guard<std::mutex> lock(queueMutex);
                if (stage1Queue.empty()) continue;
                data = stage1Queue.front();
                stage1Queue.pop();
            }
            auto result = stage2Process(data);
            saveResult(result);
        }
    });
}
```

---

## Combining Timer & ThreadPool

### Parallel Performance Analysis

```cpp
void benchmarkParallelization(const std::vector<Real>& data) {
    Timer timer;
    
    // Sequential
    timer.Start();
    for (auto& x : data) {
        heavyComputation(x);
    }
    timer.MarkTime("Sequential");
    
    // Parallel with 2 threads
    ThreadPool pool2(2);
    for (auto& x : data) {
        pool2.enqueue([&x]() { heavyComputation(x); });
    }
    pool2.wait_for_tasks();
    timer.MarkTime("Parallel (2 threads)");
    
    // Parallel with 4 threads
    ThreadPool pool4(4);
    for (auto& x : data) {
        pool4.enqueue([&x]() { heavyComputation(x); });
    }
    pool4.wait_for_tasks();
    timer.MarkTime("Parallel (4 threads)");
    
    timer.Print();
    
    // Calculate speedup
    double seqTime = timer.GetIntervalTime(0);
    double par2Time = timer.GetIntervalTime(1);
    double par4Time = timer.GetIntervalTime(2);
    
    std::cout << "Speedup (2 threads): " << seqTime / par2Time << "x\n";
    std::cout << "Speedup (4 threads): " << seqTime / par4Time << "x\n";
}
```

### Optimal Thread Count Detection

```cpp
size_t findOptimalThreadCount(std::function<void()> task, 
                              int iterations = 100) {
    std::vector<std::pair<size_t, double>> results;
    
    for (size_t threads = 1; threads <= std::thread::hardware_concurrency(); 
         threads *= 2) 
    {
        Timer timer;
        timer.Start();
        
        ThreadPool pool(threads);
        for (int i = 0; i < iterations; ++i) {
            pool.enqueue(task);
        }
        pool.wait_for_tasks();
        
        timer.MarkTime();
        double time = timer.GetIntervalTime(0);
        results.emplace_back(threads, time);
        
        std::cout << threads << " threads: " << time << "s\n";
    }
    
    // Find fastest
    auto fastest = std::min_element(results.begin(), results.end(),
        [](auto& a, auto& b) { return a.second < b.second; });
    
    return fastest->first;
}
```

---

## See Also

- **[Visualizer](Visualizer.md)**: Visualizing performance results
- **[Serializer](Serializer.md)**: Saving profiling data
- **[ConsolePrinter](ConsolePrinter.md)**: Formatted performance reports

---

## Summary

**Timer** provides:
- ✅ Nanosecond-precision timing
- ✅ Multiple checkpoints with labels
- ✅ Interval and cumulative measurements
- ✅ Formatted output

**ThreadPool** provides:
- ✅ Efficient parallel task execution
- ✅ Thread reuse (no creation overhead)
- ✅ FIFO task queue
- ✅ Thread-safe operations
- ✅ RAII cleanup

**Use cases:**
- Algorithm performance comparison
- Parallel numerical computations
- Monte Carlo simulations
- Large-scale data processing
- Performance profiling
