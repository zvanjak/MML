///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        PrecisionTestFramework.h                                            ///
///  Description: Unified framework for precision testing of numerical algorithms     ///
///               Provides standardized result collection, statistics, and output     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_PRECISION_TEST_FRAMEWORK_H
#define MML_PRECISION_TEST_FRAMEWORK_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <functional>

namespace MML::PrecisionTesting
{
    /********************************************************************************************************************/
    /********                                    PRECISION TEST RESULT                                           ********/
    /********************************************************************************************************************/
    
    /**
     * @brief Single precision test result capturing all relevant metrics
     */
    struct PrecisionTestResult
    {
        std::string algorithm_name;     ///< Name of the algorithm being tested
        std::string test_function;      ///< Name of the test function
        std::string parameters;         ///< Algorithm parameters (e.g., "h=1e-6", "order=4")
        
        double exact_value;             ///< Known exact/analytical value
        double computed_value;          ///< Value computed by the algorithm
        double absolute_error;          ///< |computed - exact|
        double relative_error;          ///< |computed - exact| / |exact|
        
        int iterations;                 ///< Number of iterations (if applicable)
        double time_ms;                 ///< Computation time in milliseconds
        
        bool converged;                 ///< Whether the algorithm converged
        std::string notes;              ///< Additional notes or warnings
        
        PrecisionTestResult() 
            : exact_value(0), computed_value(0), absolute_error(0), relative_error(0),
              iterations(0), time_ms(0), converged(true) {}
        
        PrecisionTestResult(const std::string& algo, const std::string& func,
                           double exact, double computed)
            : algorithm_name(algo), test_function(func),
              exact_value(exact), computed_value(computed),
              iterations(0), time_ms(0), converged(true)
        {
            absolute_error = std::abs(computed - exact);
            relative_error = (exact != 0.0) ? absolute_error / std::abs(exact) : absolute_error;
        }
        
        /// @brief Get error order of magnitude (e.g., -8 means error ~ 1e-8)
        int getErrorOrder() const {
            if (absolute_error == 0.0) return -16;  // Machine precision
            return static_cast<int>(std::floor(std::log10(absolute_error)));
        }
    };

    /********************************************************************************************************************/
    /********                                    STATISTICS CALCULATOR                                           ********/
    /********************************************************************************************************************/
    
    /**
     * @brief Statistics for a collection of test results
     */
    struct PrecisionStatistics
    {
        double min_abs_error;
        double max_abs_error;
        double avg_abs_error;
        double stddev_abs_error;
        
        double min_rel_error;
        double max_rel_error;
        double avg_rel_error;
        double stddev_rel_error;
        
        int total_tests;
        int converged_tests;
        double avg_time_ms;
        
        int best_error_order;   ///< Best (most negative) error order
        int worst_error_order;  ///< Worst (least negative) error order
        
        PrecisionStatistics() 
            : min_abs_error(0), max_abs_error(0), avg_abs_error(0), stddev_abs_error(0),
              min_rel_error(0), max_rel_error(0), avg_rel_error(0), stddev_rel_error(0),
              total_tests(0), converged_tests(0), avg_time_ms(0),
              best_error_order(0), worst_error_order(0) {}
    };

    /********************************************************************************************************************/
    /********                                    PRECISION TEST SUITE                                            ********/
    /********************************************************************************************************************/
    
    /**
     * @brief Collection of precision test results with analysis and output capabilities
     */
    class PrecisionTestSuite
    {
    private:
        std::string _suite_name;
        std::string _description;
        std::vector<PrecisionTestResult> _results;
        
    public:
        PrecisionTestSuite(const std::string& name, const std::string& desc = "")
            : _suite_name(name), _description(desc) {}
        
        /// @brief Add a single test result
        void addResult(const PrecisionTestResult& result) {
            _results.push_back(result);
        }
        
        /// @brief Add result with timing
        void addResult(const std::string& algo, const std::string& func,
                      double exact, double computed, double time_ms = 0.0) {
            PrecisionTestResult r(algo, func, exact, computed);
            r.time_ms = time_ms;
            _results.push_back(r);
        }
        
        /// @brief Get all results
        const std::vector<PrecisionTestResult>& getResults() const { return _results; }
        
        /// @brief Get number of results
        size_t size() const { return _results.size(); }
        
        /// @brief Clear all results
        void clear() { _results.clear(); }
        
        //----------------------------------------------------------------------------------------------------------
        // STATISTICS
        //----------------------------------------------------------------------------------------------------------
        
        /// @brief Calculate statistics for all results
        PrecisionStatistics calculateStatistics() const {
            PrecisionStatistics stats;
            if (_results.empty()) return stats;
            
            stats.total_tests = static_cast<int>(_results.size());
            
            // Initialize with first result
            stats.min_abs_error = stats.max_abs_error = _results[0].absolute_error;
            stats.min_rel_error = stats.max_rel_error = _results[0].relative_error;
            stats.best_error_order = stats.worst_error_order = _results[0].getErrorOrder();
            
            double sum_abs = 0, sum_rel = 0, sum_time = 0;
            
            for (const auto& r : _results) {
                stats.min_abs_error = std::min(stats.min_abs_error, r.absolute_error);
                stats.max_abs_error = std::max(stats.max_abs_error, r.absolute_error);
                stats.min_rel_error = std::min(stats.min_rel_error, r.relative_error);
                stats.max_rel_error = std::max(stats.max_rel_error, r.relative_error);
                
                stats.best_error_order = std::min(stats.best_error_order, r.getErrorOrder());
                stats.worst_error_order = std::max(stats.worst_error_order, r.getErrorOrder());
                
                sum_abs += r.absolute_error;
                sum_rel += r.relative_error;
                sum_time += r.time_ms;
                
                if (r.converged) stats.converged_tests++;
            }
            
            stats.avg_abs_error = sum_abs / stats.total_tests;
            stats.avg_rel_error = sum_rel / stats.total_tests;
            stats.avg_time_ms = sum_time / stats.total_tests;
            
            // Calculate standard deviations
            double sum_sq_abs = 0, sum_sq_rel = 0;
            for (const auto& r : _results) {
                sum_sq_abs += (r.absolute_error - stats.avg_abs_error) * (r.absolute_error - stats.avg_abs_error);
                sum_sq_rel += (r.relative_error - stats.avg_rel_error) * (r.relative_error - stats.avg_rel_error);
            }
            stats.stddev_abs_error = std::sqrt(sum_sq_abs / stats.total_tests);
            stats.stddev_rel_error = std::sqrt(sum_sq_rel / stats.total_tests);
            
            return stats;
        }
        
        /// @brief Calculate statistics filtered by algorithm name
        PrecisionStatistics calculateStatistics(const std::string& algorithm) const {
            PrecisionTestSuite filtered(_suite_name + " [" + algorithm + "]");
            for (const auto& r : _results) {
                if (r.algorithm_name == algorithm) {
                    filtered.addResult(r);
                }
            }
            return filtered.calculateStatistics();
        }
        
        //----------------------------------------------------------------------------------------------------------
        // CONSOLE OUTPUT
        //----------------------------------------------------------------------------------------------------------
        
        /// @brief Print header with suite name
        void printHeader(std::ostream& os = std::cout) const {
            os << "\n" << std::string(100, '=') << "\n";
            os << "  PRECISION TEST SUITE: " << _suite_name << "\n";
            if (!_description.empty()) {
                os << "  " << _description << "\n";
            }
            os << std::string(100, '=') << "\n\n";
        }
        
        /// @brief Print detailed results table
        void printDetailedTable(std::ostream& os = std::cout) const {
            if (_results.empty()) {
                os << "  No results to display.\n";
                return;
            }
            
            // Header
            os << std::left << std::setw(20) << "Algorithm"
               << std::setw(20) << "Function"
               << std::right << std::setw(14) << "Exact"
               << std::setw(14) << "Computed"
               << std::setw(12) << "Abs Error"
               << std::setw(12) << "Rel Error"
               << std::setw(8) << "Order"
               << "\n";
            os << std::string(100, '-') << "\n";
            
            // Results
            for (const auto& r : _results) {
                os << std::left << std::setw(20) << r.algorithm_name.substr(0, 19)
                   << std::setw(20) << r.test_function.substr(0, 19)
                   << std::right << std::scientific << std::setprecision(6)
                   << std::setw(14) << r.exact_value
                   << std::setw(14) << r.computed_value
                   << std::setw(12) << r.absolute_error
                   << std::setw(12) << r.relative_error
                   << std::fixed << std::setw(8) << r.getErrorOrder()
                   << "\n";
            }
            os << std::string(100, '-') << "\n";
        }
        
        /// @brief Print compact error order table (algorithms vs functions)
        void printErrorOrderMatrix(std::ostream& os = std::cout) const {
            if (_results.empty()) return;
            
            // Collect unique algorithms and functions
            std::vector<std::string> algorithms, functions;
            for (const auto& r : _results) {
                if (std::find(algorithms.begin(), algorithms.end(), r.algorithm_name) == algorithms.end())
                    algorithms.push_back(r.algorithm_name);
                if (std::find(functions.begin(), functions.end(), r.test_function) == functions.end())
                    functions.push_back(r.test_function);
            }
            
            os << "\nERROR ORDER MATRIX (log10 of absolute error):\n";
            os << std::string(80, '-') << "\n";
            
            // Header row with function names
            os << std::left << std::setw(15) << "Algorithm";
            for (const auto& func : functions) {
                os << std::setw(10) << func.substr(0, 9);
            }
            os << std::setw(10) << "AVG" << "\n";
            os << std::string(80, '-') << "\n";
            
            // Data rows
            for (const auto& algo : algorithms) {
                os << std::left << std::setw(15) << algo.substr(0, 14);
                double sum_order = 0;
                int count = 0;
                
                for (const auto& func : functions) {
                    // Find matching result
                    bool found = false;
                    for (const auto& r : _results) {
                        if (r.algorithm_name == algo && r.test_function == func) {
                            int order = r.getErrorOrder();
                            os << std::right << std::setw(10) << order;
                            sum_order += order;
                            count++;
                            found = true;
                            break;
                        }
                    }
                    if (!found) os << std::setw(10) << "-";
                }
                
                // Average order for this algorithm
                if (count > 0) {
                    os << std::fixed << std::setprecision(1) << std::setw(10) << (sum_order / count);
                }
                os << "\n";
            }
            os << std::string(80, '-') << "\n";
        }
        
        /// @brief Print summary statistics
        void printSummary(std::ostream& os = std::cout) const {
            PrecisionStatistics stats = calculateStatistics();
            
            os << "\nSUMMARY STATISTICS:\n";
            os << std::string(50, '-') << "\n";
            os << "  Total tests:      " << stats.total_tests << "\n";
            os << "  Converged:        " << stats.converged_tests << " (" 
               << std::fixed << std::setprecision(1) 
               << (100.0 * stats.converged_tests / stats.total_tests) << "%)\n";
            os << std::scientific << std::setprecision(2);
            os << "  Absolute error:   min=" << stats.min_abs_error 
               << "  max=" << stats.max_abs_error 
               << "  avg=" << stats.avg_abs_error << "\n";
            os << "  Relative error:   min=" << stats.min_rel_error 
               << "  max=" << stats.max_rel_error 
               << "  avg=" << stats.avg_rel_error << "\n";
            os << "  Error order:      best=" << stats.best_error_order 
               << "  worst=" << stats.worst_error_order << "\n";
            if (stats.avg_time_ms > 0) {
                os << std::fixed << std::setprecision(3);
                os << "  Avg time:         " << stats.avg_time_ms << " ms\n";
            }
            os << std::string(50, '-') << "\n";
        }
        
        //----------------------------------------------------------------------------------------------------------
        // MARKDOWN OUTPUT (for documentation)
        //----------------------------------------------------------------------------------------------------------
        
        /// @brief Generate markdown table
        std::string toMarkdown() const {
            std::ostringstream md;
            
            md << "## " << _suite_name << "\n\n";
            if (!_description.empty()) {
                md << _description << "\n\n";
            }
            
            // Results table
            md << "| Algorithm | Function | Exact | Computed | Abs Error | Rel Error | Order |\n";
            md << "|-----------|----------|-------|----------|-----------|-----------|-------|\n";
            
            for (const auto& r : _results) {
                md << "| " << r.algorithm_name 
                   << " | " << r.test_function
                   << " | " << std::scientific << std::setprecision(4) << r.exact_value
                   << " | " << r.computed_value
                   << " | " << r.absolute_error
                   << " | " << r.relative_error
                   << " | " << r.getErrorOrder()
                   << " |\n";
            }
            
            // Statistics
            PrecisionStatistics stats = calculateStatistics();
            md << "\n### Statistics\n\n";
            md << "- **Total tests:** " << stats.total_tests << "\n";
            md << "- **Best error order:** 10^" << stats.best_error_order << "\n";
            md << "- **Worst error order:** 10^" << stats.worst_error_order << "\n";
            md << "- **Average absolute error:** " << std::scientific << stats.avg_abs_error << "\n";
            
            return md.str();
        }
        
        /// @brief Generate error order matrix in markdown
        std::string toMarkdownMatrix() const {
            if (_results.empty()) return "";
            
            std::ostringstream md;
            
            // Collect unique algorithms and functions
            std::vector<std::string> algorithms, functions;
            for (const auto& r : _results) {
                if (std::find(algorithms.begin(), algorithms.end(), r.algorithm_name) == algorithms.end())
                    algorithms.push_back(r.algorithm_name);
                if (std::find(functions.begin(), functions.end(), r.test_function) == functions.end())
                    functions.push_back(r.test_function);
            }
            
            md << "### Error Order Matrix\n\n";
            md << "| Algorithm |";
            for (const auto& func : functions) {
                md << " " << func << " |";
            }
            md << "\n|-----------|";
            for (size_t i = 0; i < functions.size(); i++) {
                md << "------|";
            }
            md << "\n";
            
            for (const auto& algo : algorithms) {
                md << "| " << algo << " |";
                for (const auto& func : functions) {
                    bool found = false;
                    for (const auto& r : _results) {
                        if (r.algorithm_name == algo && r.test_function == func) {
                            md << " " << r.getErrorOrder() << " |";
                            found = true;
                            break;
                        }
                    }
                    if (!found) md << " - |";
                }
                md << "\n";
            }
            
            return md.str();
        }
        
        //----------------------------------------------------------------------------------------------------------
        // FILE EXPORT
        //----------------------------------------------------------------------------------------------------------
        
        /// @brief Export results to CSV file
        bool exportCSV(const std::string& filename) const {
            std::ofstream ofs(filename);
            if (!ofs) return false;
            
            // Header
            ofs << "algorithm,function,parameters,exact,computed,abs_error,rel_error,error_order,iterations,time_ms,converged,notes\n";
            
            // Data
            for (const auto& r : _results) {
                ofs << "\"" << r.algorithm_name << "\","
                    << "\"" << r.test_function << "\","
                    << "\"" << r.parameters << "\","
                    << std::scientific << std::setprecision(15)
                    << r.exact_value << ","
                    << r.computed_value << ","
                    << r.absolute_error << ","
                    << r.relative_error << ","
                    << r.getErrorOrder() << ","
                    << r.iterations << ","
                    << std::fixed << std::setprecision(3) << r.time_ms << ","
                    << (r.converged ? "true" : "false") << ","
                    << "\"" << r.notes << "\"\n";
            }
            
            return true;
        }
        
        /// @brief Export markdown to file
        bool exportMarkdown(const std::string& filename) const {
            std::ofstream ofs(filename);
            if (!ofs) return false;
            
            ofs << toMarkdown();
            ofs << "\n" << toMarkdownMatrix();
            
            return true;
        }
    };

    /********************************************************************************************************************/
    /********                                    TIMING UTILITIES                                                ********/
    /********************************************************************************************************************/
    
    /**
     * @brief Simple timer for measuring algorithm execution time
     */
    class Timer
    {
    private:
        std::chrono::high_resolution_clock::time_point _start;
        
    public:
        Timer() : _start(std::chrono::high_resolution_clock::now()) {}
        
        void reset() {
            _start = std::chrono::high_resolution_clock::now();
        }
        
        double elapsedMs() const {
            auto end = std::chrono::high_resolution_clock::now();
            return std::chrono::duration<double, std::milli>(end - _start).count();
        }
        
        double elapsedUs() const {
            auto end = std::chrono::high_resolution_clock::now();
            return std::chrono::duration<double, std::micro>(end - _start).count();
        }
    };

    /********************************************************************************************************************/
    /********                                    TEST RUNNER HELPERS                                             ********/
    /********************************************************************************************************************/
    
    /**
     * @brief Run a test function and capture result with timing
     * 
     * @tparam Func Callable that returns computed value
     * @param suite Target test suite
     * @param algo Algorithm name
     * @param func Function name
     * @param exact Known exact value
     * @param compute Callable to compute the value
     */
    template<typename Func>
    void runTest(PrecisionTestSuite& suite, 
                 const std::string& algo, 
                 const std::string& func,
                 double exact,
                 Func&& compute)
    {
        Timer timer;
        double computed = compute();
        double time_ms = timer.elapsedMs();
        
        PrecisionTestResult result(algo, func, exact, computed);
        result.time_ms = time_ms;
        suite.addResult(result);
    }
    
    /**
     * @brief Run multiple tests with varying parameters
     * 
     * @tparam ParamType Type of the varying parameter
     * @tparam Func Callable that takes parameter and returns computed value
     */
    template<typename ParamType, typename Func>
    void runParameterSweep(PrecisionTestSuite& suite,
                           const std::string& algo,
                           const std::string& func,
                           double exact,
                           const std::vector<ParamType>& params,
                           Func&& compute,
                           std::function<std::string(ParamType)> paramToString = nullptr)
    {
        for (const auto& p : params) {
            Timer timer;
            double computed = compute(p);
            double time_ms = timer.elapsedMs();
            
            PrecisionTestResult result(algo, func, exact, computed);
            result.time_ms = time_ms;
            if (paramToString) {
                result.parameters = paramToString(p);
            }
            suite.addResult(result);
        }
    }

    /********************************************************************************************************************/
    /********                                    PRECISION RATING                                                ********/
    /********************************************************************************************************************/
    
    /**
     * @brief Rate precision quality for user-friendly output
     */
    enum class PrecisionRating {
        Excellent,  ///< Error < 1e-12
        Good,       ///< Error < 1e-8
        Acceptable, ///< Error < 1e-4
        Poor,       ///< Error < 1e-2
        Failed      ///< Error >= 1e-2
    };
    
    inline PrecisionRating ratePrecision(double abs_error) {
        if (abs_error < 1e-12) return PrecisionRating::Excellent;
        if (abs_error < 1e-8)  return PrecisionRating::Good;
        if (abs_error < 1e-4)  return PrecisionRating::Acceptable;
        if (abs_error < 1e-2)  return PrecisionRating::Poor;
        return PrecisionRating::Failed;
    }
    
    inline const char* ratingToString(PrecisionRating rating) {
        switch (rating) {
            case PrecisionRating::Excellent:  return "EXCELLENT";
            case PrecisionRating::Good:       return "GOOD";
            case PrecisionRating::Acceptable: return "ACCEPTABLE";
            case PrecisionRating::Poor:       return "POOR";
            case PrecisionRating::Failed:     return "FAILED";
        }
        return "UNKNOWN";
    }
    
    inline const char* ratingToEmoji(PrecisionRating rating) {
        switch (rating) {
            case PrecisionRating::Excellent:  return "🟢";
            case PrecisionRating::Good:       return "🟢";
            case PrecisionRating::Acceptable: return "🟡";
            case PrecisionRating::Poor:       return "🟠";
            case PrecisionRating::Failed:     return "🔴";
        }
        return "⚪";
    }

} // namespace MML::PrecisionTesting

#endif // MML_PRECISION_TEST_FRAMEWORK_H
