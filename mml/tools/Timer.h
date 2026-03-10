///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Timer.h                                                             ///
///  Description: High-resolution timer for performance measurement                   ///
///               RAII-based timing with millisecond and microsecond precision        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_TIMER_H
#define MML_TIMER_H

#include "MMLBase.h"

#include <chrono>
#include <string>
#include <vector>
#include <functional>

namespace MML
{
	/// @brief Multi-mark timer for detailed performance profiling
	/// @details Timer allows recording multiple named time marks during execution,
	/// then analyzing intervals between marks. Use for complex profiling scenarios.
	/// 
	/// **Usage**:
	/// @code
	/// Timer t;
	/// t.Start();
	/// // ... phase 1 ...
	/// t.MarkTime("Phase 1");
	/// // ... phase 2 ...
	/// t.MarkTime("Phase 2");
	/// t.Print();  // Shows intervals for each phase
	/// @endcode
	/// 
	/// @see ScopedTimer for simple RAII-based single-operation timing
	class Timer
	{
	private:
		std::chrono::steady_clock::time_point _startTime;
		std::vector<std::chrono::steady_clock::time_point> _markedTimes;
		std::vector<std::string> _markedNames;

	public:
		Timer() {}
		~Timer() {}
		
		Timer(const Timer&) = delete;
		Timer& operator=(const Timer&) = delete;

		/// @brief Start (or restart) the timer, clearing all previous marks.
		/// Use this at the beginning of a new timing session.
		void Start()
		{
			_markedTimes.clear();
			_markedNames.clear();
			_startTime = std::chrono::steady_clock::now();
		}

		/// @brief Reset the timer completely (same as Start).
		/// Provided for semantic clarity.
		void Reset() { Start(); }

		/// @brief Get the number of recorded marks.
		size_t GetMarkCount() const noexcept { return _markedTimes.size(); }
		
		/// @brief Record a time mark with default name
		void MarkTime()
		{
			_markedNames.push_back("Unnamed Mark");
			_markedTimes.push_back(std::chrono::steady_clock::now());
		}
		
		/// @brief Record a time mark with the given name
		/// @param name Descriptive name for this timing point
		void MarkTime(const std::string& name)
		{
			_markedNames.push_back(name);
			_markedTimes.push_back(std::chrono::steady_clock::now());
		}
		
		/// @brief Get the time interval for a specific mark
		/// @param index Mark index (0-based). Index 0 returns time from Start() to first mark.
		///              Other indices return time from previous mark to this mark.
		/// @return Interval time in seconds
		/// @throws std::out_of_range if index is out of range
		/// @throws std::runtime_error if no marks have been recorded
		double GetIntervalTime(size_t index = 0) const
		{
			if (_markedTimes.empty())
				throw std::runtime_error("No marked times available. Use MarkTime() before accessing intervals.");
			if (index >= _markedTimes.size())
				throw std::out_of_range("Index out of range for marked times.");

			if (index == 0)
				return std::chrono::duration<double>(_markedTimes[0] - _startTime).count();
			else
				return std::chrono::duration<double>(_markedTimes[index] - _markedTimes[index - 1]).count();
		}
		
		/// @brief Get the elapsed time from Start() to a specific mark
		/// @param index Mark index (0-based)
		/// @return Time in seconds from Start() to the specified mark
		/// @throws std::out_of_range if index is out of range
		/// @throws std::runtime_error if no marks have been recorded
		double GetMarkTimeFromStart(size_t index = 0) const
		{
			if (_markedTimes.empty())
				throw std::runtime_error("No marked times available. Use MarkTime() before accessing intervals.");
			if (index >= _markedTimes.size())
				throw std::out_of_range("Index out of range for marked times.");
			
			return std::chrono::duration<double>(_markedTimes[index] - _startTime).count();
		}
		
		/// @brief Get the total elapsed time from Start() to the last mark
		/// @return Total time in seconds
		/// @throws std::runtime_error if no marks have been recorded
		double GetTotalTime() const
		{
			if (_markedTimes.empty())
				throw std::runtime_error("No marked times available. Use MarkTime() before accessing total time.");
			return std::chrono::duration<double>(_markedTimes.back() - _startTime).count();
		}

		/// @brief Print all timing results to stdout
		void Print()
		{
			if (_markedTimes.empty())
			{
				std::cout << "No marked times available." << std::endl;
				return;
			}
			std::cout << "Timer Results:" << std::endl;
			for (size_t i = 0; i < _markedTimes.size(); ++i)
			{
				std::cout << "Mark " << i << " (" << std::setw(25) << _markedNames[i] << "): "
									<< GetIntervalTime(i) << " seconds" << std::endl;
			}
			std::cout << "Total Time: " << GetTotalTime() << " seconds" << std::endl;
		}
	};

	/// @brief RAII-based timer that automatically measures scope duration
	/// @details ScopedTimer starts timing on construction and reports elapsed time on destruction.
	/// Perfect for timing individual operations without manual Start()/MarkTime() calls.
	/// 
	/// **Basic Usage** (prints to stdout):
	/// @code
	/// {
	///     ScopedTimer timer("Matrix multiplication");
	///     // ... expensive operation ...
	/// }  // Prints: "Matrix multiplication: 1.234 seconds"
	/// @endcode
	/// 
	/// **Silent mode** (stores result only):
	/// @code
	/// double elapsed;
	/// {
	///     ScopedTimer timer("Operation", elapsed);
	///     // ... operation ...
	/// }
	/// std::cout << "Took " << elapsed << " seconds\n";
	/// @endcode
	/// 
	/// **With callback**:
	/// @code
	/// {
	///     ScopedTimer timer("Query", [](const std::string& name, double secs) {
	///         logger.info("{} completed in {} ms", name, secs * 1000);
	///     });
	///     // ... operation ...
	/// }
	/// @endcode
	/// 
	/// @note Non-copyable, non-movable. Designed for stack allocation only.
	/// @see Timer for multi-mark timing scenarios
	class ScopedTimer
	{
	public:
		/// Callback signature: (operation_name, elapsed_seconds)
		using Callback = std::function<void(const std::string&, double)>;

	private:
		std::string _name;
		std::chrono::steady_clock::time_point _startTime;
		double* _resultPtr;
		Callback _callback;
		bool _silent;

	public:
		/// @brief Construct a ScopedTimer that prints to stdout on destruction
		/// @param name Descriptive name for the timed operation
		explicit ScopedTimer(const std::string& name)
			: _name(name)
			, _startTime(std::chrono::steady_clock::now())
			, _resultPtr(nullptr)
			, _callback(nullptr)
			, _silent(false)
		{}

		/// @brief Construct a ScopedTimer that stores the result in a variable
		/// @param name Descriptive name for the timed operation
		/// @param result Reference to a double that will receive the elapsed time
		/// @param silent If true, does not print to stdout (default: true)
		ScopedTimer(const std::string& name, double& result, bool silent = true)
			: _name(name)
			, _startTime(std::chrono::steady_clock::now())
			, _resultPtr(&result)
			, _callback(nullptr)
			, _silent(silent)
		{}

		/// @brief Construct a ScopedTimer with a custom callback
		/// @param name Descriptive name for the timed operation
		/// @param callback Function called with (name, elapsed_seconds) on destruction
		/// @param silent If true, does not print to stdout (default: true)
		ScopedTimer(const std::string& name, Callback callback, bool silent = true)
			: _name(name)
			, _startTime(std::chrono::steady_clock::now())
			, _resultPtr(nullptr)
			, _callback(std::move(callback))
			, _silent(silent)
		{}

		/// @brief Destructor - reports elapsed time
		~ScopedTimer()
		{
			double elapsed = std::chrono::duration<double>(
				std::chrono::steady_clock::now() - _startTime).count();
			
			// Store result if requested
			if (_resultPtr)
				*_resultPtr = elapsed;
			
			// Invoke callback if provided
			if (_callback)
			{
				try {
					_callback(_name, elapsed);
				}
				catch (...) {
					// Ignore exceptions in destructor
				}
			}
			
			// Print unless silent
			if (!_silent)
				std::cout << _name << ": " << elapsed << " seconds" << std::endl;
		}

		// Non-copyable and non-movable
		ScopedTimer(const ScopedTimer&) = delete;
		ScopedTimer& operator=(const ScopedTimer&) = delete;
		ScopedTimer(ScopedTimer&&) = delete;
		ScopedTimer& operator=(ScopedTimer&&) = delete;

		/// @brief Get the elapsed time so far without stopping the timer
		/// @return Elapsed time in seconds since construction
		double elapsed() const noexcept
		{
			return std::chrono::duration<double>(
				std::chrono::steady_clock::now() - _startTime).count();
		}

		/// @brief Get the name of this timer
		/// @return The operation name
		const std::string& name() const noexcept { return _name; }
	};
}
#endif 