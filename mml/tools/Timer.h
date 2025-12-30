///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Timer.h                                                             ///
///  Description: High-resolution timer for performance measurement                   ///
///               RAII-based timing with millisecond and microsecond precision        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_TIMER_H
#define MML_TIMER_H

#include "MMLBase.h"

#include <chrono>
#include <string>
#include <vector>

namespace MML
{
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

		void Start()
		{
			std::chrono::steady_clock::time_point t(std::chrono::steady_clock::now());
			_startTime = t;
		}
		void MarkTime()
		{
			_markedNames.push_back("Unnamed Mark");
			_markedTimes.push_back(std::chrono::steady_clock::now());
		}
		void MarkTime(const std::string& name)
		{
			_markedNames.push_back(name);
			_markedTimes.push_back(std::chrono::steady_clock::now());
		}
		double GetIntervalTime(int index = 0) const
		{
			// handling exceptions for index out of range and empty marked times
			if (index < 0 || index >= _markedTimes.size())
				throw std::out_of_range("Index out of range for marked times.");
			if( _markedTimes.empty() && index >= 0)
				throw std::runtime_error("No marked times available. Use MarkTime() before accessing intervals.");

			if (index == 0)
				return std::chrono::duration<double>(_markedTimes[0] - _startTime).count();
			else
				return std::chrono::duration<double>(_markedTimes[index] - _markedTimes[index - 1]).count();
		}
		double GetMarkTimeFromStart(int index = 0) const
		{
			if (index < 0 || index >= _markedTimes.size())
				throw std::out_of_range("Index out of range for marked times.");
			if( _markedTimes.empty() && index >= 0)
				throw std::runtime_error("No marked times available. Use MarkTime() before accessing intervals.");
			
			return std::chrono::duration<double>(_markedTimes[index] - _startTime).count();
		}
		double GetTotalTime() const
		{
			if (_markedTimes.empty())
				throw std::runtime_error("No marked times available. Use MarkTime() before accessing total time.");
			return std::chrono::duration<double>(_markedTimes.back() - _startTime).count();
		}

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
}
#endif 