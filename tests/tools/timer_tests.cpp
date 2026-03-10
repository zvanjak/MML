///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        timer_tests.cpp                                                     ///
///  Description: Unit tests for Timer                                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include <thread>
#include <chrono>

#include "MMLBase.h"
#include "tools/Timer.h"

using namespace MML;

#define TEST_PRECISION_INFO() \
	INFO("Test precision: " << Constants::Eps << " (" << (sizeof(Real) == 8 ? "double" : "float") << ")")

namespace MML::Tests::Tools::TimerTests {

	/////////////////////////////////////////////////////////////////////////////////////
	///                         BASIC TIMER OPERATIONS                                ///
	/////////////////////////////////////////////////////////////////////////////////////

	// NOTE: [!mayfail] - Timing tests can be flaky under heavy system load
	TEST_CASE("Timer - Basic start and total time", "[timer][basic][!mayfail]") {
		TEST_PRECISION_INFO();

		SECTION("Timer measures elapsed time correctly") {
			Timer timer;
			timer.Start();
			
			// Sleep for a longer period to reduce flakiness
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
			timer.MarkTime("End");  // Must mark time before GetTotalTime()
			
			Real elapsed = timer.GetTotalTime();
			
			// Should be at least 50ms (generous margin for system overhead under load)
			REQUIRE(elapsed >= 0.05);
			// Should be less than 2s (very generous upper bound for heavily loaded systems)
			REQUIRE(elapsed < 2.0);
		}

		SECTION("Timer can be restarted and clears marks") {
			Timer timer;
			timer.Start();
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
			timer.MarkTime("First End");
			
			Real firstMeasure = timer.GetTotalTime();
			REQUIRE(timer.GetMarkCount() == 1);
			
			timer.Start();  // Restart clears marks
			REQUIRE(timer.GetMarkCount() == 0);  // Marks cleared
			
			std::this_thread::sleep_for(std::chrono::milliseconds(50));
			timer.MarkTime("Second End");
			
			REQUIRE(timer.GetMarkCount() == 1);
			Real secondMeasure = timer.GetTotalTime();
			
			// Second measurement (50ms) should be less than first (100ms)
			REQUIRE(secondMeasure < firstMeasure);
			// Second measurement should be at least 25ms (generous lower bound)
			REQUIRE(secondMeasure >= 0.025);
			REQUIRE(secondMeasure < 1.0);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         MARK TIME FUNCTIONALITY                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	// NOTE: [!mayfail] - Timing tests can be flaky under heavy system load
	TEST_CASE("Timer - MarkTime and interval timing", "[timer][marks][!mayfail]") {
		TEST_PRECISION_INFO();

		SECTION("MarkTime records named intervals") {
			Timer timer;
			timer.Start();
			
			std::this_thread::sleep_for(std::chrono::milliseconds(75));
			timer.MarkTime("First");
			
			std::this_thread::sleep_for(std::chrono::milliseconds(75));
			timer.MarkTime("Second");
			
			// GetIntervalTime returns interval between marks
			Real firstInterval = timer.GetIntervalTime(0);
			Real secondInterval = timer.GetIntervalTime(1);
			
			// Both intervals should be at least 40ms (generous lower bound)
			REQUIRE(firstInterval >= 0.04);
			REQUIRE(firstInterval < 1.0);
			REQUIRE(secondInterval >= 0.04);
			REQUIRE(secondInterval < 1.0);
		}

		SECTION("GetMarkTimeFromStart returns cumulative time") {
			Timer timer;
			timer.Start();
			
			std::this_thread::sleep_for(std::chrono::milliseconds(75));
			timer.MarkTime("First");
			
			std::this_thread::sleep_for(std::chrono::milliseconds(75));
			timer.MarkTime("Second");
			
			Real firstFromStart = timer.GetMarkTimeFromStart(0);
			Real secondFromStart = timer.GetMarkTimeFromStart(1);
			
			// Second mark should be further from start than first
			REQUIRE(secondFromStart > firstFromStart);
			// First mark ~75ms, second mark ~150ms (generous lower bounds)
			REQUIRE(firstFromStart >= 0.04);
			REQUIRE(secondFromStart >= 0.08);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         START/RESET SEMANTICS                                 ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Timer - Start/Reset semantics", "[timer][semantics]") {
		TEST_PRECISION_INFO();

		SECTION("Start() clears previous marks - contract test") {
			// This test explicitly validates the fixed behavior:
			// Start() must clear marks to avoid confusing totals after restart
			Timer timer;
			timer.Start();
			timer.MarkTime("OldMark1");
			timer.MarkTime("OldMark2");
			REQUIRE(timer.GetMarkCount() == 2);
			
			// Restart should clear all marks
			timer.Start();
			REQUIRE(timer.GetMarkCount() == 0);
			
			// New marks work correctly
			timer.MarkTime("NewMark");
			REQUIRE(timer.GetMarkCount() == 1);
			Real elapsed = timer.GetTotalTime();
			REQUIRE(elapsed >= 0.0);
			REQUIRE(elapsed < 0.01);  // Should be very small (no sleep)
		}

		SECTION("Reset() is equivalent to Start()") {
			Timer timer;
			timer.Start();
			timer.MarkTime("Mark1");
			timer.MarkTime("Mark2");
			REQUIRE(timer.GetMarkCount() == 2);
			
			timer.Reset();
			REQUIRE(timer.GetMarkCount() == 0);
		}

		SECTION("GetTotalTime after restart uses new start time") {
			Timer timer;
			timer.Start();
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
			timer.MarkTime("First");
			Real firstTotal = timer.GetTotalTime();
			REQUIRE(firstTotal >= 0.05);
			
			// Restart and add a quick mark
			timer.Start();
			timer.MarkTime("AfterRestart");
			Real secondTotal = timer.GetTotalTime();
			
			// Second total should be near zero (no sleep after restart)
			REQUIRE(secondTotal >= 0.0);
			REQUIRE(secondTotal < 0.1);  // Generous upper bound
			REQUIRE(secondTotal < firstTotal);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         TIMING MONOTONICITY                                   ///
	/////////////////////////////////////////////////////////////////////////////////////

	// NOTE: [!mayfail] - Timing tests can be flaky under heavy system load
	TEST_CASE("Timer - Monotonicity", "[timer][monotonicity][!mayfail]") {
		TEST_PRECISION_INFO();

		SECTION("Mark times from start increase monotonically") {
			Timer timer;
			timer.Start();
			
			// Add marks and collect times (longer sleep for reliability)
			for (int i = 0; i < 5; ++i) {
				std::this_thread::sleep_for(std::chrono::milliseconds(30));
				timer.MarkTime("Sample" + std::to_string(i));
			}
			
			// Collect samples from marks
			std::vector<Real> samples;
			for (int i = 0; i < 5; ++i) {
				samples.push_back(timer.GetMarkTimeFromStart(i));
			}
			
			// Each sample should be greater than the previous
			for (size_t i = 1; i < samples.size(); ++i) {
				REQUIRE(samples[i] >= samples[i-1]);
			}
		}

		SECTION("Mark times from start are monotonically increasing") {
			Timer timer;
			timer.Start();
			
			for (int i = 0; i < 5; ++i) {
				std::this_thread::sleep_for(std::chrono::milliseconds(30));
				timer.MarkTime("Mark" + std::to_string(i));
			}
			
			for (int i = 1; i < 5; ++i) {
				REQUIRE(timer.GetMarkTimeFromStart(i) >= timer.GetMarkTimeFromStart(i-1));
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         EDGE CASES                                            ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Timer - Edge cases", "[timer][edge]") {
		TEST_PRECISION_INFO();

		SECTION("Timer works with zero sleep") {
			Timer timer;
			timer.Start();
			timer.MarkTime("Immediate");  // Mark immediately
			Real elapsed = timer.GetTotalTime();
			
			// Should be a very small non-negative number
			REQUIRE(elapsed >= 0.0);
			REQUIRE(elapsed < 0.01);  // Less than 10ms for immediate call
		}

		SECTION("Multiple marks can be recorded") {
			Timer timer;
			timer.Start();
			
			constexpr int NUM_MARKS = 10;
			for (int i = 0; i < NUM_MARKS; ++i) {
				timer.MarkTime("Mark" + std::to_string(i));
			}
			
			// All marks should be accessible
			for (int i = 0; i < NUM_MARKS; ++i) {
				Real time = timer.GetMarkTimeFromStart(i);
				REQUIRE(time >= 0.0);
			}
		}
	}

} // namespace MML::Tests::Tools::TimerTests

