/**
 * @file visualizer_tests.cpp
 * @brief Infrastructure tests for Visualizer.h
 * 
 * These tests verify the VisualizerResult struct functionality.
 * 
 * Note: The Visualizer class methods all require:
 * - File system operations (creating data files)
 * - External visualization applications
 * - GUI interaction
 * 
 * Therefore, only the VisualizerResult struct is tested here.
 * Full integration testing would require the visualization tools.
 */

#include <catch2/catch_all.hpp>
#include <string>

#include "tools/Visualizer.h"

using namespace MML;

namespace MML::Tests::Tools::VisualizerTests {

///////////////////////////////////////////////////////////////////////////////
//                        VISUALIZER RESULT - SUCCESS                        //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("VisualizerResult::Success - Default", "[Visualizer][Result]") {
    auto result = VisualizerResult::Success();
    
    REQUIRE(result.success == true);
    REQUIRE(result.exitCode == 0);
    REQUIRE(result.errorMessage.empty());
    REQUIRE(result.dataFilePath.empty());
}

TEST_CASE("VisualizerResult::Success - With data path", "[Visualizer][Result]") {
    auto result = VisualizerResult::Success("/path/to/data.mml");
    
    REQUIRE(result.success == true);
    REQUIRE(result.exitCode == 0);
    REQUIRE(result.dataFilePath == "/path/to/data.mml");
    REQUIRE(result.errorMessage.empty());
}

TEST_CASE("VisualizerResult::Success - With Windows path", "[Visualizer][Result]") {
    auto result = VisualizerResult::Success("C:\\Users\\data\\file.mml");
    
    REQUIRE(result.success == true);
    REQUIRE(result.dataFilePath == "C:\\Users\\data\\file.mml");
}

TEST_CASE("VisualizerResult::Success - With relative path", "[Visualizer][Result]") {
    auto result = VisualizerResult::Success("./results/output.dat");
    
    REQUIRE(result.success == true);
    REQUIRE(result.dataFilePath == "./results/output.dat");
}

///////////////////////////////////////////////////////////////////////////////
//                        VISUALIZER RESULT - FAILURE                        //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("VisualizerResult::Failure - Message only", "[Visualizer][Result]") {
    auto result = VisualizerResult::Failure("Something went wrong");
    
    REQUIRE(result.success == false);
    REQUIRE(result.exitCode == -1);  // Default error code
    REQUIRE(result.errorMessage == "Something went wrong");
    REQUIRE(result.dataFilePath.empty());
}

TEST_CASE("VisualizerResult::Failure - With data path", "[Visualizer][Result]") {
    auto result = VisualizerResult::Failure("Error occurred", "/path/to/partial.mml");
    
    REQUIRE(result.success == false);
    REQUIRE(result.exitCode == -1);
    REQUIRE(result.errorMessage == "Error occurred");
    REQUIRE(result.dataFilePath == "/path/to/partial.mml");
}

TEST_CASE("VisualizerResult::Failure - With exit code", "[Visualizer][Result]") {
    auto result = VisualizerResult::Failure("Process failed", "/path/to/file.mml", 42);
    
    REQUIRE(result.success == false);
    REQUIRE(result.exitCode == 42);
    REQUIRE(result.errorMessage == "Process failed");
    REQUIRE(result.dataFilePath == "/path/to/file.mml");
}

TEST_CASE("VisualizerResult::Failure - With negative exit code", "[Visualizer][Result]") {
    auto result = VisualizerResult::Failure("Crashed", "", -9);
    
    REQUIRE(result.success == false);
    REQUIRE(result.exitCode == -9);
}

TEST_CASE("VisualizerResult::Failure - Empty message", "[Visualizer][Result]") {
    auto result = VisualizerResult::Failure("");
    
    REQUIRE(result.success == false);
    REQUIRE(result.errorMessage.empty());
}

///////////////////////////////////////////////////////////////////////////////
//                      BOOL CONVERSION OPERATOR                             //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("VisualizerResult - Success converts to true", "[Visualizer][Result]") {
    auto result = VisualizerResult::Success();
    
    REQUIRE(static_cast<bool>(result) == true);
    
    // Test in if statement
    if (result) {
        SUCCEED("Success correctly converts to true");
    } else {
        FAIL("Success should convert to true");
    }
}

TEST_CASE("VisualizerResult - Failure converts to false", "[Visualizer][Result]") {
    auto result = VisualizerResult::Failure("error");
    
    REQUIRE(static_cast<bool>(result) == false);
    
    // Test in if statement
    if (result) {
        FAIL("Failure should convert to false");
    } else {
        SUCCEED("Failure correctly converts to false");
    }
}

TEST_CASE("VisualizerResult - Bool in NOT expressions", "[Visualizer][Result]") {
    auto success = VisualizerResult::Success();
    auto failure = VisualizerResult::Failure("fail");
    
    REQUIRE(!failure == true);
    REQUIRE(!success == false);
}

TEST_CASE("VisualizerResult - Bool in AND expressions", "[Visualizer][Result]") {
    auto s1 = VisualizerResult::Success();
    auto s2 = VisualizerResult::Success();
    auto f = VisualizerResult::Failure("fail");
    
    REQUIRE((s1 && s2) == true);
    REQUIRE((s1 && f) == false);
    REQUIRE((f && s1) == false);
}

TEST_CASE("VisualizerResult - Bool in OR expressions", "[Visualizer][Result]") {
    auto s = VisualizerResult::Success();
    auto f1 = VisualizerResult::Failure("fail1");
    auto f2 = VisualizerResult::Failure("fail2");
    
    REQUIRE((f1 || s) == true);
    REQUIRE((s || f1) == true);
    REQUIRE((f1 || f2) == false);
}

TEST_CASE("VisualizerResult - Ternary operator", "[Visualizer][Result]") {
    auto success = VisualizerResult::Success();
    auto failure = VisualizerResult::Failure("fail");
    
    std::string msg1 = success ? "ok" : "error";
    std::string msg2 = failure ? "ok" : "error";
    
    REQUIRE(msg1 == "ok");
    REQUIRE(msg2 == "error");
}

///////////////////////////////////////////////////////////////////////////////
//                       STRUCT MEMBER ACCESS                                //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("VisualizerResult - Aggregate initialization", "[Visualizer][Result]") {
    VisualizerResult result{true, false, 0, "", "/data/file.mml"};
    
    REQUIRE(result.success == true);
    REQUIRE(result.exitCode == 0);
    REQUIRE(result.errorMessage == "");
    REQUIRE(result.dataFilePath == "/data/file.mml");
}

TEST_CASE("VisualizerResult - Failure aggregate initialization", "[Visualizer][Result]") {
    VisualizerResult result{false, false, 127, "Command not found", "/usr/bin/viz"};
    
    REQUIRE(result.success == false);
    REQUIRE(result.exitCode == 127);
    REQUIRE(result.errorMessage == "Command not found");
    REQUIRE(result.dataFilePath == "/usr/bin/viz");
}

TEST_CASE("VisualizerResult - Members are modifiable", "[Visualizer][Result]") {
    VisualizerResult result = VisualizerResult::Success();
    
    result.success = false;
    result.exitCode = 99;
    result.errorMessage = "Modified";
    result.dataFilePath = "/new/path.mml";
    
    REQUIRE(result.success == false);
    REQUIRE(result.exitCode == 99);
    REQUIRE(result.errorMessage == "Modified");
    REQUIRE(result.dataFilePath == "/new/path.mml");
}

///////////////////////////////////////////////////////////////////////////////
//                       STRING CONTENT TESTS                                //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("VisualizerResult - Long error message", "[Visualizer][Result]") {
    std::string longError = "This is a very long error message that describes "
                            "in great detail what went wrong during the visualization "
                            "process, including potentially stack traces, file paths, "
                            "and other diagnostic information that would be helpful "
                            "for debugging purposes.";
    
    auto result = VisualizerResult::Failure(longError);
    
    REQUIRE(result.errorMessage == longError);
    REQUIRE(result.errorMessage.length() > 200);
}

TEST_CASE("VisualizerResult - Unicode in paths", "[Visualizer][Result]") {
    std::string unicodePath = "/path/to/données.mml";
    auto result = VisualizerResult::Success(unicodePath);
    
    REQUIRE(result.dataFilePath == unicodePath);
}

TEST_CASE("VisualizerResult - Special characters in error message", "[Visualizer][Result]") {
    std::string specialMsg = "Error: File \"test.txt\" not found!\n\tCheck path & retry.";
    auto result = VisualizerResult::Failure(specialMsg);
    
    REQUIRE(result.errorMessage == specialMsg);
    REQUIRE(result.errorMessage.find('\n') != std::string::npos);
    REQUIRE(result.errorMessage.find('\t') != std::string::npos);
}

///////////////////////////////////////////////////////////////////////////////
//                       COPY AND MOVE SEMANTICS                             //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("VisualizerResult - Copy construction", "[Visualizer][Result]") {
    auto original = VisualizerResult::Failure("original error", "/original/path.mml", 42);
    VisualizerResult copy = original;
    
    REQUIRE(copy.success == original.success);
    REQUIRE(copy.exitCode == original.exitCode);
    REQUIRE(copy.errorMessage == original.errorMessage);
    REQUIRE(copy.dataFilePath == original.dataFilePath);
    
    // Modify copy, original should be unchanged
    copy.errorMessage = "modified";
    REQUIRE(original.errorMessage == "original error");
}

TEST_CASE("VisualizerResult - Copy assignment", "[Visualizer][Result]") {
    auto original = VisualizerResult::Success("/path/to/file.mml");
    VisualizerResult copy = VisualizerResult::Failure("dummy");
    
    copy = original;
    
    REQUIRE(copy.success == true);
    REQUIRE(copy.dataFilePath == "/path/to/file.mml");
}

TEST_CASE("VisualizerResult - Move construction", "[Visualizer][Result]") {
    auto original = VisualizerResult::Success("/path/to/file.mml");
    VisualizerResult moved = std::move(original);
    
    REQUIRE(moved.success == true);
    REQUIRE(moved.dataFilePath == "/path/to/file.mml");
}

TEST_CASE("VisualizerResult - Move assignment", "[Visualizer][Result]") {
    auto original = VisualizerResult::Failure("error", "/path.mml", 10);
    VisualizerResult target = VisualizerResult::Success();
    
    target = std::move(original);
    
    REQUIRE(target.success == false);
    REQUIRE(target.exitCode == 10);
}

///////////////////////////////////////////////////////////////////////////////
//                       EDGE CASES                                          //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("VisualizerResult - Default construction", "[Visualizer][Result]") {
    VisualizerResult result{};
    
    REQUIRE(result.success == false);  // bool default-initializes to false
    REQUIRE(result.exitCode == 0);     // int default-initializes to 0
    REQUIRE(result.errorMessage.empty());
    REQUIRE(result.dataFilePath.empty());
}

TEST_CASE("VisualizerResult - Exit code zero with failure", "[Visualizer][Result]") {
    // Rare case: process exits with 0 but we mark as failure
    VisualizerResult result{false, false, 0, "Logic error", ""};
    
    REQUIRE(result.success == false);
    REQUIRE(result.exitCode == 0);
    REQUIRE(static_cast<bool>(result) == false);  // success field determines bool
}

TEST_CASE("VisualizerResult - Non-zero exit with success", "[Visualizer][Result]") {
    // Rare case: marking as success despite non-zero exit (shouldn't happen normally)
    VisualizerResult result{true, false, 1, "", ""};
    
    REQUIRE(result.success == true);
    REQUIRE(result.exitCode == 1);
    REQUIRE(static_cast<bool>(result) == true);  // success field determines bool
}

} // namespace MML::Tests::Tools::VisualizerTests

