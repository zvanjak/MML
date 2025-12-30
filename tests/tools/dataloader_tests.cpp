///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        dataloader_tests.cpp                                                ///
///  Description: Unit tests for DataLoader                                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include <cmath>
#include <fstream>
#include <string>

#include "MMLBase.h"
#include "tools/DataLoader.h"

using namespace MML;
using namespace MML::Data;

#define TEST_PRECISION_INFO() \
	INFO("Test precision: " << Constants::Eps << " (" << (sizeof(Real) == 8 ? "double" : "float") << ")")

namespace MML::Tests::Tools::DataLoaderTests {

	/////////////////////////////////////////////////////////////////////////////////////
	///                         MISSING VALUE DETECTION                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Missing value detection", "[dataloader][missing]") {
		TEST_PRECISION_INFO();

		SECTION("Common missing value patterns") {
			REQUIRE(IsMissingValue(""));
			REQUIRE(IsMissingValue("NA"));
			REQUIRE(IsMissingValue("na"));
			REQUIRE(IsMissingValue("N/A"));
			REQUIRE(IsMissingValue("n/a"));
			REQUIRE(IsMissingValue("NaN"));
			REQUIRE(IsMissingValue("nan"));
			REQUIRE(IsMissingValue("null"));
			REQUIRE(IsMissingValue("NULL"));
			REQUIRE(IsMissingValue("none"));
			REQUIRE(IsMissingValue("None"));
			REQUIRE(IsMissingValue("."));
			REQUIRE(IsMissingValue("?"));
			REQUIRE(IsMissingValue("-"));
		}

		SECTION("Non-missing values") {
			REQUIRE_FALSE(IsMissingValue("0"));
			REQUIRE_FALSE(IsMissingValue("1"));
			REQUIRE_FALSE(IsMissingValue("hello"));
			REQUIRE_FALSE(IsMissingValue("3.14"));
			REQUIRE_FALSE(IsMissingValue("true"));
			REQUIRE_FALSE(IsMissingValue("false"));
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         TYPE INFERENCE                                        ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Type inference", "[dataloader][inference]") {
		TEST_PRECISION_INFO();

		SECTION("Integer detection") {
			std::vector<std::string> values = {"1", "2", "3", "42", "-10", "0"};
			REQUIRE(InferColumnType(values) == ColumnType::INT);
		}

		SECTION("Real detection") {
			std::vector<std::string> values = {"1.5", "2.0", "3.14159", "-0.5", "1e10"};
			REQUIRE(InferColumnType(values) == ColumnType::REAL);
		}

		SECTION("Mixed int and real -> REAL") {
			std::vector<std::string> values = {"1", "2.5", "3", "4.0", "5"};
			REQUIRE(InferColumnType(values) == ColumnType::REAL);
		}

		SECTION("Boolean detection") {
			std::vector<std::string> values = {"true", "false", "True", "FALSE", "yes", "no"};
			REQUIRE(InferColumnType(values) == ColumnType::BOOL);
		}

		SECTION("Date detection") {
			std::vector<std::string> values = {"2024-01-15", "2024-02-20", "2024-12-31"};
			REQUIRE(InferColumnType(values) == ColumnType::DATE);
		}

		SECTION("Time detection") {
			std::vector<std::string> values = {"10:30:00", "14:45:30", "08:00:00"};
			REQUIRE(InferColumnType(values) == ColumnType::TIME);
		}

		SECTION("DateTime detection") {
			std::vector<std::string> values = {"2024-01-15T10:30:00", "2024-02-20 14:45:30"};
			REQUIRE(InferColumnType(values) == ColumnType::DATETIME);
		}

		SECTION("String fallback") {
			std::vector<std::string> values = {"hello", "world", "foo", "bar"};
			REQUIRE(InferColumnType(values) == ColumnType::STRING);
		}

		SECTION("With missing values") {
			std::vector<std::string> values = {"1", "NA", "3", "null", "5"};
			REQUIRE(InferColumnType(values) == ColumnType::INT);
		}

		SECTION("Empty input") {
			std::vector<std::string> values;
			REQUIRE(InferColumnType(values) == ColumnType::STRING);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         STRING UTILITIES                                      ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - String utilities", "[dataloader][utility]") {
		TEST_PRECISION_INFO();

		SECTION("Simple split") {
			auto result = SplitLine("a,b,c", ',');
			REQUIRE(result.size() == 3);
			REQUIRE(result[0] == "a");
			REQUIRE(result[1] == "b");
			REQUIRE(result[2] == "c");
		}

		SECTION("Split with quotes") {
			auto result = SplitLine("\"hello,world\",b,c", ',');
			REQUIRE(result.size() == 3);
			REQUIRE(result[0] == "hello,world");
			REQUIRE(result[1] == "b");
			REQUIRE(result[2] == "c");
		}

		SECTION("Split with escaped quotes") {
			auto result = SplitLine("\"he said \"\"hello\"\"\",b", ',');
			REQUIRE(result.size() == 2);
			REQUIRE(result[0] == "he said \"hello\"");
		}

		SECTION("Empty fields") {
			auto result = SplitLine("a,,c", ',');
			REQUIRE(result.size() == 3);
			REQUIRE(result[0] == "a");
			REQUIRE(result[1] == "");
			REQUIRE(result[2] == "c");
		}

		SECTION("Trim whitespace") {
			REQUIRE(Trim("  hello  ") == "hello");
			REQUIRE(Trim("\t\ntest\r\n") == "test");
			REQUIRE(Trim("   ") == "");
			REQUIRE(Trim("") == "");
		}

		SECTION("Remove BOM") {
			std::string withBOM = "\xEF\xBB\xBFhello";
			REQUIRE(RemoveBOM(withBOM) == "hello");

			std::string noBOM = "hello";
			REQUIRE(RemoveBOM(noBOM) == "hello");
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         CSV LOADING FROM STRING                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - CSV from string", "[dataloader][csv]") {
		TEST_PRECISION_INFO();

		SECTION("Simple CSV with header") {
			std::string csv = "name,age,score\nAlice,30,95.5\nBob,25,87.3\nCharlie,35,92.1";

			Dataset ds = LoadFromString(csv, DataFormat::CSV, true);

			REQUIRE(ds.NumColumns() == 3);
			REQUIRE(ds.NumRows() == 3);
			REQUIRE(ds.HasColumn("name"));
			REQUIRE(ds.HasColumn("age"));
			REQUIRE(ds.HasColumn("score"));

			REQUIRE(ds["name"].type == ColumnType::STRING);
			REQUIRE(ds["age"].type == ColumnType::INT);
			REQUIRE(ds["score"].type == ColumnType::REAL);

			REQUIRE(ds["name"].stringData[0] == "Alice");
			REQUIRE(ds["name"].stringData[1] == "Bob");
			REQUIRE(ds["age"].intData[0] == 30);
			REQUIRE(ds["age"].intData[1] == 25);
			REQUIRE(ds["score"].realData[0] == Catch::Approx(95.5));
		}

		SECTION("CSV without header") {
			std::string csv = "1,2,3\n4,5,6\n7,8,9";

			Dataset ds = LoadFromString(csv, DataFormat::CSV, false);

			REQUIRE(ds.NumColumns() == 3);
			REQUIRE(ds.NumRows() == 3);
			REQUIRE(ds["Column0"].intData[0] == 1);
			REQUIRE(ds["Column2"].intData[2] == 9);
		}

		SECTION("CSV with missing values") {
			std::string csv = "a,b,c\n1,NA,3\n4,,6\n7,8,null";

			Dataset ds = LoadFromString(csv, DataFormat::CSV, true);

			REQUIRE(ds.NumRows() == 3);
			REQUIRE(ds["a"].IsMissing(0) == false);
			REQUIRE(ds["b"].IsMissing(0) == true);  // NA
			REQUIRE(ds["b"].IsMissing(1) == true);  // empty
			REQUIRE(ds["c"].IsMissing(2) == true);  // null
		}

		SECTION("CSV with quoted fields") {
			std::string csv = "name,description\n\"Smith, John\",\"A \"\"great\"\" person\"";

			Dataset ds = LoadFromString(csv, DataFormat::CSV, true);

			REQUIRE(ds.NumRows() == 1);
			REQUIRE(ds["name"].stringData[0] == "Smith, John");
			REQUIRE(ds["description"].stringData[0] == "A \"great\" person");
		}

		SECTION("Boolean column") {
			std::string csv = "flag\ntrue\nfalse\nyes\nno";

			Dataset ds = LoadFromString(csv, DataFormat::CSV, true);

			REQUIRE(ds["flag"].type == ColumnType::BOOL);
			REQUIRE(ds["flag"].boolData[0] == true);
			REQUIRE(ds["flag"].boolData[1] == false);
			REQUIRE(ds["flag"].boolData[2] == true);
			REQUIRE(ds["flag"].boolData[3] == false);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         TSV LOADING                                           ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - TSV from string", "[dataloader][tsv]") {
		TEST_PRECISION_INFO();

		SECTION("Simple TSV") {
			std::string tsv = "x\ty\tz\n1.0\t2.0\t3.0\n4.0\t5.0\t6.0";

			Dataset ds = LoadFromString(tsv, DataFormat::TSV, true);

			REQUIRE(ds.NumColumns() == 3);
			REQUIRE(ds.NumRows() == 2);
			REQUIRE(ds["x"].realData[0] == Catch::Approx(1.0));
			REQUIRE(ds["z"].realData[1] == Catch::Approx(6.0));
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         JSON LOADING                                          ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - JSON from string", "[dataloader][json]") {
		TEST_PRECISION_INFO();

		SECTION("Simple JSON array") {
			std::string json = R"([
				{"name": "Alice", "age": 30, "active": true},
				{"name": "Bob", "age": 25, "active": false},
				{"name": "Charlie", "age": 35, "active": true}
			])";

			Dataset ds = LoadFromString(json, DataFormat::JSON);

			REQUIRE(ds.NumColumns() == 3);
			REQUIRE(ds.NumRows() == 3);
			REQUIRE(ds.HasColumn("name"));
			REQUIRE(ds.HasColumn("age"));
			REQUIRE(ds.HasColumn("active"));

			REQUIRE(ds["name"].stringData[0] == "Alice");
			REQUIRE(ds["age"].intData[1] == 25);
			REQUIRE(ds["active"].boolData[2] == true);
		}

		SECTION("JSON with null values") {
			std::string json = R"([
				{"value": 1},
				{"value": null},
				{"value": 3}
			])";

			Dataset ds = LoadFromString(json, DataFormat::JSON);

			REQUIRE(ds.NumRows() == 3);
			REQUIRE(ds["value"].IsMissing(1) == true);
		}

		SECTION("JSON with mixed numeric") {
			std::string json = R"([
				{"x": 1},
				{"x": 2.5},
				{"x": 3}
			])";

			Dataset ds = LoadFromString(json, DataFormat::JSON);

			REQUIRE(ds["x"].type == ColumnType::REAL);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         DATASET ACCESS METHODS                                ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Dataset access", "[dataloader][dataset]") {
		TEST_PRECISION_INFO();

		std::string csv = "name,value,flag\nA,1.5,true\nB,2.5,false\nC,3.5,true\nD,4.5,false\nE,5.5,true";
		Dataset ds = LoadFromString(csv, DataFormat::CSV, true);

		SECTION("Column access by name") {
			REQUIRE(ds["name"].name == "name");
			REQUIRE(ds["value"].realData[0] == Catch::Approx(1.5));
			REQUIRE_THROWS_AS(ds["nonexistent"], DataError);
		}

		SECTION("Column access by index") {
			REQUIRE(ds[0].name == "name");
			REQUIRE(ds[1].name == "value");
			REQUIRE_THROWS_AS(ds[100], DataError);
		}

		SECTION("GetRealColumn") {
			Vector<Real> values = ds.GetRealColumn("value");
			REQUIRE(values.size() == 5);
			REQUIRE(values[0] == Catch::Approx(1.5));
			REQUIRE(values[4] == Catch::Approx(5.5));
		}

		SECTION("GetStringColumn") {
			auto names = ds.GetStringColumn("name");
			REQUIRE(names.size() == 5);
			REQUIRE(names[0] == "A");
			REQUIRE(names[4] == "E");
		}

		SECTION("GetBoolColumn") {
			std::vector<bool> flags = ds.GetBoolColumn("flag");
			REQUIRE(flags.size() == 5);
			REQUIRE(flags[0] == true);
			REQUIRE(flags[1] == false);
		}

		SECTION("Column names") {
			auto names = ds.GetColumnNames();
			REQUIRE(names.size() == 3);
			REQUIRE(names[0] == "name");
			REQUIRE(names[1] == "value");
			REQUIRE(names[2] == "flag");
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         HEAD/TAIL/SAMPLE                                       ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Head/Tail/Sample", "[dataloader][display]") {
		TEST_PRECISION_INFO();

		std::string csv = "id\n1\n2\n3\n4\n5\n6\n7\n8\n9\n10";
		Dataset ds = LoadFromString(csv, DataFormat::CSV, true);

		SECTION("Head default") {
			std::string head = ds.Head();
			REQUIRE(head.find("1") != std::string::npos);
			REQUIRE(head.find("5") != std::string::npos);
		}

		SECTION("Tail default") {
			std::string tail = ds.Tail();
			REQUIRE(tail.find("10") != std::string::npos);
			REQUIRE(tail.find("6") != std::string::npos);
		}

		SECTION("Sample") {
			std::string sample = ds.Sample(3);
			// Just verify it returns something
			REQUIRE_FALSE(sample.empty());
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         FILTER AND SELECT                                      ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Filter and Select", "[dataloader][transform]") {
		TEST_PRECISION_INFO();

		std::string csv = "name,value,category\nA,10,X\nB,20,Y\nC,30,X\nD,40,Y\nE,50,X";
		Dataset ds = LoadFromString(csv, DataFormat::CSV, true);

		SECTION("FilterRows") {
			Dataset filtered = ds.FilterRows("category", [](const std::string& v) {
				return v == "X";
			});

			REQUIRE(filtered.NumRows() == 3);
			REQUIRE(filtered["name"].stringData[0] == "A");
			REQUIRE(filtered["name"].stringData[1] == "C");
			REQUIRE(filtered["name"].stringData[2] == "E");
		}

		SECTION("SelectColumns") {
			Dataset selected = ds.SelectColumns({"name", "value"});

			REQUIRE(selected.NumColumns() == 2);
			REQUIRE(selected.HasColumn("name"));
			REQUIRE(selected.HasColumn("value"));
			REQUIRE_FALSE(selected.HasColumn("category"));
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         PRINT SUMMARY                                         ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Print summary", "[dataloader][summary]") {
		TEST_PRECISION_INFO();

		std::string csv = "name,score,passed\nAlice,95.5,true\nBob,87.3,true\nCharlie,NA,false";
		Dataset ds = LoadFromString(csv, DataFormat::CSV, true);

		std::string summary = ds.PrintSummary();

		// Check summary contains expected info
		REQUIRE(summary.find("3 rows") != std::string::npos);
		REQUIRE(summary.find("3 columns") != std::string::npos);
		REQUIRE(summary.find("name") != std::string::npos);
		REQUIRE(summary.find("score") != std::string::npos);
		REQUIRE(summary.find("Missing") != std::string::npos);
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         FORMAT DETECTION                                       ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Format detection", "[dataloader][format]") {
		REQUIRE(DetectFormat("data.csv") == DataFormat::CSV);
		REQUIRE(DetectFormat("data.CSV") == DataFormat::CSV);
		REQUIRE(DetectFormat("data.tsv") == DataFormat::TSV);
		REQUIRE(DetectFormat("data.tab") == DataFormat::TSV);
		REQUIRE(DetectFormat("data.json") == DataFormat::JSON);
		REQUIRE(DetectFormat("data.JSON") == DataFormat::JSON);
		REQUIRE(DetectFormat("data.txt") == DataFormat::CSV);  // Default
		REQUIRE(DetectFormat("noextension") == DataFormat::CSV);  // Default
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         ERROR HANDLING                                         ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Error handling", "[dataloader][error]") {
		TEST_PRECISION_INFO();

		SECTION("File not found") {
			REQUIRE_THROWS_AS(LoadCSV("nonexistent_file.csv"), DataError);
			REQUIRE_THROWS_AS(LoadJSON("nonexistent_file.json"), DataError);
		}

		SECTION("Empty content") {
			REQUIRE_THROWS_AS(LoadFromString("", DataFormat::CSV), DataError);
		}

		SECTION("Invalid column access") {
			std::string csv = "a,b\n1,2";
			Dataset ds = LoadFromString(csv, DataFormat::CSV);
			REQUIRE_THROWS_AS(ds["nonexistent"], DataError);
			REQUIRE_THROWS_AS(ds.GetRealColumn("nonexistent"), DataError);  // column doesn't exist
			
			// GetRealColumn allows INT->REAL conversion (no throw)
			REQUIRE_NOTHROW(ds.GetRealColumn("a"));
			REQUIRE(ds.GetRealColumn("a")[0] == Catch::Approx(1.0));
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         REAL-WORLD DATA PATTERNS                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DataLoader - Real-world patterns", "[dataloader][realworld]") {
		TEST_PRECISION_INFO();

		SECTION("Scientific notation") {
			std::string csv = "value\n1.5e10\n2.3e-5\n-1.2E+3";
			Dataset ds = LoadFromString(csv, DataFormat::CSV);

			REQUIRE(ds["value"].type == ColumnType::REAL);
			REQUIRE(ds["value"].realData[0] == Catch::Approx(1.5e10));
			REQUIRE(ds["value"].realData[1] == Catch::Approx(2.3e-5));
		}

		SECTION("Dates in various formats") {
			std::string csv = "date\n2024-01-15\n2024/02/20\n2024-12-31";
			Dataset ds = LoadFromString(csv, DataFormat::CSV);

			REQUIRE(ds["date"].type == ColumnType::DATE);
		}

		SECTION("Large dataset simulation") {
			// Generate larger CSV
			std::ostringstream oss;
			oss << "id,x,y\n";
			for (int i = 0; i < 1000; ++i) {
				oss << i << "," << (i * 0.1) << "," << (i * i) << "\n";
			}

			Dataset ds = LoadFromString(oss.str(), DataFormat::CSV);

			REQUIRE(ds.NumRows() == 1000);
			REQUIRE(ds["id"].type == ColumnType::INT);
			REQUIRE(ds["x"].type == ColumnType::REAL);
			REQUIRE(ds["y"].type == ColumnType::INT);
		}
	}

}  // namespace MML::Tests::Tools::DataLoaderTests
