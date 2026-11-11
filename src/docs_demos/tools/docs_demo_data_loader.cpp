///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: DataLoader - Multi-format Dataset Loading
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml/MMLBase.h>
#include <mml/tools/DataLoader.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <filesystem>
#include <fstream>

using namespace MML;
using namespace MML::Data;

namespace MML::docs_demos::data_loader
{
    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              CSV LOADING                                         ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_csv_loading()
    {
        std::cout << "\n=== CSV Loading ===\n";

        // Sample CSV content
        std::string csvContent = 
            "Name,Age,Score,Active\n"
            "Alice,25,95.5,true\n"
            "Bob,30,88.2,false\n"
            "Charlie,22,91.8,true\n"
            "Diana,28,87.3,true\n";

        // Load from string
        Dataset data = LoadFromCSVString(csvContent, ',', true);

        std::cout << "Loaded " << data.NumColumns() << " columns, " 
                  << data.NumRows() << " rows\n\n";

        // Print column info
        for (size_t i = 0; i < data.NumColumns(); ++i) {
            const auto& col = data[i];
            std::cout << "Column '" << col.name << "' (" << ColumnTypeToString(col.type) << "): ";
            
            // Print first few values
            size_t printCount = std::min(size_t(3), col.Size());
            for (size_t j = 0; j < printCount; ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << col.GetAsString(j);
            }
            if (col.Size() > 3) std::cout << ", ...";
            std::cout << "\n";
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              TSV LOADING                                         ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_tsv_loading()
    {
        std::cout << "\n=== TSV Loading ===\n";

        // Sample TSV content (tab-separated)
        std::string tsvContent = 
            "x\ty\tz\n"
            "0.0\t0.0\t1.0\n"
            "1.0\t0.5\t2.0\n"
            "2.0\t1.0\t3.0\n"
            "3.0\t1.5\t4.0\n";

        // Load from string with tab delimiter
        Dataset data = LoadFromCSVString(tsvContent, '\t', true);

        std::cout << "Loaded " << data.NumColumns() << " columns\n";
        
        // Access numeric columns directly
        for (size_t c = 0; c < data.NumColumns(); ++c) {
            const auto& col = data[c];
            std::cout << col.name << ": [";
            for (size_t i = 0; i < col.Size(); ++i) {
                if (i > 0) std::cout << ", ";
                if (col.type == ColumnType::REAL) {
                    std::cout << col.realData[i];
                } else if (col.type == ColumnType::INT) {
                    std::cout << col.intData[i];
                } else {
                    std::cout << col.GetAsString(i);
                }
            }
            std::cout << "]\n";
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              JSON LOADING                                        ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_json_loading()
    {
        std::cout << "\n=== JSON Loading ===\n";

        // Sample JSON content (array of objects)
        std::string jsonContent = R"([
            {"name": "Alpha", "value": 1.5, "count": 10},
            {"name": "Beta", "value": 2.7, "count": 20},
            {"name": "Gamma", "value": 3.2, "count": 15}
        ])";

        Dataset data = LoadFromJSONString(jsonContent);

        std::cout << "Loaded " << data.NumColumns() << " columns from JSON\n";
        
        for (size_t c = 0; c < data.NumColumns(); ++c) {
            const auto& col = data[c];
            std::cout << col.name << " (" << ColumnTypeToString(col.type) << "): ";
            
            for (size_t i = 0; i < col.Size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << col.GetAsString(i);
            }
            std::cout << "\n";
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              FORMAT DETECTION                                    ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_format_detection()
    {
        std::cout << "\n=== Format Detection ===\n";

        std::vector<std::string> filenames = {
            "data.csv",
            "measurements.tsv",
            "results.tab",
            "config.json",
            "unknown_file"
        };

        for (const auto& filename : filenames) {
            DataFormat format = DetectFormat(filename);
            std::string formatStr;
            switch (format) {
                case DataFormat::CSV:  formatStr = "CSV"; break;
                case DataFormat::TSV:  formatStr = "TSV"; break;
                case DataFormat::JSON: formatStr = "JSON"; break;
                case DataFormat::Auto: formatStr = "Auto"; break;
            }
            std::cout << filename << " -> " << formatStr << "\n";
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              TYPE INFERENCE                                      ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_type_inference()
    {
        std::cout << "\n=== Type Inference ===\n";

        std::vector<std::string> values = {
            "123",
            "45.67",
            "true",
            "false",
            "Hello",
            "-999",
            "3.14159",
            "yes",
            "no"
        };

        for (const auto& val : values) {
            // InferColumnType takes a vector of values to infer the column type
            std::vector<std::string> singleValue = { val };
            ColumnType type = InferColumnType(singleValue);
            std::cout << "\"" << val << "\" -> " << ColumnTypeToString(type) << "\n";
        }
    }

    void demo_numeric_import_contract()
    {
        std::cout << "\n=== Numeric Import Contract ===\n";
        const std::string csv = "value\n1\n2\n3\n4\n5\n6\n7\n8\n9\n10.75\n";
        Dataset data = LoadFromCSVString(csv);
        if (data["value"].type != ColumnType::REAL || data["value"].realData[9] != Real(10.75) ||
            data["value"].CountMissing() != 0)
            throw DataError("Numeric import demo lost the fractional value");
        std::cout << data["value"].realData[9] << '\n';

        std::string invalid = csv;
        invalid.replace(invalid.find("10.75"), 5, "10oops");
        bool rejected = false;
        try {
            LoadFromCSVString(invalid);
        } catch (const DataError& error) {
            rejected = true;
            const std::string message = error.what();
            if (message.find("10oops") == std::string::npos || message.find("data row 10") == std::string::npos)
                throw;
            std::cout << "Rejected invalid input: " << message << '\n';
        }
        if (!rejected)
            throw DataError("Numeric import demo accepted trailing junk");
    }

    void demo_csv_record_contract()
    {
        std::cout << "\n=== CSV Record Contract ===\n";
        const std::string csv = "name,value\r\n\"hello\r\nworld\",1\r\n";
        Dataset data = LoadFromCSVString(csv);
        std::cout << data.NumRows() << '\n';
        std::cout << data["name"].stringData[0] << '\n';
        if (data.NumRows() != 1 || data["name"].stringData[0] != "hello\r\nworld" ||
            data["value"].intData[0] != 1)
            throw DataError("CSV record demo lost a multiline field or its original CRLF");

        bool rejected = false;
        try {
            LoadFromCSVString("name,value\n\"closed\"junk,1\n");
        } catch (const DataError& error) {
            rejected = true;
            if (std::string(error.what()).find("CSV line 2, field 1") == std::string::npos)
                throw;
            std::cout << "Rejected malformed quotes: " << error.what() << '\n';
        }
        if (!rejected)
            throw DataError("CSV record demo accepted trailing junk after a quoted field");
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              SAFE LOADING (WITH RESULT)                          ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_safe_loading()
    {
        std::cout << "\n=== Safe Loading (Error Handling) ===\n";

        const auto directory = std::filesystem::temp_directory_path() /
            ("mml_data_loader_demo_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        if (!std::filesystem::create_directory(directory))
            throw FileIOError("Cannot create DataLoader demo directory", directory.string());
        struct Cleanup {
            std::filesystem::path directory;
            ~Cleanup() {
                std::error_code error;
                std::filesystem::remove_all(directory, error);
                if (error)
                    std::cerr << "Cannot clean DataLoader demo directory: " << error.message() << '\n';
            }
        } cleanup{directory};
        auto write = [&directory](const std::string& name, const std::string& contents) {
            const std::string filename = (directory / name).string();
            std::ofstream output(filename, std::ios::binary);
            if (!output)
                throw FileIOError("Cannot create DataLoader demo fixture", filename);
            output << contents;
            const IOResult written = FinishOutputFile(output, filename, IOResult::Success(filename));
            if (!written)
                throw FileIOError(written.message, written.filename);
            return filename;
        };
        auto checkSample = [](const Dataset& data) {
            if (data.NumRows() != 2 || data.GetColumnNames() != std::vector<std::string>{"a", "b", "c"})
                throw DataError("DataLoader demo loaded an unexpected dataset shape or schema");
            for (size_t column = 0; column < data.NumColumns(); ++column) {
                if (data[column].type != ColumnType::INT || data[column].CountMissing() != 0)
                    throw DataError("DataLoader demo lost integer column types or values");
                for (size_t row = 0; row < data.NumRows(); ++row)
                    if (data[column].GetReal(row) != static_cast<Real>(1 + row * 3 + column))
                        throw DataError("DataLoader demo loaded incorrect values");
            }
        };

        const std::string goodCSV = "a,b,c\n1,2,3\n4,5,6\n";
        const Dataset parsed = LoadFromCSVString(goodCSV);
        checkSample(parsed);
        std::cout << "In-memory CSV: verified 2 rows and 3 typed columns\n";

        const std::string csvFilename = write("sample.csv", goodCSV);
        const std::string tsvFilename = write("sample.tsv", "a\tb\tc\n1\t2\t3\n4\t5\t6\n");
        const std::string tabFilename = write("sample.tab", "a\tb\tc\n1\t2\t3\n4\t5\t6\n");
        const std::string jsonFilename = write("sample.json", R"([{"a":1,"b":2,"c":3},{"a":4,"b":5,"c":6}])");
        Dataset loaded;
        const IOResult result = LoadCSVSafe(csvFilename, loaded, true, ',');
        if (!result)
            throw FileIOError(result.message, result.filename);
        if (result.error != IOError::OK || result.filename != csvFilename)
            throw DataError("Successful CSV load lost its status or filename");
        checkSample(loaded);
        std::cout << "CSV file load: " << IOErrorName(result.error) << '\n';

        checkSample(LoadCSV(csvFilename));
        checkSample(LoadTSV(tsvFilename));
        checkSample(LoadJSON(jsonFilename));
        for (const auto& filename : {csvFilename, tsvFilename, tabFilename, jsonFilename}) {
            checkSample(LoadData(filename));
            Dataset automatic;
            const IOResult detected = LoadData(filename, automatic);
            if (!detected || detected.filename != filename)
                throw FileIOError(detected.message, filename);
            checkSample(automatic);
        }
        Dataset tsv, json;
        const IOResult tsvResult = LoadTSVSafe(tsvFilename, tsv);
        const IOResult jsonResult = LoadJSONSafe(jsonFilename, json);
        if (!tsvResult)
            throw FileIOError(tsvResult.message, tsvResult.filename);
        if (!jsonResult)
            throw FileIOError(jsonResult.message, jsonResult.filename);
        checkSample(tsv);
        checkSample(json);
        checkSample(LoadFromString(goodCSV, DataFormat::CSV));
        checkSample(LoadFromString("a\tb\tc\n1\t2\t3\n4\t5\t6\n", DataFormat::TSV));
        checkSample(LoadFromString(R"([{"a":1,"b":2,"c":3},{"a":4,"b":5,"c":6}])", DataFormat::JSON));
        std::cout << "CSV, TSV, TAB and JSON file/string entry points verified\n";

        const std::string noHeaderFilename = write("no_header.csv", "1.25;alpha\n2.5;beta\n");
        CSVLoadOptions options;
        options.hasHeader = false;
        options.delimiter = ';';
        options.inferTypes = false;
        options.strictSchema = true;
        Dataset strings = LoadCSV(noHeaderFilename, options);
        if (strings.NumRows() != 2 || strings.GetColumnNames() != std::vector<std::string>{"Column0", "Column1"} ||
            strings[0].type != ColumnType::STRING || strings[1].type != ColumnType::STRING ||
            strings[0].stringData != std::vector<std::string>{"1.25", "2.5"} ||
            strings[1].stringData != std::vector<std::string>{"alpha", "beta"})
            throw DataError("CSVLoadOptions demo did not preserve headerless string columns");
        strings.columns[0].name = "x";
        strings.columns[1].name = "label";
        if (!strings.HasColumn("x") || strings["label"].GetAsString(1) != "beta")
            throw DataError("DataLoader demo column naming failed");
        Dataset explicitOptions;
        const IOResult configured = LoadCSV(noHeaderFilename, explicitOptions, options);
        if (!configured || explicitOptions[0].stringData != strings[0].stringData)
            throw DataError("Result-based CSVLoadOptions demo failed");
        std::cout << "CSVLoadOptions: delimiter, header, string storage and column naming verified\n";

        const std::string originalName = loaded.name;
        const std::string missingFilename = (directory / "missing.csv").string();
        const IOResult missing = LoadCSVSafe(missingFilename, loaded);
        if (missing || missing.error != IOError::FILE_NOT_OPENED ||
            missing.filename != missingFilename || missing.message.empty() || loaded.name != originalName)
            throw DataError("Missing CSV did not report FILE_NOT_OPENED or preserve the destination");
        checkSample(loaded);
        std::cout << "Expected missing-file failure: " << IOErrorName(missing.error) << ": " << missing.message << '\n';

        const std::string malformedFilename = write("malformed.csv", "a,b,c\n\"closed\"junk,2,3\n");
        const IOResult malformed = LoadCSVSafe(malformedFilename, loaded);
        if (malformed || malformed.error != IOError::MALFORMED_INPUT || malformed.filename != malformedFilename ||
            malformed.offset != 2 || malformed.offset_kind != IOOffsetKind::Line || loaded.name != originalName)
            throw DataError("Malformed CSV did not preserve its destination and line diagnostic");
        checkSample(loaded);
        std::cout << "Expected malformed-input failure at line " << malformed.offset << ": " << malformed.message << '\n';
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              COLUMN OPERATIONS                                   ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_column_operations()
    {
        std::cout << "\n=== Column Operations ===\n";

        // Create a dataset
        std::string csv = "x,y,label\n1.0,2.0,A\n2.0,4.0,B\n3.0,6.0,C\n4.0,8.0,D\n5.0,10.0,E\n";
        Dataset data = LoadFromCSVString(csv, ',', true);

        const DataColumn& xColumn = data["x"];
        const DataColumn& labelColumn = data["label"];
        if (data.NumRows() != 5 || data.NumColumns() != 3 ||
            xColumn.type != ColumnType::REAL || labelColumn.type != ColumnType::STRING)
            throw DataError("Column access demo loaded incorrect types or dimensions");
        for (size_t row = 0; row < data.NumRows(); ++row) {
            if (xColumn.IsMissing(row) || xColumn.GetReal(row) != static_cast<Real>(row + 1))
                throw DataError("Column access demo lost a numeric value");
            std::cout << labelColumn.GetAsString(row) << ": " << xColumn.GetReal(row) << '\n';
        }

        // Get column by name
        Vector<Real> xVec = data.GetRealColumn("x");
        Vector<Real> yVec = data.GetRealColumn("y");
        if (xVec != Vector<Real>{1, 2, 3, 4, 5} || yVec != Vector<Real>{2, 4, 6, 8, 10} ||
            data.GetStringColumn("label") != std::vector<std::string>{"A", "B", "C", "D", "E"})
            throw DataError("Dataset column extraction demo failed");

        std::cout << "Column 'x' has " << xVec.size() << " elements\n";
        std::cout << "Column 'y' has " << yVec.size() << " elements\n";

        // Compute statistics
        double sumX = 0.0, sumY = 0.0;
        for (int i = 0; i < xVec.size(); ++i) {
            sumX += xVec[i];
            sumY += yVec[i];
        }

        std::cout << "Sum of x: " << sumX << "\n";
        std::cout << "Sum of y: " << sumY << "\n";
        std::cout << "Mean x: " << sumX / xVec.size() << "\n";
        std::cout << "Mean y: " << sumY / yVec.size() << "\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              DATASET SUMMARY                                     ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_dataset_summary()
    {
        std::cout << "\n=== Dataset Summary ===\n";

        std::string csv = 
            "id,name,score,passed\n"
            "1,Alice,95.5,true\n"
            "2,Bob,88.2,true\n"
            "3,Charlie,72.1,false\n"
            "4,Diana,91.3,true\n";

        Dataset data = LoadFromCSVString(csv, ',', true);
        data.name = "Student Records";

        // Print dataset summary
        std::cout << data.PrintSummary() << "\n";
    }

    void Run()
    {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "       MML Documentation Demo: DataLoader\n";
        std::cout << std::string(70, '=') << "\n";

        demo_csv_loading();
        demo_tsv_loading();
        demo_json_loading();
        demo_format_detection();
        demo_type_inference();
        demo_numeric_import_contract();
        demo_csv_record_contract();
        demo_safe_loading();
        demo_column_operations();
        demo_dataset_summary();

        std::cout << "\n" << std::string(70, '=') << "\n";
    }
} // namespace

void Docs_Demo_DataLoader() {
    MML::docs_demos::data_loader::Run();
}

void Docs_Demo_DataLoaderUnicode() {
    const Dataset data = LoadFromJSONString(
        R"([{"name":"Caf\u00e9","active":"TrUe"},{"name":"N\u00c4N","active":"FaLsE"}])");
    if (data["name"].stringData[0] != "Caf\xC3\xA9" ||
        data["name"].stringData[1] != "N\xC3\x84" "N" ||
        data["name"].CountMissing() != 0 ||
        !data["active"].boolData[0] || data["active"].boolData[1])
        throw DataError("DataLoader Unicode and ASCII token contract failed");
    std::cout << data["name"].stringData[0] << '\n';
}
