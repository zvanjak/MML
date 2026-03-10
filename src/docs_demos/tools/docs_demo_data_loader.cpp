///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: DataLoader - Multi-format Dataset Loading
///////////////////////////////////////////////////////////////////////////////////////////

#include "MMLBase.h"
#include "tools/DataLoader.h"

#include <iostream>
#include <sstream>
#include <iomanip>

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

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              SAFE LOADING (WITH RESULT)                          ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_safe_loading()
    {
        std::cout << "\n=== Safe Loading (Error Handling) ===\n";

        // Well-formed data
        std::string goodCSV = "a,b,c\n1,2,3\n4,5,6\n";
        LoadResult result = LoadCSVSafe(goodCSV, true, true);  // fromString=true
        
        if (result.success) {
            std::cout << "Good CSV: Loaded successfully (" 
                      << result.data.NumColumns() << " columns)\n";
        }

        // Example of handling potential file load
        std::cout << "\nNote: Use LoadCSVSafe/LoadJSONSafe for file I/O ";
        std::cout << "with proper error handling.\n";
        std::cout << "The LoadResult struct contains:\n";
        std::cout << "  - success: bool\n";
        std::cout << "  - data: Dataset\n";
        std::cout << "  - errorMessage: std::string\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              COLUMN OPERATIONS                                   ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_column_operations()
    {
        std::cout << "\n=== Column Operations ===\n";

        // Create a dataset
        std::string csv = "x,y\n1.0,2.0\n2.0,4.0\n3.0,6.0\n4.0,8.0\n5.0,10.0\n";
        Dataset data = LoadFromCSVString(csv, ',', true);

        // Get column by name
        Vector<Real> xVec = data.GetRealColumn("x");
        Vector<Real> yVec = data.GetRealColumn("y");

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
        demo_safe_loading();
        demo_column_operations();
        demo_dataset_summary();

        std::cout << "\n" << std::string(70, '=') << "\n";
    }
} // namespace

void Docs_Demo_DataLoader() {
    MML::docs_demos::data_loader::Run();
}
