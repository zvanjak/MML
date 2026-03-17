///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixPrintFormat.h                                                 ///
///  Description: Formatting options for matrix console output                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_PRINT_FORMAT_H
#define MML_MATRIX_PRINT_FORMAT_H

#include <string>

namespace MML
{
  // Matrix printing format options
  struct MatrixPrintFormat {
    int width = 10;               // Field width for each element
    int precision = 3;            // Number of decimal places
    bool scientific = false;      // Use scientific notation (e.g., 1.23e+05)
    bool fixed = true;            // Use fixed-point notation
    bool showHeader = true;       // Show "Rows: N Cols: M" header
    bool showBrackets = true;     // Show [ ] brackets around rows
    bool compactMode = false;     // Single line for small matrices
    std::string delimiter = ", "; // Delimiter between elements

    // Predefined formats
    static MatrixPrintFormat Default() {
      return MatrixPrintFormat{10, 3, false, true, true, true, false, ", "};
    }

    static MatrixPrintFormat Compact() {
      return MatrixPrintFormat{8, 2, false, true, false, true, true, ", "};
    }

    static MatrixPrintFormat Scientific() {
      return MatrixPrintFormat{12, 6, true, false, true, true, false, ", "};
    }

    static MatrixPrintFormat HighPrecision() {
      return MatrixPrintFormat{15, 10, false, true, true, true, false, ", "};
    }

    static MatrixPrintFormat NoDelimiter() {
      return MatrixPrintFormat{10, 3, false, true, true, true, false, " "};
    }
  };

} // namespace MML

#endif // MML_MATRIX_PRINT_FORMAT_H
