///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ConsolePrinter.h                                                    ///
///  Description: Console output utilities for vectors, matrices, functions           ///
///               Formatted printing with precision and alignment control             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CONSOLE_PRINTER_H
#define MML_CONSOLE_PRINTER_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <memory>

#include "../base/Vector/Vector.h"

namespace MML {

	///////////////////////////       UTF-8 HELPERS       ///////////////////////////

	/// @brief Namespace for UTF-8 safe string operations in console printing
	namespace Utf8 {
		/// @brief Repeat a string (possibly multi-byte) a given number of times.
		///
		/// This function correctly handles UTF-8 multi-byte characters like box-drawing
		/// characters (─, ═, ━, etc.) which are 3 bytes each.
		///
		/// @param str The string to repeat (can be multi-byte UTF-8 character)
		/// @param count Number of times to repeat
		/// @return Concatenated string
		inline std::string repeatString(const std::string& str, size_t count) {
			if (str.empty() || count == 0) return "";
			std::string result;
			result.reserve(str.size() * count);
			for (size_t i = 0; i < count; ++i) {
				result += str;
			}
			return result;
		}

		/// @brief Check if a string represents exactly one display character.
		///
		/// A display character can be:
		/// - 1 byte: ASCII (0x00-0x7F)
		/// - 2 bytes: UTF-8 (0xC0-0xDF lead byte)
		/// - 3 bytes: UTF-8 (0xE0-0xEF lead byte) - most box-drawing chars
		/// - 4 bytes: UTF-8 (0xF0-0xF7 lead byte)
		///
		/// @param str The string to check
		/// @return true if the string is exactly one UTF-8 display character
		inline bool isSingleDisplayChar(const std::string& str) {
			if (str.empty()) return false;
			unsigned char lead = static_cast<unsigned char>(str[0]);
			size_t expectedLen = 1;
			if ((lead & 0x80) == 0x00) expectedLen = 1;       // ASCII
			else if ((lead & 0xE0) == 0xC0) expectedLen = 2;  // 2-byte
			else if ((lead & 0xF0) == 0xE0) expectedLen = 3;  // 3-byte (box-drawing)
			else if ((lead & 0xF8) == 0xF0) expectedLen = 4;  // 4-byte
			else return false; // Invalid lead byte
			return str.size() == expectedLen;
		}
	} // namespace Utf8

	///////////////////////////       CSV HELPERS       ///////////////////////////

	/// @brief Namespace for RFC 4180 compliant CSV utilities
	namespace Csv {
		/// @brief Check if a string needs quoting per RFC 4180
		///
		/// A field needs quoting if it contains:
		/// - Comma (,) - field delimiter
		/// - Double quote (") - must be escaped
		/// - Newline (\n or \r) - record delimiter
		///
		/// @param str The string to check
		/// @return true if the string needs to be quoted
		inline bool needsQuoting(const std::string& str) {
			return str.find_first_of(",\"\n\r") != std::string::npos;
		}

		/// @brief Escape a string for RFC 4180 compliant CSV output
		///
		/// Per RFC 4180:
		/// - Fields containing comma, quote, or newline must be quoted
		/// - Double quotes inside fields are escaped by doubling them
		/// - Example: 'He said "Hi"' becomes '"He said ""Hi"""'
		///
		/// @param str The string to escape
		/// @return The escaped string (quoted if necessary)
		inline std::string escape(const std::string& str) {
			if (!needsQuoting(str)) {
				return str;
			}

			std::string result;
			result.reserve(str.size() + 2);  // At least for the surrounding quotes
			result += '"';

			for (char c : str) {
				if (c == '"') {
					result += "\"\"";	// Escape quote by doubling
				} else {
					result += c;
				}
			}

			result += '"';
			return result;
		}

		/// @brief Parse an RFC 4180 escaped CSV field back to original value
		///
		/// This reverses the escape() operation:
		/// - Removes surrounding quotes if present
		/// - Converts doubled quotes back to single quotes
		///
		/// @param str The escaped string to unescape
		/// @return The original unescaped string
		inline std::string unescape(const std::string& str) {
			if (str.empty()) return str;

			// Check if quoted
			if (str.front() != '"' || str.back() != '"' || str.size() < 2) {
				return str;  // Not quoted, return as-is
			}

			std::string result;
			result.reserve(str.size() - 2);

			// Skip opening and closing quotes
			for (size_t i = 1; i < str.size() - 1; ++i) {
				if (str[i] == '"' && i + 1 < str.size() - 1 && str[i + 1] == '"') {
					result += '"';
					++i;	// Skip the second quote
				} else {
					result += str[i];
				}
			}

			return result;
		}
	} // namespace Csv

	///////////////////////////       MODERN API       ///////////////////////////

	// Format type enumeration
	enum class FormatType {
		General,	 // Default formatting
		Fixed,		 // std::fixed
		Scientific,	 // std::scientific
		Hexadecimal, // std::hex
		Integer		 // No decimal point
	};

	// Alignment enumeration
	enum class Alignment { Left, Right, Center };

	// Border style enumeration
	enum class BorderStyle {
		None,	  // No borders
		Simple,	  // ASCII: + - |
		Markdown, // | - |
		Rounded,  // Unicode: ╭─╮│╰─╯
		Double,	  // Unicode: ╔═╗║╚═╝
		Bold	  // Unicode: ┏━┓┃┗━┛
	};

	// Export format enumeration
	enum class ExportFormat {
		Console,  // Pretty-printed to console
		CSV,	  // Comma-separated values
		TSV,	  // Tab-separated values
		Markdown, // Markdown table
		LaTeX,	  // LaTeX tabular
		HTML	  // HTML table
	};

	///////////////////////////       RAII Stream Guard       ///////////////////////////

	// RAII Stream State Guard - preserves stream state across operations
	class StreamStateGuard {
	private:
		std::ostream& m_stream;
		std::ios_base::fmtflags m_flags;
		std::streamsize m_precision;
		std::streamsize m_width;

	public:
		explicit StreamStateGuard(std::ostream& os)
			: m_stream(os)
			, m_flags(os.flags())
			, m_precision(os.precision())
			, m_width(os.width()) {}

		~StreamStateGuard() {
			m_stream.flags(m_flags);
			m_stream.precision(m_precision);
			m_stream.width(m_width);
		}

		// Non-copyable, non-movable
		StreamStateGuard(const StreamStateGuard&) = delete;
		StreamStateGuard& operator=(const StreamStateGuard&) = delete;
		StreamStateGuard(StreamStateGuard&&) = delete;
		StreamStateGuard& operator=(StreamStateGuard&&) = delete;
	};

	///////////////////////////       Column Format       ///////////////////////////

	// Column Format Class with Builder Pattern
	class ColumnFormat {
	private:
		std::string m_name;
		int m_width;
		int m_precision;
		FormatType m_formatType;
		Alignment m_alignment;
		int m_calculatedWidth; // -1 means use m_width

	public:
		static constexpr int AUTO_WIDTH = -1;

		// Constructor
		explicit ColumnFormat(std::string name)
			: m_name(std::move(name))
			, m_width(12)
			, m_precision(6)
			, m_formatType(FormatType::General)
			, m_alignment(Alignment::Right)
			, m_calculatedWidth(-1) {}

		// Builder pattern methods
		ColumnFormat& width(int w) {
			m_width = w;
			return *this;
		}

		ColumnFormat& precision(int p) {
			m_precision = p;
			return *this;
		}

		ColumnFormat& format(FormatType fmt) {
			m_formatType = fmt;
			return *this;
		}

		ColumnFormat& align(Alignment a) {
			m_alignment = a;
			return *this;
		}

		ColumnFormat& autoWidth() {
			m_width = AUTO_WIDTH;
			return *this;
		}

		// Accessors
		const std::string& name() const { return m_name; }
		int width() const { return m_calculatedWidth >= 0 ? m_calculatedWidth : m_width; }
		int rawWidth() const { return m_width; }
		int precision() const { return m_precision; }
		FormatType formatType() const { return m_formatType; }
		Alignment alignment() const { return m_alignment; }

		// Internal - for auto-width calculation
		void setCalculatedWidth(int w) { m_calculatedWidth = w; }

		// Helper to format a value according to this column's settings
		template<typename T>
		std::string formatValue(const T& value) const {
			std::ostringstream oss;

			// Apply format type
			switch (m_formatType) {
			case FormatType::Fixed:
				oss << std::fixed;
				break;
			case FormatType::Scientific:
				oss << std::scientific;
				break;
			case FormatType::Hexadecimal:
				oss << std::hex;
				break;
			case FormatType::Integer:
				// No decimal formatting
				break;
			case FormatType::General:
			default:
				// Use default formatting
				break;
			}

			// Apply precision
			if (m_formatType != FormatType::Integer && m_formatType != FormatType::Hexadecimal) {
				oss << std::setprecision(m_precision);
			}

			// Output value
			oss << value;
			return oss.str();
		}

		// Helper to format with alignment and width
		std::string formatAligned(const std::string& content) const {
			int w = width();
			if (content.length() >= static_cast<size_t>(w)) {
				return content;
			}

			int padding = w - static_cast<int>(content.length());

			switch (m_alignment) {
			case Alignment::Left:
				return content + std::string(padding, ' ');
			case Alignment::Right:
				return std::string(padding, ' ') + content;
			case Alignment::Center: {
				int leftPad = padding / 2;
				int rightPad = padding - leftPad;
				return std::string(leftPad, ' ') + content + std::string(rightPad, ' ');
			}
			default:
				return content;
			}
		}
	};

	///////////////////////////       Table Style       ///////////////////////////

	// Table Style Class
	class TableStyle {
	private:
		BorderStyle m_border;
		bool m_showHeader;
		bool m_showRowSeparators;
		bool m_compactMode;

	public:
		TableStyle()
			: m_border(BorderStyle::Simple)
			, m_showHeader(true)
			, m_showRowSeparators(false)
			, m_compactMode(false) {}

		// Builder pattern methods
		TableStyle& border(BorderStyle style) {
			m_border = style;
			return *this;
		}

		TableStyle& header(bool show) {
			m_showHeader = show;
			return *this;
		}

		TableStyle& rowSeparators(bool show) {
			m_showRowSeparators = show;
			return *this;
		}

		TableStyle& compact(bool enable) {
			m_compactMode = enable;
			return *this;
		}

		// Accessors
		BorderStyle borderStyle() const { return m_border; }
		bool showHeader() const { return m_showHeader; }
		bool showRowSeparators() const { return m_showRowSeparators; }
		bool compactMode() const { return m_compactMode; }

		// Border characters
		struct BorderChars {
			std::string topLeft, topRight, bottomLeft, bottomRight;
			std::string horizontal, vertical;
			std::string leftT, rightT, topT, bottomT, cross;
		};

		BorderChars getChars() const {
			switch (m_border) {
			case BorderStyle::None:
				return {"", "", "", "", "", "", "", "", "", "", ""};
			case BorderStyle::Simple:
				return {"+", "+", "+", "+", "-", "|", "+", "+", "+", "+", "+"};
			case BorderStyle::Markdown:
				return {"|", "|", "|", "|", "-", "|", "|", "|", "|", "|", "|"};
			case BorderStyle::Rounded:
				return {"╭", "╮", "╰", "╯", "─", "│", "├", "┤", "┬", "┴", "┼"};
			case BorderStyle::Double:
				return {"╔", "╗", "╚", "╝", "═", "║", "╠", "╣", "╦", "╩", "╬"};
			case BorderStyle::Bold:
				return {"┏", "┓", "┗", "┛", "━", "┃", "┣", "┫", "┳", "┻", "╋"};
			default:
				return {"+", "+", "+", "+", "-", "|", "+", "+", "+", "+", "+"};
			}
		}
	};

	///////////////////////////       Table Printer       ///////////////////////////

	// Modern TablePrinter with full features
	template<typename RowTag, typename CellValue>
	class TablePrinter {
	private:
		std::string m_tagColumnName;
		ColumnFormat m_tagFormat;
		std::vector<ColumnFormat> m_columnFormats;
		std::vector<RowTag> m_rowTags;
		std::vector<std::vector<CellValue>> m_data;
		TableStyle m_style;

	public:
		// Construction with simple column names
		TablePrinter(std::string tagName, std::vector<std::string> columnNames)
			: m_tagColumnName(std::move(tagName))
			, m_tagFormat(m_tagColumnName) {
			for (const auto& name : columnNames) {
				m_columnFormats.emplace_back(name);
			}
		}

		// Construction with explicit tag format
		TablePrinter(ColumnFormat tagFormat, std::vector<ColumnFormat> columnFormats)
			: m_tagColumnName(tagFormat.name())
			, m_tagFormat(std::move(tagFormat))
			, m_columnFormats(std::move(columnFormats)) {}

		// Builder pattern for configuration
		TablePrinter& style(const TableStyle& s) {
			m_style = s;
			return *this;
		}

		TablePrinter& tagFormat(const ColumnFormat& fmt) {
			m_tagFormat = fmt;
			return *this;
		}

		TablePrinter& columnFormat(size_t col, const ColumnFormat& fmt) {
			if (col < m_columnFormats.size()) {
				m_columnFormats[col] = fmt;
			}
			return *this;
		}

		// Data manipulation
		void addRow(RowTag tag, std::vector<CellValue> values) {
			if (values.size() != m_columnFormats.size()) {
				throw VectorDimensionError("Number of values does not match number of columns", 
				                           static_cast<int>(m_columnFormats.size()), static_cast<int>(values.size()));
			}
			m_rowTags.push_back(std::move(tag));
			m_data.push_back(std::move(values));
		}

		void clear() {
			m_rowTags.clear();
			m_data.clear();
		}

		void reserve(size_t rows) {
			m_rowTags.reserve(rows);
			m_data.reserve(rows);
		}

		// Query
		size_t rowCount() const { return m_data.size(); }
		size_t columnCount() const { return m_columnFormats.size(); }

		// Output to console (default)
		void print(std::ostream& os = std::cout) const { exportTo(os, ExportFormat::Console); }

		// Print to file
		void printToFile(const std::string& filename) const {
			std::ofstream file(filename);
			if (!file) {
				throw FileIOError("Cannot open file: " + filename, filename);
			}
			print(file);
		}

		// Convert to string
		std::string toString() const {
			std::ostringstream oss;
			print(oss);
			return oss.str();
		}

		// Export to different formats
		void exportTo(std::ostream& os, ExportFormat format) const {
			StreamStateGuard guard(os);

			switch (format) {
			case ExportFormat::Console:
				exportConsole(os);
				break;
			case ExportFormat::CSV:
				exportCSV(os);
				break;
			case ExportFormat::TSV:
				exportTSV(os);
				break;
			case ExportFormat::Markdown:
				exportMarkdown(os);
				break;
			case ExportFormat::LaTeX:
				exportLaTeX(os);
				break;
			case ExportFormat::HTML:
				exportHTML(os);
				break;
			}
		}

		void exportToFile(const std::string& filename, ExportFormat format) const {
			std::ofstream file(filename);
			if (!file) {
				throw FileIOError("Cannot open file: " + filename, filename);
			}
			exportTo(file, format);
		}

	private:
		// Console export with borders and styling
		void exportConsole(std::ostream& os) const {
			auto bc = m_style.getChars();
			bool hasBorders = m_style.borderStyle() != BorderStyle::None;

			// Calculate auto-widths if needed
			calculateAutoWidths();

			// Print top border
			if (hasBorders && m_style.showHeader()) {
				printHorizontalLine(os, bc.topLeft, bc.horizontal, bc.topT, bc.topRight);
			}

			// Print header
			if (m_style.showHeader()) {
				if (hasBorders)
					os << bc.vertical << " ";
				os << m_tagFormat.formatAligned(m_tagFormat.name());

				for (const auto& fmt : m_columnFormats) {
					if (hasBorders)
						os << " " << bc.vertical << " ";
					else
						os << "  ";
					os << fmt.formatAligned(fmt.name());
				}

				if (hasBorders)
					os << " " << bc.vertical;
				os << "\n";

				// Header separator
				if (hasBorders) {
					printHorizontalLine(os, bc.leftT, bc.horizontal, bc.cross, bc.rightT);
				}
			}

			// Print data rows
			for (size_t row = 0; row < m_data.size(); ++row) {
				if (hasBorders)
					os << bc.vertical << " ";
				os << m_tagFormat.formatAligned(m_tagFormat.formatValue(m_rowTags[row]));

				for (size_t col = 0; col < m_columnFormats.size(); ++col) {
					if (hasBorders)
						os << " " << bc.vertical << " ";
					else
						os << "  ";
					os << m_columnFormats[col].formatAligned(m_columnFormats[col].formatValue(m_data[row][col]));
				}

				if (hasBorders)
					os << " " << bc.vertical;
				os << "\n";

				// Row separator
				if (hasBorders && m_style.showRowSeparators() && row < m_data.size() - 1) {
					printHorizontalLine(os, bc.leftT, bc.horizontal, bc.cross, bc.rightT);
				}
			}

			// Print bottom border
			if (hasBorders) {
				printHorizontalLine(os, bc.bottomLeft, bc.horizontal, bc.bottomT, bc.bottomRight);
			}
		}

		// CSV export (RFC 4180 compliant)
		void exportCSV(std::ostream& os) const {
			// Header - escape column names in case they contain special characters
			os << Csv::escape(m_tagFormat.name());
			for (const auto& fmt : m_columnFormats) {
				os << "," << Csv::escape(fmt.name());
			}
			os << "\n";

			// Data - escape row tags and values
			for (size_t row = 0; row < m_data.size(); ++row) {
				// Row tag (convert to string and escape)
				std::ostringstream tagStr;
				tagStr << m_rowTags[row];
				os << Csv::escape(tagStr.str());

				// Values (convert to string and escape)
				for (const auto& value : m_data[row]) {
					std::ostringstream valStr;
					valStr << value;
					os << "," << Csv::escape(valStr.str());
				}
				os << "\n";
			}
		}

		// TSV export
		void exportTSV(std::ostream& os) const {
			// Header
			os << m_tagFormat.name();
			for (const auto& fmt : m_columnFormats) {
				os << "\t" << fmt.name();
			}
			os << "\n";

			// Data
			for (size_t row = 0; row < m_data.size(); ++row) {
				os << m_rowTags[row];
				for (const auto& value : m_data[row]) {
					os << "\t" << value;
				}
				os << "\n";
			}
		}

		// Markdown export
		void exportMarkdown(std::ostream& os) const {
			// Header
			os << "| " << m_tagFormat.name();
			for (const auto& fmt : m_columnFormats) {
				os << " | " << fmt.name();
			}
			os << " |\n";

			// Separator
			os << "|" << std::string(m_tagFormat.name().length() + 2, '-');
			for (const auto& fmt : m_columnFormats) {
				os << "|" << std::string(fmt.name().length() + 2, '-');
			}
			os << "|\n";

			// Data
			for (size_t row = 0; row < m_data.size(); ++row) {
				os << "| " << m_rowTags[row];
				for (const auto& value : m_data[row]) {
					os << " | " << value;
				}
				os << " |\n";
			}
		}

		// LaTeX export
		void exportLaTeX(std::ostream& os) const {
			// Begin tabular
			os << "\\begin{tabular}{|";
			for (size_t i = 0; i <= m_columnFormats.size(); ++i) {
				os << "c|";
			}
			os << "}\n\\hline\n";

			// Header
			os << m_tagFormat.name();
			for (const auto& fmt : m_columnFormats) {
				os << " & " << fmt.name();
			}
			os << " \\\\\n\\hline\n";

			// Data
			for (size_t row = 0; row < m_data.size(); ++row) {
				os << m_rowTags[row];
				for (const auto& value : m_data[row]) {
					os << " & " << value;
				}
				os << " \\\\\n";
			}

			// End tabular
			os << "\\hline\n\\end{tabular}\n";
		}

		// HTML export
		void exportHTML(std::ostream& os) const {
			os << "<table border=\"1\">\n";

			// Header
			os << "  <thead>\n    <tr>\n";
			os << "      <th>" << m_tagFormat.name() << "</th>\n";
			for (const auto& fmt : m_columnFormats) {
				os << "      <th>" << fmt.name() << "</th>\n";
			}
			os << "    </tr>\n  </thead>\n";

			// Body
			os << "  <tbody>\n";
			for (size_t row = 0; row < m_data.size(); ++row) {
				os << "    <tr>\n";
				os << "      <td>" << m_rowTags[row] << "</td>\n";
				for (const auto& value : m_data[row]) {
					os << "      <td>" << value << "</td>\n";
				}
				os << "    </tr>\n";
			}
			os << "  </tbody>\n";

			os << "</table>\n";
		}

		// Helper: print horizontal line for borders (UTF-8 safe)
		void printHorizontalLine(std::ostream& os, const std::string& left, const std::string& horiz, const std::string& cross,
								 const std::string& right) const {
			os << left;
			os << Utf8::repeatString(horiz, static_cast<size_t>(m_tagFormat.width() + 2));

			for (size_t i = 0; i < m_columnFormats.size(); ++i) {
				os << cross;
				os << Utf8::repeatString(horiz, static_cast<size_t>(m_columnFormats[i].width() + 2));
			}

			os << right << "\n";
		}

		// Helper: calculate auto-widths
		void calculateAutoWidths() const {
			// Mutable to allow calculation during const operations
			auto& tagFmt = const_cast<ColumnFormat&>(m_tagFormat);
			auto& colFmts = const_cast<std::vector<ColumnFormat>&>(m_columnFormats);

			// Tag column
			if (tagFmt.rawWidth() == ColumnFormat::AUTO_WIDTH) {
				int maxWidth = static_cast<int>(tagFmt.name().length());
				for (const auto& tag : m_rowTags) {
					std::string formatted = tagFmt.formatValue(tag);
					maxWidth = std::max(maxWidth, static_cast<int>(formatted.length()));
				}
				tagFmt.setCalculatedWidth(maxWidth);
			}

			// Value columns
			for (size_t col = 0; col < colFmts.size(); ++col) {
				if (colFmts[col].rawWidth() == ColumnFormat::AUTO_WIDTH) {
					int maxWidth = static_cast<int>(colFmts[col].name().length());
					for (const auto& row : m_data) {
						std::string formatted = colFmts[col].formatValue(row[col]);
						maxWidth = std::max(maxWidth, static_cast<int>(formatted.length()));
					}
					colFmts[col].setCalculatedWidth(maxWidth);
				}
			}
		}
	};

	///////////////////////////       Vector Table Printer       ///////////////////////////

	// Specialized printer for vectors
	class VectorTablePrinter {
	private:
		std::vector<ColumnFormat> m_formats;
		std::vector<const Vector<Real>*> m_vectors;
		TableStyle m_style;
		size_t m_maxLength;

	public:
		VectorTablePrinter()
			: m_maxLength(0) {}

		// Constructor that accepts formats and vectors
		VectorTablePrinter(std::vector<ColumnFormat> formats, std::vector<Vector<Real>*> vectors)
			: m_maxLength(0) {
			for (size_t i = 0; i < vectors.size(); ++i) {
				if (i < formats.size()) {
					addVector(formats[i], *vectors[i]);
				} else {
					addVector("Column" + std::to_string(i), *vectors[i]);
				}
			}
		}

		void addVector(const std::string& name, const Vector<Real>& vec) {
			m_formats.emplace_back(name);
			m_vectors.push_back(&vec);
			m_maxLength = std::max(m_maxLength, static_cast<size_t>(vec.size()));
		}

		void addVector(const ColumnFormat& format, const Vector<Real>& vec) {
			m_formats.push_back(format);
			m_vectors.push_back(&vec);
			m_maxLength = std::max(m_maxLength, static_cast<size_t>(vec.size()));
		}

		VectorTablePrinter& style(const TableStyle& s) {
			m_style = s;
			return *this;
		}

		void print(std::ostream& os = std::cout) const {
			StreamStateGuard guard(os);

			auto bc = m_style.getChars();
			bool hasBorders = m_style.borderStyle() != BorderStyle::None;

			// Print top border
			if (hasBorders && m_style.showHeader()) {
				printHorizontalLine(os, bc.topLeft, bc.horizontal, bc.topT, bc.topRight);
			}

			// Print header
			if (m_style.showHeader()) {
				if (hasBorders)
					os << bc.vertical << " ";

				for (size_t i = 0; i < m_formats.size(); ++i) {
					if (i > 0) {
						if (hasBorders)
							os << " " << bc.vertical << " ";
						else
							os << "  ";
					}
					os << m_formats[i].formatAligned(m_formats[i].name());
				}

				if (hasBorders)
					os << " " << bc.vertical;
				os << "\n";

				// Header separator
				if (hasBorders) {
					printHorizontalLine(os, bc.leftT, bc.horizontal, bc.cross, bc.rightT);
				}
			}

			// Print data rows
			for (size_t row = 0; row < m_maxLength; ++row) {
				if (hasBorders)
					os << bc.vertical << " ";

				for (size_t col = 0; col < m_vectors.size(); ++col) {
					if (col > 0) {
						if (hasBorders)
							os << " " << bc.vertical << " ";
						else
							os << "  ";
					}

					if (row < static_cast<size_t>(m_vectors[col]->size())) {
						Real value = (*m_vectors[col])[row];
						os << m_formats[col].formatAligned(m_formats[col].formatValue(value));
					} else {
						os << m_formats[col].formatAligned("");
					}
				}

				if (hasBorders)
					os << " " << bc.vertical;
				os << "\n";
			}

			// Print bottom border
			if (hasBorders) {
				printHorizontalLine(os, bc.bottomLeft, bc.horizontal, bc.bottomT, bc.bottomRight);
			}
		}

		std::string toString() const {
			std::ostringstream oss;
			print(oss);
			return oss.str();
		}

	private:
		// Helper: print horizontal line for borders (UTF-8 safe)
		void printHorizontalLine(std::ostream& os, const std::string& left, const std::string& horiz, const std::string& cross,
								 const std::string& right) const {
			os << left;

			for (size_t i = 0; i < m_formats.size(); ++i) {
				if (i > 0)
					os << cross;
				os << Utf8::repeatString(horiz, static_cast<size_t>(m_formats[i].width() + 2));
			}

			os << right << "\n";
		}
	};

} // namespace MML

#endif
