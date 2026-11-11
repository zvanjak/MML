///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixIO.h                                                          ///
///  Description: Matrix text/CSV I/O and persistence format dispatch                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_IO_H
#define MML_MATRIX_IO_H

#include <mml/MMLBase.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/tools/CsvUtils.h>
#include <mml/tools/persistence/MatrixBinary.h>
#include <mml/tools/persistence/MatrixJSON.h>
#include <mml/tools/TextIOFormat.h>

#include <fstream>
#include <locale>
#include <sstream>
#include <vector>
#include <string>

namespace MML
{
namespace Persistence
{
	template<class Type>
	IOResult Save(std::ostream& out, const Matrix<Type>& mat, Format format, const SaveOptions& options);

	///////////////////////////////////////////////////////////////////////////////////////////
	///                               Text File I/O                                         ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Load matrix from text file (space-separated values)
	 * @tparam Type Matrix element type
	 * @param filename Path to input file
	 * @param outMat Matrix to populate
	 * @return IOResult with success/failure info
	 * 
	 * File format:
	 * - First line: rows cols
	 * - Following lines: space-separated matrix elements (row by row)
	 */
	template<class Type>
	IOResult LoadMatrixFromFile(const std::string& filename, Matrix<Type>& outMat) {
		std::ifstream file(filename);
		TextIO::ApplyClassicFormat(file);
		if (!file.is_open()) {
			return IOResult{false, IOError::FILE_NOT_OPENED, "Could not open file " + filename + " for reading", filename}.WithFile(filename);
		}

		int rows = 0, cols = 0;
		file >> rows >> cols;
		if (file.fail() || rows <= 0 || cols <= 0 || rows > 100000 || cols > 100000) {
			return InputResult(file, SerializeFailure(IOError::INVALID_FORMAT, "Invalid matrix dimensions"), filename);
		}
		outMat.resize(rows, cols);
		
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				file >> outMat(i, j);
		
		if (file.fail()) {
			return InputResult(file, SerializeFailure(IOError::READ_FAILED, "Stream error while reading matrix data"), filename);
		}
		
		file.close();
		return (IOResult::Success(filename)).WithFile(filename);
	}

	/**
	 * @brief Save matrix to text file
	 * @tparam Type Matrix element type
	 * @param mat Matrix to save
	 * @param filename Path to output file
	 * @return IOResult with success/failure info
	 */
	template<class Type>
	IOResult SaveMatrixToFile(const Matrix<Type>& mat, const std::string& filename) {
		const auto validation = Detail::ValidateSaveFilename(filename);
		if (!validation.success) return validation;
		std::ofstream file(filename);
		TextIO::ApplyClassicFormat(file);
		if (!file.is_open()) {
			return IOResult{false, IOError::FILE_NOT_OPENED, "Could not create file " + filename + " for writing", filename}.WithFile(filename);
		}
		
		return FinishOutputFile(file, filename, Save(file, mat, Format::Text, SaveOptions{}));
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                               CSV File I/O                                          ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Load matrix from CSV file
	 * @tparam Type Matrix element type
	 * @param filename Path to CSV file
	 * @param outMat Matrix to populate
	 * @return IOResult with success/failure info
	 */
	template<class Type>
	IOResult LoadMatrixFromCSV(const std::string& filename, Matrix<Type>& outMat) {
		std::ifstream file(filename);
		TextIO::ApplyClassicFormat(file);
		if (!file.is_open())
			return IOResult{false, IOError::FILE_NOT_OPENED, "Could not open CSV file " + filename, filename}.WithFile(filename);

		std::vector<std::vector<Type>> data;
		std::string line;
		std::int64_t lineNumber = 0;

		while (std::getline(file, line)) {
			++lineNumber;
			std::vector<Type> row;
			for (const std::string& cell : CsvUtils::SplitLine(line, ',')) {
				std::stringstream cellStream(cell);
				cellStream.imbue(std::locale::classic());
				Type value;
				cellStream >> value;
				if (cellStream.fail())
					return SerializeFailure(IOError::INVALID_FORMAT, "Unparsable CSV cell: '" + cell + "'")
						.WithFile(filename)
						.At(lineNumber, IOOffsetKind::Line);
				row.push_back(value);
			}
			if (!row.empty() && !data.empty() && row.size() != data.front().size())
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Inconsistent CSV matrix row width")
					.WithFile(filename)
					.At(lineNumber, IOOffsetKind::Line);
			if (!row.empty())
				data.push_back(row);
		}
		if (file.bad() || !file.eof())
			return InputResult(file, SerializeFailure(IOError::READ_FAILED, "CSV matrix stream read failed"), filename);

		file.close();
		if (data.empty())
			return SerializeFailure(IOError::READ_FAILED, "No matrix data found").WithFile(filename).At(1, IOOffsetKind::Line);

		outMat = Matrix<Type>(data);
		return (IOResult::Success(filename)).WithFile(filename);
	}

	/**
	 * @brief Save matrix to CSV file
	 * @tparam Type Matrix element type
	 * @param mat Matrix to save
	 * @param filename Path to CSV file
	 * @return IOResult with success/failure info
	 */
	template<class Type>
	IOResult SaveMatrixToCSV(const Matrix<Type>& mat, const std::string& filename) {
		const auto validation = Detail::ValidateSaveFilename(filename);
		if (!validation.success) return validation;
		std::ofstream file(filename);
		TextIO::ApplyClassicFormat(file);
		if (!file.is_open())
			return IOResult{false, IOError::FILE_NOT_OPENED, "Could not create CSV file " + filename, filename}.WithFile(filename);

		return FinishOutputFile(file, filename, Save(file, mat, Format::CSV, SaveOptions{}));
	}

	template<class Type>
	IOResult Save(std::ostream& out, const Matrix<Type>& mat, Format format) {
		return Save(out, mat, format, SaveOptions{});
	}

	template<class Type>
	IOResult Save(std::ostream& out, const Matrix<Type>& mat, Format format, const SaveOptions& options) {
		if (format == Format::Binary) {
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown)
				return SaveBinary(out, mat, options);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Unsupported Matrix binary scalar");
		}
		if (format == Format::JSON) {
			if constexpr (std::is_floating_point_v<Type>)
				return SaveJson(out, mat, options);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Matrix JSON requires real floating-point scalars");
		}
		if (format != Format::Text && format != Format::CSV)
			return SerializeFailure(IOError::INVALID_FORMAT, "Stream serialization requires an explicit format");
		if (options.precision < 0)
			return SerializeFailure(IOError::INVALID_PARAMETERS, "Text precision must be zero (automatic) or positive");
		const int precision = options.precision == 0 ? TextIO::RoundTripPrecision<Type>() : options.precision;
		TextIO::ScopedTextFormat guard(out, precision);
		switch (format)
		{
		case Format::Text:
			out << mat.rows() << " " << mat.cols() << std::endl;
			for (int i = 0; i < mat.rows(); ++i) {
				for (int j = 0; j < mat.cols(); ++j)
					out << mat(i, j) << " ";
				out << std::endl;
			}
			break;
		case Format::CSV:
			for (int i = 0; i < mat.rows(); ++i) {
				for (int j = 0; j < mat.cols(); ++j) {
					if constexpr (is_std_complex<Type>::value) {
						std::ostringstream cell;
						TextIO::ApplyClassicFormat(cell, precision);
						cell << mat(i, j);
						out << CsvUtils::Escape(cell.str());
					} else {
						out << mat(i, j);
					}
					if (j < mat.cols() - 1) out << ",";
				}
				out << "\n";
			}
			break;
		default:
			return {false, IOError::INVALID_FORMAT, "Stream serialization requires an explicit format"};
		}

		if (out.fail())
			return {false, IOError::WRITE_FAILED, "Stream error while writing matrix"};
		return {true, IOError::OK, "Success"};
	}

	template<class Type>
	IOResult Load(std::istream& in, Matrix<Type>& outMat, Format format) {
		if (format == Format::Binary) {
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown)
				return LoadBinary(in, outMat, LoadOptions{});
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Unsupported Matrix binary scalar");
		}
		if (format == Format::JSON) {
			if constexpr (std::is_floating_point_v<Type>)
				return LoadJson(in, outMat, LoadOptions{});
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Matrix JSON requires real floating-point scalars");
		}
		TextIO::ScopedTextFormat guard(in);
		switch (format)
		{
		case Format::Text:
		{
			int rows, cols;
			in >> rows >> cols;
			if (in.fail() || rows <= 0 || cols <= 0 || rows > 100000 || cols > 100000)
				return {false, IOError::INVALID_FORMAT, "Invalid matrix dimensions in stream"};

			outMat.resize(rows, cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					in >> outMat(i, j);
			break;
		}
		case Format::CSV:
		{
			std::vector<std::vector<Type>> data;
			std::string line;
			while (std::getline(in, line)) {
				std::vector<Type> row;
				for (const auto& cell : CsvUtils::SplitLine(line, ',')) {
					std::stringstream cellStream(cell);
					cellStream.imbue(std::locale::classic());
					Type value;
					cellStream >> value;
					if (cellStream.fail())
						return {false, IOError::INVALID_FORMAT, "Unparsable CSV cell: '" + cell + "'"};
					row.push_back(value);
				}
				if (!row.empty() && !data.empty() && row.size() != data.front().size())
					return SerializeFailure(IOError::SCHEMA_MISMATCH, "Inconsistent CSV matrix row width");
				if (!row.empty())
					data.push_back(row);
			}

			if (in.bad() || !in.eof())
				return SerializeFailure(IOError::READ_FAILED, "CSV matrix stream read failed");
			if (data.empty())
				return {false, IOError::READ_FAILED, "No matrix data found in CSV stream"};
			outMat = Matrix<Type>(data);
			in.clear(in.rdstate() & ~std::ios::failbit);
			break;
		}
		case Format::Auto:
		default:
			return {false, IOError::INVALID_FORMAT, "Stream deserialization requires an explicit format"};
		}

		if (in.fail())
			return {false, IOError::READ_FAILED, "Stream error while reading matrix"};
		return {true, IOError::OK, "Success"};
	}

	template<class Type>
	IOResult Load(std::istream& in, Matrix<Type>& outMat, Format format, const LoadOptions& options) {
		if (format == Format::Binary) {
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown)
				return LoadBinary(in, outMat, options);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Unsupported Matrix binary scalar");
		}
		if (format == Format::JSON) {
			if constexpr (std::is_floating_point_v<Type>)
				return LoadJson(in, outMat, options);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Matrix JSON requires real floating-point scalars");
		}
		return Load(in, outMat, format);
	}

	template<class Type>
	IOResult Save(const Matrix<Type>& mat, const std::string& filename, Format format = Format::Auto) {
		const auto validation = Detail::ValidateSaveFilename(filename);
		if (!validation.success) return validation;
		const Format resolvedFormat = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		switch (resolvedFormat)
		{
		case Format::Text:
			return (SaveMatrixToFile(mat, filename)).WithFile(filename);
		case Format::CSV:
			return (SaveMatrixToCSV(mat, filename)).WithFile(filename);
		case Format::Binary:
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown)
				return (SaveBinary(mat, filename)).WithFile(filename);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Unsupported Matrix binary scalar").WithFile(filename);
		case Format::JSON:
			if constexpr (std::is_floating_point_v<Type>)
				return (SaveJson(mat, filename, SaveOptions{})).WithFile(filename);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Matrix JSON requires real floating-point scalars").WithFile(filename);
		case Format::Auto:
		default:
			return IOResult{false, IOError::INVALID_FORMAT, "Could not infer serialization format from " + filename}.WithFile(filename);
		}
	}

	template<class Type>
	IOResult Save(const Matrix<Type>& mat, const std::string& filename, const SaveOptions& options, Format format = Format::Auto) {
		const auto validation = Detail::ValidateSaveFilename(filename);
		if (!validation.success) return validation;
		const Format resolvedFormat = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		if (resolvedFormat == Format::JSON) {
			if constexpr (std::is_floating_point_v<Type>)
				return (SaveJson(mat, filename, options)).WithFile(filename);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Matrix JSON requires real floating-point scalars").WithFile(filename);
		}
		if (resolvedFormat == Format::Auto)
			return IOResult{false, IOError::INVALID_FORMAT, "Could not infer serialization format from " + filename}.WithFile(filename);

		if (resolvedFormat == Format::Binary) {
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown)
				return (SaveBinary(mat, filename, options)).WithFile(filename);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Unsupported Matrix binary scalar").WithFile(filename);
		}

		if (resolvedFormat != Format::Text && resolvedFormat != Format::CSV)
			return SerializeFailure(IOError::INVALID_FORMAT, "Matrix persistence requires Text, CSV, JSON or Binary format").WithFile(filename);
		if (options.precision < 0)
			return SerializeFailure(IOError::INVALID_PARAMETERS, "Text precision must be zero (automatic) or positive").WithFile(filename);
		const std::ios::openmode mode = std::ios::openmode{};
		std::ofstream file(filename, mode);
		TextIO::ApplyClassicFormat(file);
		if (!file.is_open())
			return IOResult{false, IOError::FILE_NOT_OPENED, "Could not create file " + filename + " for writing", filename}.WithFile(filename);

		return FinishOutputFile(file, filename, Save(file, mat, resolvedFormat, options));
	}

	template<class Type>
	IOResult Load(const std::string& filename, Matrix<Type>& outMat, Format format = Format::Auto) {
		const Format resolvedFormat = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		switch (resolvedFormat)
		{
		case Format::Text:
			return (LoadMatrixFromFile(filename, outMat)).WithFile(filename);
		case Format::CSV:
			return (LoadMatrixFromCSV(filename, outMat)).WithFile(filename);
		case Format::Binary:
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown)
				return (LoadBinary(filename, outMat)).WithFile(filename);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Unsupported Matrix binary scalar").WithFile(filename);
		case Format::JSON:
			if constexpr (std::is_floating_point_v<Type>)
				return (LoadJson(filename, outMat, LoadOptions{})).WithFile(filename);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Matrix JSON requires real floating-point scalars").WithFile(filename);
		case Format::Auto:
		default:
			return IOResult{false, IOError::INVALID_FORMAT, "Could not infer serialization format from " + filename}.WithFile(filename);
		}
	}

	template<class Type>
	IOResult Load(const std::string& filename, Matrix<Type>& outMat, const LoadOptions& options, Format format = Format::Auto) {
		const Format resolvedFormat = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		if (resolvedFormat == Format::JSON) {
			if constexpr (std::is_floating_point_v<Type>)
				return (LoadJson(filename, outMat, options)).WithFile(filename);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Matrix JSON requires real floating-point scalars").WithFile(filename);
		}
		if (resolvedFormat == Format::Binary) {
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown)
				return (LoadBinary(filename, outMat, options)).WithFile(filename);
			else
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Unsupported Matrix binary scalar").WithFile(filename);
		}
		return (Load(filename, outMat, resolvedFormat)).WithFile(filename);
	}

} // namespace Persistence
} // namespace MML

#endif // MML_MATRIX_IO_H
