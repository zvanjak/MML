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

#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace MML
{
namespace Persistence
{
	///////////////////////////////////////////////////////////////////////////////////////////
	///                               Text File I/O                                         ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Load matrix from text file (space-separated values)
	 * @tparam Type Matrix element type
	 * @param filename Path to input file
	 * @param outMat Matrix to populate
	 * @return SerializeResult with success/failure info
	 * 
	 * File format:
	 * - First line: rows cols
	 * - Following lines: space-separated matrix elements (row by row)
	 */
	template<class Type>
	SerializeResult LoadMatrixFromFile(const std::string& filename, Matrix<Type>& outMat)
	{
		std::ifstream file(filename);
		if (!file.is_open()) {
			return {false, SerializeError::FILE_NOT_OPENED, "Could not open file " + filename + " for reading"};
		}
		
		int rows, cols;
		file >> rows >> cols;
		if (file.fail() || rows <= 0 || cols <= 0 || rows > 100000 || cols > 100000) {
			return {false, SerializeError::INVALID_FORMAT, "Invalid dimensions (" + std::to_string(rows) + "x" + std::to_string(cols) + ") in " + filename};
		}
		outMat.Resize(rows, cols);
		
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				file >> outMat(i, j);
		
		if (file.fail()) {
			return {false, SerializeError::READ_FAILED, "Stream error while reading matrix data from " + filename};
		}
		
		file.close();
		return {true, SerializeError::OK, "Success"};
	}

	/**
	 * @brief Save matrix to text file
	 * @tparam Type Matrix element type
	 * @param mat Matrix to save
	 * @param filename Path to output file
	 * @return SerializeResult with success/failure info
	 */
	template<class Type>
	SerializeResult SaveMatrixToFile(const Matrix<Type>& mat, const std::string& filename)
	{
		std::ofstream file(filename);
		if (!file.is_open()) {
			return {false, SerializeError::FILE_NOT_OPENED, "Could not create file " + filename + " for writing"};
		}
		
		file << mat.rows() << " " << mat.cols() << std::endl;
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j)
				file << mat(i, j) << " ";
			file << std::endl;
		}
		
		if (file.fail()) {
			return {false, SerializeError::WRITE_FAILED, "Stream error while writing matrix to " + filename};
		}
		
		file.close();
		return {true, SerializeError::OK, "Success"};
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                               CSV File I/O                                          ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Load matrix from CSV file
	 * @tparam Type Matrix element type
	 * @param filename Path to CSV file
	 * @param outMat Matrix to populate
	 * @return SerializeResult with success/failure info
	 */
	template<class Type>
	SerializeResult LoadMatrixFromCSV(const std::string& filename, Matrix<Type>& outMat)
	{
		std::ifstream file(filename);
		if (!file.is_open())
			return {false, SerializeError::FILE_NOT_OPENED, "Could not open CSV file " + filename};
		
		std::vector<std::vector<Type>> data;
		std::string line;
		
		while (std::getline(file, line)) {
			std::vector<Type> row;
			for (const std::string& cell : CsvUtils::SplitLine(line, ',')) {
				std::stringstream cellStream(cell);
				Type value;
				cellStream >> value;
				row.push_back(value);
			}
			if (!row.empty())
				data.push_back(row);
		}
		
		file.close();
		if (data.empty())
			return {false, SerializeError::READ_FAILED, "No data found in CSV file " + filename};
		
		outMat = Matrix<Type>(data);
		return {true, SerializeError::OK, "Success"};
	}

	/**
	 * @brief Save matrix to CSV file
	 * @tparam Type Matrix element type
	 * @param mat Matrix to save
	 * @param filename Path to CSV file
	 * @return SerializeResult with success/failure info
	 */
	template<class Type>
	SerializeResult SaveMatrixToCSV(const Matrix<Type>& mat, const std::string& filename)
	{
		std::ofstream file(filename);
		if (!file.is_open())
			return {false, SerializeError::FILE_NOT_OPENED, "Could not create CSV file " + filename};
		
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				file << mat(i, j);
				if (j < mat.cols() - 1) file << ",";
			}
			file << "\n";
		}
		
		file.close();
		return {true, SerializeError::OK, "Success"};
	}

	template<class Type>
	SerializeResult Save(std::ostream& out, const Matrix<Type>& mat, Format format)
	{
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
					out << mat(i, j);
					if (j < mat.cols() - 1) out << ",";
				}
				out << "\n";
			}
			break;
		case Format::Binary:
			return SaveBinary(out, mat, SaveOptions{});
		case Format::JSON:
			return SaveJson(out, mat, SaveOptions{});
		case Format::Auto:
		default:
			return { false, SerializeError::INVALID_FORMAT, "Stream serialization requires an explicit format" };
		}

		if (out.fail())
			return { false, SerializeError::WRITE_FAILED, "Stream error while writing matrix" };
		return { true, SerializeError::OK, "Success" };
	}

	template<class Type>
	SerializeResult Save(std::ostream& out, const Matrix<Type>& mat, Format format, const SaveOptions& options)
	{
		const std::streamsize oldPrecision = out.precision();
		out << std::setprecision(options.precision);
		SerializeResult result = Save(out, mat, format);
		out.precision(oldPrecision);
		return result;
	}

	template<class Type>
	SerializeResult Load(std::istream& in, Matrix<Type>& outMat, Format format)
	{
		switch (format)
		{
		case Format::Text:
		{
			int rows, cols;
			in >> rows >> cols;
			if (in.fail() || rows <= 0 || cols <= 0 || rows > 100000 || cols > 100000)
				return { false, SerializeError::INVALID_FORMAT, "Invalid matrix dimensions in stream" };

			outMat.Resize(rows, cols);
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
				std::stringstream ss(line);
				std::string cell;
				while (std::getline(ss, cell, ',')) {
					std::stringstream cellStream(cell);
					Type value;
					cellStream >> value;
					row.push_back(value);
				}
				if (!row.empty())
					data.push_back(row);
			}

			if (data.empty())
				return { false, SerializeError::READ_FAILED, "No matrix data found in CSV stream" };
			outMat = Matrix<Type>(data);
			break;
		}
		case Format::Binary:
			return LoadBinary(in, outMat, LoadOptions{});
		case Format::JSON:
			return LoadJson(in, outMat, LoadOptions{});
		case Format::Auto:
		default:
			return { false, SerializeError::INVALID_FORMAT, "Stream deserialization requires an explicit format" };
		}

		if (in.fail())
			return { false, SerializeError::READ_FAILED, "Stream error while reading matrix" };
		return { true, SerializeError::OK, "Success" };
	}

	template<class Type>
	SerializeResult Load(std::istream& in, Matrix<Type>& outMat, Format format, const LoadOptions& options)
	{
		(void)options;
		return Load(in, outMat, format);
	}

	template<class Type>
	SerializeResult Save(const Matrix<Type>& mat, const std::string& filename, Format format = Format::Auto)
	{
		const Format resolvedFormat = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		switch (resolvedFormat)
		{
		case Format::Text:
			return SaveMatrixToFile(mat, filename);
		case Format::CSV:
			return SaveMatrixToCSV(mat, filename);
		case Format::Binary:
			return SaveBinary(mat, filename);
		case Format::JSON:
			return SaveJson(mat, filename, SaveOptions{});
		case Format::Auto:
		default:
			return { false, SerializeError::INVALID_FORMAT, "Could not infer serialization format from " + filename };
		}
	}

	template<class Type>
	SerializeResult Save(const Matrix<Type>& mat, const std::string& filename, const SaveOptions& options, Format format = Format::Auto)
	{
		const Format resolvedFormat = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		if (resolvedFormat == Format::JSON)
			return SaveJson(mat, filename, options);
		if (resolvedFormat == Format::Auto)
			return { false, SerializeError::INVALID_FORMAT, "Could not infer serialization format from " + filename };

		if (resolvedFormat == Format::Binary)
			return SaveBinary(mat, filename, options);

		const std::ios::openmode mode = std::ios::openmode{};
		std::ofstream file(filename, mode);
		if (!file.is_open())
			return { false, SerializeError::FILE_NOT_OPENED, "Could not create file " + filename + " for writing" };

		return Save(file, mat, resolvedFormat, options);
	}

	template<class Type>
	SerializeResult Load(const std::string& filename, Matrix<Type>& outMat, Format format = Format::Auto)
	{
		const Format resolvedFormat = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		switch (resolvedFormat)
		{
		case Format::Text:
			return LoadMatrixFromFile(filename, outMat);
		case Format::CSV:
			return LoadMatrixFromCSV(filename, outMat);
		case Format::Binary:
			return LoadBinary(filename, outMat);
		case Format::JSON:
			return LoadJson(filename, outMat, LoadOptions{});
		case Format::Auto:
		default:
			return { false, SerializeError::INVALID_FORMAT, "Could not infer serialization format from " + filename };
		}
	}

	template<class Type>
	SerializeResult Load(const std::string& filename, Matrix<Type>& outMat, const LoadOptions& options, Format format = Format::Auto)
	{
		const Format resolvedFormat = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		if (resolvedFormat == Format::JSON)
			return LoadJson(filename, outMat, options);
		return Load(filename, outMat, resolvedFormat);
	}

} // namespace Persistence
} // namespace MML

#endif // MML_MATRIX_IO_H
