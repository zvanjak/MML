///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixIO.h                                                          ///
///  Description: Matrix file I/O operations - text, CSV, and binary file support    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_IO_H
#define MML_MATRIX_IO_H

#include "mml/MMLBase.h"
#include "mml/base/Matrix/Matrix.h"
#include "mml/tools/serializer/SerializerBase.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace MML
{
namespace Serializer
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

	///////////////////////////////////////////////////////////////////////////////////////////
	///                               Binary File I/O                                       ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Save matrix to binary file (compact, exact precision)
	 * @tparam Type Matrix element type
	 * @param mat Matrix to save
	 * @param filename Path to output file
	 * @return SerializeResult with success/failure info
	 * 
	 * Binary format:
	 * - 4 bytes: magic number 0x4D4D4C4D ("MMLM" for MML Matrix)
	 * - 4 bytes: version (currently 1)
	 * - 4 bytes: rows (int32)
	 * - 4 bytes: cols (int32)
	 * - 4 bytes: element size in bytes (sizeof(Type))
	 * - 4 bytes: reserved (0)
	 * - rows*cols*sizeof(Type) bytes: data in row-major order
	 */
	template<class Type>
	SerializeResult SaveMatrixToBinary(const Matrix<Type>& mat, const std::string& filename)
	{
		std::ofstream file(filename, std::ios::binary);
		if (!file.is_open()) {
			return {false, SerializeError::FILE_NOT_OPENED, "Could not create binary file " + filename};
		}
		
		// Write header
		const uint32_t magic = BinaryFormat::MAGIC_MATRIX;
		const uint32_t version = BinaryFormat::VERSION_MATRIX;
		const uint32_t rows = static_cast<uint32_t>(mat.rows());
		const uint32_t cols = static_cast<uint32_t>(mat.cols());
		const uint32_t elemSize = static_cast<uint32_t>(sizeof(Type));
		const uint32_t reserved = 0;
		
		file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
		file.write(reinterpret_cast<const char*>(&version), sizeof(version));
		file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
		file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
		file.write(reinterpret_cast<const char*>(&elemSize), sizeof(elemSize));
		file.write(reinterpret_cast<const char*>(&reserved), sizeof(reserved));
		
		// Write data - use data() for contiguous access
		file.write(reinterpret_cast<const char*>(mat.data()), 
		           static_cast<std::streamsize>(rows) * cols * sizeof(Type));
		
		if (file.fail()) {
			return {false, SerializeError::WRITE_FAILED, "Write error while saving binary matrix to " + filename};
		}
		
		file.close();
		return {true, SerializeError::OK, "Success"};
	}

	/**
	 * @brief Load matrix from binary file
	 * @tparam Type Matrix element type
	 * @param filename Path to input file
	 * @param outMat Matrix to populate
	 * @return SerializeResult with success/failure info
	 */
	template<class Type>
	SerializeResult LoadMatrixFromBinary(const std::string& filename, Matrix<Type>& outMat)
	{
		std::ifstream file(filename, std::ios::binary);
		if (!file.is_open()) {
			return {false, SerializeError::FILE_NOT_OPENED, "Could not open binary file " + filename};
		}
		
		// Read and verify header
		uint32_t magic, version, rows, cols, elemSize, reserved;
		
		file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
		if (magic != BinaryFormat::MAGIC_MATRIX) {
			return {false, SerializeError::INVALID_FORMAT, "Invalid magic number in " + filename};
		}
		
		file.read(reinterpret_cast<char*>(&version), sizeof(version));
		if (version != BinaryFormat::VERSION_MATRIX) {
			return {false, SerializeError::INVALID_FORMAT, "Unsupported version " + std::to_string(version) + " in " + filename};
		}
		
		file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
		file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
		file.read(reinterpret_cast<char*>(&elemSize), sizeof(elemSize));
		file.read(reinterpret_cast<char*>(&reserved), sizeof(reserved));
		
		if (elemSize != sizeof(Type)) {
			return {false, SerializeError::INVALID_FORMAT, "Element size mismatch in " + filename + " (file: " + std::to_string(elemSize) + ", expected: " + std::to_string(sizeof(Type)) + ")"};
		}
		if (rows == 0 || cols == 0 || rows > 100000 || cols > 100000) {
			return {false, SerializeError::INVALID_FORMAT, "Invalid dimensions (" + std::to_string(rows) + "x" + std::to_string(cols) + ") in " + filename};
		}
		
		// Read data
		outMat.Resize(static_cast<int>(rows), static_cast<int>(cols));
		file.read(reinterpret_cast<char*>(outMat.data()), 
		          static_cast<std::streamsize>(rows) * cols * sizeof(Type));
		
		if (file.fail()) {
			return {false, SerializeError::READ_FAILED, "Read error while loading binary matrix from " + filename};
		}
		
		file.close();
		return {true, SerializeError::OK, "Success"};
	}

	/**
	 * @brief Get estimated binary file size in bytes
	 * @tparam Type Matrix element type
	 * @param rows Number of rows
	 * @param cols Number of columns
	 * @return Estimated file size in bytes
	 */
	template<class Type>
	size_t EstimateMatrixBinaryFileSize(int rows, int cols)
	{
		return 24 + static_cast<size_t>(rows) * static_cast<size_t>(cols) * sizeof(Type);
	}

} // namespace Serializer
} // namespace MML

#endif // MML_MATRIX_IO_H
