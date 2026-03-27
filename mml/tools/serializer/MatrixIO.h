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

#include <iostream>
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
	 * @return true if successful
	 * 
	 * File format:
	 * - First line: rows cols
	 * - Following lines: space-separated matrix elements (row by row)
	 */
	template<class Type>
	bool LoadMatrixFromFile(const std::string& filename, Matrix<Type>& outMat)
	{
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Error: could not open file " << filename << " for reading." << std::endl;
			return false;
		}
		
		int rows, cols;
		file >> rows >> cols;
		if (file.fail() || rows <= 0 || cols <= 0 || rows > 100000 || cols > 100000) {
			std::cerr << "Error: invalid dimensions (" << rows << "x" << cols << ") in " << filename << std::endl;
			return false;
		}
		outMat.Resize(rows, cols);
		
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				file >> outMat(i, j);
		
		file.close();
		return true;
	}

	/**
	 * @brief Save matrix to text file
	 * @tparam Type Matrix element type
	 * @param mat Matrix to save
	 * @param filename Path to output file
	 * @return true if successful
	 */
	template<class Type>
	bool SaveMatrixToFile(const Matrix<Type>& mat, const std::string& filename)
	{
		std::ofstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Error: could not create file " << filename << " for writing." << std::endl;
			return false;
		}
		
		file << mat.rows() << " " << mat.cols() << std::endl;
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j)
				file << mat(i, j) << " ";
			file << std::endl;
		}
		
		file.close();
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                               CSV File I/O                                          ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Load matrix from CSV file
	 * @tparam Type Matrix element type
	 * @param filename Path to CSV file
	 * @param outMat Matrix to populate
	 * @return true if successful
	 */
	template<class Type>
	bool LoadMatrixFromCSV(const std::string& filename, Matrix<Type>& outMat)
	{
		std::ifstream file(filename);
		if (!file.is_open()) return false;
		
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
		if (data.empty()) return false;
		
		outMat = Matrix<Type>(data);
		return true;
	}

	/**
	 * @brief Save matrix to CSV file
	 * @tparam Type Matrix element type
	 * @param mat Matrix to save
	 * @param filename Path to CSV file
	 * @return true if successful
	 */
	template<class Type>
	bool SaveMatrixToCSV(const Matrix<Type>& mat, const std::string& filename)
	{
		std::ofstream file(filename);
		if (!file.is_open()) return false;
		
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				file << mat(i, j);
				if (j < mat.cols() - 1) file << ",";
			}
			file << "\n";
		}
		
		file.close();
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                               Binary File I/O                                       ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Save matrix to binary file (compact, exact precision)
	 * @tparam Type Matrix element type
	 * @param mat Matrix to save
	 * @param filename Path to output file
	 * @return true if successful
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
	bool SaveMatrixToBinary(const Matrix<Type>& mat, const std::string& filename)
	{
		std::ofstream file(filename, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "Error: could not create binary file " << filename << std::endl;
			return false;
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
		
		file.close();
		return true;
	}

	/**
	 * @brief Load matrix from binary file
	 * @tparam Type Matrix element type
	 * @param filename Path to input file
	 * @param outMat Matrix to populate
	 * @return true if successful
	 */
	template<class Type>
	bool LoadMatrixFromBinary(const std::string& filename, Matrix<Type>& outMat)
	{
		std::ifstream file(filename, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "Error: could not open binary file " << filename << std::endl;
			return false;
		}
		
		// Read and verify header
		uint32_t magic, version, rows, cols, elemSize, reserved;
		
		file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
		if (magic != BinaryFormat::MAGIC_MATRIX) {
			std::cerr << "Error: invalid magic number in " << filename << std::endl;
			return false;
		}
		
		file.read(reinterpret_cast<char*>(&version), sizeof(version));
		if (version != BinaryFormat::VERSION_MATRIX) {
			std::cerr << "Error: unsupported version " << version << " in " << filename << std::endl;
			return false;
		}
		
		file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
		file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
		file.read(reinterpret_cast<char*>(&elemSize), sizeof(elemSize));
		file.read(reinterpret_cast<char*>(&reserved), sizeof(reserved));
		
		if (elemSize != sizeof(Type)) {
			std::cerr << "Error: element size mismatch in " << filename 
			          << " (file: " << elemSize << ", expected: " << sizeof(Type) << ")" << std::endl;
			return false;
		}
		if (rows == 0 || cols == 0 || rows > 100000 || cols > 100000) {
			std::cerr << "Error: invalid dimensions (" << rows << "x" << cols << ") in " << filename << std::endl;
			return false;
		}
		
		// Read data
		outMat.Resize(static_cast<int>(rows), static_cast<int>(cols));
		file.read(reinterpret_cast<char*>(outMat.data()), 
		          static_cast<std::streamsize>(rows) * cols * sizeof(Type));
		
		file.close();
		return true;
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
