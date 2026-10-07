///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        persistence/MatrixBinary.h                                          ///
///  Description: Canonical binary persistence for dynamic matrices                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_MATRIX_BINARY_H
#define MML_PERSISTENCE_MATRIX_BINARY_H

#include <mml/base/Matrix/Matrix.h>
#include <mml/tools/persistence/BinaryBase.h>

#include <fstream>
#include <limits>
#include <type_traits>

namespace MML
{
namespace Persistence
{
	template<class Type>
	SerializeResult SaveBinary(std::ostream& out, const Matrix<Type>& mat, const SaveOptions& = {})
	{
		static_assert(Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown,
		              "Matrix binary persistence supports float, double, long double, complex<float>, and complex<double>");

		const uint64_t rows = static_cast<uint64_t>(mat.rows());
		const uint64_t cols = static_cast<uint64_t>(mat.cols());
		const uint64_t scalarBytes = Detail::BinaryScalarByteSize<Type>();
		if (cols != 0 && rows > (std::numeric_limits<uint64_t>::max() - 16) / cols / scalarBytes)
			return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary payload size overflows uint64");

		BinaryEnvelopeHeader header;
		header.object_kind = BinaryObjectKind::Matrix;
		header.object_schema_version = 1;
		header.scalar_type = Detail::BinaryScalarTypeFor<Type>();
		header.scalar_byte_size = static_cast<uint32_t>(scalarBytes);
		header.payload_byte_count = 16 + rows * cols * scalarBytes;

		SerializeResult result = WriteBinaryEnvelopeHeader(out, header);
		if (!result.success) return result;
		WriteUInt64LE(out, rows);
		WriteUInt64LE(out, cols);
		for (int row = 0; row < mat.rows(); ++row)
			for (int col = 0; col < mat.cols(); ++col)
			{
				result = Detail::WriteBinaryScalarLE(out, mat(row, col));
				if (!result.success) return result;
			}

		return out.fail() ? SerializeFailure(SerializeError::WRITE_FAILED, "Failed to write Matrix binary payload") : SerializeSuccess();
	}

	template<class Type>
	SerializeResult LoadBinary(std::istream& in, Matrix<Type>& mat, const LoadOptions& options = {})
	{
		static_assert(Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown,
		              "Matrix binary persistence supports float, double, long double, complex<float>, and complex<double>");

		BinaryEnvelopeReadResult envelope = ReadBinaryEnvelopeHeader(in, options);
		if (!envelope.result.success) return envelope.result;
		if (envelope.header.object_kind != BinaryObjectKind::Matrix)
			return SerializeFailure(SerializeError::TYPE_MISMATCH, "Binary payload is not a Matrix");
		SerializeResult result = Detail::ValidateBinaryScalarType<Type>(envelope.header.scalar_type, envelope.header.scalar_byte_size);
		if (!result.success) return result;
		if (envelope.header.payload_byte_count < 16)
			return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix binary payload is too small for shape header");

		char shapeBytes[16]{};
		in.read(shapeBytes, 16);
		if (in.gcount() != 16)
			return SerializeFailure(SerializeError::TRUNCATED_INPUT, "Truncated Matrix binary shape header");
		const uint64_t rows = ReadUInt64LE(shapeBytes);
		const uint64_t cols = ReadUInt64LE(shapeBytes + 8);
		if (rows > static_cast<uint64_t>(std::numeric_limits<int>::max()) || cols > static_cast<uint64_t>(std::numeric_limits<int>::max()))
			return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary dimensions exceed supported range");
		const uint64_t count = rows * cols;
		if (cols != 0 && count / cols != rows)
			return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary shape overflows element count");
		const uint64_t scalarBytes = Detail::BinaryScalarByteSize<Type>();
		if (count > std::numeric_limits<uint64_t>::max() / scalarBytes)
			return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary byte count overflows uint64");
		if (envelope.header.payload_byte_count != 16 + count * scalarBytes)
			return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix binary payload size does not match shape");
		if (count * scalarBytes > options.max_allocation_bytes)
			return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary data exceeds configured allocation limit");

		mat.Resize(static_cast<int>(rows), static_cast<int>(cols));
		for (int row = 0; row < mat.rows(); ++row)
			for (int col = 0; col < mat.cols(); ++col)
			{
				Type value{};
				result = Detail::ReadBinaryScalarLE(in, value);
				if (!result.success) return result;
				mat(row, col) = value;
			}
		return SerializeSuccess();
	}

	template<class Type>
	SerializeResult SaveBinary(const Matrix<Type>& mat, const std::string& filename, const SaveOptions& options = {})
	{
		std::ofstream file(filename, std::ios::binary);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not create binary file " + filename);
		return SaveBinary(file, mat, options);
	}

	template<class Type>
	SerializeResult LoadBinary(const std::string& filename, Matrix<Type>& mat, const LoadOptions& options = {})
	{
		std::ifstream file(filename, std::ios::binary);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not open binary file " + filename);
		return LoadBinary(file, mat, options);
	}

	template<class Type>
	SerializeResult SaveMatrixToBinary(const Matrix<Type>& mat, const std::string& filename)
	{
		return SaveBinary(mat, filename);
	}

	template<class Type>
	SerializeResult LoadMatrixFromBinary(const std::string& filename, Matrix<Type>& mat)
	{
		return LoadBinary(filename, mat);
	}

	template<class Type>
	std::size_t EstimateMatrixBinaryFileSize(int rows, int cols)
	{
		return BinaryEnvelope::HEADER_SIZE + 16 + static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols) * Detail::BinaryScalarByteSize<Type>();
	}

} // namespace Persistence
} // namespace MML

#endif // MML_PERSISTENCE_MATRIX_BINARY_H
