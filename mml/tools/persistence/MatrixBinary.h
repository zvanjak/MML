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
	namespace Detail
	{
		inline constexpr uint32_t MATRIX_BINARY_SCHEMA_VERSION = 1;

		template<class Type>
		IOResult BuildMatrixBinaryHeader(uint64_t rows, uint64_t cols, BinaryEnvelopeHeader& header) {
			const uint64_t scalarBytes = BinaryScalarByteSize<Type>();
			if (cols != 0 && rows > (std::numeric_limits<uint64_t>::max() - 16) / cols / scalarBytes)
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary payload size overflows uint64");
			header.object_kind = BinaryObjectKind::Matrix;
			header.magic = BinaryMagic(header.object_kind);
			header.object_schema_version = MATRIX_BINARY_SCHEMA_VERSION;
			header.scalar_type = BinaryScalarTypeFor<Type>();
			header.scalar_byte_size = static_cast<uint32_t>(scalarBytes);
			header.payload_byte_count = 16 + rows * cols * scalarBytes;
			return ValidateBinaryEnvelopeHeader(header);
		}
	}

	template<class Type>
	IOResult SaveBinary(std::ostream& out, const Matrix<Type>& mat, const SaveOptions& = {}) {
		static_assert(Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown,
		              "Matrix binary persistence supports float, double, long double, complex<float>, and complex<double>");

		const uint64_t rows = static_cast<uint64_t>(mat.rows());
		const uint64_t cols = static_cast<uint64_t>(mat.cols());
		BinaryEnvelopeHeader header;
		IOResult result = Detail::BuildMatrixBinaryHeader<Type>(rows, cols, header);
		if (!result.success) return result;
		result = WriteBinaryEnvelopeHeader(out, header);
		if (!result.success) return result;
		WriteUInt64LE(out, rows);
		WriteUInt64LE(out, cols);
		for (int row = 0; row < mat.rows(); ++row)
			for (int col = 0; col < mat.cols(); ++col)
			{
				result = Detail::WriteBinaryScalarLE(out, mat(row, col));
				if (!result.success) return result;
			}

		return out.fail() ? SerializeFailure(IOError::WRITE_FAILED, "Failed to write Matrix binary payload") : SerializeSuccess();
	}

	template<class Type>
	IOResult LoadBinary(std::istream& in, Matrix<Type>& mat, const LoadOptions& options = {}) {
		static_assert(Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown,
		              "Matrix binary persistence supports float, double, long double, complex<float>, and complex<double>");

		BinaryEnvelopeReadResult envelope = ReadBinaryEnvelopeHeader(in, options);
		if (!envelope.result.success) return envelope.result;
		if (envelope.header.object_kind != BinaryObjectKind::Matrix)
			return SerializeFailure(IOError::TYPE_MISMATCH, "Binary payload is not a Matrix");
		if (envelope.header.object_schema_version != Detail::MATRIX_BINARY_SCHEMA_VERSION)
			return SerializeFailure(IOError::UNSUPPORTED_VERSION, "Unsupported Matrix binary object schema version");
		IOResult result = Detail::ValidateBinaryScalarType<Type>(envelope.header.scalar_type, envelope.header.scalar_byte_size);
		if (!result.success) return result;
		if (envelope.header.payload_byte_count < 16)
			return SerializeFailure(IOError::SCHEMA_MISMATCH, "Matrix binary payload is too small for shape header");

		char shapeBytes[16]{};
		in.read(shapeBytes, 16);
		if (in.gcount() != 16)
			return SerializeFailure(IOError::TRUNCATED_INPUT, "Truncated Matrix binary shape header");
		const uint64_t rows = ReadUInt64LE(shapeBytes);
		const uint64_t cols = ReadUInt64LE(shapeBytes + 8);
		if (rows > static_cast<uint64_t>(std::numeric_limits<int>::max()) || cols > static_cast<uint64_t>(std::numeric_limits<int>::max()))
			return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary dimensions exceed supported range");
		const uint64_t count = rows * cols;
		if (cols != 0 && count / cols != rows)
			return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary shape overflows element count");
		const uint64_t scalarBytes = Detail::BinaryScalarByteSize<Type>();
		if (count > std::numeric_limits<uint64_t>::max() / scalarBytes)
			return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary byte count overflows uint64");
		if (envelope.header.payload_byte_count != 16 + count * scalarBytes)
			return SerializeFailure(IOError::SCHEMA_MISMATCH, "Matrix binary payload size does not match shape");
		if (count * scalarBytes > options.max_allocation_bytes)
			return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Matrix binary data exceeds configured allocation limit");

		mat.resize(static_cast<int>(rows), static_cast<int>(cols));
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
	IOResult SaveBinary(const Matrix<Type>& mat, const std::string& filename, const SaveOptions& options = {}) {
		auto validation = Detail::ValidateSaveFilename(filename);
		if (!validation.success) return validation;
		BinaryEnvelopeHeader header;
		validation = Detail::BuildMatrixBinaryHeader<Type>(static_cast<uint64_t>(mat.rows()), static_cast<uint64_t>(mat.cols()), header);
		if (!validation.success) return validation.WithFile(filename);
		std::ofstream file(filename, std::ios::binary);
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not create binary file " + filename).WithFile(filename);
		return FinishOutputFile(file, filename, SaveBinary(file, mat, options));
	}

	template<class Type>
	IOResult LoadBinary(const std::string& filename, Matrix<Type>& mat, const LoadOptions& options = {}) {
		std::ifstream file(filename, std::ios::binary);
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not open binary file " + filename).WithFile(filename);
		return InputResult(file, LoadBinary(file, mat, options), filename);
	}

	template<class Type>
	IOResult SaveMatrixToBinary(const Matrix<Type>& mat, const std::string& filename) {
		return (SaveBinary(mat, filename)).WithFile(filename);
	}

	template<class Type>
	IOResult LoadMatrixFromBinary(const std::string& filename, Matrix<Type>& mat) {
		return (LoadBinary(filename, mat)).WithFile(filename);
	}

	template<class Type>
	std::size_t EstimateMatrixBinaryFileSize(int rows, int cols)
	{
		return BinaryEnvelope::HEADER_SIZE + 16 + static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols) * Detail::BinaryScalarByteSize<Type>();
	}

} // namespace Persistence
} // namespace MML

#endif // MML_PERSISTENCE_MATRIX_BINARY_H
