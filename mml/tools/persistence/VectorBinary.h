///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        persistence/VectorBinary.h                                          ///
///  Description: Canonical binary persistence for real and complex vectors           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_VECTOR_BINARY_H
#define MML_PERSISTENCE_VECTOR_BINARY_H

#include <mml/base/Vector/Vector.h>
#include <mml/tools/persistence/BinaryBase.h>

#include <complex>
#include <fstream>
#include <limits>
#include <type_traits>
#include <vector>

namespace MML
{
namespace Persistence
{
	namespace Detail
	{
		template<class Type, class SizeFunction, class ValueFunction>
		IOResult SaveVectorBinary(std::ostream& out, SizeFunction size, ValueFunction value) {
			const uint64_t count = static_cast<uint64_t>(size());
			const uint64_t scalarBytes = BinaryScalarByteSize<Type>();
			if (count > (std::numeric_limits<uint64_t>::max() - 8) / scalarBytes)
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector binary payload size overflows uint64");

			BinaryEnvelopeHeader header;
			header.object_kind = BinaryObjectKind::Vector;
			header.object_schema_version = 1;
			header.scalar_type = BinaryScalarTypeFor<Type>();
			header.scalar_byte_size = static_cast<uint32_t>(scalarBytes);
			header.payload_byte_count = 8 + count * scalarBytes;

			IOResult result = WriteBinaryEnvelopeHeader(out, header);
			if (!result.success) return result;
			WriteUInt64LE(out, count);
			for (uint64_t index = 0; index < count; ++index)
			{
				result = WriteBinaryScalarLE(out, value(index));
				if (!result.success) return result;
			}
			return out.fail() ? SerializeFailure(IOError::WRITE_FAILED, "Failed to write Vector binary payload") : SerializeSuccess();
		}

		template<class Type, class ResizeFunction, class AssignFunction>
		IOResult LoadVectorBinary(std::istream& in, ResizeFunction resize, AssignFunction assign, const LoadOptions& options) {
			BinaryEnvelopeReadResult envelope = ReadBinaryEnvelopeHeader(in, options);
			if (!envelope.result.success) return envelope.result;
			if (envelope.header.object_kind != BinaryObjectKind::Vector)
				return SerializeFailure(IOError::TYPE_MISMATCH, "Binary payload is not a Vector");
			IOResult result = ValidateBinaryScalarType<Type>(envelope.header.scalar_type, envelope.header.scalar_byte_size);
			if (!result.success) return result;
			if (envelope.header.payload_byte_count < 8)
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Vector binary payload is too small for shape header");

			char countBytes[8]{};
			in.read(countBytes, 8);
			if (in.gcount() != 8)
				return SerializeFailure(IOError::TRUNCATED_INPUT, "Truncated Vector binary shape header");
			const uint64_t count = ReadUInt64LE(countBytes);
			if (count > static_cast<uint64_t>(std::numeric_limits<int>::max()))
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector binary count exceeds supported dimension range");
			const uint64_t scalarBytes = BinaryScalarByteSize<Type>();
			if (count > std::numeric_limits<uint64_t>::max() / scalarBytes)
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector binary byte count overflows uint64");
			if (envelope.header.payload_byte_count != 8 + count * scalarBytes)
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Vector binary payload size does not match element count");
			if (count * scalarBytes > options.max_allocation_bytes)
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector binary data exceeds configured allocation limit");

			resize(static_cast<int>(count));
			for (uint64_t index = 0; index < count; ++index)
			{
				Type value{};
				result = ReadBinaryScalarLE(in, value);
				if (!result.success) return result;
				assign(static_cast<int>(index), value);
			}
			return SerializeSuccess();
		}
	}

	template<class Type>
	IOResult SaveBinary(std::ostream& out, const Vector<Type>& vec, const SaveOptions& = {}) {
		static_assert(Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown,
		              "Vector binary persistence supports float, double, long double, complex<float>, and complex<double>");
		return Detail::SaveVectorBinary<Type>(out, [&] { return vec.size(); }, [&](uint64_t index) { return vec[static_cast<int>(index)]; });
	}

	template<class Type>
	IOResult LoadBinary(std::istream& in, Vector<Type>& vec, const LoadOptions& options = {}) {
		static_assert(Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown,
		              "Vector binary persistence supports float, double, long double, complex<float>, and complex<double>");
		return Detail::LoadVectorBinary<Type>(in, [&](int count) { vec.resize(count); },
		                                      [&](int index, Type value) { vec[index] = value; }, options);
	}

	template<class Type>
	IOResult SaveBinary(std::ostream& out, const std::vector<std::complex<Type>>& vec, const SaveOptions& = {}) {
		static_assert(std::is_same_v<Type, float> || std::is_same_v<Type, double>,
		              "Complex Vector binary persistence currently supports float and double components");
		using ComplexType = std::complex<Type>;
		return Detail::SaveVectorBinary<ComplexType>(out, [&] { return vec.size(); }, [&](uint64_t index) { return vec[static_cast<std::size_t>(index)]; });
	}

	template<class Type>
	IOResult LoadBinary(std::istream& in, std::vector<std::complex<Type>>& vec, const LoadOptions& options = {}) {
		static_assert(std::is_same_v<Type, float> || std::is_same_v<Type, double>,
		              "Complex Vector binary persistence currently supports float and double components");
		using ComplexType = std::complex<Type>;
		return Detail::LoadVectorBinary<ComplexType>(in, [&](int count) { vec.resize(static_cast<std::size_t>(count)); },
		                                             [&](int index, ComplexType value) { vec[static_cast<std::size_t>(index)] = value; }, options);
	}

	template<class VectorType>
	IOResult SaveVectorBinaryFile(const VectorType& vec, const std::string& filename, const SaveOptions& options = {}) {
		std::ofstream file(filename, std::ios::binary);
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not create binary file " + filename).WithFile(filename);
		return FinishOutputFile(file, filename, SaveBinary(file, vec, options));
	}

	template<class VectorType>
	IOResult LoadVectorBinaryFile(const std::string& filename, VectorType& vec, const LoadOptions& options = {}) {
		std::ifstream file(filename, std::ios::binary);
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not open binary file " + filename).WithFile(filename);
		return InputResult(file, LoadBinary(file, vec, options), filename);
	}

	template<class Type>
	IOResult SaveBinary(const Vector<Type>& vec, const std::string& filename, const SaveOptions& options = {}) {
		return (SaveVectorBinaryFile(vec, filename, options)).WithFile(filename);
	}

	template<class Type>
	IOResult LoadBinary(const std::string& filename, Vector<Type>& vec, const LoadOptions& options = {}) {
		return (LoadVectorBinaryFile(filename, vec, options)).WithFile(filename);
	}

	template<class Type>
	IOResult SaveBinary(const std::vector<std::complex<Type>>& vec, const std::string& filename, const SaveOptions& options = {}) {
		return (SaveVectorBinaryFile(vec, filename, options)).WithFile(filename);
	}

	template<class Type>
	IOResult LoadBinary(const std::string& filename, std::vector<std::complex<Type>>& vec, const LoadOptions& options = {}) {
		return (LoadVectorBinaryFile(filename, vec, options)).WithFile(filename);
	}

	template<class Type>
	IOResult SaveComplexVector(const std::vector<std::complex<Type>>& vec, const std::string& filename) {
		return (SaveBinary(vec, filename)).WithFile(filename);
	}

	template<class Type>
	IOResult LoadComplexVector(const std::string& filename, std::vector<std::complex<Type>>& vec) {
		return (LoadBinary(filename, vec)).WithFile(filename);
	}

	template<class Type>
	IOResult SaveRealAsComplex(const std::vector<Type>& vec, const std::string& filename) {
		using ComponentType = std::conditional_t<std::is_same_v<Type, long double>, double, Type>;
		std::vector<std::complex<ComponentType>> complexVec;
		complexVec.reserve(vec.size());
		for (const Type& value : vec)
			complexVec.emplace_back(static_cast<ComponentType>(value), ComponentType(0));
		return (SaveBinary(complexVec, filename)).WithFile(filename);
	}

	template<class Type>
	IOResult SaveVectorAsComplex(const Vector<Type>& vec, const std::string& filename) {
		using ComponentType = std::conditional_t<std::is_same_v<Type, long double>, double, Type>;
		std::vector<std::complex<ComponentType>> complexVec;
		complexVec.reserve(static_cast<std::size_t>(vec.size()));
		for (int index = 0; index < vec.size(); ++index)
			complexVec.emplace_back(static_cast<ComponentType>(vec[index]), ComponentType(0));
		return (SaveBinary(complexVec, filename)).WithFile(filename);
	}

	template<class Type>
	IOResult LoadComplexAsReal(const std::string& filename, std::vector<Type>& vec) {
		using ComponentType = std::conditional_t<std::is_same_v<Type, long double>, double, Type>;
		std::vector<std::complex<ComponentType>> complexVec;
		IOResult result = LoadBinary(filename, complexVec);
		if (!result.success)
			return (result).WithFile(filename);
		vec.resize(complexVec.size());
		for (std::size_t index = 0; index < complexVec.size(); ++index)
			vec[index] = static_cast<Type>(complexVec[index].real());
		return (SerializeSuccess()).WithFile(filename);
	}

	template<class Type>
	IOResult LoadComplexAsVector(const std::string& filename, Vector<Type>& vec) {
		using ComponentType = std::conditional_t<std::is_same_v<Type, long double>, double, Type>;
		std::vector<std::complex<ComponentType>> complexVec;
		IOResult result = LoadBinary(filename, complexVec);
		if (!result.success)
			return (result).WithFile(filename);
		vec.resize(static_cast<int>(complexVec.size()));
		for (int index = 0; index < vec.size(); ++index)
			vec[index] = static_cast<Type>(complexVec[static_cast<std::size_t>(index)].real());
		return (SerializeSuccess()).WithFile(filename);
	}

} // namespace Persistence
} // namespace MML

#endif // MML_PERSISTENCE_VECTOR_BINARY_H
