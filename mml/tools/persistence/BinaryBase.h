///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        persistence/BinaryBase.h                                            ///
///  Description: Common binary envelope and scalar encoding helpers                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_BINARY_BASE_H
#define MML_PERSISTENCE_BINARY_BASE_H

#include <mml/tools/persistence/PersistenceBase.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <istream>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>

namespace MML
{
namespace Persistence
{
	enum class BinaryObjectKind : uint32_t
	{
		Unknown = 0,
		Vector = 1,
		VectorN = 2,
		Matrix = 3,
		MatrixNM = 4,
		SparseCOO = 5,
		SparseCSR = 6,
		SparseCSC = 7,
		Tensor = 8,
		ODESolution = 9,
		SampledData = 10
	};

	enum class BinaryScalarType : uint32_t
	{
		Unknown = 0,
		Float32 = 1,
		Float64 = 2,
		LongDouble = 3,
		ComplexFloat32 = 4,
		ComplexFloat64 = 5
	};

	enum class BinarySidecarPolicy : uint32_t
	{
		None = 0,
		Optional = 1,
		Required = 2
	};

	struct BinaryEnvelopeHeader
	{
		std::array<char, 8> magic = { 'M', 'M', 'L', '_', 'M', 'A', 'T', 'X' };
		uint16_t container_version = 1;
		uint16_t header_size = 48;
		BinaryObjectKind object_kind = BinaryObjectKind::Matrix;
		uint32_t object_schema_version = 1;
		BinaryScalarType scalar_type = BinaryScalarType::Float64;
		uint32_t scalar_byte_size = 8;
		uint32_t endian_marker = 0x01020304;
		uint64_t payload_byte_count = 0;
		uint32_t flags = 0;
		BinarySidecarPolicy sidecar_policy = BinarySidecarPolicy::None;
	};

	struct BinaryEnvelopeReadResult
	{
		SerializeResult result;
		BinaryEnvelopeHeader header;
	};

	namespace BinaryEnvelope
	{
		static inline constexpr uint16_t VERSION = 1;
		static inline constexpr uint16_t HEADER_SIZE = 48;
		static inline constexpr uint32_t ENDIAN_MARKER = 0x01020304;
	}

	inline const char* BinaryObjectCode(BinaryObjectKind kind)
	{
		switch (kind)
		{
		case BinaryObjectKind::Vector: return "VECT";
		case BinaryObjectKind::VectorN: return "VECN";
		case BinaryObjectKind::Matrix: return "MATX";
		case BinaryObjectKind::MatrixNM: return "MATN";
		case BinaryObjectKind::SparseCOO: return "COO_";
		case BinaryObjectKind::SparseCSR: return "CSR_";
		case BinaryObjectKind::SparseCSC: return "CSC_";
		case BinaryObjectKind::Tensor: return "TENS";
		case BinaryObjectKind::ODESolution: return "ODES";
		case BinaryObjectKind::SampledData: return "SMPL";
		default: return "????";
		}
	}

	inline BinaryObjectKind BinaryObjectKindFromCode(const char* code)
	{
		const std::string value(code, 4);
		if (value == "VECT") return BinaryObjectKind::Vector;
		if (value == "VECN") return BinaryObjectKind::VectorN;
		if (value == "MATX") return BinaryObjectKind::Matrix;
		if (value == "MATN") return BinaryObjectKind::MatrixNM;
		if (value == "COO_") return BinaryObjectKind::SparseCOO;
		if (value == "CSR_") return BinaryObjectKind::SparseCSR;
		if (value == "CSC_") return BinaryObjectKind::SparseCSC;
		if (value == "TENS") return BinaryObjectKind::Tensor;
		if (value == "ODES") return BinaryObjectKind::ODESolution;
		if (value == "SMPL") return BinaryObjectKind::SampledData;
		return BinaryObjectKind::Unknown;
	}

	inline std::array<char, 8> BinaryMagic(BinaryObjectKind kind)
	{
		std::array<char, 8> magic = { 'M', 'M', 'L', '_', '?', '?', '?', '?' };
		const char* code = BinaryObjectCode(kind);
		for (int i = 0; i < 4; ++i)
			magic[4 + i] = code[i];
		return magic;
	}

	constexpr uint32_t ExpectedScalarByteSize(BinaryScalarType scalarType)
	{
		switch (scalarType)
		{
		case BinaryScalarType::Float32: return 4;
		case BinaryScalarType::Float64: return 8;
		case BinaryScalarType::LongDouble: return 24;
		case BinaryScalarType::ComplexFloat32: return 8;
		case BinaryScalarType::ComplexFloat64: return 16;
		default: return 0;
		}
	}

	inline void WriteUInt16LE(std::ostream& out, uint16_t value)
	{
		const char bytes[2] = {
			static_cast<char>(value & 0xFF),
			static_cast<char>((value >> 8) & 0xFF)
		};
		out.write(bytes, 2);
	}

	inline void WriteUInt32LE(std::ostream& out, uint32_t value)
	{
		const char bytes[4] = {
			static_cast<char>(value & 0xFF),
			static_cast<char>((value >> 8) & 0xFF),
			static_cast<char>((value >> 16) & 0xFF),
			static_cast<char>((value >> 24) & 0xFF)
		};
		out.write(bytes, 4);
	}

	inline void WriteUInt64LE(std::ostream& out, uint64_t value)
	{
		for (int i = 0; i < 8; ++i)
		{
			const char byte = static_cast<char>((value >> (8 * i)) & 0xFF);
			out.write(&byte, 1);
		}
	}

	inline uint16_t ReadUInt16LE(const char* bytes)
	{
		return static_cast<uint16_t>(static_cast<unsigned char>(bytes[0])) |
		       static_cast<uint16_t>(static_cast<unsigned char>(bytes[1]) << 8);
	}

	inline uint32_t ReadUInt32LE(const char* bytes)
	{
		return static_cast<uint32_t>(static_cast<unsigned char>(bytes[0])) |
		       (static_cast<uint32_t>(static_cast<unsigned char>(bytes[1])) << 8) |
		       (static_cast<uint32_t>(static_cast<unsigned char>(bytes[2])) << 16) |
		       (static_cast<uint32_t>(static_cast<unsigned char>(bytes[3])) << 24);
	}

	inline uint64_t ReadUInt64LE(const char* bytes)
	{
		uint64_t value = 0;
		for (int i = 0; i < 8; ++i)
			value |= static_cast<uint64_t>(static_cast<unsigned char>(bytes[i])) << (8 * i);
		return value;
	}

	inline SerializeResult ValidateBinaryEnvelopeHeader(const BinaryEnvelopeHeader& header, const LoadOptions& options = {})
	{
		if (header.magic[0] != 'M' || header.magic[1] != 'M' || header.magic[2] != 'L' || header.magic[3] != '_')
			return SerializeFailure(SerializeError::INVALID_FORMAT, "Binary magic must start with MML_");

		const BinaryObjectKind magicKind = BinaryObjectKindFromCode(header.magic.data() + 4);
		if (magicKind == BinaryObjectKind::Unknown)
			return SerializeFailure(SerializeError::TYPE_MISMATCH, "Unknown binary object code in magic");
		if (magicKind != header.object_kind)
			return SerializeFailure(SerializeError::TYPE_MISMATCH, "Binary magic object code and object-kind field disagree");
		if (header.container_version != BinaryEnvelope::VERSION)
			return SerializeFailure(SerializeError::UNSUPPORTED_VERSION, "Unsupported binary envelope version");
		if (header.header_size < BinaryEnvelope::HEADER_SIZE)
			return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Binary envelope header size is too small");
		if (header.endian_marker != BinaryEnvelope::ENDIAN_MARKER)
			return SerializeFailure(SerializeError::ENDIAN_MISMATCH, "Binary envelope endian marker is invalid");
		if (header.flags != 0)
			return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Binary envelope v1 does not allow nonzero flags");
		if (static_cast<uint32_t>(header.sidecar_policy) > static_cast<uint32_t>(BinarySidecarPolicy::Required))
			return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Binary envelope sidecar policy is invalid");
		if (header.object_schema_version == 0)
			return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Binary object schema version must be positive");

		const uint32_t expectedScalarBytes = ExpectedScalarByteSize(header.scalar_type);
		if (expectedScalarBytes == 0)
			return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Binary scalar type is unsupported");
		if (header.scalar_byte_size != expectedScalarBytes)
			return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Binary scalar byte size does not match scalar type");
		if (header.payload_byte_count > options.max_allocation_bytes)
			return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Binary payload exceeds configured allocation limit");
		return SerializeSuccess();
	}

	inline SerializeResult WriteBinaryEnvelopeHeader(std::ostream& out, BinaryEnvelopeHeader header)
	{
		header.magic = BinaryMagic(header.object_kind);
		header.container_version = BinaryEnvelope::VERSION;
		header.header_size = BinaryEnvelope::HEADER_SIZE;
		header.endian_marker = BinaryEnvelope::ENDIAN_MARKER;

		SerializeResult validation = ValidateBinaryEnvelopeHeader(header);
		if (!validation.success)
			return validation;

		out.write(header.magic.data(), static_cast<std::streamsize>(header.magic.size()));
		WriteUInt16LE(out, header.container_version);
		WriteUInt16LE(out, header.header_size);
		WriteUInt32LE(out, static_cast<uint32_t>(header.object_kind));
		WriteUInt32LE(out, header.object_schema_version);
		WriteUInt32LE(out, static_cast<uint32_t>(header.scalar_type));
		WriteUInt32LE(out, header.scalar_byte_size);
		WriteUInt32LE(out, header.endian_marker);
		WriteUInt64LE(out, header.payload_byte_count);
		WriteUInt32LE(out, header.flags);
		WriteUInt32LE(out, static_cast<uint32_t>(header.sidecar_policy));

		if (out.fail())
			return SerializeFailure(SerializeError::WRITE_FAILED, "Failed to write binary envelope header");
		return SerializeSuccess();
	}

	inline BinaryEnvelopeReadResult ReadBinaryEnvelopeHeader(std::istream& in, const LoadOptions& options = {})
	{
		BinaryEnvelopeReadResult read;
		std::array<char, BinaryEnvelope::HEADER_SIZE> bytes{};
		in.read(bytes.data(), static_cast<std::streamsize>(bytes.size()));
		if (in.gcount() != static_cast<std::streamsize>(bytes.size()))
		{
			read.result = SerializeFailure(SerializeError::TRUNCATED_INPUT, "Could not read complete binary envelope header");
			return read;
		}

		for (std::size_t i = 0; i < read.header.magic.size(); ++i)
			read.header.magic[i] = bytes[i];
		read.header.container_version = ReadUInt16LE(bytes.data() + 8);
		read.header.header_size = ReadUInt16LE(bytes.data() + 10);
		read.header.object_kind = static_cast<BinaryObjectKind>(ReadUInt32LE(bytes.data() + 12));
		read.header.object_schema_version = ReadUInt32LE(bytes.data() + 16);
		read.header.scalar_type = static_cast<BinaryScalarType>(ReadUInt32LE(bytes.data() + 20));
		read.header.scalar_byte_size = ReadUInt32LE(bytes.data() + 24);
		read.header.endian_marker = ReadUInt32LE(bytes.data() + 28);
		read.header.payload_byte_count = ReadUInt64LE(bytes.data() + 32);
		read.header.flags = ReadUInt32LE(bytes.data() + 40);
		read.header.sidecar_policy = static_cast<BinarySidecarPolicy>(ReadUInt32LE(bytes.data() + 44));
		read.result = ValidateBinaryEnvelopeHeader(read.header, options);
		return read;
	}

	namespace Detail
	{
		template<class Type>
		constexpr BinaryScalarType BinaryScalarTypeFor()
		{
			if constexpr (std::is_same_v<Type, float>) return BinaryScalarType::Float32;
			else if constexpr (std::is_same_v<Type, double>) return BinaryScalarType::Float64;
			else if constexpr (std::is_same_v<Type, long double>) return BinaryScalarType::LongDouble;
			else if constexpr (std::is_same_v<Type, std::complex<float>>) return BinaryScalarType::ComplexFloat32;
			else if constexpr (std::is_same_v<Type, std::complex<double>>) return BinaryScalarType::ComplexFloat64;
			else return BinaryScalarType::Unknown;
		}

		template<class Type>
		constexpr uint32_t BinaryScalarByteSize()
		{
			return ExpectedScalarByteSize(BinaryScalarTypeFor<Type>());
		}

		template<class Type>
		SerializeResult WriteFloatingScalarLE(std::ostream& out, Type value)
		{
			static_assert(std::is_same_v<Type, float> || std::is_same_v<Type, double> || std::is_same_v<Type, long double>);
			if constexpr (std::is_same_v<Type, float>)
			{
				uint32_t bits = 0;
				std::memcpy(&bits, &value, sizeof(bits));
				WriteUInt32LE(out, bits);
			}
			else if constexpr (std::is_same_v<Type, double>)
			{
				uint64_t bits = 0;
				std::memcpy(&bits, &value, sizeof(bits));
				WriteUInt64LE(out, bits);
			}
			else
			{
				static_assert(std::numeric_limits<long double>::radix == 2);
				static_assert(std::numeric_limits<long double>::digits <= 128);

				const uint8_t valueClass = std::isnan(value) ? 2 : (std::isinf(value) ? 1 : 0);
				const uint8_t sign = std::signbit(value) ? 1 : 0;
				const uint16_t precision = static_cast<uint16_t>(std::numeric_limits<long double>::digits);
				int exponent = 0;
				uint64_t significandLow = 0;
				uint64_t significandHigh = 0;
				if (valueClass == 0 && value != 0)
				{
					long double fraction = std::frexp(std::fabs(value), &exponent);
					for (int bitIndex = precision - 1; bitIndex >= 0; --bitIndex)
					{
						fraction = std::ldexp(fraction, 1);
						if (fraction >= 1)
						{
							if (bitIndex < 64) significandLow |= uint64_t{1} << bitIndex;
							else significandHigh |= uint64_t{1} << (bitIndex - 64);
							fraction -= 1;
						}
					}
				}

				out.put(static_cast<char>(valueClass));
				out.put(static_cast<char>(sign));
				WriteUInt16LE(out, precision);
				WriteUInt32LE(out, static_cast<uint32_t>(static_cast<int32_t>(exponent)));
				WriteUInt64LE(out, significandLow);
				WriteUInt64LE(out, significandHigh);
			}
			return out.fail() ? SerializeFailure(SerializeError::WRITE_FAILED, "Failed to write binary scalar") : SerializeSuccess();
		}

		template<class Type>
		SerializeResult ReadFloatingScalarLE(std::istream& in, Type& value)
		{
			static_assert(std::is_same_v<Type, float> || std::is_same_v<Type, double> || std::is_same_v<Type, long double>);
			constexpr std::size_t byteCount = std::is_same_v<Type, long double> ? 24 : sizeof(Type);
			char bytes[byteCount]{};
			in.read(bytes, byteCount);
			if (in.gcount() != static_cast<std::streamsize>(byteCount))
				return SerializeFailure(SerializeError::TRUNCATED_INPUT, "Truncated binary scalar");
			if constexpr (std::is_same_v<Type, float>)
			{
				uint32_t bits = ReadUInt32LE(bytes);
				std::memcpy(&value, &bits, sizeof(value));
			}
			else if constexpr (std::is_same_v<Type, double>)
			{
				uint64_t bits = ReadUInt64LE(bytes);
				std::memcpy(&value, &bits, sizeof(value));
			}
			else
			{
				static_assert(std::numeric_limits<long double>::radix == 2);
				const uint8_t valueClass = static_cast<uint8_t>(bytes[0]);
				const uint8_t sign = static_cast<uint8_t>(bytes[1]);
				const uint16_t precision = ReadUInt16LE(bytes + 2);
				const int32_t exponent = static_cast<int32_t>(ReadUInt32LE(bytes + 4));
				const uint64_t significandLow = ReadUInt64LE(bytes + 8);
				const uint64_t significandHigh = ReadUInt64LE(bytes + 16);
				if (valueClass > 2 || sign > 1 || precision == 0 || precision > 128)
					return SerializeFailure(SerializeError::INVALID_FORMAT, "Invalid canonical long-double scalar");
				if (valueClass == 0 && precision > std::numeric_limits<long double>::digits)
					return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Stored long double exceeds target precision");

				if (valueClass == 1)
					value = std::numeric_limits<long double>::infinity();
				else if (valueClass == 2)
					value = std::numeric_limits<long double>::quiet_NaN();
				else
				{
					value = std::ldexp(static_cast<long double>(significandHigh), 64);
					value += static_cast<long double>(significandLow);
					value = std::ldexp(value, exponent - precision);
				}
				if (sign != 0) value = -value;
			}
			return SerializeSuccess();
		}

		template<class Type>
		SerializeResult WriteBinaryScalarLE(std::ostream& out, const Type& value)
		{
			static_assert(BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown,
			              "Binary persistence supports float, double, long double, complex<float>, and complex<double>");
			if constexpr (std::is_floating_point_v<Type>)
				return WriteFloatingScalarLE(out, value);
			else
			{
				SerializeResult result = WriteFloatingScalarLE(out, value.real());
				return result.success ? WriteFloatingScalarLE(out, value.imag()) : result;
			}
		}

		template<class Type>
		SerializeResult ReadBinaryScalarLE(std::istream& in, Type& value)
		{
			static_assert(BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown,
			              "Binary persistence supports float, double, long double, complex<float>, and complex<double>");
			if constexpr (std::is_floating_point_v<Type>)
				return ReadFloatingScalarLE(in, value);
			else
			{
				typename Type::value_type real{};
				typename Type::value_type imag{};
				SerializeResult result = ReadFloatingScalarLE(in, real);
				if (!result.success) return result;
				result = ReadFloatingScalarLE(in, imag);
				if (result.success) value = Type(real, imag);
				return result;
			}
		}

		template<class Type>
		SerializeResult ValidateBinaryScalarType(BinaryScalarType scalarType, uint32_t scalarByteSize)
		{
			if (scalarType != BinaryScalarTypeFor<Type>() || scalarByteSize != BinaryScalarByteSize<Type>())
				return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Binary scalar metadata does not match target type");
			return SerializeSuccess();
		}
	}

} // namespace Persistence
} // namespace MML

#endif // MML_PERSISTENCE_BINARY_BASE_H
