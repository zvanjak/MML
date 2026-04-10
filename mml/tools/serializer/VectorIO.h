///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorIO.h                                                          ///
///  Description: Binary serialization for complex vectors                            ///
///               Save/load std::vector<std::complex<T>> to/from binary files         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
///
/// File format (.mmlc) - Version 2:
///   Header (20 bytes):
///     - 4 bytes: Magic number (0x4D4D4C43 = "MMLC")
///     - 4 bytes: Format version (2)
///     - 4 bytes: Count (number of elements)
///     - 4 bytes: Element size (sizeof(T), typically 8 for double)
///     - 4 bytes: Endianness marker (0x01020304 = native byte order)
///   Data:
///     - count * sizeof(T) bytes: Real parts (contiguous)
///     - count * sizeof(T) bytes: Imaginary parts (contiguous)
///
/// Version 1 files (16-byte header, no endianness) are still readable.
///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_VECTOR_IO_H
#define MML_VECTOR_IO_H

#include "mml/MMLBase.h"
#include "mml/tools/serializer/SerializerBase.h"

#include <vector>
#include <complex>
#include <string>
#include <fstream>

namespace MML
{
namespace Serializer
{

namespace Detail
{
	/// Endianness marker: written as 0x01020304 in native byte order.
	/// On read, if the value is 0x04030201 the file was written on opposite endianness.
	static inline constexpr uint32_t ENDIAN_MARKER = 0x01020304;
	static inline constexpr uint32_t ENDIAN_MARKER_SWAPPED = 0x04030201;
	static inline constexpr uint32_t VERSION_VECTOR_COMPLEX_V2 = 2;
} // namespace Detail

/// @brief Save complex vector to binary file (v2 format with endianness marker)
/// @tparam T Scalar type (typically Real = double)
/// @param vec Vector of complex values
/// @param filename Output file path
/// @return SerializeResult with success/failure info
template<typename T>
SerializeResult SaveComplexVector(const std::vector<std::complex<T>>& vec, 
                                  const std::string& filename)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return {false, SerializeError::FILE_NOT_OPENED, "Could not create binary file " + filename};
    }
    
    // Write header (v2 with endianness marker)
    const uint32_t magic = BinaryFormat::MAGIC_VECTOR_COMPLEX;
    const uint32_t version = Detail::VERSION_VECTOR_COMPLEX_V2;
    const uint32_t count = static_cast<uint32_t>(vec.size());
    const uint32_t elemSize = static_cast<uint32_t>(sizeof(T));
    const uint32_t endianMarker = Detail::ENDIAN_MARKER;
    
    file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));
    file.write(reinterpret_cast<const char*>(&count), sizeof(count));
    file.write(reinterpret_cast<const char*>(&elemSize), sizeof(elemSize));
    file.write(reinterpret_cast<const char*>(&endianMarker), sizeof(endianMarker));
    
    // Write real parts (contiguous for cache efficiency)
    for (const auto& v : vec) {
        T realPart = v.real();
        file.write(reinterpret_cast<const char*>(&realPart), sizeof(T));
    }
    
    // Write imaginary parts (contiguous)
    for (const auto& v : vec) {
        T imagPart = v.imag();
        file.write(reinterpret_cast<const char*>(&imagPart), sizeof(T));
    }
    
    if (file.fail()) {
        return {false, SerializeError::WRITE_FAILED, "Write error while saving complex vector to " + filename};
    }
    
    file.close();
    return {true, SerializeError::OK, "Success"};
}

/// @brief Load complex vector from binary file (supports v1 and v2 formats)
/// @tparam T Scalar type (typically Real = double)
/// @param filename Input file path
/// @param vec Vector to populate with loaded values
/// @return SerializeResult with success/failure info
template<typename T>
SerializeResult LoadComplexVector(const std::string& filename,
                                  std::vector<std::complex<T>>& vec)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return {false, SerializeError::FILE_NOT_OPENED, "Could not open binary file " + filename};
    }
    
    // Read and verify header
    uint32_t magic, version, count, elemSize;
    
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (magic != BinaryFormat::MAGIC_VECTOR_COMPLEX) {
        return {false, SerializeError::INVALID_FORMAT, "Invalid magic number in " + filename};
    }
    
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != BinaryFormat::VERSION_VECTOR_COMPLEX && version != Detail::VERSION_VECTOR_COMPLEX_V2) {
        return {false, SerializeError::INVALID_FORMAT, "Unsupported version " + std::to_string(version) + " in " + filename};
    }
    
    file.read(reinterpret_cast<char*>(&count), sizeof(count));
    file.read(reinterpret_cast<char*>(&elemSize), sizeof(elemSize));
    
    if (elemSize != sizeof(T)) {
        return {false, SerializeError::INVALID_FORMAT, "Element size mismatch in " + filename + " (file: " + std::to_string(elemSize) + ", expected: " + std::to_string(sizeof(T)) + ")"};
    }

    // V2: read and check endianness marker
    if (version == Detail::VERSION_VECTOR_COMPLEX_V2) {
        uint32_t endianMarker;
        file.read(reinterpret_cast<char*>(&endianMarker), sizeof(endianMarker));
        if (endianMarker == Detail::ENDIAN_MARKER_SWAPPED) {
            return {false, SerializeError::INVALID_FORMAT, "Endianness mismatch: file " + filename + " was written on a machine with different byte order"};
        }
        if (endianMarker != Detail::ENDIAN_MARKER) {
            return {false, SerializeError::INVALID_FORMAT, "Corrupted endianness marker in " + filename};
        }
    }

    // Validate count against actual file size
    auto currentPos = file.tellg();
    file.seekg(0, std::ios::end);
    auto fileSize = file.tellg();
    file.seekg(currentPos);

    auto expectedDataSize = static_cast<std::streamoff>(count) * elemSize * 2;
    if (fileSize - currentPos < expectedDataSize) {
        return {false, SerializeError::INVALID_FORMAT, "File " + filename + " is too small for declared count " + std::to_string(count)};
    }

    // Read real parts
    std::vector<T> realParts(count);
    for (uint32_t i = 0; i < count; ++i) {
        file.read(reinterpret_cast<char*>(&realParts[i]), sizeof(T));
    }
    
    // Read imaginary parts
    std::vector<T> imagParts(count);
    for (uint32_t i = 0; i < count; ++i) {
        file.read(reinterpret_cast<char*>(&imagParts[i]), sizeof(T));
    }
    
    if (file.fail()) {
        return {false, SerializeError::READ_FAILED, "Read error while loading complex vector from " + filename};
    }
    
    // Construct complex values
    vec.resize(count);
    for (uint32_t i = 0; i < count; ++i) {
        vec[i] = std::complex<T>(realParts[i], imagParts[i]);
    }
    
    file.close();
    return {true, SerializeError::OK, "Success"};
}

/// @brief Save real vector as complex (zero imaginary parts)
/// @details Convenience for storing eigenvalues from symmetric matrices
/// @tparam T Scalar type
/// @param vec Vector of real values
/// @param filename Output file path
/// @return SerializeResult with success/failure info
template<typename T>
SerializeResult SaveRealAsComplex(const std::vector<T>& vec, const std::string& filename)
{
    std::vector<std::complex<T>> complexVec;
    complexVec.reserve(vec.size());
    for (const auto& v : vec) {
        complexVec.emplace_back(v, T(0));
    }
    return SaveComplexVector(complexVec, filename);
}

/// @brief Save MML Vector<T> as complex (zero imaginary parts)
/// @details For symmetric matrix eigenvalues stored as Vector<Real>
/// @tparam T Scalar type
/// @param vec MML Vector of real values
/// @param filename Output file path
/// @return SerializeResult with success/failure info
template<typename T>
SerializeResult SaveVectorAsComplex(const Vector<T>& vec, const std::string& filename)
{
    std::vector<std::complex<T>> complexVec;
    complexVec.reserve(vec.size());
    for (int i = 0; i < vec.size(); ++i) {
        complexVec.emplace_back(vec[i], T(0));
    }
    return SaveComplexVector(complexVec, filename);
}

/// @brief Load complex vector as real values (extract real parts only)
/// @details Use when you know values are real (e.g., symmetric matrix eigenvalues)
/// @tparam T Scalar type
/// @param filename Input file path
/// @param vec Vector to populate with real parts
/// @return SerializeResult with success/failure info
template<typename T>
SerializeResult LoadComplexAsReal(const std::string& filename, std::vector<T>& vec)
{
    std::vector<std::complex<T>> complexVec;
    SerializeResult result = LoadComplexVector(filename, complexVec);
    if (!result) {
        return result;
    }
    
    vec.resize(complexVec.size());
    for (size_t i = 0; i < complexVec.size(); ++i) {
        vec[i] = complexVec[i].real();
    }
    return {true, SerializeError::OK, "Success"};
}

/// @brief Load complex vector into MML Vector (real parts only)
/// @details For symmetric matrix eigenvalues as Vector<Real>
/// @tparam T Scalar type
/// @param filename Input file path
/// @param vec MML Vector to populate
/// @return SerializeResult with success/failure info
template<typename T>
SerializeResult LoadComplexAsVector(const std::string& filename, Vector<T>& vec)
{
    std::vector<std::complex<T>> complexVec;
    SerializeResult result = LoadComplexVector(filename, complexVec);
    if (!result) {
        return result;
    }
    
    vec = Vector<T>(static_cast<int>(complexVec.size()));
    for (size_t i = 0; i < complexVec.size(); ++i) {
        vec[static_cast<int>(i)] = complexVec[i].real();
    }
    return {true, SerializeError::OK, "Success"};
}

} // namespace Serializer
} // namespace MML

#endif // MML_VECTOR_IO_H
