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
/// File format (.mmlc):
///   Header (16 bytes):
///     - 4 bytes: Magic number (0x4D4D4C43 = "MMLC")
///     - 4 bytes: Format version (currently 1)
///     - 4 bytes: Count (number of elements)
///     - 4 bytes: Element size (sizeof(T), typically 8 for double)
///   Data:
///     - count * sizeof(T) bytes: Real parts (contiguous)
///     - count * sizeof(T) bytes: Imaginary parts (contiguous)
///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_VECTOR_IO_H
#define MML_VECTOR_IO_H

#include "mml/MMLBase.h"

#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <iostream>

namespace MML
{
namespace Serializer
{

/// @brief Save complex vector to binary file
/// @tparam T Scalar type (typically Real = double)
/// @param vec Vector of complex values
/// @param filename Output file path
/// @return true if successful
template<typename T>
bool SaveComplexVector(const std::vector<std::complex<T>>& vec, 
                       const std::string& filename)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: could not create binary file " << filename << std::endl;
        return false;
    }
    
    // Write header
    const uint32_t magic = BinaryFormat::MAGIC_VECTOR_COMPLEX;
    const uint32_t version = BinaryFormat::VERSION_VECTOR_COMPLEX;
    const uint32_t count = static_cast<uint32_t>(vec.size());
    const uint32_t elemSize = static_cast<uint32_t>(sizeof(T));
    
    file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));
    file.write(reinterpret_cast<const char*>(&count), sizeof(count));
    file.write(reinterpret_cast<const char*>(&elemSize), sizeof(elemSize));
    
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
    
    file.close();
    return true;
}

/// @brief Load complex vector from binary file
/// @tparam T Scalar type (typically Real = double)
/// @param filename Input file path
/// @param vec Vector to populate with loaded values
/// @return true if successful
template<typename T>
bool LoadComplexVector(const std::string& filename,
                       std::vector<std::complex<T>>& vec)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: could not open binary file " << filename << std::endl;
        return false;
    }
    
    // Read and verify header
    uint32_t magic, version, count, elemSize;
    
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (magic != BinaryFormat::MAGIC_VECTOR_COMPLEX) {
        std::cerr << "Error: invalid magic number in " << filename << std::endl;
        return false;
    }
    
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != BinaryFormat::VERSION_VECTOR_COMPLEX) {
        std::cerr << "Error: unsupported version " << version << " in " << filename << std::endl;
        return false;
    }
    
    file.read(reinterpret_cast<char*>(&count), sizeof(count));
    file.read(reinterpret_cast<char*>(&elemSize), sizeof(elemSize));
    
    if (elemSize != sizeof(T)) {
        std::cerr << "Error: element size mismatch in " << filename 
                  << " (file: " << elemSize << ", expected: " << sizeof(T) << ")" << std::endl;
        return false;
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
    
    // Construct complex values
    vec.resize(count);
    for (uint32_t i = 0; i < count; ++i) {
        vec[i] = std::complex<T>(realParts[i], imagParts[i]);
    }
    
    file.close();
    return true;
}

/// @brief Save real vector as complex (zero imaginary parts)
/// @details Convenience for storing eigenvalues from symmetric matrices
/// @tparam T Scalar type
/// @param vec Vector of real values
/// @param filename Output file path
/// @return true if successful
template<typename T>
bool SaveRealAsComplex(const std::vector<T>& vec, const std::string& filename)
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
/// @return true if successful
template<typename T>
bool SaveVectorAsComplex(const Vector<T>& vec, const std::string& filename)
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
/// @return true if successful
template<typename T>
bool LoadComplexAsReal(const std::string& filename, std::vector<T>& vec)
{
    std::vector<std::complex<T>> complexVec;
    if (!LoadComplexVector(filename, complexVec)) {
        return false;
    }
    
    vec.resize(complexVec.size());
    for (size_t i = 0; i < complexVec.size(); ++i) {
        vec[i] = complexVec[i].real();
    }
    return true;
}

/// @brief Load complex vector into MML Vector (real parts only)
/// @details For symmetric matrix eigenvalues as Vector<Real>
/// @tparam T Scalar type
/// @param filename Input file path
/// @param vec MML Vector to populate
/// @return true if successful
template<typename T>
bool LoadComplexAsVector(const std::string& filename, Vector<T>& vec)
{
    std::vector<std::complex<T>> complexVec;
    if (!LoadComplexVector(filename, complexVec)) {
        return false;
    }
    
    vec = Vector<T>(static_cast<int>(complexVec.size()));
    for (size_t i = 0; i < complexVec.size(); ++i) {
        vec[static_cast<int>(i)] = complexVec[i].real();
    }
    return true;
}

} // namespace Serializer
} // namespace MML

#endif // MML_VECTOR_IO_H
