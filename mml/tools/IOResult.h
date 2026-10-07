///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        tools/IOResult.h                                                    ///
///  Description: Shared result and error types for MML data I/O                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_IO_RESULT_H
#define MML_IO_RESULT_H

#include <string>

namespace MML
{
	enum class SerializeError
	{
		OK,
		FILE_NOT_OPENED,
		INVALID_PARAMETERS,
		WRITE_FAILED,
		READ_FAILED,
		INVALID_FORMAT,
		UNSUPPORTED_VERSION,
		TYPE_MISMATCH,
		ENDIAN_MISMATCH,
		ALLOCATION_LIMIT_EXCEEDED,
		SCHEMA_MISMATCH,
		MALFORMED_INPUT,
		UNSUPPORTED_SCALAR,
		TRUNCATED_INPUT
	};

	inline const char* SerializeErrorName(SerializeError error)
	{
		switch (error)
		{
		case SerializeError::OK: return "OK";
		case SerializeError::FILE_NOT_OPENED: return "FILE_NOT_OPENED";
		case SerializeError::INVALID_PARAMETERS: return "INVALID_PARAMETERS";
		case SerializeError::WRITE_FAILED: return "WRITE_FAILED";
		case SerializeError::READ_FAILED: return "READ_FAILED";
		case SerializeError::INVALID_FORMAT: return "INVALID_FORMAT";
		case SerializeError::UNSUPPORTED_VERSION: return "UNSUPPORTED_VERSION";
		case SerializeError::TYPE_MISMATCH: return "TYPE_MISMATCH";
		case SerializeError::ENDIAN_MISMATCH: return "ENDIAN_MISMATCH";
		case SerializeError::ALLOCATION_LIMIT_EXCEEDED: return "ALLOCATION_LIMIT_EXCEEDED";
		case SerializeError::SCHEMA_MISMATCH: return "SCHEMA_MISMATCH";
		case SerializeError::MALFORMED_INPUT: return "MALFORMED_INPUT";
		case SerializeError::UNSUPPORTED_SCALAR: return "UNSUPPORTED_SCALAR";
		case SerializeError::TRUNCATED_INPUT: return "TRUNCATED_INPUT";
		default: return "UNKNOWN_SERIALIZE_ERROR";
		}
	}

	struct SerializeResult
	{
		bool success;
		SerializeError error;
		std::string message;
	};

	inline SerializeResult SerializeSuccess(std::string message = "Success")
	{
		return { true, SerializeError::OK, message };
	}

	inline SerializeResult SerializeFailure(SerializeError error, std::string message)
	{
		return { false, error, std::string(SerializeErrorName(error)) + ": " + message };
	}
}

#endif // MML_IO_RESULT_H
