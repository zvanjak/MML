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
#include <cstdint>
#include <fstream>
#include <utility>

namespace MML
{
	enum class IOError {
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
		TRUNCATED_INPUT,
		PROCESS_FAILED,
		TIMEOUT
	};

	using SerializeError [[deprecated("use IOError")]] = IOError;

	inline const char* IOErrorName(IOError error) {
		switch (error)
		{
		case IOError::OK:
			return "OK";
		case IOError::FILE_NOT_OPENED:
			return "FILE_NOT_OPENED";
		case IOError::INVALID_PARAMETERS:
			return "INVALID_PARAMETERS";
		case IOError::WRITE_FAILED:
			return "WRITE_FAILED";
		case IOError::READ_FAILED:
			return "READ_FAILED";
		case IOError::INVALID_FORMAT:
			return "INVALID_FORMAT";
		case IOError::UNSUPPORTED_VERSION:
			return "UNSUPPORTED_VERSION";
		case IOError::TYPE_MISMATCH:
			return "TYPE_MISMATCH";
		case IOError::ENDIAN_MISMATCH:
			return "ENDIAN_MISMATCH";
		case IOError::ALLOCATION_LIMIT_EXCEEDED:
			return "ALLOCATION_LIMIT_EXCEEDED";
		case IOError::SCHEMA_MISMATCH:
			return "SCHEMA_MISMATCH";
		case IOError::MALFORMED_INPUT:
			return "MALFORMED_INPUT";
		case IOError::UNSUPPORTED_SCALAR:
			return "UNSUPPORTED_SCALAR";
		case IOError::TRUNCATED_INPUT:
			return "TRUNCATED_INPUT";
		case IOError::PROCESS_FAILED:
			return "PROCESS_FAILED";
		case IOError::TIMEOUT:
			return "TIMEOUT";
		default: return "UNKNOWN_SERIALIZE_ERROR";
		}
	}

	enum class IOOffsetKind { Unknown, Byte, Line };

	struct [[nodiscard]] IOResult {
		bool success = false;
		IOError error = IOError::INVALID_PARAMETERS;
		std::string message{};
		std::string filename{};
		std::int64_t offset = -1;
		IOOffsetKind offset_kind = IOOffsetKind::Unknown;
		bool timedOut = false;
		int exitCode = 0;
		bool dataOnly = false;

		explicit operator bool() const noexcept { return success; }
		bool IsSuccess() const noexcept { return success; }
		bool isDataOnly() const noexcept { return dataOnly; }

		IOResult WithFile(const std::string& path) const {
			IOResult result = *this;
			result.filename = path;
			return result;
		}

		IOResult At(std::int64_t position, IOOffsetKind kind = IOOffsetKind::Byte) const {
			IOResult result = *this;
			if (position >= 0) {
				result.offset = position;
				result.offset_kind = kind;
			}
			return result;
		}

		static IOResult Success(const std::string& path = "") { return {true, IOError::OK, "", path}; }
		static IOResult Failure(const std::string& message, const std::string& path = "", int code = -1,
														IOError error = IOError::PROCESS_FAILED) {
			IOResult result{false, error, message, path};
			result.exitCode = code;
			return result;
		}
		static IOResult Timeout(const std::string& message, const std::string& path = "") {
			IOResult result = Failure(message, path);
			result.error = IOError::TIMEOUT;
			result.timedOut = true;
			return result;
		}
		static IOResult DataSaved(const std::string& path) {
			IOResult result = Success(path);
			result.message = "Data saved (no visualizer launched)";
			result.dataOnly = true;
			return result;
		}
	};

	using SerializeResult [[deprecated("use IOResult")]] = IOResult;
	using VisualizerResult [[deprecated("use IOResult")]] = IOResult;

	inline IOResult SerializeSuccess(std::string message = "Success") { return {true, IOError::OK, std::move(message)}; }

	inline IOResult SerializeFailure(IOError error, std::string message) {
		return {false, error, std::string(IOErrorName(error)) + ": " + message};
	}

	[[deprecated("use IOErrorName")]] inline const char* SerializeErrorName(IOError error) { return IOErrorName(error); }

	inline IOResult InputResult(std::istream& input, IOResult result, const std::string& filename = "") {
		result.filename = filename;
		if (!result.success && result.offset < 0) {
			const auto position = input.rdbuf()->pubseekoff(0, std::ios::cur, std::ios::in);
			result = result.At(static_cast<std::int64_t>(position));
		}
		return result;
	}

	inline IOResult FinishOutputFile(std::ofstream& file, const std::string& filename, IOResult result) {
		result.filename = filename;
		const auto position = file.rdbuf()->pubseekoff(0, std::ios::cur, std::ios::out);
		file.flush();
		const bool writeFailed = !file;
		file.close();
		if (result.success && (writeFailed || file.fail()))
			result = SerializeFailure(IOError::WRITE_FAILED, "Writing or closing file failed").WithFile(filename);
		if (!result.success && result.offset < 0)
			result = result.At(static_cast<std::int64_t>(position));
		return result;
	}
}

#endif // MML_IO_RESULT_H
