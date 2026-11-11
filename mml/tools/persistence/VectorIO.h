///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        tools/persistence/VectorIO.h                                         ///
///  Description: Vector text persistence and transactional file/stream dispatch        ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_VECTOR_IO_H
#define MML_PERSISTENCE_VECTOR_IO_H

#include <mml/base/Vector/Vector.h>
#include <mml/tools/persistence/VectorBinary.h>
#include <mml/tools/persistence/VectorJSON.h>
#include <mml/tools/TextIOFormat.h>

#include <charconv>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>

namespace MML::Persistence {
	namespace Detail {
		template<class Type>
		constexpr bool VectorTextScalarSupported() {
			if constexpr (is_std_complex<Type>::value)
				return std::is_floating_point_v<typename Type::value_type>;
			else
				return std::is_floating_point_v<Type> || (std::is_integral_v<Type> && !std::is_same_v<Type, bool>);
		}

		template<class Type>
		IOResult ValidateVectorFormat(Format format) {
			switch (format) {
				case Format::Text:
					if constexpr (VectorTextScalarSupported<Type>()) return SerializeSuccess();
					else return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Vector text requires integer, real or complex floating-point scalars");
				case Format::JSON:
					if constexpr (std::is_floating_point_v<Type>) return SerializeSuccess();
					else return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Vector JSON requires real floating-point scalars");
				case Format::Binary:
					if constexpr (BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown) return SerializeSuccess();
					else return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Unsupported Vector binary scalar");
				default:
					return SerializeFailure(IOError::INVALID_FORMAT, "Vector persistence requires Text, JSON or Binary format");
			}
		}

		class VectorTextReader {
			std::istream& _input;
			std::int64_t _line = 1;
			std::int64_t _tokenLine = 1;

			static bool IsSpace(int c) {
				return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\v' || c == '\f';
			}
			void Consume() {
				if (_input.get() == '\n')
					++_line;
			}
		public:
			explicit VectorTextReader(std::istream& input) : _input(input) {}
			std::int64_t getTokenLine() const { return _tokenLine; }

			IOResult SkipWhitespace() {
				while (!_input.eof()) {
					const int c = _input.peek();
					if (c == std::char_traits<char>::eof() || !IsSpace(c))
						break;
					Consume();
				}
				if (_input.bad() || (_input.fail() && !_input.eof()))
					return SerializeFailure(IOError::READ_FAILED, "Vector text stream read failed").At(_line, IOOffsetKind::Line);
				return SerializeSuccess();
			}

			IOResult ReadToken(std::string& token, bool component = false) {
				const auto skipped = SkipWhitespace();
				if (!skipped) return skipped;
				_tokenLine = _line;
				if (_input.eof())
					return SerializeFailure(IOError::TRUNCATED_INPUT, "Missing Vector text value").At(_line, IOOffsetKind::Line);
				token.clear();
				while (!_input.eof()) {
					const int c = _input.peek();
					if (c == std::char_traits<char>::eof() || IsSpace(c) || (component && (c == ',' || c == ')')))
						break;
					token += static_cast<char>(c);
					Consume();
				}
				if (_input.bad())
					return SerializeFailure(IOError::READ_FAILED, "Vector text stream read failed").At(_line, IOOffsetKind::Line);
				if (token.empty())
					return SerializeFailure(IOError::MALFORMED_INPUT, "Empty Vector text value").At(_line, IOOffsetKind::Line);
				return SerializeSuccess();
			}

			IOResult Expect(char expected) {
				const auto skipped = SkipWhitespace();
				if (!skipped) return skipped;
				if (_input.eof())
					return SerializeFailure(IOError::TRUNCATED_INPUT, "Incomplete complex Vector value").At(_line, IOOffsetKind::Line);
				if (_input.peek() != expected)
					return SerializeFailure(IOError::MALFORMED_INPUT,
						"Expected '" + std::string(1, expected) + "' in complex Vector value").At(_line, IOOffsetKind::Line);
				Consume();
				return SerializeSuccess();
			}

			IOResult Finish(bool strict) {
				if (!strict) return SerializeSuccess();
				const auto skipped = SkipWhitespace();
				if (!skipped) return skipped;
				if (!_input.eof())
					return SerializeFailure(IOError::MALFORMED_INPUT, "Unexpected trailing Vector text data").At(_line, IOOffsetKind::Line);
				return SerializeSuccess();
			}
		};

		template<class Type>
		bool ParseVectorNumber(std::string_view token, Type& value) {
			if (!token.empty() && token.front() == '+') {
				token.remove_prefix(1);
				if (token.empty() || token.front() == '-' || token.front() == '+') return false;
			}
			if (token.empty()) return false;
			if constexpr (std::is_floating_point_v<Type>) {
				if (token == "nan") { value = std::numeric_limits<Type>::quiet_NaN(); return true; }
				if (token == "inf") { value = std::numeric_limits<Type>::infinity(); return true; }
				if (token == "-inf") { value = -std::numeric_limits<Type>::infinity(); return true; }
				if constexpr (std::is_same_v<Type, long double>) {
					// Streams may accept hex floats; retain the decimal-only from_chars grammar.
					return token.find_first_of("xXpP") == std::string_view::npos &&
						TextIO::ParseFloatClassic(std::string(token), value) && std::isfinite(value);
				} else {
					const auto parsed = std::from_chars(token.data(), token.data() + token.size(), value, std::chars_format::general);
					return parsed.ec == std::errc{} && parsed.ptr == token.data() + token.size() && std::isfinite(value);
				}
			} else {
				const auto parsed = std::from_chars(token.data(), token.data() + token.size(), value);
				return parsed.ec == std::errc{} && parsed.ptr == token.data() + token.size();
			}
		}

		template<class Type>
		void WriteVectorTextValue(std::ostream& out, const Type& value) {
			if constexpr (is_std_complex<Type>::value) {
				out << '(';
				WriteVectorTextValue(out, value.real());
				out << ',';
				WriteVectorTextValue(out, value.imag());
				out << ')';
			} else if constexpr (std::is_floating_point_v<Type>) {
				if (std::isnan(value)) out << "nan";
				else if (std::isinf(value)) out << (std::signbit(value) ? "-inf" : "inf");
				else out << value;
			} else if constexpr (std::is_signed_v<Type>)
				out << static_cast<std::intmax_t>(value);
			else
				out << static_cast<std::uintmax_t>(value);
		}

		template<class Type>
		IOResult ReadVectorTextValue(VectorTextReader& reader, Type& value) {
			std::string token;
			if constexpr (is_std_complex<Type>::value) {
				using Component = typename Type::value_type;
				auto result = reader.Expect('(');
				if (!result) return result;
				Component real{}, imaginary{};
				result = reader.ReadToken(token, true);
				if (!result) return result;
				if (!ParseVectorNumber(token, real))
					return SerializeFailure(IOError::MALFORMED_INPUT, "Invalid complex real component: " + token).At(reader.getTokenLine(), IOOffsetKind::Line);
				result = reader.Expect(',');
				if (!result) return result;
				result = reader.ReadToken(token, true);
				if (!result) return result;
				if (!ParseVectorNumber(token, imaginary))
					return SerializeFailure(IOError::MALFORMED_INPUT, "Invalid complex imaginary component: " + token).At(reader.getTokenLine(), IOOffsetKind::Line);
				result = reader.Expect(')');
				if (!result) return result;
				value = Type(real, imaginary);
			} else {
				const auto result = reader.ReadToken(token);
				if (!result) return result;
				if (!ParseVectorNumber(token, value))
					return SerializeFailure(IOError::MALFORMED_INPUT, "Invalid or out-of-range Vector value: " + token).At(reader.getTokenLine(), IOOffsetKind::Line);
			}
			return SerializeSuccess();
		}

		template<class Type>
		IOResult ReadVectorText(std::istream& in, Vector<Type>& candidate, const LoadOptions& options) {
			TextIO::ScopedTextFormat format(in);
			VectorTextReader reader(in);
			std::string token;
			auto result = reader.ReadToken(token);
			if (!result) return result;
			std::string_view countText(token);
			if (!countText.empty() && countText.front() == '+')
				countText.remove_prefix(1);
			std::uint64_t count{};
			const auto parsed = std::from_chars(countText.data(), countText.data() + countText.size(), count);
			if (parsed.ec == std::errc::result_out_of_range)
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector text count exceeds supported capacity").At(reader.getTokenLine(), IOOffsetKind::Line);
			if (countText.empty() || parsed.ec != std::errc{} || parsed.ptr != countText.data() + countText.size())
				return SerializeFailure(IOError::MALFORMED_INPUT, "Vector text count must be a nonnegative integer").At(reader.getTokenLine(), IOOffsetKind::Line);
			if (count > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
			    count > options.max_allocation_bytes / sizeof(Type))
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector text count exceeds dimension or allocation limit").At(reader.getTokenLine(), IOOffsetKind::Line);
			candidate.resize(static_cast<int>(count));
			for (int i = 0; i < candidate.size(); ++i) {
				result = ReadVectorTextValue(reader, candidate[i]);
				if (!result) return result;
			}
			return reader.Finish(options.strict_schema);
		}
	}

	template<class Type>
	IOResult Save(std::ostream& out, const Vector<Type>& vec, Format format, const SaveOptions& options = {}) {
		const auto validation = Detail::ValidateVectorFormat<Type>(format);
		if (!validation) return validation;
		if (format == Format::Text) {
			if constexpr (Detail::VectorTextScalarSupported<Type>()) {
				if (options.precision < 0)
					return SerializeFailure(IOError::INVALID_PARAMETERS, "Vector text precision must be zero or positive");
				TextIO::ScopedTextFormat scoped(out, options.precision == 0 ? TextIO::RoundTripPrecision<Type>() : options.precision);
				out << vec.size() << '\n';
				for (int i = 0; i < vec.size(); ++i) {
					Detail::WriteVectorTextValue(out, vec[i]);
					out << '\n';
				}
				return out ? SerializeSuccess() : SerializeFailure(IOError::WRITE_FAILED, "Vector text stream write failed");
			}
		} else if (format == Format::JSON) {
			if constexpr (std::is_floating_point_v<Type>) return SaveJson(out, vec, options);
		} else if (format == Format::Binary) {
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown) return SaveBinary(out, vec, options);
		}
		return SerializeFailure(IOError::INVALID_FORMAT, "Unhandled Vector persistence format");
	}

	template<class Type>
	IOResult Load(std::istream& in, Vector<Type>& vec, Format format, const LoadOptions& options = {}) {
		const auto validation = Detail::ValidateVectorFormat<Type>(format);
		if (!validation) return validation;
		if (in.fail())
			return SerializeFailure(IOError::READ_FAILED, "Vector input stream is already failed");
		try {
			Vector<Type> candidate;
			IOResult result;
			if (format == Format::Text) {
				if constexpr (Detail::VectorTextScalarSupported<Type>()) result = Detail::ReadVectorText(in, candidate, options);
			} else if (format == Format::JSON) {
				if constexpr (std::is_floating_point_v<Type>) result = LoadJson(in, candidate, options);
			} else if (format == Format::Binary) {
				if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown) result = LoadBinary(in, candidate, options);
			}
			if (!result) return result;
			if (in.bad())
				return SerializeFailure(IOError::READ_FAILED, "Vector input stream read failed");
			vec = std::move(candidate);
			return result;
		} catch (const std::bad_alloc&) {
			return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector persistence allocation failed");
		} catch (const std::length_error&) {
			return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector persistence allocation exceeds container capacity");
		}
	}

	template<class Type>
	IOResult Save(const Vector<Type>& vec, const std::string& filename, const SaveOptions& options, Format format = Format::Auto) {
		const Format resolved = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		const auto validation = Detail::ValidateVectorFormat<Type>(resolved);
		if (!validation) return validation.WithFile(filename);
		const auto pathValidation = Detail::ValidateSaveFilename(filename);
		if (!pathValidation.success) return pathValidation;
		if (resolved == Format::JSON) {
			if constexpr (std::is_floating_point_v<Type>) return SaveJson(vec, filename, options);
		}
		if (resolved == Format::Binary) {
			if constexpr (Detail::BinaryScalarTypeFor<Type>() != BinaryScalarType::Unknown) return SaveBinary(vec, filename, options);
		}
		if (resolved != Format::Binary && options.precision < 0)
			return SerializeFailure(IOError::INVALID_PARAMETERS, "Vector text precision must be zero or positive").WithFile(filename);
		std::ofstream file(filename, std::ios::out | (resolved == Format::Binary ? std::ios::binary : std::ios::openmode{}));
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not create vector file").WithFile(filename);
		return FinishOutputFile(file, filename, Save(file, vec, resolved, options));
	}

	template<class Type>
	IOResult Save(const Vector<Type>& vec, const std::string& filename, Format format = Format::Auto) {
		return Save(vec, filename, SaveOptions{}, format);
	}

	template<class Type>
	IOResult Load(const std::string& filename, Vector<Type>& vec, const LoadOptions& options, Format format = Format::Auto) {
		const Format resolved = format == Format::Auto ? DetectFormatFromPath(filename) : format;
		const auto validation = Detail::ValidateVectorFormat<Type>(resolved);
		if (!validation) return validation.WithFile(filename);
		if (filename.empty() || filename.find('\0') != std::string::npos)
			return SerializeFailure(IOError::INVALID_PARAMETERS, "Invalid vector filename").WithFile(filename);
		std::ifstream file(filename, std::ios::in | (resolved == Format::Binary ? std::ios::binary : std::ios::openmode{}));
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not open vector file").WithFile(filename);
		return InputResult(file, Load(file, vec, resolved, options), filename);
	}

	template<class Type>
	IOResult Load(const std::string& filename, Vector<Type>& vec, Format format = Format::Auto) {
		return Load(filename, vec, LoadOptions{}, format);
	}
}

#endif // MML_PERSISTENCE_VECTOR_IO_H
