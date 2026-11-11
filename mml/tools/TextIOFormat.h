///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        tools/TextIOFormat.h                                                ///
///  Description: Locale-independent numeric text I/O helpers shared by every         ///
///               MML text reader and writer (serializer, persistence, data loader)   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_TOOLS_TEXT_IO_FORMAT_H
#define MML_TOOLS_TEXT_IO_FORMAT_H

#include <mml/MMLBase.h>

#include <ios>
#include <complex>
#include <limits>
#include <locale>
#include <sstream>
#include <string>

namespace MML::TextIO
{
	template<class Type>
	constexpr int RoundTripPrecision() {
		if constexpr (is_std_complex<Type>::value)
			return std::numeric_limits<typename Type::value_type>::max_digits10;
		else
			return std::numeric_limits<Type>::max_digits10;
	}

	/// @brief Apply canonical text formatting without leaking it to a caller's stream.
	class ScopedTextFormat {
		std::ios& _stream;
		std::locale _locale;
		std::ios::fmtflags _flags;
		std::streamsize _precision;
		std::streamsize _width;
		char _fill;
	public:
		explicit ScopedTextFormat(std::ios& stream, int precision = RoundTripPrecision<Real>())
			: _stream(stream), _locale(stream.getloc()), _flags(stream.flags()),
			  _precision(stream.precision()), _width(stream.width()), _fill(stream.fill()) {
			stream.imbue(std::locale::classic());
			stream.flags(std::ios::skipws | std::ios::dec);
			stream.precision(precision);
			stream.width(0);
			stream.fill(' ');
		}
		~ScopedTextFormat() {
			_stream.imbue(_locale);
			_stream.flags(_flags);
			_stream.precision(_precision);
			_stream.width(_width);
			_stream.fill(_fill);
		}
		ScopedTextFormat(const ScopedTextFormat&) = delete;
		ScopedTextFormat& operator=(const ScopedTextFormat&) = delete;
	};

	/// @brief Put a text stream into the one format MML files are defined in: the classic "C"
	/// locale (a '.' decimal point, no digit grouping, regardless of the process/global locale)
	/// and, by default, round-trip precision for Real.
	///
	/// Every MML text writer and reader must call this on its stream. A process that has called
	/// std::locale::global(std::locale("de_DE")) otherwise writes "3,14" and reads "3.14" as 3
	/// (analysis v4, TOOLS T1); std::cout-style default precision (6 digits) otherwise truncates
	/// doubles to ~1e-6 relative error on the way to disk (TOOLS T2).
	inline void ApplyClassicFormat(std::ios& stream, int precision = std::numeric_limits<Real>::max_digits10)
	{
		stream.imbue(std::locale::classic());
		stream.precision(precision);
	}

	/// @brief Locale-independent parse of a whole string as a floating-point number.
	/// @return true if the entire string (modulo surrounding whitespace) is a valid, in-range number.
	/// @note std::stod / strtod honour the C locale and would reject "3.14" under a comma-decimal
	///       locale; this always parses the canonical '.' form.
	template<typename FloatType>
	inline bool ParseFloatClassic(const std::string& text, FloatType& out)
	{
		std::istringstream in(text);
		in.imbue(std::locale::classic());
		FloatType value{};
		in >> value;
		if (in.fail())
			return false;
		in >> std::ws;
		if (!in.eof())
			return false;
		out = value;
		return true;
	}

	inline bool ParseReal(const std::string& text, Real& out) { return ParseFloatClassic<Real>(text, out); }
}

#endif // MML_TOOLS_TEXT_IO_FORMAT_H
