#if !defined MML_PRINT_FORMATTING_H
#define MML_PRINT_FORMATTING_H

#include <mml/MMLBase.h>
#include <mml/base/Matrix/MatrixPrintFormat.h>

#include <iomanip>
#include <cmath>
#include <ios>
#include <locale>
#include <ostream>
#include <sstream>
#include <string>

namespace MML::PrintFormatting {
	class StreamState {
		std::ostream& _stream;
		std::ios::fmtflags _flags;
		std::streamsize _precision, _width;
		char _fill;
	public:
		explicit StreamState(std::ostream& stream)
			: _stream(stream), _flags(stream.flags()), _precision(stream.precision()),
			  _width(stream.width()), _fill(stream.fill()) {
			stream.width(0);
		}
		~StreamState() {
			_stream.flags(_flags);
			_stream.precision(_precision);
			_stream.width(_width);
			_stream.fill(_fill);
		}
		StreamState(const StreamState&) = delete;
		StreamState& operator=(const StreamState&) = delete;
	};

	enum class Notation { Inherit, General, Fixed, Scientific, Hexadecimal, Integer };

	inline void ApplyNotation(std::ostream& out, Notation notation) {
		switch (notation) {
			case Notation::General: out << std::defaultfloat; break;
			case Notation::Fixed: out << std::fixed; break;
			case Notation::Scientific: out << std::scientific; break;
			case Notation::Hexadecimal: out << std::hex; break;
			case Notation::Inherit:
			case Notation::Integer: break;
		}
	}

	template<class Type, class Threshold = Real>
	void Value(std::ostream& out, const Type& value, int width, int precision,
	           Threshold threshold = 0, bool replaceZero = false) {
		out << std::setw(width) << std::setprecision(precision);
		if constexpr (std::is_arithmetic_v<Type> || is_std_complex<Type>::value) {
			bool zero = false;
			if (threshold > 0 || replaceZero) {
				if constexpr (std::is_integral_v<Type>)
					zero = std::abs(static_cast<Real>(value)) <= threshold;
				else
					zero = Abs(value) <= threshold;
			}
			if (zero) {
				if (replaceZero) out << Real(0);
				else out << Type{0};
				return;
			}
		}
		out << value;
	}

	template<class Getter>
	std::ostream& Vector(std::ostream& out, int count, Getter get, int width, int precision,
	                     Real threshold = 0, bool replaceZero = false) {
		StreamState state(out);
		out << '[';
		for (int i = 0; i < count; ++i) {
			if (i) out << ", ";
			Value(out, get(i), width, precision, threshold, replaceZero);
		}
		out << ']';
		return out;
	}

	template<class Getter>
	std::ostream& Column(std::ostream& out, int count, Getter get, int width, int precision) {
		StreamState state(out);
		for (int i = 0; i < count; ++i) {
			Value(out, get(i), width, precision);
			out << '\n';
		}
		return out;
	}

	template<class Getter>
	void Matrix(std::ostream& out, int rows, int cols, Getter get, const MatrixPrintFormat& format,
	            const char* kind = "") {
		StreamState state(out);
		if (rows == 0 || cols == 0) {
			if (format.showHeader)
				out << "Rows: " << rows << " Cols: " << cols << " - Empty matrix\n";
			return;
		}
		if (format.showHeader)
			out << "Rows: " << rows << " Cols: " << cols << kind << '\n';
		ApplyNotation(out, format.scientific ? Notation::Scientific : format.fixed ? Notation::Fixed : Notation::Inherit);
		if (format.compactMode && rows <= 3 && cols <= 3) {
			if (format.showBrackets) out << '[';
			for (int i = 0; i < rows; ++i) {
				if (i) out << "; ";
				for (int j = 0; j < cols; ++j) {
					if (j) out << format.delimiter;
					Value(out, get(i, j), format.width, format.precision, format.zeroThreshold);
				}
			}
			if (format.showBrackets) out << ']';
			out << '\n';
		} else {
			for (int i = 0; i < rows; ++i) {
				if (format.showBrackets) out << "[ ";
				for (int j = 0; j < cols; ++j) {
					if (j) out << format.delimiter;
					Value(out, get(i, j), format.width, format.precision, format.zeroThreshold);
				}
				if (format.showBrackets) out << " ]";
				if (i < rows - 1) out << '\n';
			}
		}
	}

	template<class Type>
	std::string Text(const Type& value, int precision, Notation notation) {
		std::ostringstream out;
		ApplyNotation(out, notation);
		if (notation == Notation::Integer) {
			if constexpr (std::is_floating_point_v<Type>)
				// Round before formatting, without narrowing to any integer storage type.
				out << std::fixed << std::setprecision(0) << std::round(value);
			else if constexpr (std::is_integral_v<Type>)
				out << std::dec << +value;
			else
				throw ArgumentError("PrintFormatting::Text: Integer requires an integral or floating-point scalar");
		} else if (notation == Notation::Hexadecimal) {
			if constexpr (std::is_integral_v<Type>)
				out << +value;
			else
				throw ArgumentError("PrintFormatting::Text: Hexadecimal requires an integral scalar");
		} else
			Value(out, value, 0, precision);
		return out.str();
	}
}

#endif // MML_PRINT_FORMATTING_H
