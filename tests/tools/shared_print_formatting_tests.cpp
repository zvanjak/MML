#include <catch2/catch_all.hpp>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/MatrixTriDiag.h>
#include <mml/tools/ConsolePrinter.h>

#include <array>
#include <functional>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <streambuf>

using namespace MML;

TEST_CASE("Dynamic fixed and structured matrices share every format option", "[PrintFormatting][Matrix]")
{
	const Matrix<Real> dynamic(2, 2, {1, 2, 2, 3});
	const MatrixNM<Real, 2, 2> fixed{1, 2, 2, 3};
	const MatrixSym<Real> symmetric(2, {1, 2, 3});
	const TridiagonalMatrix<Real> tridiagonal(2, {1, 2, 2, 3});
	for (bool compact : {false, true})
		for (bool brackets : {false, true})
			for (bool scientific : {false, true}) {
				MatrixPrintFormat format;
				format.showHeader = false;
				format.compactMode = compact;
				format.showBrackets = brackets;
				format.scientific = scientific;
				format.width = 0;
				format.precision = 2;
				format.delimiter = "|";
				format.zeroThreshold = 2;
				std::ostringstream a, b, c, d;
				dynamic.Print(a, format);
				fixed.Print(b, format);
				symmetric.Print(c, format);
				tridiagonal.Print(d, format);
				REQUIRE(a.str() == b.str());
				REQUIRE(a.str() == c.str());
				REQUIRE(a.str() == d.str());
			}
	MatrixPrintFormat compact = MatrixPrintFormat::Compact();
	compact.width = 0;
	compact.zeroThreshold = 2;
	std::ostringstream output;
	dynamic.Print(output, compact);
	REQUIRE(output.str() == "[0.00, 0.00; 0.00, 3.00]\n");
}

TEST_CASE("Shared printing restores caller formatting on success and stream exceptions", "[PrintFormatting][StreamState]")
{
	class ThrowingBuffer : public std::streambuf {
		std::streamsize xsputn(const char*, std::streamsize) override { throw std::runtime_error("print failed"); }
		int_type overflow(int_type) override { throw std::runtime_error("print failed"); }
	};
	const Matrix<Real> matrix(1, 2, {1, 2});
	const Vector<Real> vector{1, 2};
	const VectorN<Real, 2> fixed{1, 2};
	const MatrixNM<Real, 1, 2> fixedMatrix{1, 2};
	const MatrixSym<Real> symmetric(2, {1, 2, 3});
	const TridiagonalMatrix<Real> tridiagonal(2, {1, 2, 2, 3});
	const std::array<std::function<void(std::ostream&)>, 7> printers{{
		[&](std::ostream& out) { matrix.Print(out); },
		[&](std::ostream& out) { vector.Print(out, 5, 3); },
		[&](std::ostream& out) { fixed.Print(out, 5, 3); },
		[&](std::ostream& out) { fixedMatrix.Print(out); },
		[&](std::ostream& out) { symmetric.Print(out); },
		[&](std::ostream& out) { tridiagonal.Print(out); },
		[&](std::ostream& out) { vector.PrintCol(out, 5, 3); }
	}};
	for (const auto& print : printers)
		for (bool fail : {false, true}) {
			ThrowingBuffer buffer;
			std::ostringstream successful;
			std::ostream failing(&buffer);
			std::ostream& out = fail ? failing : successful;
			out << std::scientific << std::showpos;
			out.precision(7); out.width(13); out.fill('@');
			out.exceptions(std::ios::badbit);
			const auto flags = out.flags();
			const auto locale = out.getloc();
			if (fail) REQUIRE_THROWS(print(out));
			else print(out);
			REQUIRE(out.flags() == flags);
			REQUIRE(out.precision() == 7);
			REQUIRE(out.width() == 13);
			REQUIRE(out.fill() == '@');
			REQUIRE(out.getloc() == locale);
		}
}

TEST_CASE("Shared vector printing preserves layouts and nonfinite values", "[PrintFormatting][Vector]")
{
	Vector<Real> dynamic{1, 2};
	VectorN<Real, 2> fixed{1, 2};
	std::ostringstream a, b;
	dynamic.Print(a, 0, 3);
	fixed.Print(b, 0, 3);
	REQUIRE(a.str() == "[1, 2]");
	REQUIRE(a.str() == b.str());
	Vector<Real> nonfinite{std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()};
	std::ostringstream result;
	nonfinite.Print(result, 0, 3, Real(1));
	REQUIRE(result.str().find("nan") != std::string::npos);
	REQUIRE(result.str().find("inf") != std::string::npos);
	const auto column = ColumnFormat("value").precision(3).format(FormatType::Fixed);
	REQUIRE(column.formatValue(Real(1) / 4) == "0.250");
	REQUIRE(ColumnFormat("value").format(FormatType::Hexadecimal).formatValue(255) == "ff");
}
