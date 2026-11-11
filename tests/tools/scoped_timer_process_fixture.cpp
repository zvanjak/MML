#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/Timer.h>
#endif

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <streambuf>
#include <string>

namespace {
	enum class Fault { None, Write, ShortWrite, Flush, FlushState, UnknownWrite };

	class OutputBuffer : public std::streambuf {
	public:
		Fault fault = Fault::None;
		std::string text;
	protected:
		std::streamsize xsputn(const char* data, std::streamsize count) override {
			if (fault == Fault::Write) throw std::runtime_error("fixture write failure");
			if (fault == Fault::UnknownWrite) throw 123;
			if (fault == Fault::ShortWrite) return 0;
			text.append(data, static_cast<std::size_t>(count));
			return count;
		}
		int_type overflow(int_type value) override {
			if (traits_type::eq_int_type(value, traits_type::eof())) return traits_type::not_eof(value);
			const char character = traits_type::to_char_type(value);
			return xsputn(&character, 1) == 1 ? value : traits_type::eof();
		}
		int sync() override {
			if (fault == Fault::Flush) throw std::runtime_error("fixture flush failure");
			return fault == Fault::FlushState ? -1 : 0;
		}
	};

	struct RestoreOutput {
		std::streambuf* buffer = std::cout.rdbuf();
		std::ios::iostate state = std::cout.rdstate();
		std::ios::iostate exceptions = std::cout.exceptions();
		~RestoreOutput() {
			std::cout.exceptions(std::ios::goodbit);
			std::cout.rdbuf(buffer);
			std::cout.clear(state);
			std::cout.exceptions(exceptions);
		}
	};
}

int main(int argc, char** argv) {
	if (argc != 3) return 91;
	std::set_terminate([] { std::_Exit(86); });
#ifdef _WIN32
	FILE* redirected = nullptr;
	if (freopen_s(&redirected, argv[2], "w", stderr) != 0) return 92;
#else
	if (!std::freopen(argv[2], "w", stderr)) return 92;
#endif
	const std::string mode = argv[1];
	OutputBuffer buffer;
	RestoreOutput restore;
	std::cout.rdbuf(&buffer);
	if (mode == "short-write") buffer.fault = Fault::ShortWrite;
	else if (mode == "flush-state") buffer.fault = Fault::FlushState;
	else if (mode == "throw-flush") buffer.fault = Fault::Flush;
	else if (mode == "throw-unknown") buffer.fault = Fault::UnknownWrite;
	else if (mode == "throw-write" || mode == "callback-writes" || mode == "both" ||
		mode == "unwind-output" || mode == "unwind-both") buffer.fault = Fault::Write;
	else if (mode != "bad-state" && mode != "silent-bad" && mode != "callback-std" &&
		mode != "callback-unknown" && mode != "callback-print") return 91;
	if (buffer.fault == Fault::Write || buffer.fault == Fault::UnknownWrite || buffer.fault == Fault::Flush)
		std::cout.exceptions(std::ios::badbit | std::ios::failbit);
	if (mode == "bad-state" || mode == "silent-bad") std::cout.setstate(std::ios::badbit);
	const auto configuredMask = std::cout.exceptions();
	const auto configuredState = std::cout.rdstate();

	const bool callback = mode.find("callback") != std::string::npos || mode == "both" || mode == "unwind-both";
	const bool unwinding = mode.find("unwind") != std::string::npos;
	int calls = 0;
	double observed = -1;
	bool primaryPreserved = false;
	const auto before = std::chrono::steady_clock::now();
	if (callback) {
		auto action = [&](const std::string& name, double elapsed) {
			++calls;
			observed = elapsed;
			if (name != "fixture") throw std::logic_error("wrong callback name");
			if (mode == "callback-unknown") throw 123;
			if (mode == "callback-writes") std::cout << "callback output";
			throw std::runtime_error("fixture callback failure");
		};
		try {
			MML::ScopedTimer timer("fixture", action, mode != "both" && mode != "unwind-both" && mode != "callback-print");
			if (unwinding) throw std::logic_error("unrelated operation failure");
		} catch (const std::logic_error& error) {
			primaryPreserved = std::string(error.what()) == "unrelated operation failure";
		}
	} else {
		try {
			MML::ScopedTimer timer("fixture", observed, mode == "silent-bad");
			if (unwinding) throw std::logic_error("unrelated operation failure");
		} catch (const std::logic_error& error) {
			primaryPreserved = std::string(error.what()) == "unrelated operation failure";
		}
	}
	const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - before).count();
	if (!std::isfinite(observed) || observed < 0 || observed > total || calls != (callback ? 1 : 0) ||
		primaryPreserved != unwinding) return 93;
	if (mode == "callback-print" && buffer.text.find("fixture: ") == std::string::npos) return 94;
	if (std::cout.exceptions() != configuredMask) return 96;
	if ((buffer.fault != Fault::None || mode == "bad-state") && !std::cout.fail()) return 97;
	if (mode == "silent-bad" && std::cout.rdstate() != configuredState) return 98;
	if (std::fflush(stderr) != 0) return 95;
	return 0;
}
