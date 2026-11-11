#include <filesystem>
#include <fstream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <thread>
#ifndef _WIN32
#include <csignal>
#include <unistd.h>
#endif

int main(int argc, char** argv)
{
#ifndef _WIN32
	if (argc > 1 && std::string(argv[1]) == "-W") {
		int firstArgument = 1;
		while (firstArgument < argc && std::string(argv[firstArgument]) != "--args") ++firstArgument;
		if (firstArgument == argc) return 94;
		argc -= firstArgument;
		argv += firstArgument;
	}
	if (argc > 2 && (std::string(argv[1]) == "--sleep-report" ||
		std::string(argv[1]) == "--terminate-report" || std::string(argv[1]) == "--exit-after")) {
		if (std::signal(SIGTERM, std::string(argv[1]) == "--sleep-report" ? SIG_IGN : SIG_DFL) == SIG_ERR) return 95;
		std::ofstream report(argv[2]);
		report << getpid() << '\n';
		report.close();
		if (!report) return 96;
		std::this_thread::sleep_for(std::chrono::milliseconds(argc > 3 ? std::stoi(argv[3]) : 5000));
		return 0;
	}
	if (argc > 1 && std::string(argv[1]) == "--signal") {
		if (std::signal(SIGTERM, SIG_DFL) == SIG_ERR) return 95;
		return std::raise(SIGTERM);
	}
	if (argc > 2 && std::string(argv[1]) == "--environment") {
		std::ofstream report(argv[2]);
		for (const auto* key : {"LD_LIBRARY_PATH", "QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH", "PATH"}) {
			const auto* value = std::getenv(key);
			report << (value ? value : "<unset>") << '\n';
		}
		report.close();
		return report ? 0 : 96;
	}
#endif
	if (argc < 2)
		return 91;
	if (std::string(argv[1]) == "--sleep") {
		std::this_thread::sleep_for(std::chrono::seconds(5));
		return 0;
	}
	int exitCode = 0;
	std::ifstream input(argv[1]);
	if (!(input >> exitCode)) {
		input.clear();
		input.seekg(0);
		std::string token;
		std::getline(input, token);
		if (token != "MML_REAL_FUNCTION" && token != "MML_VECTOR_FIELD_2D_CARTESIAN" &&
			token != "MML_VECTOR_FIELD_3D_CARTESIAN")
			return 92;
		exitCode = 0;
	}
	std::ofstream report(std::string(argv[1]) + ".argv");
	report << argc - 1 << '\n' << std::filesystem::current_path().string() << '\n';
	for (int i = 1; i < argc; ++i)
		report << argv[i] << '\n';
	report.close();
	if (!report)
		return 93;
	return exitCode;
}
