#include <filesystem>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>

int main(int argc, char** argv)
{
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
		if (token != "MML_REAL_FUNCTION")
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
