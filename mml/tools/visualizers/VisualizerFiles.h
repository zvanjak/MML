///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        tools/visualizers/VisualizerFiles.h                                  ///
///  Description: Visualizer input resolution and contained output-file paths          ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VISUALIZER_FILES_H
#define MML_VISUALIZER_FILES_H

#include <mml/tools/IOResult.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace MML::detail::visualizer {
	inline std::string SanitizeFilename(const std::string& filename) {
		const std::string basename = std::filesystem::path(filename).filename().string();
		std::string result;
		result.reserve(basename.size());
		for (unsigned char c : basename)
			result += std::isalnum(c) || c == '_' || c == '-' || c == '.' ? static_cast<char>(c) : '_';
		const auto first = result.find_first_not_of('.');
		return first == std::string::npos ? "output" : result.substr(first);
	}

	inline IOResult MakeOutputPath(const std::string& filename, const std::string& resultsRoot,
	                               const std::string& suffix = {}) {
		if (resultsRoot.empty() || filename.find('\0') != std::string::npos ||
		    resultsRoot.find('\0') != std::string::npos || suffix.find('\0') != std::string::npos)
			return IOResult::Failure("Invalid visualizer output path", filename, -1, IOError::INVALID_PARAMETERS);

		const std::filesystem::path directory(resultsRoot);
		std::error_code error;
		std::filesystem::create_directories(directory, error);
		if (error)
			return IOResult::Failure("Cannot create visualizer results directory: " + error.message(),
			                         directory.string(), error.value(), IOError::FILE_NOT_OPENED);
		const auto canonicalDirectory = std::filesystem::weakly_canonical(directory, error);
		if (error)
			return IOResult::Failure("Cannot resolve visualizer results directory: " + error.message(),
			                         directory.string(), error.value(), IOError::FILE_NOT_OPENED);
		const auto output = directory / (SanitizeFilename(filename) + suffix);
		const auto canonicalOutput = std::filesystem::weakly_canonical(output, error);
		if (error)
			return IOResult::Failure("Cannot resolve visualizer output path: " + error.message(),
			                         output.string(), error.value(), IOError::FILE_NOT_OPENED);
		const auto mismatch = std::mismatch(canonicalDirectory.begin(), canonicalDirectory.end(),
		                                    canonicalOutput.begin(), canonicalOutput.end());
		if (mismatch.first != canonicalDirectory.end())
			return IOResult::Failure("Visualizer output path escapes results directory", output.string(),
			                         -1, IOError::INVALID_PARAMETERS);
		return IOResult::Success(canonicalOutput.string());
	}

	inline IOResult ResolveInputFile(const std::string& filename, const std::string& resultsRoot) {
		if (filename.empty() || filename.find('\0') != std::string::npos)
			return IOResult::Failure("Invalid visualizer input filename", filename, -1, IOError::INVALID_PARAMETERS);

		std::string candidate = filename;
		for (int attempt = 0; attempt < 2; ++attempt) {
			if (attempt == 1) {
				const auto fallback = MakeOutputPath(filename, resultsRoot);
				if (!fallback)
					return fallback;
				candidate = fallback.filename;
			}
			std::error_code error;
			const auto status = std::filesystem::status(candidate, error);
			if (error && error != std::errc::no_such_file_or_directory && error != std::errc::not_a_directory)
				return IOResult::Failure("Cannot inspect visualizer input: " + error.message(), candidate,
				                         error.value(), IOError::FILE_NOT_OPENED);
			if (!std::filesystem::exists(status))
				continue;
			if (!std::filesystem::is_regular_file(status))
				return IOResult::Failure("Visualizer input is not a regular file", candidate, -1, IOError::FILE_NOT_OPENED);
			const auto canonical = std::filesystem::weakly_canonical(candidate, error);
			if (error)
				return IOResult::Failure("Cannot resolve visualizer input: " + error.message(), candidate,
				                         error.value(), IOError::FILE_NOT_OPENED);
			return IOResult::Success(canonical.string());
		}
		return IOResult::Failure("Data file not found: " + filename, filename, -1, IOError::FILE_NOT_OPENED);
	}

	inline IOResult ResolveInputFiles(const std::vector<std::string>& filenames, const std::string& resultsRoot,
	                                  std::vector<std::string>& resolved) {
		if (filenames.empty())
			return IOResult::Failure("No visualizer input files specified", {}, -1, IOError::INVALID_PARAMETERS);
		std::vector<std::string> paths;
		paths.reserve(filenames.size());
		for (const auto& filename : filenames) {
			const auto result = ResolveInputFile(filename, resultsRoot);
			if (!result)
				return result;
			paths.push_back(result.filename);
		}
		resolved = std::move(paths);
		return IOResult::Success(resolved.front());
	}
}

#endif // MML_VISUALIZER_FILES_H
