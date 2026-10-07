///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        tools/persistence/PersistenceBase.h                                 ///
///  Description: Common types and options for durable object persistence             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_PERSISTENCE_BASE_H
#define MML_PERSISTENCE_PERSISTENCE_BASE_H

#include <mml/tools/IOResult.h>

#include <cstddef>
#include <map>
#include <string>

namespace MML
{
namespace Persistence
{
	enum class Format
	{
		Auto,
		Text,
		CSV,
		JSON,
		Binary
	};

	struct Metadata
	{
		std::string title;
		std::string description;
		std::string units;
		std::string coordinate_system;
		std::map<std::string, std::string> tags;
	};

	struct SaveOptions
	{
		int precision = 17;
		bool include_metadata = true;
		bool pretty_json = true;
		bool allow_lossy = false;
		Metadata metadata;
	};

	struct LoadOptions
	{
		bool strict_schema = true;
		bool allow_scalar_conversion = false;
		std::size_t max_allocation_bytes = 1ull << 30;
	};

	inline char ToLowerAscii(char c)
	{
		return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c;
	}

	inline std::string FileExtensionLower(const std::string& path)
	{
		const size_t slashPos = path.find_last_of("/\\");
		const size_t dotPos = path.find_last_of('.');
		if (dotPos == std::string::npos || (slashPos != std::string::npos && dotPos < slashPos))
			return "";

		std::string extension = path.substr(dotPos + 1);
		for (char& character : extension)
			character = ToLowerAscii(character);
		return extension;
	}

	inline Format DetectFormatFromPath(const std::string& path)
	{
		const std::string extension = FileExtensionLower(path);
		if (extension == "mml" || extension == "txt" || extension == "dat")
			return Format::Text;
		if (extension == "csv")
			return Format::CSV;
		if (extension == "mmlj")
			return Format::JSON;
		if (extension == "mmlb" || extension == "bin")
			return Format::Binary;
		return Format::Auto;
	}
}
}

#endif // MML_PERSISTENCE_PERSISTENCE_BASE_H
