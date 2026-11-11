#if !defined MML_EXAMPLES_VISUALIZATION_RESULT_H
#define MML_EXAMPLES_VISUALIZATION_RESULT_H

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/tools/IOResult.h>
#endif

#include <iostream>
#include <string_view>

namespace MML::Examples {
	inline void ReportVisualizationResult(const IOResult& result, std::string_view operation,
	                                      std::ostream& errors = std::cerr) {
		if (result.success)
			return;
		errors << operation << " failed [" << IOErrorName(result.error) << "]: " << result.message << '\n';
		if (!result.filename.empty())
			errors << "  File: " << result.filename << '\n';
		if (result.exitCode != 0)
			errors << "  Exit code: " << result.exitCode << '\n';
		if (result.timedOut)
			errors << "  Timed out: yes\n";
	}
}

#endif // MML_EXAMPLES_VISUALIZATION_RESULT_H
