#include <mml/MMLVisualizators.h>

#ifndef MML_VISUALIZER_CONFIG_H
#error The legacy configuration header must forward to the canonical visualizer configuration.
#endif

#if defined(MML_VISUALIZER_FACADE_H) || defined(MML_VISUALIZER_PROCESS_H) || \
	defined(MML_VISUALIZER_FILES_H) || defined(MML_SERIALIZER_BASE_H) || \
	defined(_WINDOWS_) || defined(_INC_WINDOWS)
#error Configuration-only includes must not load serialization, the facade, or platform process headers.
#endif

#include <type_traits>

static_assert(std::is_same_v<decltype(MML::GetVisualizerBackend()), MML::VisualizerBackend>);
