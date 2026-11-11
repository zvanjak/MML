#include <mml/base/Vector/Vector.h>

#if defined(MML_PERSISTENCE_VECTOR_IO_H) || defined(MML_IO_RESULT_H)
#error The mathematical vector must not depend on persistence tools.
#endif

#include <string>

template<class Type>
concept HasVectorFileSave = requires(const Type& value, const std::string& filename) {
	Type::SaveToFile(value, filename);
};

template<class Type>
concept HasVectorFileLoad = requires(Type& value, const std::string& filename) {
	Type::LoadFromFile(filename, value);
};

static_assert(!HasVectorFileSave<MML::Vector<MML::Real>>);
static_assert(!HasVectorFileLoad<MML::Vector<MML::Real>>);
