#include <mml/MMLExceptions.h>
#include <type_traits>

template<class Type>
concept RootExposesDomainException = requires { sizeof(Type); };

static_assert(RootExposesDomainException<MML::FourierError>);
static_assert(RootExposesDomainException<MML::StatisticsError>);
static_assert(RootExposesDomainException<MML::CurveFittingError>);
static_assert(std::is_base_of_v<MML::MMLException, MML::FourierError>);
static_assert(std::is_base_of_v<MML::MMLException, MML::StatisticsError>);
static_assert(std::is_base_of_v<MML::MMLException, MML::CurveFittingError>);
