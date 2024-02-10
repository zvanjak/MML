#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Intervals.h"
#endif

using namespace MML;

class TestFuncClass
{
    std::string _name;
    const IInterval &_interval;

public:
    TestFuncClass(std::string inName, const IInterval &inInterval) : _name(inName), _interval(inInterval) { }


};

void Demo_Intervals()
{
    OpenClosedInterval int1(0.0, 1.0);

}