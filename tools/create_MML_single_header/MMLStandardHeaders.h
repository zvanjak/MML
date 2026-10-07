/////////////////////////////////////////////////////////////////////////////////////
//  __ __ __ __ _                                                                  //
// |  |  |  |  | |    Minimal Math Library for Modern C++                          //
// | | | | | | | |__  version 2.0                                                  //
// |_| |_|_| |_|____| https://github.com/zvanjak/mml                               //
//                                                                                 //
// Copyright: 2023 - 2026, Zvonimir Vanjak                                         //
//                                                                                 //
// License: MIT License (see LICENSE.md)                                           //
//                                                                                 //
// This is a single-header version of MML - all headers combined into one file.    //
// For the multi-header version, see the mml/ directory.                           //
/////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SINGLE_HEADER
#define MML_SINGLE_HEADER

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

// Standard library headers
#include <stdexcept>
#include <initializer_list>
#include <memory>
#include <functional>
#include <type_traits>
#include <concepts>

#include <string>
#include <vector>
#include <array>
#include <list>
#include <tuple>
#include <utility>
#include <map>
#include <set>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <filesystem>

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <limits>
#include <complex>
#include <numbers>
#include <optional>
#include <regex>
#include <cstdint>

#include <random>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <queue>
#include <stack>


#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#else
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#endif
