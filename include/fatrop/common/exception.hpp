//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_common_exception_hpp__
#define __fatrop_common_exception_hpp__

#include <exception>
#include <iostream>
#include <string>

// Define assertion macro
#define fatrop_assert(condition)                                                                   \
    do                                                                                             \
    {                                                                                              \
        if (!(condition))                                                                          \
        {                                                                                          \
            std::cerr << "Assertion failed: " << #condition << "at: " << __FILE__ << ":"           \
                      << __LINE__ << "\n";                                                         \
            std::abort();                                                                          \
        }                                                                                          \
    } while (0)

// Define debug assertion macro
#ifdef NDEBUG
#define fatrop_dbg_assert(condition) ((void)0) // No-op in release builds
#else
#define fatrop_dbg_assert(condition) fatrop_assert(condition)
#endif

#endif // __fatrop_common_exception_hpp__
