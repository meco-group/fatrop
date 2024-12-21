//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_common_exception_hpp__
#define __fatrop_common_exception_hpp__

#include <exception>
#include <iostream>
#include <string>
#include <sstream>

// Custom exception class
class FatropException : public std::exception {
public:
    FatropException(const std::string& condition, const char* file, int line)
        : message_("Assertion failed: " + condition + " at: " + file + ":" + std::to_string(line)) {}

    const char* what() const noexcept override {
        return message_.c_str();
    }

private:
    std::string message_;
};

// Define assertion macro
#define fatrop_assert(condition)                                                                   \
    do                                                                                             \
    {                                                                                              \
        if (!(condition))                                                                          \
        {                                                                                          \
            std::ostringstream oss;                                                                \
            oss << "Assertion failed: " << #condition << " at: " << __FILE__ << ":"                \
                << __LINE__ << "\n";                                                               \
            std::cerr << oss.str();                                                                \
            throw FatropException(#condition, __FILE__, __LINE__);                                 \
        }                                                                                          \
    } while (0)

// Define debug assertion macro
#ifdef NDEBUG
#define fatrop_dbg_assert(condition) ((void)0) // No-op in release builds
#else
#define fatrop_dbg_assert(condition) fatrop_assert(condition)
#endif

#endif // __fatrop_common_exception_hpp__
