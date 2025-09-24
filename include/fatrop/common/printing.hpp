//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_common_printing_hpp__
#define __fatrop_common_printing_hpp__

#include <iostream>
#include <memory>

namespace fatrop
{
    class OutputStreamManager
    {
    public:
        // Delete copy constructor and assignment operator
        OutputStreamManager(const OutputStreamManager &) = delete;
        OutputStreamManager &operator=(const OutputStreamManager &) = delete;

        // Get the singleton instance
        static OutputStreamManager &get_instance()
        {
            static OutputStreamManager instance;
            return instance;
        }

        // Get the current printing stream
        static std::ostream &get_stream() { return *get_instance().stream_; }

        // Set a new printing stream
        static void set_stream(std::unique_ptr<std::ostream> stream)
        {
            if (get_instance().owns_stream_)
            {
                delete get_instance().stream_;
            }
            get_instance().stream_ = stream.release();
            get_instance().owns_stream_ = true;
        }

    private:
        OutputStreamManager() : stream_(&std::cout), owns_stream_(false) {}

        std::ostream *stream_;
        bool owns_stream_;
    };

    enum class PrintLevel
    {
        None = 0,
        Iterations = 5, // consistent with Ipopt and legacy fatrop
        Debug = 6,
        Diagnostic = 7,
    };
    class NullBuffer : public std::streambuf
    {
    public:
        int overflow(int c) override { return c; }
    };

    class NullStream : public std::ostream
    {
    public:
        NullStream() : std::ostream(&null_buffer_) {}

    private:
        NullBuffer null_buffer_;
    };

    class PrintLevelManager
    {
    public:
        PrintLevelManager(PrintLevel print_level) : this_print_level_(print_level) {};
        template <typename T> std::ostream &operator<<(const T &value)
        {
            if (global_print_level_ >= this_print_level_)
                return OutputStreamManager::get_stream() << value;
            return null_stream_;
        };

        static void set_print_level(const PrintLevel& print_level) { global_print_level_ = print_level; }
        static void set_print_level(const int& print_level) {global_print_level_ = (PrintLevel) print_level;}

    private:
        const PrintLevel this_print_level_;
        inline static PrintLevel global_print_level_ = PrintLevel::Iterations;
        inline static NullStream null_stream_;
    };

#define PRINT_ITERATIONS PrintLevelManager(PrintLevel::Iterations)
#define PRINT_DEBUG PrintLevelManager(PrintLevel::Debug)
#define PRINT_DIAGNOSTIC PrintLevelManager(PrintLevel::Diagnostic)

} // namespace fatrop

#endif // __fatrop_common_printing_hpp__
