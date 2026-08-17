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
        ~OutputStreamManager()
        {
            if (owns_stream_)
            {
                delete stream_;
                stream_ = nullptr;
                owns_stream_ = false;
            }
        }

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
        Error = 1,
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
        static void set_print_level(const PrintLevel& print_level) { global_print_level_ = print_level; }
        static void set_print_level(const int& print_level) {global_print_level_ = (PrintLevel) print_level;}

        // Query whether output at the given level would actually be printed. Use this to
        // skip forming output (argument evaluation and stream formatting) altogether.
        static bool is_enabled(const PrintLevel &print_level)
        {
            return global_print_level_ >= print_level;
        }

    private:
        inline static PrintLevel global_print_level_ = PrintLevel::Iterations;
    };

#ifndef FATROP_VERSION
#define FATROP_VERSION "unknown"
#endif

    /**
     * @brief Prints the fatrop copyright/version banner.
     *
     * The banner is printed at most once per process. It can be suppressed
     * with the "suppress_banner" option.
     */
    class Banner
    {
    public:
        // Print the banner if it has not been printed or suppressed.
        static void print_once()
        {
            if (printed_ || suppressed_)
                return;
            printed_ = true;
            OutputStreamManager::get_stream()
                << "\n"
                << "**************************************************************************\n"
                << "   This program contains Fatrop " FATROP_VERSION ", a nonlinear optimization solver\n"
                << "                    for optimal control and robotics.\n"
                << "                 Copyright (c) Lander Vanroye, KU Leuven.\n"
                << "         Fatrop is dual licensed under BSD-2-Clause and EPL-2.0.\n"
                << "     For more information visit https://github.com/meco-group/fatrop\n"
                << "**************************************************************************\n\n";
        }
        // Suppress (or re-enable) printing of the banner.
        static void set_suppress(const bool &suppress) { suppressed_ = suppress; }

    private:
        inline static bool printed_ = false;
        inline static bool suppressed_ = false;
    };

    // Helper making the whole print expression void so it can be the false branch of the
    // conditional in FATROP_PRINT. operator& binds looser than operator<<, so the full
    // insertion chain is swallowed as one operand.
    class OStreamVoidify
    {
    public:
        void operator&(std::ostream &) {}
    };

// When the level is disabled the conditional short-circuits: the insertion chain after the
// macro, including its arguments, is never evaluated.
#define FATROP_PRINT(level)                                                                        \
    (!fatrop::PrintLevelManager::is_enabled(level))                                                \
        ? (void)0                                                                                  \
        : fatrop::OStreamVoidify() & fatrop::OutputStreamManager::get_stream()

#define PRINT_ERROR FATROP_PRINT(fatrop::PrintLevel::Error)
#define PRINT_ITERATIONS FATROP_PRINT(fatrop::PrintLevel::Iterations)
#define PRINT_DEBUG FATROP_PRINT(fatrop::PrintLevel::Debug)
#define PRINT_DIAGNOSTIC FATROP_PRINT(fatrop::PrintLevel::Diagnostic)

} // namespace fatrop

#endif // __fatrop_common_printing_hpp__
