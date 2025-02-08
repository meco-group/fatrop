//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_common_printing_hpp__
#define __fatrop_common_printing_hpp__

#include <iostream>
#include <memory>

namespace fatrop
{
class Printing
{
public:
    // Delete copy constructor and assignment operator
    Printing(const Printing&) = delete;
    Printing& operator=(const Printing&) = delete;

    // Get the singleton instance
    static Printing& get_instance()
    {
        static Printing instance;
        return instance;
    }

    // Get the current printing stream
    static std::ostream& get_stream()
    {
        return *get_instance().stream_;
    }

    // Set a new printing stream
    static void set_stream(std::unique_ptr<std::ostream> stream)
    {
        if (get_instance().owns_stream_) {
            delete get_instance().stream_;
        }
        get_instance().stream_ = stream.release();
        get_instance().owns_stream_ = true;
    }

private:
    Printing() : stream_(&std::cout), owns_stream_(false) {}
    
    std::ostream* stream_;
    bool owns_stream_;
};

// Global stream accessor
inline std::ostream& f_out = Printing::get_stream();

} // namespace fatrop

#endif // __fatrop_common_printing_hpp__
