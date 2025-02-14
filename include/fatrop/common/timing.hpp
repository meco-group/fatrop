
//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_common_timing_hpp__
#define __fatrop_common_timing_hpp__

extern "C"
{
#include <blasfeo.h> // we use the blasfeo timer internally
}
#include "fatrop/common/printing.hpp"
#define fatrop_timer blasfeo_timer
#define fatrop_tic blasfeo_tic
#define fatrop_toc blasfeo_toc
namespace fatrop
{
    struct Timer
    {
        void start() { fatrop_tic(&timer); }
        double stop()
        {
            el_time += fatrop_toc(&timer);
            double ret = el_time;
            el_time = 0.;
            return ret;
        }
        double pause()
        {
            el_time += fatrop_toc(&timer);
            return el_time;
        }
        double elapsed() const { return el_time; }
        void reset() { el_time = 0.; }
        blasfeo_timer timer;
        double el_time = 0.;
        friend std::ostream &operator<<(std::ostream &os, const Timer &el_time)
        {
            os << el_time;
            return os;
        }
    };
}

#endif // __fatrop_common_timing_hpp__
