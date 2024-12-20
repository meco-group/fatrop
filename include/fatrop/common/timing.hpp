
//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_common_timing_hpp__
#define __fatrop_common_timing_hpp__

extern "C"
{
#include <blasfeo/include/blasfeo_timing.h>
}
#define fatrop_timer blasfeo_timer
#define fatrop_tic blasfeo_tic
#define fatrop_toc blasfeo_toc

#endif // __fatrop_common_timing_hpp__
