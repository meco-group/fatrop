#ifndef __fatrop_context_context_included__
#define __fatrop_context_context_included__

#if defined(FATROP_CONTEXT_EIGEN)
    #include "eigen.hpp"
#else
    #include "blasfeo.hpp"

    // forward declaration of blasfeo types and functions
    extern "C"
    {
        struct VEC;
        struct MAT;
    }
#endif

#endif //__fatrop_context_context_included__
