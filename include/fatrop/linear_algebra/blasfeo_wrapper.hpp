#ifndef __fatrop_linear_algebra_blasfeo_wrapper_hpp__
#define __fatrop_linear_algebra_blasfeo_wrapper_hpp__

#include "fatrop/common/exception.hpp"
extern "C"
{
#include <blasfeo.h>
}

namespace fatrop
{
    static inline void blasfeo_daxpby_debug(int m, double alpha, blasfeo_dvec *x, int xi,
                                            double beta, blasfeo_dvec *y, int yi, blasfeo_dvec *z,
                                            int zi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_daxpby(m, alpha, x, xi, beta, y, yi, z, zi);
    }

    static inline void blasfeo_daxpy_debug(int m, double alpha, blasfeo_dvec *x, int xi,
                                           blasfeo_dvec *y, int yi, blasfeo_dvec *z, int zi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_daxpy(m, alpha, x, xi, y, yi, z, zi);
    }

    static inline void blasfeo_dvecse_debug(int m, double alpha, blasfeo_dvec *x, int xi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && "Vector index must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m &&
                          "Vector index plus size must not exceed vector dimension");
        blasfeo_dvecse(m, alpha, x, xi);
    }

    static inline void blasfeo_dvecsc_debug(int m, double alpha, blasfeo_dvec *x, int xi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && "Vector index must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m &&
                          "Vector index plus size must not exceed vector dimension");
        blasfeo_dvecsc(m, alpha, x, xi);
    }

    static inline void blasfeo_dveccpsc_debug(int m, double alpha, blasfeo_dvec *x, int xi,
                                              blasfeo_dvec *y, int yi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_dveccpsc(m, alpha, x, xi, y, yi);
    }

    static inline void blasfeo_dveccp_debug(int m, blasfeo_dvec *x, int xi, blasfeo_dvec *y, int yi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_dveccp(m, x, xi, y, yi);
    }

    static inline void blasfeo_dvecmul_debug(int m, blasfeo_dvec *x, int xi, blasfeo_dvec *y,
                                             int yi, blasfeo_dvec *z, int zi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_dvecmul(m, x, xi, y, yi, z, zi);
    }
}

#endif // __fatrop_linear_algebra_blasfeo_wrapper_hpp__
