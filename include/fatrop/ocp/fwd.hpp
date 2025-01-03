//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_fwd_hpp__
#define __fatrop_ocp_fwd_hpp__

namespace fatrop
{
    class OcpType;
    struct OcpDims;
    template <typename T> struct ProblemInfo;
    template <> struct ProblemInfo<OcpType>;
    template <typename T> struct Jacobian;
    template <> struct Jacobian<OcpType>;
    template <typename T> struct Hessian;
    template <> struct Hessian<OcpType>;
} // namespace fatrop

#endif // __fatrop_ocp_fwd_hpp__