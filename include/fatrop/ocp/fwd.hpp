//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_fwd_hpp__
#define __fatrop_ocp_fwd_hpp__

namespace fatrop
{
    class OcpType;
    template <typename T> struct ProblemDims;
    template <> struct ProblemDims<OcpType>;
    template <typename T> struct ProblemInfo;
    template <> struct ProblemInfo<OcpType>;
    template <typename T> struct Jacobian;
    template <> struct Jacobian<OcpType>;
    template <typename T> struct Hessian;
    template <> struct Hessian<OcpType>;
    template <typename T> struct PdSolverOrig;
    template <> class PdSolverOrig<OcpType>;
    template <typename ProblemType> class AugSystemSolver;
    template <> class AugSystemSolver<OcpType>;
    template <typename ProblemType> class PdSystemOrig;
    template <> class PdSystemOrig<OcpType>;
    template <typename ProblemType> class PdSystemResto;
    template <> class PdSystemResto<OcpType>;
    template <typename T> struct PdSolverResto;
    template <> class PdSolverResto<OcpType>;
} // namespace fatrop

#endif // __fatrop_ocp_fwd_hpp__