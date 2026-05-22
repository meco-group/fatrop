//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_fwd_hpp__
#define __fatrop_dense_fwd_hpp__

namespace fatrop
{
    struct DenseType;
    template <typename T> struct ProblemDims;
    template <> struct ProblemDims<DenseType>;
    template <typename T> struct ProblemInfo;
    template <> struct ProblemInfo<DenseType>;
    template <typename T> struct Jacobian;
    template <> struct Jacobian<DenseType>;
    template <typename T> struct Hessian;
    template <> struct Hessian<DenseType>;
    template <typename T> struct PdSolverOrig;
    template <> class PdSolverOrig<DenseType>;
    template <typename ProblemType> class AugSystemSolver;
    template <> class AugSystemSolver<DenseType>;
    template <typename ProblemType> class PdSystemOrig;
    template <> class PdSystemOrig<DenseType>;
    template <typename ProblemType> class PdSystemResto;
    template <> class PdSystemResto<DenseType>;
    template <typename T> struct PdSolverResto;
    template <> class PdSolverResto<DenseType>;
} // namespace fatrop

#endif // __fatrop_dense_fwd_hpp__
