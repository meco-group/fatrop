//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_fwd_hpp__
#define __fatrop_graph_fwd_hpp__

namespace fatrop
{
    // Block-Cholesky linear-solver side (tagged with GraphType).
    struct GraphType;
    class BlockSparsityPattern;
    class BlockPdMatrix;
    class BlockEliminationOrder;
    class BlockSymbolicFactorization;
    class BlockCholeskySolver;

    // GraphProblem-side forward declarations (the NLP-level interface that
    // re-uses the block solver above).
    struct GraphProblem;
    template <typename T> struct ProblemDims;
    template <> struct ProblemDims<GraphProblem>;
    template <typename T> struct ProblemInfo;
    template <> struct ProblemInfo<GraphProblem>;
    template <typename T> struct Jacobian;
    template <> struct Jacobian<GraphProblem>;
    template <typename T> struct Hessian;
    template <> struct Hessian<GraphProblem>;
    template <typename T> class PdSolverOrig;
    template <> class PdSolverOrig<GraphProblem>;
    template <typename T> class PdSolverResto;
    template <> class PdSolverResto<GraphProblem>;
    template <typename T> class PdSystemResto;
    template <> class PdSystemResto<GraphProblem>;
    template <typename T> class AugSystemSolver;
    template <> class AugSystemSolver<GraphProblem>;
    class GraphAbstract;
    class NlpGraph;
} // namespace fatrop

#endif // __fatrop_graph_fwd_hpp__
