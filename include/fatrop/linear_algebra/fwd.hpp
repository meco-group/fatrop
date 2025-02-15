//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_fwd_hpp__
#define __fatrop_linear_algebra_fwd_hpp__

namespace fatrop
{
    // Base vector template class
    template <typename Derived> class VecReal;

    // Vector block and operations
    template <typename Derived> class VecRealBlock;
    template <typename Dep1> class VecRealAbs;
    template <typename Dep1> class VecRealLog;
    template <typename Dep1> class VecRealExp;
    template <typename Dep1> class VecRealSin;
    template <typename Dep1> class VecRealCos;
    template <typename Dep1> class VecRealSqrt;
    template <typename IfElseOp, typename Dep1, typename Dep2> class VecRealIfElse;
    template <typename Dep1> class VecRealTimesScalar;
    template <typename Dep1, typename Dep2> class VecRealPlusVecReal;
    template <typename Dep1> class ScalarDivVecReal;
    template <typename Dep1, typename Dep2> class VecRealMinusVecReal;
    template <typename Dep1, typename Dep2> class VecRealTimesVecReal;
    template <typename Dep1, typename Dep2> class VecRealDivVecReal;
    template <typename Dep1, typename Dep2> class VecRealMin;
    template <typename Dep1, typename Dep2> class VecRealMax;
    template <typename Dep1, typename Dep2> class VecRealConcat;

    // Non-template vector classes
    class VecRealScalar;
    class VecRealView;
    class VecRealAllocated;

    // Base class for vector operation specializations
    template <typename Derived> class VecOperationSpecialization;

    // Matrix-related forward declarations
    template <typename Derived> class MatReal;

    template <typename Derived> class MatrixBlock;

    class MatRealAllocated;
    template <typename Derived> class MatRealView1D;
    template <typename Dep1> class MatRealTranspose;

    class MatRealRowView;
    class MatRealColView;
    class MatRealDiagonalView;
    class MatRealView;

    template<typename LsType> class LinearSystem;
    template<typename Derived, typename LsType> class LinearSolver;

} // namespace fatrop

#endif // __fatrop_fwd_hpp__
