//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_fwd_hpp__
#define __fatrop_linear_algebra_fwd_hpp__

namespace fatrop
{
    // Base vector template class
    template <typename Derived> class Vec;

    // Vector block and operations
    template <typename Derived> class VecBlock;
    template <typename Dep1> class VecAbs;
    template <typename Dep1> class VecLog;
    template <typename Dep1> class VecExp;
    template <typename Dep1> class VecSin;
    template <typename Dep1> class VecCos;
    template <typename IfElseOp, typename Dep1, typename Dep2> class VecIfElse;
    template <typename Dep1> class VecTimesScalar;
    template <typename Dep1, typename Dep2> class VecPlusVec;
    template <typename Dep1> class ScalarDivVec;
    template <typename Dep1, typename Dep2> class VecMinusVec;
    template <typename Dep1, typename Dep2> class VecTimesVec;
    template <typename Dep1, typename Dep2> class VecDivVec;

    // Non-template vector classes
    class VecScalar;
    class VecNumeric;
    class VecAllocated;

    // Base class for vector operation specializations
    template <typename Derived> class VecOperationSpecialization;

    // Matrix-related forward declarations
    template <typename Derived> class Mat;

    template <typename Derived> class MatrixBlock;

    class MatrixAllocated;
    template <typename Derived> class MatView1D;

    class MatRowView;
    class MatColView;
    class MatDiagonalView;
    class MatrixNumeric;

} // namespace fatrop

#endif // __fatrop_fwd_hpp__
