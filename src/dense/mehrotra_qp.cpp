//
// Copyright (C) 2026 Lander Vanroye, KU Leuven
//
// Explicit instantiation of the Mehrotra predictor-corrector QP solver for
// the DenseType problem class.
//
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/dense/pd_system_orig.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/qp/fatrop_qp.hxx"

namespace fatrop
{
    template class MehrotraQpAlgorithm<DenseType>;
    template class MehrotraQpBuilder<DenseType>;
} // namespace fatrop
