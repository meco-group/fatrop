//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/linear_system.hpp"

#include "fatrop/graph/block_pd_matrix.hpp"
#include "fatrop/linear_algebra/blasfeo_operations.hpp"

using namespace fatrop;

LinearSystem<GraphType>::LinearSystem(const BlockPdMatrix &matrix, VecRealView &rhs)
    : matrix_(matrix), rhs_(rhs), m_(matrix.sparsity().total_size())
{
}

void LinearSystem<GraphType>::get_rhs(VecRealView &out) { veccp(m_, rhs_, 0, out, 0); }

void LinearSystem<GraphType>::set_rhs(const VecRealView &in) { veccp(m_, in, 0, rhs_, 0); }

void LinearSystem<GraphType>::apply_on_right(const VecRealView &x, Scalar alpha,
                                             const VecRealView &y, VecRealView &out)
{
    matrix_.apply_on_right(x, alpha, y, out);
}
