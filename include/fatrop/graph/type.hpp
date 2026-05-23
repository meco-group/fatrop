//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_type_hpp__
#define __fatrop_graph_type_hpp__

/**
 * @file type.hpp
 * @brief Defines the GraphType tag for block-sparse symmetric positive
 *        definite linear systems.
 *
 * A graph-typed linear system has the form @c A x = -b where @c A is a
 * symmetric positive definite matrix with a block-sparse structure described
 * by an undirected adjacency graph. Each vertex of the graph corresponds to a
 * variable block of a user-prescribed size; each edge corresponds to a dense
 * off-diagonal block (and its symmetric counterpart). The matrix is the
 * natural generalisation of a scalar SPD system to one whose unknowns are
 * vectors rather than scalars.
 *
 * The corresponding linear solver applies a block-supernodal Cholesky
 * factorisation using a (block-)minimum-degree elimination order.
 *
 * @see fatrop::BlockSparsityPattern
 * @see fatrop::BlockPdMatrix
 * @see fatrop::BlockCholeskySolver
 */

namespace fatrop
{
    struct GraphType
    {
    };
} // namespace fatrop

#endif // __fatrop_graph_type_hpp__
