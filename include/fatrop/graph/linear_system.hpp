//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_linear_system_hpp__
#define __fatrop_graph_linear_system_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/type.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/linear_system.hpp"
#include "fatrop/linear_algebra/vector.hpp"

namespace fatrop
{
    /**
     * @brief Linear system @c M x = -b backed by a @ref BlockPdMatrix.
     *
     * Mirrors the @c LinearSystem interface used by the OCP and dense
     * formulations: @c get_rhs / @c set_rhs read and write the stored
     * right-hand side @c b, while @c apply_on_right computes
     * @c out = M x + alpha y. The right-hand side is held by reference (a
     * @c VecRealView) so iterative refinement can mutate it in place.
     */
    template <> class LinearSystem<GraphType>
    {
        friend class BlockCholeskySolver;

    public:
        LinearSystem(const BlockPdMatrix &matrix, VecRealView &rhs);

        Index m() const { return m_; }

        void get_rhs(VecRealView &out);
        void set_rhs(const VecRealView &in);
        void apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView &y,
                            VecRealView &out);

        const BlockPdMatrix &matrix() const { return matrix_; }

    private:
        const BlockPdMatrix &matrix_;
        VecRealView &rhs_;
        const Index m_;
    };
} // namespace fatrop

#endif // __fatrop_graph_linear_system_hpp__
