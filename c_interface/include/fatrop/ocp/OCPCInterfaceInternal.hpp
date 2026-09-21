#ifndef __fatrop_ocp_solver_ocp_c_interface_internal_hpp__
#define __fatrop_ocp_solver_ocp_c_interface_internal_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/nlp/nlp.hpp"
namespace fatrop
{
#include "fatrop/ocp/OCPCInterface.h"
}
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/problem_info.hpp"
namespace fatrop
{
    class FatropOcpCMapping : public Nlp<OcpType>
    {
    public:
        FatropOcpCInterface *ocp;
        FatropOcpCDims s;
        FatropOcpCMapping(FatropOcpCInterface *ocp);
        const NlpDims &nlp_dims() const override;
        const ProblemDims<OcpType> &problem_dims() const override;
        Index eval_lag_hess(const ProblemInfo<OcpType> &info, const Scalar objective_scale,
                            const VecRealView &primal_x, const VecRealView &primal_s,
                            const VecRealView &lam, Hessian<OcpType> &hess) override;
        Index eval_constr_jac(const ProblemInfo<OcpType> &info, const VecRealView &primal_x,
                              const VecRealView &primal_s, Jacobian<OcpType> &jac) override;
        Index eval_constraint_violation(const ProblemInfo<OcpType> &info,
                                        const VecRealView &primal_x, const VecRealView &primal_s,
                                        VecRealView &res) override;
        Index eval_objective_gradient(const ProblemInfo<OcpType> &info,
                                      const Scalar objective_scale, const VecRealView &primal_x,
                                      const VecRealView &primal_s, VecRealView &grad_x,
                                      VecRealView &grad_s) override;
        Index eval_objective(const ProblemInfo<OcpType> &info, const Scalar objective_scale,
                             const VecRealView &primal_x, const VecRealView &primal_s,
                             Scalar &res) override;
        Index get_bounds(const ProblemInfo<OcpType> &info, VecRealView &lower_bounds,
                         VecRealView &upper_bounds) override;
        Index get_initial_primal(const ProblemInfo<OcpType> &info, VecRealView &primal_x) override;
        void get_primal_damping(const ProblemInfo<OcpType> &info, VecRealView &damping) override;
        void apply_jacobian_s_transpose(const ProblemInfo<OcpType> &info,
                                        const VecRealView &multipliers, const Scalar alpha,
                                        const VecRealView &y, VecRealView &out) override;

    private:
        ProblemDims<OcpType> ocp_dims_;
        NlpDims nlp_dims_;
        Index K_;
        std::vector<MAT> matrix_buffer_[3];
    };
}
#endif // __fatrop_ocp_solver_ocp_c_interface_internal_hpp__