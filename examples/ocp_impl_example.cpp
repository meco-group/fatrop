#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/ocp/ocp_abstract.hpp"
using namespace fatrop;

// example problem 2D point mass
// states: [x, y, vx, vy]
// inputs: [fx, fy]
// dynamics: [xk+1 = xk + dt*vxk, yk+1 = yk + dt*vyk, vxk+1 = vxk + dt*fx/m, vyk+1 = vyk + dt*fy/m]
// cost: fx^2 + fy^2
// constraints:
//  at k = 0: x = 0, y = 0, vx = 0, vy = 0
//  at k = K: x = 1, y = 1, vx = 0, vy = 0
class OcpTest : public OcpAbstract
{
public:
    virtual Index get_nx(const Index k) const { return 4; }

    virtual Index get_nu(const Index k) const
    {
        if (k == 0)
        {
            return 2;
        }
        else if (k == K_ - 1)
        {
            return 0;
        }
        else
        {
            return 2;
        }
    }
    virtual Index get_ng(const Index k) const
    {
        if (k == 0)
        {
            return 4;
        }
        else if (k == K_ - 1)
        {
            return 4;
        }
        else
        {
            return 0;
        }
    };

    virtual Index get_ng_ineq(const Index k) const { return 0; };
    virtual Index get_horizon_length() const { return K_; };
    virtual Index eval_BAbt(const Scalar *states_kp1, const Scalar *inputs_k,
                            const Scalar *states_k, MAT *res, const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        // Matrix B
        // [  0,    0  ]
        // [  0,    0  ]
        // [ dt/m,  0  ]
        // [  0,   dt/m ]

        // Matrix A
        // [ 1,  0,  dt,  0  ]
        // [ 0,  1,   0,  dt ]
        // [ 0,  0,   1,  0  ]
        // [ 0,  0,   0,  1  ]

        blasfeo_matel_wrap(res, 0, 2) = dt_ / m_;
        blasfeo_matel_wrap(res, 1, 3) = dt_ / m_;

        blasfeo_diare_wrap(4, 1.0, res, 0, 0);
        blasfeo_matel_wrap(res, 2, 0) = dt_;
        blasfeo_matel_wrap(res, 3, 1) = dt_;
        return 0;
    }
    virtual Index eval_RSQrqt(const Scalar *objective_scale, const Scalar *inputs_k,
                              const Scalar *states_k, const Scalar *lam_dyn_k,
                              const Scalar *lam_eq_k, const Scalar *lam_eq_ineq_k, MAT *res,
                              const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        // Matrix R
        // [ 2,  0 ]
        // [ 0,  2 ]
        if (k < K_ - 1)
        {
            blasfeo_diare_wrap(2, 2.0, res, 0, 0);
        }
        return 0;
    };
    virtual Index eval_Ggt(const Scalar *inputs_k, const Scalar *states_k, MAT *res, const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        if (k == 0)
            blasfeo_diare_wrap(4, 1.0, res, 0, 0);
        if (k == K_ - 1)
            blasfeo_diare_wrap(4, 1.0, res, 0, 0);
        return 0;
    }
    virtual Index eval_Ggt_ineq(const Scalar *inputs_k, const Scalar *states_k, MAT *res,
                                const Index k)
    {
        return 0;
    };
    virtual Index eval_b(const Scalar *states_kp1, const Scalar *inputs_k, const Scalar *states_k,
                         Scalar *res, const Index k)
    {
        res[0] = states_kp1[0] - states_k[0] - dt_ * states_k[2];
        res[1] = states_kp1[1] - states_k[1] - dt_ * states_k[3];
        res[2] = states_kp1[2] - states_k[2] - dt_ * inputs_k[0] / m_;
        res[3] = states_kp1[3] - states_k[3] - dt_ * inputs_k[1] / m_;
        return 0;
    }

    virtual Index eval_g(const Scalar *inputs_k, const Scalar *states_k, Scalar *res, const Index k)
    {
        if (k == 0)
        {
            res[0] = states_k[0];
            res[1] = states_k[1];
            res[2] = states_k[2];
            res[3] = states_k[3];
        }
        else if (k == K_ - 1)
        {
            res[0] = states_k[0] - 1.0;
            res[1] = states_k[1] - 1.0;
            res[2] = states_k[2];
            res[3] = states_k[3];
        }
        return 0;
    };
    virtual Index eval_gineq(const Scalar *inputs_k, const Scalar *states_k, Scalar *res,
                             const Index k)
    {
        return 0;
    };
    virtual Index eval_rq(const Scalar *objective_scale, const Scalar *inputs_k,
                          const Scalar *states_k, Scalar *res, const Index k)
    {
        if (k == K_ - 1)
        {
            res[0] = 0.;
            res[1] = 0.;
            res[2] = 0.;
            res[3] = 0.;
        }
        else
        {
            res[0] = objective_scale[0] * inputs_k[0];
            res[1] = objective_scale[0] * inputs_k[1];
            res[2] = 0.;
            res[3] = 0.;
            res[4] = 0.;
            res[5] = 0.;
            res[6] = 0.;
        }
        return 0;
    }
    virtual Index eval_L(const Scalar *objective_scale, const Scalar *inputs_k,
                         const Scalar *states_k, Scalar *res, const Index k)
    {
        if (k == K_ - 1)
        {
            *res = 0.;
        }
        else
        {
            *res = objective_scale[0] * (inputs_k[0] * inputs_k[0] + inputs_k[1] * inputs_k[1]);
        }
        return 0;
    }
    virtual Index get_bounds(Scalar *lower, Scalar *upper, const Index k) const { return 0; }
    virtual Index get_initial_xk(Scalar *xk, const Index k) const { return 0; };
    virtual Index get_initial_uk(Scalar *uk, const Index k) const { return 0; };
    virtual ~OcpTest() = default;

private:
    const Index K_ = 100;
    const Scalar m_ = 1.0;
    const Scalar dt_ = 0.05;
};

int main() { return 0; }