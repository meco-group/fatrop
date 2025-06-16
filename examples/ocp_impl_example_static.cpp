#include <limits>
#include <memory>
#include <fatrop/fatrop.hpp>
#include "fatrop/ocp/nlp_ocp.hxx"
using namespace fatrop;

/**
 * This is an example on how to use a specialization for the OcpAbstractTpl class, rather than using
 * the dynamic polymorphism by using the OcpAbstractDynamic tag. Using static polymorphism reduces
 * the overhead and enables compiler optimizations for the NlpOcpTpl class.
 *
 */

// example problem 2D point mass with a small nonlinearity in the dynamics
// states: [x, y, vx, vy]
// inputs: [fx, fy]
// dynamics: [xk+1 = xk + dt*vxk, yk+1 = yk + dt*vyk, vxk+1 = vxk + dt*fx/m +0.5*dt*fy**2/m, vyk+1
// = vyk + dt*fy/m] cost: fx^2 + fy^2
// constraints:
//  at k = 0: x = 0, y = 0, vx = 0, vy = 0
//  at k = K: x = 1, y = 1, vx = 0, vy = 0

struct MyOcpTag
{
};

typedef fatrop::OcpAbstractTpl<MyOcpTag> MyOcp;

template <> class fatrop::OcpAbstractTpl<MyOcpTag>
{
public:
    Index get_nx(const Index k) const { return 4; }

    Index get_nu(const Index k) const
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
    Index get_ng(const Index k) const
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

    Index get_ng_ineq(const Index k) const { return k == K_ - 1 ? 0 : 2; };
    Index get_horizon_length() const { return K_; };
    Index eval_BAbt(const Scalar *states_kp1, const Scalar *inputs_k, const Scalar *states_k,
                    MAT *res, const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        // Matrix B
        // [  0,    0        ]
        // [  0,    0        ]
        // [ dt/m,  dt/m*fy  ]
        // [  0,    dt/m     ]

        // Matrix A
        // [ 1,  0,  dt,  0  ]
        // [ 0,  1,   0,  dt ]
        // [ 0,  0,   1,  0  ]
        // [ 0,  0,   0,  1  ]

        blasfeo_matel_wrap(res, 0, 2) = dt_ / m_;
        blasfeo_matel_wrap(res, 1, 3) = dt_ / m_;
        blasfeo_matel_wrap(res, 1, 2) = dt_ / m_ * (inputs_k[1] - 1.);

        blasfeo_diare_wrap(4, 1.0, res, 2, 0);
        blasfeo_matel_wrap(res, 4, 0) = dt_;
        blasfeo_matel_wrap(res, 5, 1) = dt_;
        return 0;
    }
    Index eval_RSQrqt(const Scalar *objective_scale, const Scalar *inputs_k, const Scalar *states_k,
                      const Scalar *lam_dyn_k, const Scalar *lam_eq_k, const Scalar *lam_eq_ineq_k,
                      MAT *res, const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        // Matrix R
        // [ 2,  0 ]
        // [ 0,  2 ]
        if (k < K_ - 1)
        {
            blasfeo_diare_wrap(2, objective_scale[0] * 2.0, res, 0, 0);
            // add the contribution from the nonlinearity in the dynamics
            Scalar lam = lam_dyn_k[2];
            blasfeo_matel_wrap(res, 1, 1) += dt_ * lam / m_;
        }
        return 0;
    };
    Index eval_Ggt(const Scalar *inputs_k, const Scalar *states_k, MAT *res, const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        if (k == 0)
            blasfeo_diare_wrap(4, 1.0, res, 2, 0);
        if (k == K_ - 1)
            blasfeo_diare_wrap(4, 1.0, res, 0, 0);
        return 0;
    }
    Index eval_Ggt_ineq(const Scalar *inputs_k, const Scalar *states_k, MAT *res, const Index k)
    {
        if (k == K_ - 1)
            return 0;
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        blasfeo_matel_wrap(res, 0, 0) = 1.0;
        blasfeo_matel_wrap(res, 1, 1) = 1.0;
        return 0;
    };
    Index eval_b(const Scalar *states_kp1, const Scalar *inputs_k, const Scalar *states_k,
                 Scalar *res, const Index k)
    {
        res[0] = -states_kp1[0] + states_k[0] + dt_ * states_k[2];
        res[1] = -states_kp1[1] + states_k[1] + dt_ * states_k[3];
        res[2] = -states_kp1[2] + states_k[2] + dt_ * inputs_k[0] / m_ +
                 dt_ * 0.5 * (inputs_k[1] - 1.) * (inputs_k[1] - 1.) / m_;
        res[3] = -states_kp1[3] + states_k[3] + dt_ * inputs_k[1] / m_;
        return 0;
    }

    Index eval_g(const Scalar *inputs_k, const Scalar *states_k, Scalar *res, const Index k)
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
            res[1] = states_k[1] - 2.0;
            res[2] = states_k[2] - 3.0;
            res[3] = states_k[3] - 4.0;
        }
        return 0;
    };
    Index eval_gineq(const Scalar *inputs_k, const Scalar *states_k, Scalar *res, const Index k)
    {
        if (k == K_ - 1)
            return 0;
        res[0] = inputs_k[0];
        res[1] = inputs_k[1];
        return 0;
    };
    Index eval_rq(const Scalar *objective_scale, const Scalar *inputs_k, const Scalar *states_k,
                  Scalar *res, const Index k)
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
            res[0] = 2 * objective_scale[0] * inputs_k[0];
            res[1] = 2 * objective_scale[0] * inputs_k[1];
            res[2] = 0.;
            res[3] = 0.;
            res[4] = 0.;
            res[5] = 0.;
            res[6] = 0.;
        }
        return 0;
    }
    Index eval_L(const Scalar *objective_scale, const Scalar *inputs_k, const Scalar *states_k,
                 Scalar *res, const Index k)
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
    Index get_bounds(Scalar *lower, Scalar *upper, const Index k) const
    {
        if (k == K_ - 1)
            return 0;
        lower[0] = -50.;
        upper[0] = 50.;
        lower[1] = -100.;
        upper[1] = 100.;
        return 0;
    }

    Index get_initial_xk(Scalar *xk, const Index k) const
    {
        xk[0] = 0.;
        xk[1] = 0.;
        xk[2] = 0.;
        xk[3] = 0.;
        return 0;
    };
    Index get_initial_uk(Scalar *uk, const Index k) const
    {
        if (k == K_ - 1)
            return 0;
        uk[0] = 0.;
        uk[1] = 0.;
        return 0;
    };

private:
    const Index K_ = 100;
    const Scalar m_ = 1.0;
    const Scalar dt_ = 0.05;
};

int main()
{
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcpTpl<MyOcpTag>>(std::make_shared<MyOcp>()));
    std::shared_ptr<IpAlgorithm<OcpType>> ipalg = builder.build();
    IpSolverReturnFlag ret = ipalg->optimize();
    std::cout << "Return flag: " << int(ret) << std::endl;
    std::cout << "Return flag == success: " << (ret == IpSolverReturnFlag::Success) << std::endl;
    return 0;
}