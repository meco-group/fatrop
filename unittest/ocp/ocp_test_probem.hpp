#ifndef __fatrop_unittest_ocp_ocp_test_hpp__
#define __fatrop_unittest_ocp_ocp_test_hpp__

#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/ocp/ocp_abstract.hpp"
#include <limits>

// example problem 2D point mass
// states: [x, y, vx, vy]
// inputs: [fx, fy]
// dynamics: [xk+1 = xk + dt*vxk, yk+1 = yk + dt*vyk, vxk+1 = vxk + dt*fx/m, vyk+1 = vyk + dt*fy/m]
// cost: fx^2 + fy^2
// constraints:
//  at k = 0: x = 0, y = 0, vx = 0, vy = 0
//  at k = K: x = 1, y = 1, vx = 0, vy = 0
namespace fatrop
{
    namespace test
    {
        class OcpTestProblem : public OcpAbstract
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

            virtual Index get_ng_ineq(const Index k) const { return k == K_ - 1 ? 0 : 2; };
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

                blasfeo_diare_wrap(4, 1.0, res, 2, 0);
                blasfeo_matel_wrap(res, 4, 0) = dt_;
                blasfeo_matel_wrap(res, 5, 1) = dt_;
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
                    blasfeo_diare_wrap(2, *objective_scale*2.0, res, 0, 0);
                }
                return 0;
            };
            virtual Index eval_Ggt(const Scalar *inputs_k, const Scalar *states_k, MAT *res,
                                   const Index k)
            {
                // set zero
                blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
                if (k == 0)
                    blasfeo_diare_wrap(4, 1.0, res, 2, 0);
                if (k == K_ - 1)
                    blasfeo_diare_wrap(4, 1.0, res, 0, 0);
                return 0;
            }
            virtual Index eval_Ggt_ineq(const Scalar *inputs_k, const Scalar *states_k, MAT *res,
                                        const Index k)
            {
                if (k == K_ - 1)
                    return 0;
                // set zero
                blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
                blasfeo_matel_wrap(res, 0, 0) = 1.0;
                blasfeo_matel_wrap(res, 1, 1) = 1.0;
                return 0;
            };
            virtual Index eval_b(const Scalar *states_kp1, const Scalar *inputs_k,
                                 const Scalar *states_k, Scalar *res, const Index k)
            {
                res[0] = -states_kp1[0] + states_k[0] + dt_ * states_k[2] + 0.01 * k;
                res[1] = -states_kp1[1] + states_k[1] + dt_ * states_k[3] + 0.02 * k;
                res[2] = -states_kp1[2] + states_k[2] + dt_ * inputs_k[0] / m_ + 0.03 * k;
                res[3] = -states_kp1[3] + states_k[3] + dt_ * inputs_k[1] / m_ + 0.04 * k;
                return 0;
            }

            virtual Index eval_g(const Scalar *inputs_k, const Scalar *states_k, Scalar *res,
                                 const Index k)
            {
                if (k == 0)
                {
                    res[0] = states_k[0] + 0.04;
                    res[1] = states_k[1] + 0.02;
                    res[2] = states_k[2] + 0.01;
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
            virtual Index eval_gineq(const Scalar *inputs_k, const Scalar *states_k, Scalar *res,
                                     const Index k)
            {
                if (k == K_ - 1)
                    return 0;
                res[0] = inputs_k[0];
                res[1] = inputs_k[1];
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
                    res[0] = 2 * objective_scale[0] * inputs_k[0] + 0.05 * k;
                    res[1] = 2 * objective_scale[0] * inputs_k[1] + 0.06 * k;
                    res[2] = 0. +  0.07 * k;
                    res[3] = 0. +  0.08 * k;
                    res[4] = 0. +  0.09 * k;
                    res[5] = 0. +  0.10 * k;
                    res[6] = 0. +  0.11 * k;
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
                    *res = objective_scale[0] *
                           (inputs_k[0] * inputs_k[0] + inputs_k[1] * inputs_k[1]);
                }
                return 0;
            }
            virtual Index get_bounds(Scalar *lower, Scalar *upper, const Index k) const
            {
                if (k == K_ - 1)
                    return 0;
                lower[0] = -0.5 - 0.1 * k;
                upper[0] = 0.6 + 0.1 * k;
                lower[1] = -0.75 - 0.1 * k;
                upper[1] = std::numeric_limits<Scalar>::infinity();
                return 0;
            }

            virtual Index get_initial_xk(Scalar *xk, const Index k) const
            {
                xk[0] = 0.01 * k;
                xk[1] = 0.02 * k;
                xk[2] = 0.03 * k;
                xk[3] = 0.04 * k;
                return 0;
            };
            virtual Index get_initial_uk(Scalar *uk, const Index k) const
            {
                if (k == K_ - 1)
                    return 0;
                uk[0] = 0.05 * k;
                uk[1] = 0.06 * k;
                return 0;
            };
            virtual ~OcpTestProblem() = default;

        private:
            const Index K_ = 100;
            const Scalar m_ = 1.0;
            const Scalar dt_ = 0.05;
        };

    } // namespace fatrop::test
} // namespace fatrop::test
#endif // __fatrop_unittest_ocp_ocp_test_hpp__