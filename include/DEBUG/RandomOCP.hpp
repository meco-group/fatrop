#ifndef RANDOMOCPINCLUDED
#define RANDOMOCPINCLUDED
#include "OCP/BFOCP.hpp"
#include "DEBUG/FatropDebugTools.hpp"
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include <vector>
namespace fatrop
{
    class RandomOCP : public BFOCP
    {
    public:
        RandomOCP(const vector<int> &nu, const vector<int> &nx, const vector<int> &ng, const int K) : nu_(nu), nx_(nx), ng_(ng), K_(K){

                                                                                                                                 };
        int get_nxk(const int k) const override { return nx_.at(k); };
        int get_nuk(const int k) const override { return nu_.at(k); };
        int get_ngk(const int k) const override { return ng_.at(k); };
        int get_horizon_length() const override { return K_; };
        int eval_BAbtk(const double *states_kp1,
                       const double *scales_states_kp1,
                       const double *states_k,
                       const double *scales_states_k,
                       const double *inputs_k,
                       const double *scales_inputs_k,
                       const double *scales_lam,
                       MAT *res,
                       const int k) override
        {
            int nu = get_nuk(k);
            int nx = get_nxk(k);
            int nxkp1 = get_nxk(k + 1);
            FatropMatBF fmres(res);
            fmres = random_matrix(nu + nx + 1, nxkp1, seed + 1);
            return 0;
        };
        int eval_RSQrqtk(const double *objective_scale,
                         const double *states_k,
                         const double *scales_states_k,
                         const double *inputs_k,
                         const double *scales_inputs_k,
                         const double *lam_dyn_k,
                         const double *scales_lam_dyn_k,
                         const double *lam_eq_k,
                         const double *scales_lam_eq_k,
                         MAT *res,
                         const int k) override
        {
            int nu = get_nuk(k);
            int nx = get_nxk(k);
            FatropMatBF fmres(res);
            if (k < K_ - 1)
            {
                fmres = random_matrix(nu + nx + 1, nu + nx, seed + 0);
                Eig rand1 = Eig(fmres.block(0, 0, nu + nx, nu + nx));
                fmres.block(0, 0, nu + nx, nu + nx) = Eig(rand1.transpose() * rand1);
            }
            else
            {
                fmres.block(nu, nu, nx + 1, nx) = random_matrix(nx + 1, nx, seed + 6);
                Eig rand1 = Eig(fmres.block(nu, nu, nx, nx));
                fmres.block(nu, nu, nx, nx) = Eig(rand1.transpose() * rand1);
            }
            return 0;
        }
        int eval_Ggtk(const double *states_k,
                      const double *scales_states_k,
                      const double *inputs_k,
                      const double *scales_inputs_k,
                      const double *scales,
                      MAT *res,
                      const int k) override
        {
            int nu = get_nuk(k);
            int nx = get_nxk(k);
            int ng = get_ngk(k);
            FatropMatBF fmres(res);
            fmres.block(nu, 0, nx + 1, ng) = random_matrix(nx + 1, ng, seed + 4);
            return 0;
        };
        int eval_bk(
            const double *states_kp1,
            const double *scales_states_kp1,
            const double *states_k,
            const double *scales_states_k,
            const double *inputs_k,
            const double *scales_inputs_k,
            const double *scales_lam,
            double *res,
            const int k) 
        {
                assert(false); // feature not implemented yet
                return 0;
        };
        int eval_gk(
            const double *states_k,
            const double *scales_states_k,
            const double *inputs_k,
            const double *scales_inputs_k,
            const double *scales,
            double *res,
            const int k) 
        {
                assert(false); // feature not implemented yet
                return 0;
        };
        int eval_rqk(
            const double *objective_scale,
            const double *states_k,
            const double *scales_states_k,
            const double *inputs_k,
            const double *scales_inputs_k,
            const double *scales_lam_dyn_k,
            const double *scales_lam_eq_k,
            double *res,
            const int k) 
        {
                assert(false); // feature not implemented yet
                return 0;
        };

    private:
        const vector<int> nu_;
        const vector<int> nx_;
        const vector<int> ng_;
        const int K_;
        int seed = 0;
    };
}
#endif //  RANDOMOCPINCLUDED