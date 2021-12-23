#ifndef SPARSE_OCP_INCLUDED
#define SPARSE_OCP_INCLUDED
using namespace std;
#include "FatropSparse.hpp"
#include "InterfaceMUMPS.hpp"
namespace fatrop
{
#define Id Eigen::MatrixXd::Identity
    class Sparse_OCP
    {
    public:
        Sparse_OCP(OCP_dims &dims, OCP_KKT &OCP, bool Guzero = false) : dims(dims)
        {
            int K = dims.K;
            c_vec = vector<eq_sp>(dims.K, nullptr);
            // initialize variables
            for (int k = 0; k < K - 1; k++)
            {
                int nu = dims.nu.at(k);
                int nx = dims.nx.at(k);
                u_vec.push_back(KKT.get_variable(nu));
                x_vec.push_back(KKT.get_variable(nx));
            }
            {
                int nx = dims.nx.at(K - 1);
                x_vec.push_back(KKT.get_variable(nx));
            }
            for (int k = 0; k < K - 1; k++)
            {
                int nu = dims.nu.at(k);
                int nx = dims.nx.at(k);
                int ng = dims.ng.at(k);
                vector<double> rhs_dyn;
                vector<double> rhs_con;
                vector<double> grad_u;
                vector<double> grad_x;
                for (int i = 0; i < nx; i++)
                {
                    rhs_dyn.push_back(-OCP.BAbt[k].get_el(nu + nx, i));
                };
                for (int i = 0; i < ng; i++)
                {
                    rhs_con.push_back(-OCP.Ggt[k].get_el(nu + nx, i));
                };
                for (int i = 0; i < nu; i++)
                {
                    grad_u.push_back(-OCP.RSQrqt[k].get_el(nu + nx, i));
                };
                for (int i = 0; i < nx; i++)
                {
                    grad_x.push_back(-OCP.RSQrqt[k].get_el(nu + nx, nu + i));
                };
                x_vec.at(k)->set_grad(grad_x);
                u_vec.at(k)->set_grad(grad_u);
                KKT.set_hess_block(Eig(OCP.RSQrqt[k].block(0, 0, nu, nu)), u_vec.at(k), u_vec.at(k));
                KKT.set_hess_block(Eig(OCP.RSQrqt[k].block(nu, nu, nx, nx)), x_vec.at(k), x_vec.at(k));
                KKT.set_hess_block(Eig(OCP.RSQrqt[k].block(nu, 0, nx, nu)), x_vec.at(k), u_vec.at(k));
                Eig B(Eig(Eig(OCP.BAbt[k].block(0, 0, nu, nx)).transpose()));
                Eig A(Eig(Eig(OCP.BAbt[k].block(nu, 0, nx, nx)).transpose()));
                // TODO EYE IS HERE CONSIDERED AS A DENSE MATRIX
                KKT.set_equation(B * u_vec.at(k) + A * (x_vec.at(k)) + Eig(-Id(dims.nx.at(k + 1), dims.nx.at(k + 1))) * (x_vec.at(k + 1)), rhs_dyn);
                if (Guzero)
                {
                    Eig Gx(Eig(OCP.Ggt[k].block(nu, 0, nx, ng)).transpose());
                    c_vec.at(k) = (KKT.set_equation(Gx * x_vec.at(k), rhs_con));
                }
                else
                {
                    Eig Gu(Eig(OCP.Ggt[k].block(0, 0, nu, ng)).transpose());
                    Eig Gx(Eig(OCP.Ggt[k].block(nu, 0, nx, ng)).transpose());
                    c_vec.at(k) = (KKT.set_equation(Gu * u_vec.at(k) + Gx * x_vec.at(k), rhs_con));
                }
            }
            // K - 1
            {
                int nu = dims.nu.at(K - 1);
                int nx = dims.nx.at(K - 1);
                int ng = dims.ng.at(K - 1);
                vector<double> rhs_con;
                vector<double> grad_x;
                for (int i = 0; i < ng; i++)
                {
                    rhs_con.push_back(-OCP.Ggt[K - 1].get_el(nu + nx, i));
                };
                for (int i = 0; i < nx; i++)
                {
                    grad_x.push_back(-OCP.RSQrqt[K - 1].get_el(nu + nx, nu + i));
                };
                Eig Gx(Eig(OCP.Ggt[K - 1].block(nu, 0, nx, ng)).transpose());
                x_vec.at(K - 1)->set_grad(grad_x);
                KKT.set_hess_block(Eig(OCP.RSQrqt[K - 1].block(nu, nu, nx, nx)), x_vec.at(K - 1), x_vec.at(K - 1));
                c_vec.at(K - 1) = KKT.set_equation(Gx * x_vec.at(K - 1), rhs_con);
            };
        };

        void fact_solve(const fatrop_vector_bf &ux, const fatrop_vector_bf &lam)
        {

            vector<triplet> ocptripl;
            KKT.get_triplets(ocptripl);
            InterfaceMUMPS interfo(ocptripl.size(), KKT.get_size(), ocptripl);
            interfo.preprocess();
            vector<double> rhso(KKT.get_rhs());
            interfo.solve(ocptripl, rhso);

            int offs = 0;
            int offs_c = 0;
            VEC *ux_p = (VEC *)ux;
            VEC *lam_p = (VEC *)lam;
            int offs_lam = sum(dims.nu + dims.nx) - dims.nu.at(dims.K - 1);
            // put result in ux vector
            for (int k = 0; k < dims.K - 1; k++)
            {
                int nxk = dims.nx.at(k);
                int nuk = dims.nu.at(k);
                int ngk = dims.ng.at(k);
                int sol_offs_u = u_vec.at(k)->offset;
                int sol_offs_x = x_vec.at(k)->offset;
                int lam_c_offs = c_vec.at(k)->offset;
                for (int i = 0; i < nuk; i++)
                {
                    VECEL(ux_p, offs + i) = rhso.at(sol_offs_u + i);
                }
                for (int i = 0; i < nxk; i++)
                {
                    VECEL(ux_p, offs + nuk + i) = rhso.at(sol_offs_x + i);
                }
                for (int i = 0; i < ngk; i++)
                {
                    VECEL(lam_p, offs_c + i) = rhso.at(offs_lam + lam_c_offs + i);
                }
                offs += nuk + nxk;
                offs_c += ngk;
            }
            {
                int sol_offs_x = x_vec.at(dims.K - 1)->offset;
                int nxk = dims.nx.at(dims.K - 1);
                int nuk = dims.nu.at(dims.K - 1);
                int lam_c_offs = c_vec.at(dims.K - 1)->offset;
                int ngk = dims.ng.at(dims.K - 1);
                for (int i = 0; i < nxk; i++)
                {
                    VECEL(ux_p, offs + nuk + i) = rhso.at(sol_offs_x + i);
                }
                for (int i = 0; i < ngk; i++)
                {
                    VECEL(lam_p, offs_c + i) = rhso.at(offs_lam + lam_c_offs + i);
                }
            }
        }
        KKT_matrix KKT;
        OCP_dims dims;
        vector<var_sp> u_vec;
        vector<var_sp> x_vec;
        vector<eq_sp> c_vec;
    };

} // namespace fatrop

#endif // SPARSE_OCP_INCLUDED