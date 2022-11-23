#ifndef FATROPSTATSINCLUDED
#define FATROPSTATSINCLUDED
#include <iostream>
using namespace std;
namespace fatrop
{
    struct FatropStats
    {
        double compute_sd_time = 0.0;
        double duinf_time = 0.0;
        double eval_hess_time = 0.0;
        double eval_jac_time = 0.0;
        double eval_cv_time = 0.0;
        double eval_grad_time = 0.0;
        double eval_obj_time = 0.0;
        double initialization_time = 0.0;
        double time_total = 0.0;
        int eval_hess_count = 0;
        int eval_jac_count = 0;
        int eval_cv_count = 0;
        int eval_grad_count = 0;
        int eval_obj_count = 0;
        int iterations_count = 0;
        void Print()
        {
            cout << "---- stats ----" << endl;
            cout << "compute_sd:     " << compute_sd_time << " s"<< endl;
            cout << "duinf:          " << duinf_time << " s"<< endl;
            cout << "initialization: " << initialization_time << " s  count: "<< iterations_count <<endl;
            double time_FE = eval_hess_time + eval_jac_time + eval_cv_time + eval_grad_time + eval_obj_time ;
            cout << "time_FE :       " << time_FE << " s"<< endl;
            cout << "    eval hess:  " << eval_hess_time << " s  count: "<< eval_hess_count <<endl;
            cout << "    eval jac:   " << eval_jac_time << " s  count: "<< eval_jac_count <<endl;
            cout << "    eval cv:    " << eval_cv_time << " s  count: "<< eval_cv_count <<endl;
            cout << "    eval grad:  " << eval_grad_time << " s  count: "<< eval_grad_count <<endl;
            cout << "    eval obj:   " << eval_obj_time << " s  count: "<< eval_obj_count <<endl;
            cout << "rest  :       " <<time_total -  time_FE - initialization_time - duinf_time - compute_sd_time<< " s"<< endl;
            cout << "----- "<< endl;
            cout << "time_w/o_FE : " << time_total - time_FE << " s" << endl;
            cout << "time_FE :       " << time_FE << " s"<< endl;
            cout << "time_total : " << time_total << " s  iterations: " << iterations_count << endl;
        }
    };
} // namespace fatrop
#endif // FATROPSTATSINCLUDED