#ifndef FATROPITERATIONDATAINCLUDED
#define FATROPITERATIONDATAINCLUDED
#include "aux/SmartPtr.hpp"
#include <cmath>
#include <vector>
#include <iostream>
using namespace std;
namespace fatrop
{
    struct IterationData
    {
        int iter = 0;
        double mu = 0.0;
        double objective = 0.0;
        double constraint_violation = 0.0;
        double du_inf = 0.0;
        int ls = 0;
        double reg = 0.0;
        double alpha_pr = 0.0;
        double alpha_du = 0.0;
        char type = 'x';
    };
    class Journaller 
    {
    public:
        Journaller(const int maxiter)
        {
            iterationdata.reserve(maxiter + 1);
        }
        void PrintIterations()
        {
            if (print_count == 0)
            {
                printf(" it  obj                    cv                  du                  lg(mu) reg  alpha_du  alpha_pr  ls\n");
            }
            for (std::vector<double>::size_type i = print_count; i < iterationdata.size(); i++)
            {
                IterationData iterationdata_i = iterationdata.at(i);
                // if (iterationdata_i.reg == 0.0)
                // {
                //     printf("step %3d: obj %.5e, cv : %.2e, du %.2e, lg(mu), %4.1f, lg(reg)  -.-, a_d %.2e, a_p %.2e, ls %d%c \n", iterationdata_i.iter, iterationdata_i.objective, iterationdata_i.constraint_violation, iterationdata_i.du_inf, log10(iterationdata_i.mu), iterationdata_i.alpha_du, iterationdata_i.alpha_pr, iterationdata_i.ls, iterationdata_i.type);
                // }
                // else
                // {
                //     printf("step %3d: obj %.5e, cv : %.2e, du %.2e, lg(mu), %4.1f, lg(reg) %4.1f, a_d %.2e, a_p %.2e, ls %d%c \n", iterationdata_i.iter, iterationdata_i.objective, iterationdata_i.constraint_violation, iterationdata_i.du_inf, log10(iterationdata_i.mu), log10(iterationdata_i.reg), iterationdata_i.alpha_du, iterationdata_i.alpha_pr, iterationdata_i.ls, iterationdata_i.type);
                // }
                if (iterationdata_i.reg == 0.0)
                {
                    printf("%3d, %.15e, %.12e, %.12e, %4.1f,  -.-, %.2e, %.2e, %d%c \n", iterationdata_i.iter, iterationdata_i.objective, iterationdata_i.constraint_violation, iterationdata_i.du_inf, log10(iterationdata_i.mu), iterationdata_i.alpha_du, iterationdata_i.alpha_pr, iterationdata_i.ls, iterationdata_i.type);
                }
                else
                {
                    printf("%3d, %.15e, %.12e, %.12e, %4.1f, %4.1f, %.2e, %.2e, %d%c \n", iterationdata_i.iter, iterationdata_i.objective, iterationdata_i.constraint_violation, iterationdata_i.du_inf, log10(iterationdata_i.mu), log10(iterationdata_i.reg), iterationdata_i.alpha_du, iterationdata_i.alpha_pr, iterationdata_i.ls, iterationdata_i.type);
                }
            }
            print_count = iterationdata.size();
        }
        void Push()
        {
            iterationdata.push_back(it_curr);
        }
        void Reset()
        {
            iterationdata.resize(0);
            print_count = 0;
        }
        int print_count = 0;
        vector<IterationData> iterationdata;
        IterationData it_curr;
    };
} // namespace fatrop
#endif //  FATROPITERATIONDATAINCLUDED