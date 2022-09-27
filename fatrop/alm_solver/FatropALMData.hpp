#ifndef FATROPALMDATAINCLUDED
#define FATROPALMDATAINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "templates/NLPAlg.hpp"
#include <algorithm>
using namespace std;
namespace fatrop
{
    class FatropALMData
    {
    public:
        FatropALMData(int n_ineqs, int n_vars) :n_ineqs_(n_ineqs), n_vars(n_vars), inequalitymem(n_ineqs, 5), lower_bounds(inequalitymem[0]), upper_bounds(inequalitymem[1]), ineq_curr(inequalitymem[2]), auglags_Lcurr(inequalitymem[3]), auglags_Ucurr(inequalitymem[4]), initial_x(n_vars)
        {
        };
        int UpdateLags(double penalty_curr)
        {
            double* ineq_curr_p = ((VEC*) ineq_curr)->pa;
            double* lower_p = ((VEC*) lower_bounds)->pa;
            double* upper_p = ((VEC*) upper_bounds)->pa;
            double* auglagsL_p = ((VEC*) auglags_Lcurr)->pa;
            double* auglagsU_p = ((VEC*) auglags_Ucurr)->pa;
            double max_dist = 0.0;
            double lag_max = 1e2;
            for(int i = 0; i<n_ineqs_; i ++)
            {
                double dist_lowi = ineq_curr_p[i]  - lower_p[i];
                double dist_upperi = upper_p[i] - ineq_curr_p[i];
                auglagsL_p[i] = std::max(0.0, std::min(auglagsL_p[i]- penalty_curr*dist_lowi, lag_max));
                auglagsU_p[i] = std::max(0.0, std::min(auglagsU_p[i]- penalty_curr*dist_upperi, lag_max));
                max_dist = std::max(std::max(max_dist, -dist_lowi), -dist_upperi);
            }
            cout << "maximum violation " << max_dist << endl;
            return 0;
        }
        const int n_ineqs_;
        const int n_vars;
        FatropMemoryVecBF inequalitymem;
        FatropVecBF lower_bounds;
        FatropVecBF upper_bounds;
        FatropVecBF ineq_curr;
        FatropVecBF auglags_Lcurr;
        FatropVecBF auglags_Ucurr;
        vector<double> initial_x; 
    };
}
#endif // FATROPALMDATAINCLUDED