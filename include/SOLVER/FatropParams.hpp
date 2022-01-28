#ifndef FATROPPARAMSINCLUDED
#define FATROPPARAMSINCLUDED
#include "AUX/SmartPtr.hpp"
namespace fatrop
{
    class FatropParams: public RefCountedObj
    {
        public:
        int maxiter = 10;
        double smax = 100.0; 
        double lammax = 1e3;
        double tol = 1e-8;
        double mu0 = 1e-1;
        double kappa_eta = 10.0; 
        double kappa_mu = 0.2;
        double theta_mu = 1.5;
        double delta_w0 = 1e-4;
        double delta_wmin = 1e-20;
        double kappa_wmin = 1.0/3.0;
        double kappa_wplus = 8;
        double kappa_wplusem = 100;
    };

} // namespace fatrop
#endif // FatropParams