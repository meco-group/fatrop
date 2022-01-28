#ifndef FATROPPARAMSINCLUDED
#define FATROPPARAMSINCLUDED
#include "AUX/SmartPtr.hpp"
namespace fatrop
{
    class FatropParams: public RefCountedObj
    {
        public:
        int maxiter = 100;
        double smax = 100.0; 
        double lammax = 1e3;
        double tol = 1e-8;
    };

} // namespace fatrop
#endif // FatropParams