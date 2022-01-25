#ifndef FATROPPARAMSINCLUDED
#define FATROPPARAMSINCLUDED
#include "AUX/SmartPtr.hpp"
namespace fatrop
{
    class FatropParams: public RefCountedObj
    {
        public:
        int maxiter = 100;
        double smax = 0.0;
    };

} // namespace fatrop
#endif // FatropParams