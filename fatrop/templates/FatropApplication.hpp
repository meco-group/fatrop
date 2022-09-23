#ifndef FATROPAPPLICATIONINCLUDED
#define FATROPAPPLICATIONINCLUDED
using namespace std;
#include <vector>
namespace fatrop
{
    class FatropApplication
    {
        public:
        virtual void Initialize() = 0;
        virtual void Reset() = 0;
        virtual void SetBounds(const vector<double>& lower, const vector<double>& upper) = 0;
        virtual void SetInitial(const vector<double>& initial) = 0;
        virtual void GetSolution(vector<double>& sol) = 0;
        virtual int Optimize() = 0;
        // virtual void SetOption = 0;
    }; 
}
#endif // FATROPAPPLICATIONINCLUDED