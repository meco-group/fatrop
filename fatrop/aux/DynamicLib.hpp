#ifndef DYNAMICLIBINCLUDED
#define DYNAMICLIBINCLUDED
#include <dlfcn.h>
#include <string>
#include "SmartPtr.hpp"
using namespace std;
namespace fatrop
{
    class DLHandler : public RefCountedObj
    {
    public:
        DLHandler(const string &filename);
        ~DLHandler();
        void *handle;
    };

};     // namespace fatrop
#endif // DYNAMICLIBINCLUDED