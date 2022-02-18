#ifndef DYNAMICLIBINCLUDED
#define DYNAMICLIBINCLUDED
#include <dlfcn.h>
#include <string>
#include "SmartPtr.hpp"
using namespace std;
namespace fatrop
{
    class DLHandler: public RefCountedObj
    {
    public:
        DLHandler(const string &filename)
        {
            handle = dlopen(&filename[0], RTLD_LAZY);
            if (handle == 0)
            {
                printf("Cannot open f.so, error %s\n", dlerror());
            }
        }
        ~DLHandler()
        {
            dlclose(handle);
        }
        void *handle;
    };
};     // namespace fatrop
#endif // DYNAMICLIBINCLUDED