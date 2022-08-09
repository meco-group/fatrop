#ifndef DYNAMICLIBINCLUDED
#define DYNAMICLIBINCLUDED
#include <dlfcn.h>
#include <string>
#include "SmartPtr.hpp"
using namespace std;
namespace fatrop
{
    class DLHandler 
    {
    public:
        DLHandler(const string &filename);
        ~DLHandler();
        void *handle;
    };

};     // namespace fatrop
#endif // DYNAMICLIBINCLUDED