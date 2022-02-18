#include "aux/DynamicLib.hpp"
namespace fatrop
{
    DLHandler::DLHandler(const string &filename)
    {
        handle = dlopen(&filename[0], RTLD_LAZY);
        if (handle == 0)
        {
            printf("Cannot open f.so, error %s\n", dlerror());
        }
    }

    DLHandler::~DLHandler()
    {
        dlclose(handle);
    }
    void *handle;
}
