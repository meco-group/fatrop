#ifndef DYNAMICLIBINCLUDED
#define DYNAMICLIBINCLUDED
#include <dlfcn.h>
#include <string>
namespace fatrop
{
    class DLHandler 
    {
    public:
        DLHandler(const std::string &filename);
        ~DLHandler();
        void *handle;
    };

};     // namespace fatrop
#endif // DYNAMICLIBINCLUDED