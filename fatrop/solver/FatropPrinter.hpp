#define PRIORITY1 PrintPriority<1>()
#ifndef FATROPPRINTERINCLUDED
#define FATROPPRINTERINCLUDED
#include <iostream>
#include <string>
namespace fatrop
{
    static std::ostream nullstream(nullptr);
    template <int Priority>
    struct PrintPriority
    {
        PrintPriority(){};
        const int priority = Priority;
    };
    template <int priority>
    std::ostream &operator<<(std::ostream &os, const PrintPriority<priority> &p)
    {
        if (p.priority >= -1)
            return os;
        return nullstream;
    }
    class FatropPrinter
    {
    public:
        FatropPrinter(const int priority = 0, std::ostream &stream = std::cout) : priority_(priority), stream_(stream){};
        std::ostream &
        level(const int p)
        {
            if (p <= priority_)
                return stream_;
            return nullstream;
        }
        int &print_level() { return priority_; }

    private:
        int priority_;
        std::ostream& stream_;
    };

} // namespace fatrop
#endif //  FATROPITERATIONDATAINCLUDED