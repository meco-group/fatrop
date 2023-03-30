#ifndef FATROPITERATIONDATAINCLUDED
#define FATROPITERATIONDATAINCLUDED
#include <cmath>
#include <vector>
#include <iostream>
#include <solver/FatropPrinter.hpp>
namespace fatrop
{
    struct IterationData
    {
        int iter = 0;
        double mu = 0.0;
        double objective = 0.0;
        double constraint_violation = 0.0;
        double du_inf = 0.0;
        int ls = 0;
        double reg = 0.0;
        double alpha_pr = 0.0;
        double alpha_du = 0.0;
        char type = 'x';
    };
    class Journaller
    {
    public:
        Journaller(const int maxiter);
        void print_iterations();
        void push();
        void reset();
        int print_count = 0;
        std::vector<IterationData> iterationdata;
        IterationData it_curr;
    };
} // namespace fatrop
#endif //  FATROPITERATIONDATAINCLUDED