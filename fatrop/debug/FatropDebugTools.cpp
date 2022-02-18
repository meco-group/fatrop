#include "debug/FatropDebugTools.hpp"
using namespace fatrop;
namespace fatrop
{
    Eig random_matrix(const int m, const int n, const int seed)
    {
        Eig res(m, n);
        std::default_random_engine e(seed);
        std::uniform_real_distribution<> dis(0, 1); // rage 0 - 1
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
            {
                {
                    res(i, j) = dis(e);
                }
            }
        return res;
    }
}
