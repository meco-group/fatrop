#include "ocp/OCPLSRiccati.hpp"

namespace fatrop
{
    bool check_reg(const int m, MAT *sA, const int ai, const int aj)
    {
        for (int i = 0; i < m; i++)
        {
            if (MATEL(sA, ai + i, aj + i) < 1e-6)
                return false;
        }
        return true;
    }
} // namespace fatrop
