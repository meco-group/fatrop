#include "fatrop/ip_algorithm/ip_utils.hpp"
#include <limits> // for std::numeric_limits<Scalar>::epsilon()

namespace fatrop
{
    namespace internal
    {
        bool compare_le(const Scalar a, const Scalar b, const Scalar base)
        {
            Scalar eps = std::numeric_limits<Scalar>::epsilon();
            return (a - b <= 10. * eps * std::abs(base));
        }
    };
} // namespace fatrop
