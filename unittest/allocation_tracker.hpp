#ifndef FATROP_TEST_ALLOCATION_TRACKER_HPP
#define FATROP_TEST_ALLOCATION_TRACKER_HPP
#include <cstddef>

namespace fatrop::test::allocation_tracker
{
struct Counts
{
    std::size_t calls = 0;
    std::size_t bytes = 0;
};

// Test-only glibc allocation interposition. Counts malloc/calloc/realloc and
// aligned allocation calls, including C++ new and calls originating in shared
// libfatrop / BLASFEO. Construct the complete solver before start(); do not run
// optimize() or a factorization to warm it up. Assertions belong after stop().
// trace=true additionally prints at most 40 allocation backtraces to stderr.
void start(bool trace = false);
Counts stop();
}
#endif
