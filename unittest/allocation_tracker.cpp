// Test-only Linux/glibc interposition, linked into an executable with
// --export-dynamic so calls from shared libraries are counted as well.
#include "allocation_tracker.hpp"
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <execinfo.h>
#include <unistd.h>

extern "C"
{
void *__libc_malloc(std::size_t) noexcept;
void *__libc_calloc(std::size_t, std::size_t) noexcept;
void *__libc_realloc(void *, std::size_t) noexcept;
void *__libc_memalign(std::size_t, std::size_t) noexcept;
}

namespace
{
std::atomic<bool> tracking{false};
std::atomic<std::size_t> calls{0}, bytes{0};
bool tracing = false;
thread_local bool inside_trace = false;

void record(const char *operation, std::size_t size)
{
    if (!tracking.load(std::memory_order_relaxed) || inside_trace) return;
    const auto number = calls.fetch_add(1, std::memory_order_relaxed);
    bytes.fetch_add(size, std::memory_order_relaxed);
    if (!tracing || number >= 40) return;
    inside_trace = true;
    char message[128];
    const int n = std::snprintf(message, sizeof(message),
        "\nallocation %zu: %s(%zu)\n", number + 1, operation, size);
    if (n > 0) (void)::write(STDERR_FILENO, message, static_cast<std::size_t>(n));
    void *stack[24];
    const int depth = backtrace(stack, 24);
    backtrace_symbols_fd(stack, depth, STDERR_FILENO);
    inside_trace = false;
}
}

extern "C" void *malloc(std::size_t size) noexcept
{
    record("malloc", size);
    return __libc_malloc(size);
}
extern "C" void *calloc(std::size_t count, std::size_t size) noexcept
{
    record("calloc", count * size);
    return __libc_calloc(count, size);
}
extern "C" void *realloc(void *pointer, std::size_t size) noexcept
{
    if (size) record("realloc", size);
    return __libc_realloc(pointer, size);
}
extern "C" void *memalign(std::size_t alignment, std::size_t size) noexcept
{
    record("memalign", size);
    return __libc_memalign(alignment, size);
}
extern "C" void *aligned_alloc(std::size_t alignment, std::size_t size) noexcept
{
    record("aligned_alloc", size);
    return __libc_memalign(alignment, size);
}
extern "C" int posix_memalign(void **pointer, std::size_t alignment, std::size_t size) noexcept
{
    if (alignment < sizeof(void *) || (alignment & (alignment - 1))) return EINVAL;
    record("posix_memalign", size);
    void *value = __libc_memalign(alignment, size);
    if (!value) return ENOMEM;
    *pointer = value;
    return 0;
}

namespace fatrop::test::allocation_tracker
{
void start(bool trace)
{
    tracking.store(false, std::memory_order_relaxed);
    tracing = trace;
    calls.store(0, std::memory_order_relaxed);
    bytes.store(0, std::memory_order_relaxed);
    tracking.store(true, std::memory_order_relaxed);
}
Counts stop()
{
    tracking.store(false, std::memory_order_relaxed);
    return {calls.load(std::memory_order_relaxed), bytes.load(std::memory_order_relaxed)};
}
}
