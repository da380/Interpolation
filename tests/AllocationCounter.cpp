#include "AllocationCounter.hpp"

#include <atomic>
#include <cstdlib>
#include <new>

namespace {
std::atomic<std::size_t> g_allocations{0};
}

namespace InterpolationTest {

std::size_t
AllocationCount() {
    return g_allocations.load(std::memory_order_relaxed);
}

} // namespace InterpolationTest

#if INTERPOLATION_TEST_COUNTS_ALLOCATIONS

void *
operator new(std::size_t size) {
    g_allocations.fetch_add(1, std::memory_order_relaxed);
    if (void *memory = std::malloc(size == 0 ? 1 : size)) {
        return memory;
    }
    throw std::bad_alloc{};
}

void *
operator new[](std::size_t size) {
    return ::operator new(size);
}

void
operator delete(void *memory) noexcept {
    std::free(memory);
}

void
operator delete[](void *memory) noexcept {
    std::free(memory);
}

void
operator delete(void *memory, std::size_t) noexcept {
    std::free(memory);
}

void
operator delete[](void *memory, std::size_t) noexcept {
    std::free(memory);
}

#endif // INTERPOLATION_TEST_COUNTS_ALLOCATIONS
