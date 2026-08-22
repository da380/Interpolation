#ifndef INTERPOLATION_TEST_ALLOCATION_COUNTER_HPP
#define INTERPOLATION_TEST_ALLOCATION_COUNTER_HPP

#include <cstddef>

// AddressSanitizer replaces the global operator new itself, including the
// nothrow overloads, so a second replacement here means memory allocated by
// the sanitizer is released by ours and every test trips an
// alloc-dealloc-mismatch. Counting stands down under sanitizers; the
// allocation claim is still checked in the ordinary builds.
#if defined(__SANITIZE_ADDRESS__)
#define INTERPOLATION_TEST_COUNTS_ALLOCATIONS 0
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define INTERPOLATION_TEST_COUNTS_ALLOCATIONS 0
#else
#define INTERPOLATION_TEST_COUNTS_ALLOCATIONS 1
#endif
#else
#define INTERPOLATION_TEST_COUNTS_ALLOCATIONS 1
#endif

namespace InterpolationTest {

/** @brief Whether the global operator new replacement is compiled in. */
constexpr bool
AllocationCountingEnabled() {
    return INTERPOLATION_TEST_COUNTS_ALLOCATIONS != 0;
}

/**
 * @brief Number of allocations made so far through the global operator new.
 *
 * The counter lives in its own translation unit. Replacing the global
 * operator new in the same file that also constructs containers lets the
 * compiler see a standard allocation paired with our free and report a
 * mismatch, which is a false positive but a noisy one.
 */
std::size_t AllocationCount();

} // namespace InterpolationTest

#endif // INTERPOLATION_TEST_ALLOCATION_COUNTER_HPP
