#ifndef INTERPOLATION_SAMPLES_HPP
#define INTERPOLATION_SAMPLES_HPP

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <string>

#include <Interpolation/Concepts.hpp>

namespace Interpolation::Detail {

/**
 * @brief Validate sample ranges at construction.
 *
 * Preconditions on the data are checked here and reported by throwing, rather
 * than asserted. Construction is not the hot path, and an assert vanishes
 * under NDEBUG, which silently turns a caught mistake into undefined
 * behaviour in exactly the builds where it matters most.
 *
 * Abscissae must be strictly increasing. A repeated abscissa is not treated as
 * a discontinuity: every interpolator here divides by a node difference
 * somewhere, and for the barycentric Lagrange basis that division happens
 * during setup. Piecewise-continuous functions are represented by composing
 * several interpolators rather than by encoding a jump as a doubled node.
 *
 * @param x Abscissa range.
 * @param y Ordinate range.
 * @param minimumNodes Fewest nodes the caller's algorithm can use.
 * @param what Interpolator name, used in the exception message.
 * @throws std::invalid_argument if the lengths differ, there are too few
 *         nodes, or the abscissae are not strictly increasing.
 */
template <typename XRange, typename YRange>
void
ValidateSamples(const XRange &x, const YRange &y, std::size_t minimumNodes,
                const char *what) {
    const auto nx = static_cast<std::size_t>(std::ranges::size(x));
    const auto ny = static_cast<std::size_t>(std::ranges::size(y));

    if (nx != ny) {
        throw std::invalid_argument(
            std::string{what} + ": abscissa and ordinate ranges differ in " +
            "length (" + std::to_string(nx) + " and " + std::to_string(ny) +
            ")");
    }

    if (nx < minimumNodes) {
        throw std::invalid_argument(std::string{what} + ": needs at least " +
                                    std::to_string(minimumNodes) +
                                    " nodes, but was given " +
                                    std::to_string(nx));
    }

    const auto offending = std::ranges::adjacent_find(
        x, [](const auto &lhs, const auto &rhs) { return !(lhs < rhs); });
    if (offending != std::ranges::end(x)) {
        const auto position = static_cast<std::size_t>(
            std::ranges::distance(std::ranges::begin(x), offending));
        throw std::invalid_argument(
            std::string{what} +
            ": abscissae must be strictly increasing, but "
            "element " +
            std::to_string(position) + " is not less than element " +
            std::to_string(position + 1));
    }
}

/**
 * @brief Index of the segment used to evaluate a query.
 *
 * Returns `i` such that the polynomial piece on `[x[i], x[i + 1]]` is the one
 * to use. Queries outside the sample interval select the first or final
 * segment, so evaluation there continues that piece rather than failing.
 *
 * This replaces the same upper-bound-and-clamp written out separately in each
 * interpolator.
 *
 * @param x Strictly increasing abscissae, at least two of them.
 * @param query Query abscissa.
 * @return Segment index in `[0, size(x) - 1)`.
 */
template <typename XRange, typename Real>
std::size_t
LocateSegment(const XRange &x, Real query) {
    const auto n = static_cast<std::size_t>(std::ranges::size(x));
    const auto upper = std::ranges::upper_bound(x, query);

    if (upper == std::ranges::begin(x)) {
        return 0;
    }
    if (upper == std::ranges::end(x)) {
        return n - 2;
    }
    return static_cast<std::size_t>(
               std::ranges::distance(std::ranges::begin(x), upper)) -
           1;
}

} // namespace Interpolation::Detail

#endif // INTERPOLATION_SAMPLES_HPP
