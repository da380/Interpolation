#ifndef INTERPOLATION_SAMPLES_HPP
#define INTERPOLATION_SAMPLES_HPP

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <string>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Side.hpp>

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

/**
 * @brief Check that a nodal output range is the right length.
 *
 * Shared by the `EvaluateAtNodes` of every interpolator, so the message is
 * written once and reads the same wherever it comes from.
 *
 * @param given Length the caller supplied.
 * @param expected Number of nodes.
 * @param what Interpolator name, used in the exception message.
 * @throws std::invalid_argument if the lengths differ.
 */
inline void
ValidateNodeOutput(std::size_t given, std::size_t expected, const char *what) {
    if (given != expected) {
        throw std::invalid_argument(
            std::string{what} +
            "::EvaluateAtNodes: the output range has length " +
            std::to_string(given) + " but there are " +
            std::to_string(expected) + " nodes");
    }
}

/**
 * @brief Index of the segment that owns node `k`, from a chosen side.
 *
 * The counterpart of LocateSegment for a query that is known to land exactly
 * on a node, so there is nothing to search for. `Side::Right` reproduces
 * LocateSegment exactly, including its treatment of the final node, which has
 * no segment to its right and so is answered from the one to its left;
 * `Side::Left` mirrors it at the first node.
 *
 * The distinction matters only for a derivative order the interpolant does not
 * carry across a knot — the third for a cubic, the first for a linear
 * interpolant. Below that order the two sides agree to rounding.
 *
 * @param n Number of nodes, at least two.
 * @param k Node index.
 * @param side Which adjoining segment to use.
 * @return Segment index in `[0, n - 1)`.
 */
constexpr std::size_t
NodeSegment(std::size_t n, std::size_t k, Side side) {
    if (side == Side::Right) {
        return k + 1 < n ? k : n - 2;
    }
    return k > 0 ? k - 1 : 0;
}

/**
 * @brief Value or `N`th derivative of a linear piece.
 *
 * The piece spans `[0, h]` with end values `f0` and `f1`, evaluated at offset
 * `t` from its left end. Shared by the one-dimensional interpolant and by the
 * tensor-product grids, so the formula exists once.
 */
template <std::size_t N, typename Real, typename Scalar>
constexpr Scalar
LinearPiece(Real h, Real t, Scalar f0, Scalar f1) {
    if constexpr (N > 1) {
        return Scalar{};
    } else if constexpr (N == 1) {
        return (f1 - f0) / h;
    } else {
        const auto b = t / h;
        const auto a = static_cast<Real>(1) - b;
        return a * f0 + b * f1;
    }
}

/**
 * @brief Value or `N`th derivative of a cubic-spline piece.
 *
 * The piece spans `[0, h]` with end values `f0`, `f1` and end second
 * derivatives `m0`, `m1`, evaluated at offset `t` from its left end. Pieces
 * are cubic, so the third derivative is constant and higher orders vanish.
 */
template <std::size_t N, typename Real, typename Scalar>
constexpr Scalar
SplinePiece(Real h, Real t, Scalar f0, Scalar f1, Scalar m0, Scalar m1) {
    if constexpr (N > 3) {
        return Scalar{};
    } else if constexpr (N == 3) {
        return (m1 - m0) / h;
    } else {
        constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);
        const auto b = t / h;
        const auto a = static_cast<Real>(1) - b;

        if constexpr (N == 2) {
            return a * m0 + b * m1;
        } else if constexpr (N == 1) {
            return (f1 - f0) / h +
                   oneSixth * h * ((1 - 3 * a * a) * m0 + (3 * b * b - 1) * m1);
        } else {
            return a * f0 + b * f1 +
                   ((a * a * a - a) * m0 + (b * b * b - b) * m1) * h * h *
                       oneSixth;
        }
    }
}

} // namespace Interpolation::Detail

#endif // INTERPOLATION_SAMPLES_HPP
