#ifndef INTERPOLATION_SAMPLES_HPP
#define INTERPOLATION_SAMPLES_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Tridiagonal.hpp>

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

/**
 * @brief Natural-spline second derivatives at the nodes.
 *
 * Solves the same tridiagonal system CubicSpline builds, with the natural
 * condition at both ends, for a sequence of values sampled on `nodes`. It is
 * factored out here because the tensor-product grid has to solve it once per
 * row and once per column, and the caller supplies the scratch diagonals so
 * that a whole grid costs one set of buffers rather than one per line.
 *
 * @param nodes Strictly increasing abscissae.
 * @param values Sampled values, the same length as `nodes`.
 * @param curvature Output, the same length as `nodes`.
 * @param sub Scratch, the same length as `nodes`.
 * @param diag Scratch, the same length as `nodes`.
 * @param super Scratch, the same length as `nodes`.
 */
template <typename Real, typename Scalar>
void
NaturalCurvatures(std::span<const Real> nodes, std::span<const Scalar> values,
                  std::span<Scalar> curvature, std::span<Real> sub,
                  std::span<Real> diag, std::span<Real> super) {
    const auto n = nodes.size();
    assert(n >= 2);

    constexpr auto oneThird = static_cast<Real>(1) / static_cast<Real>(3);
    constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);

    std::ranges::fill(sub, Real{});
    std::ranges::fill(diag, Real{});
    std::ranges::fill(super, Real{});
    std::ranges::fill(curvature, Scalar{});

    for (std::size_t i = 1; i + 1 < n; ++i) {
        const auto hPrev = nodes[i] - nodes[i - 1];
        const auto hNext = nodes[i + 1] - nodes[i];
        sub[i] = oneSixth * hPrev;
        diag[i] = oneThird * (hPrev + hNext);
        super[i] = oneSixth * hNext;
        curvature[i] = (values[i + 1] - values[i]) / hNext -
                       (values[i] - values[i - 1]) / hPrev;
    }

    // The natural condition states directly that the end curvature is zero.
    diag[0] = 1;
    diag[n - 1] = 1;

    SolveTridiagonal<Real, Scalar>(std::span<const Real>{sub}, diag,
                                   std::span<const Real>{super}, curvature);
}

/**
 * @brief Not-a-knot second derivatives at the nodes.
 *
 * The not-a-knot condition makes the third derivative continuous across the
 * first and last interior knots, so the first two polynomial pieces are one
 * cubic and the last two are another. Unlike the natural condition it does
 * not impose anything false at the boundary, so the scheme stays fourth
 * order there and a cubic is reproduced exactly.
 *
 * It is solved in the nodal-slope formulation rather than the nodal-curvature
 * one. Written directly in curvatures the boundary row reaches outside the
 * tridiagonal band, and eliminating that entry produces a leading coefficient
 * of `h0^2 - h1^2`, which vanishes on a uniform grid. In slopes the same
 * condition gives the boundary row
 *
 * @f[
 * h_1 d_0 + (h_0 + h_1) d_1 =
 *   \frac{h_1 (2h_1 + 3h_0)\delta_0 + h_0^2 \delta_1}{h_0 + h_1},
 * @f]
 *
 * whose diagonal is @f$h_1 > 0@f$ for any spacing. The slopes are then
 * converted to curvatures so that evaluation uses the same segment formula as
 * every other spline here.
 *
 * No pivoting is used. The interior rows are diagonally dominant by a factor
 * of two; the two boundary rows are not, but the first elimination step turns
 * the second pivot into @f$h_0 + h_1@f$, and the pivots stay positive.
 *
 * @param nodes Strictly increasing abscissae, at least four of them.
 * @param values Sampled values, the same length as `nodes`.
 * @param curvature Output, the same length as `nodes`.
 * @param sub Scratch, the same length as `nodes`.
 * @param diag Scratch, the same length as `nodes`.
 * @param super Scratch, the same length as `nodes`.
 * @param slope Scratch, the same length as `nodes`.
 */
template <typename Real, typename Scalar>
void
NotAKnotCurvatures(std::span<const Real> nodes, std::span<const Scalar> values,
                   std::span<Scalar> curvature, std::span<Real> sub,
                   std::span<Real> diag, std::span<Real> super,
                   std::span<Scalar> slope) {
    const auto n = nodes.size();
    assert(n >= 4);

    const auto h = [&](std::size_t i) { return nodes[i + 1] - nodes[i]; };
    const auto secant = [&](std::size_t i) {
        return (values[i + 1] - values[i]) / h(i);
    };

    std::ranges::fill(sub, Real{});
    std::ranges::fill(diag, Real{});
    std::ranges::fill(super, Real{});
    std::ranges::fill(slope, Scalar{});

    // Interior rows: continuity of the second derivative, in slope form.
    for (std::size_t i = 1; i + 1 < n; ++i) {
        sub[i] = h(i);
        diag[i] = 2 * (h(i - 1) + h(i));
        super[i] = h(i - 1);
        slope[i] = static_cast<Real>(3) *
                   (h(i) * secant(i - 1) + h(i - 1) * secant(i));
    }

    // Left not-a-knot row.
    {
        const auto h0 = h(0);
        const auto h1 = h(1);
        diag[0] = h1;
        super[0] = h0 + h1;
        slope[0] = (h1 * (2 * h1 + 3 * h0) * secant(0) + h0 * h0 * secant(1)) /
                   (h0 + h1);
    }

    // Right not-a-knot row, the mirror of the left.
    {
        const auto hLast = h(n - 2);
        const auto hPrev = h(n - 3);
        sub[n - 1] = hLast + hPrev;
        diag[n - 1] = hPrev;
        slope[n - 1] = (hPrev * (2 * hPrev + 3 * hLast) * secant(n - 2) +
                        hLast * hLast * secant(n - 3)) /
                       (hPrev + hLast);
    }

    SolveTridiagonal<Real, Scalar>(std::span<const Real>{sub}, diag,
                                   std::span<const Real>{super}, slope);

    // Convert nodal slopes to nodal curvatures: on segment i the piece is
    // y_i + d_i t + c t^2 + e t^3, so S''(x_i) = 2c.
    for (std::size_t i = 0; i + 1 < n; ++i) {
        curvature[i] = static_cast<Real>(2) *
                       (static_cast<Real>(3) * secant(i) -
                        static_cast<Real>(2) * slope[i] - slope[i + 1]) /
                       h(i);
    }
    const auto hEnd = h(n - 2);
    curvature[n - 1] = (static_cast<Real>(2) * slope[n - 2] +
                        static_cast<Real>(4) * slope[n - 1] -
                        static_cast<Real>(6) * secant(n - 2)) /
                       hEnd;
}

} // namespace Interpolation::Detail

#endif // INTERPOLATION_SAMPLES_HPP
