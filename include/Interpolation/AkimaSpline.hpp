#ifndef INTERPOLATION_AKIMA_SPLINE_HPP
#define INTERPOLATION_AKIMA_SPLINE_HPP

#include <cmath>
#include <cstddef>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Samples.hpp>
#include <Interpolation/Side.hpp>

namespace Interpolation {

/**
 * @brief Piecewise-cubic interpolation using locally weighted secant slopes.
 *
 * Construction takes ranges, borrowing lvalues and owning rvalues. Queries
 * outside the sample interval continue the first or final cubic piece.
 *
 * @tparam XView View over real abscissae.
 * @tparam YView View over real or complex ordinates.
 */
template <typename XView, typename YView>
    requires InterpolationRanges<XView, YView> && std::ranges::view<XView> &&
             std::ranges::view<YView>
class AkimaSpline {
  public:
    /** @brief Abscissa precision. */
    using Real = std::ranges::range_value_t<XView>;
    /** @brief Ordinate type, real or complex. */
    using Scalar = std::ranges::range_value_t<YView>;

    /**
     * @brief Construct from abscissa and ordinate ranges.
     * @param x Strictly increasing abscissae, at least three.
     * @param y Ordinates, the same length as `x`.
     * @throws std::invalid_argument if the ranges differ in length, are too
     *         short, or the abscissae are not strictly increasing.
     */
    AkimaSpline(XView x, YView y) : _x{std::move(x)}, _y{std::move(y)} {
        Detail::ValidateSamples(_x, _y, 3, "AkimaSpline");
        ComputeSlopes();
    }

    /** @brief Number of interpolation nodes. */
    std::size_t Size() const {
        return static_cast<std::size_t>(std::ranges::size(_x));
    }

    /**
     * @brief Evaluate the interpolant or its `N`th derivative.
     *
     * The pieces are cubic, so derivatives of order four and above are
     * identically zero and the third is piecewise constant.
     *
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> Scalar Evaluate(Real x) const {
        if constexpr (N > 3) {
            return Scalar{};
        } else {
            const auto i = Detail::LocateSegment(_x, x);
            return Piece<N>(i, x - _x[i]);
        }
    }

    /** @brief Evaluate the interpolant; the same as `Evaluate<0>`. */
    Scalar operator()(Real x) const { return Evaluate<0>(x); }

    /**
     * @brief Evaluate the interpolant, or its `N`th derivative, at every node.
     *
     * Writes one value per node into `out`, in node order, without a segment
     * search: the segment adjoining each node is known.
     *
     * The third derivative jumps across an interior knot; `side` chooses which
     * limit is reported, and defaults to the right-hand one, so that
     * `EvaluateAtNodes<N>(out)` agrees with `Evaluate<N>(Node(k))` for every
     * `k`. Lower orders agree from both sides to rounding.
     *
     * @tparam N Derivative order; `0` is the value itself.
     * @param out Output, one per node. Overwritten.
     * @param side Which adjoining segment to answer from at a node.
     * @throws std::invalid_argument if `out` has the wrong length.
     */
    template <std::size_t N = 0>
    void EvaluateAtNodes(std::span<Scalar> out, Side side = Side::Right) const {
        const auto n = Size();
        Detail::ValidateNodeOutput(out.size(), n, "AkimaSpline");
        for (std::size_t k = 0; k < n; ++k) {
            if constexpr (N > 3) {
                out[k] = Scalar{};
            } else {
                const auto i = Detail::NodeSegment(n, k, side);
                out[k] = Piece<N>(i, _x[k] - _x[i]);
            }
        }
    }

    /**
     * @brief The nodal values or `N`th derivatives, in a fresh vector.
     *
     * The convenient form of EvaluateAtNodes. It allocates, so a caller
     * sweeping many lines should keep one buffer and call EvaluateAtNodes.
     */
    template <std::size_t N = 0>
    std::vector<Scalar> NodeValues(Side side = Side::Right) const {
        std::vector<Scalar> values(Size());
        EvaluateAtNodes<N>(std::span<Scalar>{values}, side);
        return values;
    }

    /** @brief Abscissa of node `i`. */
    Real Node(std::size_t i) const { return _x[i]; }

    /** @brief Index of the segment used to evaluate `x`. */
    std::size_t Segment(Real x) const { return Detail::LocateSegment(_x, x); }

    /**
     * @brief Integral of segment `i` from its left node over a width `t`.
     *
     * The piece is the cubic `y + s t + c t^2 + d t^3`, so the integral is
     * `y t + s t^2/2 + c t^3/3 + d t^4/4`.
     */
    Scalar SegmentIntegral(std::size_t i, Real t) const {
        const auto h = _x[i + 1] - _x[i];
        const auto two = static_cast<Scalar>(2);
        const auto three = static_cast<Scalar>(3);
        const auto c = (three * _m[i] - two * _s[i] - _s[i + 1]) / h;
        const auto d = (_s[i] + _s[i + 1] - two * _m[i]) / (h * h);
        return _y[i] * t + _s[i] * (t * t) / static_cast<Real>(2) +
               c * (t * t * t) / static_cast<Real>(3) +
               d * (t * t * t * t) / static_cast<Real>(4);
    }

  private:
    XView _x;
    YView _y;
    std::vector<Scalar> _m; // Secant slopes, one per segment.
    std::vector<Scalar> _s; // Akima-weighted slopes, one per node.

    // Value or Nth derivative of the cubic on segment i, at offset t from its
    // left node. The two quadratic and cubic coefficients follow from the end
    // values and the Akima slopes, so they are recovered here rather than
    // stored.
    template <std::size_t N> Scalar Piece(std::size_t i, Real t) const {
        const auto h = _x[i + 1] - _x[i];
        const auto two = static_cast<Scalar>(2);
        const auto three = static_cast<Scalar>(3);
        const auto c = (three * _m[i] - two * _s[i] - _s[i + 1]) / h;
        const auto d = (_s[i] + _s[i + 1] - two * _m[i]) / (h * h);

        if constexpr (N == 3) {
            return static_cast<Scalar>(6) * d;
        } else if constexpr (N == 2) {
            return two * c + static_cast<Scalar>(6) * d * t;
        } else if constexpr (N == 1) {
            return _s[i] + t * (two * c + three * d * t);
        } else {
            return _y[i] + t * (_s[i] + t * (c + d * t));
        }
    }

    void ComputeSlopes() {
        const auto n = Size();

        _m.reserve(n - 1);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            _m.push_back((_y[i + 1] - _y[i]) / (_x[i + 1] - _x[i]));
        }

        // The weights are magnitudes, so they stay real for complex ordinates.
        using Weight = decltype(std::abs(Scalar{}));
        constexpr auto oneHalf = static_cast<Weight>(0.5);

        _s.reserve(n);
        _s.push_back(_m[0]);
        _s.push_back((_m[0] + _m[1]) * oneHalf);
        for (std::size_t i = 2; i + 2 < n; ++i) {
            const auto a = std::abs(_m[i + 1] - _m[i]);
            const auto b = std::abs(_m[i - 1] - _m[i - 2]);
            const auto weightSum = a + b;
            if (weightSum == Weight{}) {
                _s.push_back((_m[i - 1] + _m[i]) * oneHalf);
            } else {
                _s.push_back((a * _m[i - 1] + b * _m[i]) / weightSum);
            }
        }
        if (n > 3) {
            _s.push_back((_m[n - 3] + _m[n - 2]) * oneHalf);
        }
        _s.push_back(_m[n - 2]);
    }
};

/// Borrow lvalue containers and own rvalue ones.
template <std::ranges::viewable_range X, std::ranges::viewable_range Y>
AkimaSpline(X &&,
            Y &&) -> AkimaSpline<std::views::all_t<X>, std::views::all_t<Y>>;

} // namespace Interpolation

#endif // INTERPOLATION_AKIMA_SPLINE_HPP
