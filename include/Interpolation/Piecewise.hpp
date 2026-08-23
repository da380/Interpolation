#ifndef INTERPOLATION_PIECEWISE_HPP
#define INTERPOLATION_PIECEWISE_HPP

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Function.hpp>
#include <Interpolation/Side.hpp>

namespace Interpolation {

/**
 * @brief A borrowed view of one piece, together with the interval it covers.
 *
 * This is what `Piecewise::Piece` hands back, and it is the point of the
 * whole arrangement: a single layer can be pulled out and used on its own —
 * handed to a quadrature or an ODE solver over that interval — while still
 * modelling `Function1D`, so it also composes with the rest of the algebra.
 *
 * It **borrows** the piece rather than copying it, and is valid only while
 * the `Piecewise` that produced it lives. That is not a compromise: a piece
 * built from rvalue samples holds an `owning_view`, which is move-only, so
 * there is no copy to hand out. It is also the same contract the
 * interpolators already have for lvalue ranges, and the view itself is
 * trivially cheap to copy and pass around.
 *
 * Evaluation is *not* clamped to the interval. The underlying interpolant
 * continues its end piece outside, which is what a solver taking a stage
 * evaluation slightly beyond a step needs. `Lower()` and `Upper()` report the
 * interval; they do not police it.
 *
 * @tparam F The piece type.
 */
template <Function1D F> class PieceView {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = typename F::Scalar;

    /** @brief Construct from a piece and the interval it covers. */
    constexpr PieceView(const F &f, Real lower, Real upper)
        : _f{&f}, _lower{lower}, _upper{upper} {}

    /**
     * @brief Evaluate the piece or its `N`th derivative.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa, not required to lie in the interval.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        return _f->template Evaluate<N>(x);
    }

    /** @brief Evaluate the piece; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

    /** @brief Left end of the interval this piece covers. */
    constexpr Real Lower() const { return _lower; }

    /** @brief Right end of the interval this piece covers. */
    constexpr Real Upper() const { return _upper; }

    /** @brief Width of the interval. */
    constexpr Real Width() const { return _upper - _lower; }

    /** @brief Whether `x` lies in the closed interval. */
    constexpr bool Contains(Real x) const { return x >= _lower && x <= _upper; }

    // Forwarded so that a borrowed piecewise-polynomial still satisfies
    // PiecewisePolynomial1D and Antidifferentiable1D, and can therefore be
    // handed to Primitive() without being copied. Each is available only when
    // the underlying function has it.

    /** @brief Number of nodes, when the borrowed function has them. */
    constexpr std::size_t Size() const
        requires requires(const F &f) { f.Size(); }
    {
        return _f->Size();
    }

    /** @brief Abscissa of node `i`, when the borrowed function has nodes. */
    constexpr Real Node(std::size_t i) const
        requires requires(const F &f, std::size_t j) { f.Node(j); }
    {
        return _f->Node(i);
    }

    /** @brief Segment containing `x`, when the borrowed function has one. */
    constexpr std::size_t Segment(Real x) const
        requires requires(const F &f, Real v) { f.Segment(v); }
    {
        return _f->Segment(x);
    }

    /** @brief Integral over part of segment `i`, when available. */
    constexpr Scalar SegmentIntegral(std::size_t i, Real t) const
        requires requires(const F &f, std::size_t j, Real v) {
            f.SegmentIntegral(j, v);
        }
    {
        return _f->SegmentIntegral(i, t);
    }

    /** @brief Closed-form antiderivative, when the borrowed function has one.
     */
    constexpr Scalar Antiderivative(Real x) const
        requires requires(const F &f, Real v) { f.Antiderivative(v); }
    {
        return _f->Antiderivative(x);
    }

    /** @brief The underlying piece, without its interval. */
    constexpr const F &Function() const { return *_f; }

  private:
    const F *_f;
    Real _lower;
    Real _upper;
};

/**
 * @brief A function defined by different pieces on contiguous intervals.
 *
 * Built for data that is genuinely piecewise continuous — a layered model,
 * say — where the value may jump at an interface. The base interpolators
 * require strictly increasing abscissae and would reject such data outright;
 * this is where the discontinuity is represented instead, explicitly, as
 * structure rather than as a coincidence in the samples.
 *
 * Three deliberate choices:
 *
 * - The pieces **tile** their interval: breakpoint `k + 1` ends piece `k` and
 *   begins piece `k + 1`. A gap would mean the function is undefined there,
 *   which is a different thing and is not modelled.
 * - Continuity at a breakpoint is **not checked**. Whether the pieces agree
 *   is the caller's business; enforcing it would only start an argument about
 *   tolerance.
 * - Evaluation is **right-continuous** by default, so piece `k` owns
 *   `[b[k], b[k+1])`. `Limits` returns both one-sided values, which at a real
 *   interface is usually what is actually wanted.
 *
 * All pieces share a type. That is not a limitation so much as a separation:
 * mixing kinds is type erasure, a different concern, and a type-erased
 * `Function1D` wrapper would itself be a `Function1D`, so `Piecewise` of it
 * would give mixed pieces without this class knowing anything about it.
 *
 * @tparam F The piece type.
 */
template <Function1D F> class Piecewise {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = typename F::Scalar;

    /**
     * @brief Construct from breakpoints and pieces.
     * @param breakpoints Strictly increasing, one longer than `pieces`.
     * @param pieces One function per interval, at least one.
     * @throws std::invalid_argument if there are no pieces, the lengths do
     *         not correspond, or the breakpoints are not strictly increasing.
     */
    Piecewise(std::vector<Real> breakpoints, std::vector<F> pieces)
        : _breakpoints{std::move(breakpoints)}, _pieces{std::move(pieces)} {
        if (_pieces.empty()) {
            throw std::invalid_argument("Piecewise: needs at least one piece");
        }
        if (_breakpoints.size() != _pieces.size() + 1) {
            throw std::invalid_argument(
                "Piecewise: expected " + std::to_string(_pieces.size() + 1) +
                " breakpoints for " + std::to_string(_pieces.size()) +
                " pieces, but was given " +
                std::to_string(_breakpoints.size()));
        }
        const auto offending = std::ranges::adjacent_find(
            _breakpoints,
            [](const auto &lhs, const auto &rhs) { return !(lhs < rhs); });
        if (offending != _breakpoints.end()) {
            const auto position = static_cast<std::size_t>(
                std::ranges::distance(_breakpoints.begin(), offending));
            throw std::invalid_argument(
                "Piecewise: breakpoints must be strictly increasing, but "
                "element " +
                std::to_string(position) + " is not less than element " +
                std::to_string(position + 1));
        }
    }

    /** @brief Number of pieces. */
    std::size_t PieceCount() const { return _pieces.size(); }

    /** @brief Breakpoint `k`, of which there are `PieceCount() + 1`. */
    Real Breakpoint(std::size_t k) const { return _breakpoints[k]; }

    /** @brief Left end of the whole domain. */
    Real Lower() const { return _breakpoints.front(); }

    /** @brief Right end of the whole domain. */
    Real Upper() const { return _breakpoints.back(); }

    /**
     * @brief Index of the piece answering for `x`.
     *
     * Outside the domain the first or last piece is named, matching how the
     * interpolators continue their end segments.
     *
     * @param x Query abscissa.
     * @param side Which piece owns a query landing exactly on a breakpoint.
     */
    std::size_t IndexOf(Real x, Side side = Side::Right) const {
        const auto begin = _breakpoints.begin();
        const auto found = side == Side::Right
                               ? std::ranges::upper_bound(_breakpoints, x)
                               : std::ranges::lower_bound(_breakpoints, x);
        const auto raw = std::ranges::distance(begin, found);
        if (raw <= 0) {
            return 0;
        }
        return std::min(static_cast<std::size_t>(raw - 1), PieceCount() - 1);
    }

    /**
     * @brief The piece on interval `k`, with the interval it belongs to.
     *
     * This is usually more useful than evaluating through the whole object:
     * a single layer can be handed to a solver that will work over exactly
     * that interval. The result borrows, so it is valid while this object
     * lives.
     */
    PieceView<F> Piece(std::size_t k) const {
        return PieceView<F>{_pieces[k], _breakpoints[k], _breakpoints[k + 1]};
    }

    /** @brief The piece answering for `x`, with its interval. */
    PieceView<F> PieceAt(Real x, Side side = Side::Right) const {
        return Piece(IndexOf(x, side));
    }

    /** @brief All pieces, each with its interval. */
    std::vector<PieceView<F>> Pieces() const {
        std::vector<PieceView<F>> all;
        all.reserve(PieceCount());
        for (std::size_t k = 0; k < PieceCount(); ++k) {
            all.push_back(Piece(k));
        }
        return all;
    }

    /**
     * @brief Evaluate the function or its `N`th derivative, right-continuous.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> Scalar Evaluate(Real x) const {
        return EvaluateFrom<N>(x, Side::Right);
    }

    /**
     * @brief Evaluate from a chosen side of a breakpoint.
     * @tparam N Derivative order.
     * @param x Query abscissa.
     * @param side Which piece answers when `x` is exactly a breakpoint.
     */
    template <std::size_t N = 0> Scalar EvaluateFrom(Real x, Side side) const {
        return _pieces[IndexOf(x, side)].template Evaluate<N>(x);
    }

    /**
     * @brief Both one-sided values at `x`.
     *
     * At a breakpoint these differ exactly when the function jumps there,
     * which for a layered model is normally the quantity of interest. Away
     * from a breakpoint they agree.
     *
     * @tparam N Derivative order.
     * @param x Query abscissa.
     * @return `{from the left, from the right}`.
     */
    template <std::size_t N = 0>
    std::pair<Scalar, Scalar> Limits(Real x) const {
        return {EvaluateFrom<N>(x, Side::Left),
                EvaluateFrom<N>(x, Side::Right)};
    }

    /** @brief Evaluate the function; the same as `Evaluate<0>`. */
    Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    std::vector<Real> _breakpoints;
    std::vector<F> _pieces;
};

/**
 * @brief Build a piecewise function from samples that repeat an abscissa.
 *
 * A doubled abscissa is one common convention for flagging a discontinuity in
 * tabulated data. The base interpolators reject it, deliberately, because
 * every one of them divides by a node difference somewhere. This is the
 * bridge: it splits the samples wherever an abscissa repeats and hands each
 * run to `build`, which returns the piece for that interval.
 *
 * @code
 * auto layered = SplitAtRepeats(radius, density, [](auto r, auto d) {
 *     return CubicSpline{std::move(r), std::move(d)};
 * });
 * auto mantle = layered.Piece(2);          // usable on its own
 * auto jump   = layered.Limits(cmbRadius); // {above, below}
 * @endcode
 *
 * @param x Non-decreasing abscissae; a repeat marks a breakpoint.
 * @param y Ordinates, the same length as `x`.
 * @param build Called with the owned abscissae and ordinates of each run.
 * @throws std::invalid_argument if the lengths differ, the abscissae
 *         decrease, or an abscissa appears more than twice.
 */
template <typename XRange, typename YRange, typename Builder>
auto
SplitAtRepeats(const XRange &x, const YRange &y, Builder build) {
    using Real = std::ranges::range_value_t<XRange>;
    using Scalar = std::ranges::range_value_t<YRange>;
    using Piece = decltype(build(std::vector<Real>{}, std::vector<Scalar>{}));

    const auto n = static_cast<std::size_t>(std::ranges::size(x));
    if (n != static_cast<std::size_t>(std::ranges::size(y))) {
        throw std::invalid_argument(
            "SplitAtRepeats: abscissa and ordinate ranges differ in length");
    }
    if (n == 0) {
        throw std::invalid_argument("SplitAtRepeats: no samples");
    }

    std::vector<Real> breakpoints;
    std::vector<Piece> pieces;
    breakpoints.push_back(x[0]);

    std::vector<Real> runX;
    std::vector<Scalar> runY;

    const auto flush = [&](Real upper) {
        pieces.push_back(build(std::move(runX), std::move(runY)));
        breakpoints.push_back(upper);
        runX.clear();
        runY.clear();
    };

    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0) {
            if (x[i] < x[i - 1]) {
                throw std::invalid_argument(
                    "SplitAtRepeats: abscissae must be non-decreasing, but "
                    "element " +
                    std::to_string(i) + " is less than element " +
                    std::to_string(i - 1));
            }
            if (x[i] == x[i - 1]) {
                if (i + 1 < n && x[i + 1] == x[i]) {
                    throw std::invalid_argument(
                        "SplitAtRepeats: abscissa at element " +
                        std::to_string(i) +
                        " appears more than twice, so the intended pieces are "
                        "ambiguous");
                }
                flush(x[i - 1]);
            }
        }
        runX.push_back(x[i]);
        runY.push_back(y[i]);
    }
    flush(x[n - 1]);

    return Piecewise<Piece>{std::move(breakpoints), std::move(pieces)};
}

} // namespace Interpolation

#endif // INTERPOLATION_PIECEWISE_HPP
