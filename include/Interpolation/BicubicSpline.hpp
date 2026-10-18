#ifndef INTERPOLATION_BICUBIC_SPLINE_HPP
#define INTERPOLATION_BICUBIC_SPLINE_HPP

#include <algorithm>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/CubicSplineSystem.hpp>
#include <Interpolation/Grid.hpp>
#include <Interpolation/Samples.hpp>

namespace Interpolation {

/**
 * @brief The endpoint conditions along one axis of a grid.
 *
 * The same pairing rules hold as in one dimension: NotAKnot and Periodic
 * constrain the whole axis, so they are used at both ends or neither.
 */
struct AxisConditions {
    /** Condition at the first node of the axis. */
    BoundaryCondition low;
    /** Condition at the final node of the axis. */
    BoundaryCondition high;
};

/**
 * @brief The prescribed first derivatives across the clamped edges of a grid.
 *
 * On a clamped edge of the first axis, `xLow` or `xHigh` holds the derivative
 * in the first variable at that edge, one per node of the second axis; on a
 * clamped edge of the second axis, `yLow` or `yHigh` holds the derivative in
 * the second variable, one per node of the first. Between those nodes the
 * edge derivative is the spline of the given values under the other axis's
 * conditions.
 *
 * The spans are read during construction and not kept, so they need not
 * outlive the spline. An edge that is not clamped ignores its span, which may
 * be left empty.
 *
 * @tparam S Value type, real or complex, matching the grid's values.
 */
template <typename S> struct EdgeSlopes {
    /** d/dx at the first node of the first axis, per second-axis node. */
    std::span<const S> xLow = {};
    /** d/dx at the final node of the first axis, per second-axis node. */
    std::span<const S> xHigh = {};
    /** d/dy at the first node of the second axis, per first-axis node. */
    std::span<const S> yLow = {};
    /** d/dy at the final node of the second axis, per first-axis node. */
    std::span<const S> yHigh = {};
};

/**
 * @brief Tensor-product bicubic spline interpolation on a rectilinear grid.
 *
 * This is the genuine tensor product rather than a spline-of-splines
 * evaluated line by line. Three arrays of curvatures are computed once at
 * construction:
 *
 * - @f$ M^x @f$, the second derivative in the first variable;
 * - @f$ M^y @f$, the second derivative in the second variable;
 * - @f$ M^{xy} @f$, the mixed fourth derivative, obtained by applying the
 *   first-axis solve to @f$ M^y @f$.
 *
 * Evaluation is then two applications of the same one-dimensional spline
 * segment formula the 1D class uses, and costs no allocation.
 *
 * Each axis takes its own pair of end conditions. One condition on all four
 * edges is the short form, NotAKnot by default; the long form takes an
 * AxisConditions per axis, and the slopes for any clamped edge. One axis may
 * be clamped, not both: clamping both would need the cross-derivative at the
 * four corners as well. A clamped first axis makes the end data of the mixed
 * solve the second-axis curvatures of the edge slopes, which is what keeps
 * the result the tensor product of the two conditions: along a clamped edge
 * the derivative across it is the spline of the given slopes, between the
 * nodes as well as at them.
 *
 * A periodic axis follows the one-dimensional convention: its final node
 * closes the period, the values and any edge slopes along it must close
 * exactly, and queries along it are reduced into one period.
 *
 * Construction takes ranges, borrowing lvalues and owning rvalues, with the
 * values flat and row-major: element `(i, j)` at `i * size(y) + j`. Queries
 * outside the grid along an axis that is not periodic continue the nearest
 * edge cell.
 *
 * @tparam XView View over the first-axis abscissae.
 * @tparam YView View over the second-axis abscissae.
 * @tparam VView View over the flat, row-major values.
 */
template <typename XView, typename YView, typename VView>
    requires RealRange<XView> && RealRange<YView> &&
             RealOrComplexRange<VView> && std::ranges::view<XView> &&
             std::ranges::view<YView> && std::ranges::view<VView>
class BicubicSpline {
  public:
    /** @brief Abscissa precision. */
    using Real = std::ranges::range_value_t<XView>;
    /** @brief Value type, real or complex. */
    using Scalar = std::ranges::range_value_t<VView>;

    /**
     * @brief Construct from two axes and a flat, row-major value range, with
     * one condition on all four edges.
     * @param x First-axis abscissae, strictly increasing.
     * @param y Second-axis abscissae, strictly increasing.
     * @param values `size(x) * size(y)` values in row-major order.
     * @param edge Condition applied on all four edges. NotAKnot is the
     *        default because it keeps the scheme fourth order at the
     *        boundary, where Natural costs an order; it needs at least four
     *        nodes on each axis, and Natural two. Periodic makes the spline
     *        periodic in both variables. Clamped is refused here, since it
     *        needs slopes; the per-axis constructor takes them.
     * @throws std::invalid_argument if an axis is too short or not strictly
     *         increasing, the value count does not match the axes, Clamped is
     *         requested, or Periodic values do not close.
     */
    BicubicSpline(XView x, YView y, VView values,
                  BoundaryCondition edge = BoundaryCondition::NotAKnot)
        : _x{std::move(x)}, _y{std::move(y)}, _v{std::move(values)} {
        if (edge == BoundaryCondition::Clamped) {
            throw std::invalid_argument(
                "BicubicSpline: Clamped would need a prescribed derivative "
                "along every edge, which this constructor does not take; use "
                "the constructor taking AxisConditions and EdgeSlopes");
        }
        Solve({edge, edge}, {edge, edge}, {});
    }

    /**
     * @brief Construct with a pair of end conditions on each axis.
     *
     * @code
     * BicubicSpline s{theta, phi, values,
     *                 AxisConditions{BoundaryCondition::Clamped,
     *                                BoundaryCondition::Clamped},
     *                 AxisConditions{BoundaryCondition::Periodic,
     *                                BoundaryCondition::Periodic},
     *                 EdgeSlopes<double>{.xLow = north, .xHigh = south}};
     * @endcode
     *
     * @param x First-axis abscissae, strictly increasing.
     * @param y Second-axis abscissae, strictly increasing.
     * @param values `size(x) * size(y)` values in row-major order.
     * @param alongX Conditions at the first and final nodes of the first
     *        axis.
     * @param alongY Conditions at the first and final nodes of the second
     *        axis.
     * @param slopes Derivatives across each clamped edge; see EdgeSlopes.
     * @throws std::invalid_argument if an axis is too short for its
     *         conditions or not strictly increasing, the value count does not
     *         match the axes, NotAKnot or Periodic is used at one end of an
     *         axis only, both axes are clamped, a clamped edge's slopes are the
     *         wrong length, or values or slopes do not close along a periodic
     *         axis.
     */
    BicubicSpline(XView x, YView y, VView values, AxisConditions alongX,
                  AxisConditions alongY, EdgeSlopes<Scalar> slopes = {})
        : _x{std::move(x)}, _y{std::move(y)}, _v{std::move(values)} {
        Solve(alongX, alongY, slopes);
    }

    /** @brief Number of nodes along the first axis. */
    std::size_t Rows() const { return _rows; }

    /** @brief Number of nodes along the second axis. */
    std::size_t Columns() const { return _columns; }

    /**
     * @brief Evaluate the interpolant or a mixed partial derivative.
     *
     * The pieces are bicubic, so an order above three in either variable is
     * identically zero.
     *
     * @tparam Nx Derivative order in the first variable.
     * @tparam Ny Derivative order in the second variable.
     * @param x First-axis query.
     * @param y Second-axis query.
     */
    template <std::size_t Nx = 0, std::size_t Ny = 0>
    Scalar Evaluate(Real x, Real y) const {
        x = Reduce(x, _x, _periodX);
        y = Reduce(y, _y, _periodY);
        const auto i = Detail::LocateSegment(_x, x);
        const auto j = Detail::LocateSegment(_y, y);

        const auto hx = _x[i + 1] - _x[i];
        const auto tx = x - _x[i];
        const auto hy = _y[j + 1] - _y[j];
        const auto ty = y - _y[j];

        // Collapse the first axis at the two bounding y-nodes: once for the
        // values, once for their second derivatives in y. Those four numbers
        // are exactly what the second-axis piece needs.
        const auto valueAtLow = Detail::SplinePiece<Nx, Real, Scalar>(
            hx, tx, V(i, j), V(i + 1, j), Mx(i, j), Mx(i + 1, j));
        const auto valueAtHigh = Detail::SplinePiece<Nx, Real, Scalar>(
            hx, tx, V(i, j + 1), V(i + 1, j + 1), Mx(i, j + 1),
            Mx(i + 1, j + 1));
        const auto curvatureAtLow = Detail::SplinePiece<Nx, Real, Scalar>(
            hx, tx, My(i, j), My(i + 1, j), Mxy(i, j), Mxy(i + 1, j));
        const auto curvatureAtHigh = Detail::SplinePiece<Nx, Real, Scalar>(
            hx, tx, My(i, j + 1), My(i + 1, j + 1), Mxy(i, j + 1),
            Mxy(i + 1, j + 1));

        return Detail::SplinePiece<Ny, Real, Scalar>(
            hy, ty, valueAtLow, valueAtHigh, curvatureAtLow, curvatureAtHigh);
    }

    /** @brief Evaluate the interpolant; the same as `Evaluate<0, 0>`. */
    Scalar operator()(Real x, Real y) const { return Evaluate<0, 0>(x, y); }

  private:
    XView _x;
    YView _y;
    VView _v;

    std::size_t _rows = 0;
    std::size_t _columns = 0;
    std::vector<Scalar> _mx;  // second derivative in the first variable
    std::vector<Scalar> _my;  // second derivative in the second variable
    std::vector<Scalar> _mxy; // mixed fourth derivative
    std::optional<Real> _periodX;
    std::optional<Real> _periodY;

    Scalar V(std::size_t i, std::size_t j) const {
        return _v[i * _columns + j];
    }
    Scalar Mx(std::size_t i, std::size_t j) const {
        return _mx[i * _columns + j];
    }
    Scalar My(std::size_t i, std::size_t j) const {
        return _my[i * _columns + j];
    }
    Scalar Mxy(std::size_t i, std::size_t j) const {
        return _mxy[i * _columns + j];
    }

    // The query a periodic axis actually answers; any other is unchanged.
    template <typename Axis>
    static Real Reduce(Real q, const Axis &axis,
                       const std::optional<Real> &period) {
        if (period) {
            return Detail::ReduceToPeriod(q, static_cast<Real>(axis[0]),
                                          *period)
                .reduced;
        }
        return q;
    }

    static bool HasClampedEnd(AxisConditions conditions) {
        return conditions.low == BoundaryCondition::Clamped ||
               conditions.high == BoundaryCondition::Clamped;
    }

    static bool IsPeriodic(AxisConditions conditions) {
        return conditions.low == BoundaryCondition::Periodic;
    }

    // Everything that can be said about the arguments before any solve, so
    // that a fault in the caller's data is reported against this class and
    // never surfaces from a solve deep in the assembly.
    void Validate(AxisConditions alongX, AxisConditions alongY,
                  const EdgeSlopes<Scalar> &slopes) const {
        const auto minimumX = Detail::SplineConditionMinimum(
            alongX.low, alongX.high, "BicubicSpline");
        const auto minimumY = Detail::SplineConditionMinimum(
            alongY.low, alongY.high, "BicubicSpline");
        if (HasClampedEnd(alongX) && HasClampedEnd(alongY)) {
            throw std::invalid_argument(
                "BicubicSpline: clamping both axes would need the "
                "cross-derivative at the four corners as well, which this "
                "class does not take; clamp one axis only");
        }
        Detail::ValidateGrid(_x, _y, _v, minimumX, minimumY, "BicubicSpline");

        const auto checkLength = [](std::span<const Scalar> edge,
                                    std::size_t expected, const char *name) {
            if (edge.size() != expected) {
                throw std::invalid_argument(
                    std::string{"BicubicSpline: the clamped edge "} + name +
                    " has " + std::to_string(edge.size()) +
                    " slopes, but there are " + std::to_string(expected) +
                    " nodes along it");
            }
        };
        if (alongX.low == BoundaryCondition::Clamped) {
            checkLength(slopes.xLow, _columns, "xLow");
        }
        if (alongX.high == BoundaryCondition::Clamped) {
            checkLength(slopes.xHigh, _columns, "xHigh");
        }
        if (alongY.low == BoundaryCondition::Clamped) {
            checkLength(slopes.yLow, _rows, "yLow");
        }
        if (alongY.high == BoundaryCondition::Clamped) {
            checkLength(slopes.yHigh, _rows, "yHigh");
        }

        // A periodic axis closes in the values and in the slopes across any
        // clamped edge running along it, exactly as in one dimension. Every
        // array derived from them then closes by construction.
        const auto checkClosure = [](bool closes, const std::string &what) {
            if (!closes) {
                throw std::invalid_argument(
                    "BicubicSpline: along a periodic axis the final node "
                    "closes the period, so " +
                    what + " must equal the first exactly");
            }
        };
        const auto edgeCloses = [](std::span<const Scalar> edge) {
            return edge.empty() || edge.back() == edge.front();
        };
        if (IsPeriodic(alongX)) {
            for (std::size_t j = 0; j < _columns; ++j) {
                checkClosure(V(_rows - 1, j) == V(0, j),
                             "the final row of values");
            }
            checkClosure(edgeCloses(slopes.yLow) && edgeCloses(slopes.yHigh),
                         "the final slope across each clamped y-edge");
        }
        if (IsPeriodic(alongY)) {
            for (std::size_t i = 0; i < _rows; ++i) {
                checkClosure(V(i, _columns - 1) == V(i, 0),
                             "the final column of values");
            }
            checkClosure(edgeCloses(slopes.xLow) && edgeCloses(slopes.xHigh),
                         "the final slope across each clamped x-edge");
        }
    }

    void Solve(AxisConditions alongX, AxisConditions alongY,
               const EdgeSlopes<Scalar> &slopes) {
        _rows = static_cast<std::size_t>(std::ranges::size(_x));
        _columns = static_cast<std::size_t>(std::ranges::size(_y));
        Validate(alongX, alongY, slopes);

        _mx.assign(_rows * _columns, Scalar{});
        _my.assign(_rows * _columns, Scalar{});
        _mxy.assign(_rows * _columns, Scalar{});

        // Copy the axes into contiguous buffers once: a view is not
        // guaranteed to be contiguous, and the systems below borrow theirs.
        std::vector<Real> xNodes(_rows);
        std::vector<Real> yNodes(_columns);
        for (std::size_t i = 0; i < _rows; ++i) {
            xNodes[i] = _x[i];
        }
        for (std::size_t j = 0; j < _columns; ++j) {
            yNodes[j] = _y[j];
        }

        // Every row shares one system and every column shares another, so the
        // two matrices are assembled and factorised once for the whole grid
        // rather than once per line. That is the tensor product paying for
        // itself: the first-axis system is used _columns times for the values
        // and again for the mixed term, and the second-axis system _rows
        // times.
        const auto xSystem = CubicSplineSystem{std::span<const Real>{xNodes},
                                               alongX.low, alongX.high};
        const auto ySystem = CubicSplineSystem{std::span<const Real>{yNodes},
                                               alongY.low, alongY.high};
        _periodX = xSystem.Period();
        _periodY = ySystem.Period();

        // The right-hand side is read while the solution is written, so the
        // two cannot share a buffer. Both are allocated once for the whole
        // grid rather than once per line.
        const auto longest = std::max(_rows, _columns);
        std::vector<Scalar> line(longest);
        std::vector<Scalar> result(longest);

        // The end derivatives a line's solve needs, read from a slope array
        // only on a clamped end; a solve ignores them anywhere else.
        const auto endData = [](BoundaryCondition condition,
                                std::span<const Scalar> edge) {
            return [=](std::size_t k) {
                return condition == BoundaryCondition::Clamped ? edge[k]
                                                               : Scalar{};
            };
        };

        const auto solveAlongX = [&](auto read, auto low, auto high,
                                     auto write) {
            const auto in = std::span<const Scalar>{line.data(), _rows};
            const auto out = std::span<Scalar>{result.data(), _rows};
            for (std::size_t j = 0; j < _columns; ++j) {
                for (std::size_t i = 0; i < _rows; ++i) {
                    line[i] = read(i, j);
                }
                xSystem.Solve(in, low(j), high(j), out);
                for (std::size_t i = 0; i < _rows; ++i) {
                    write(i, j, out[i]);
                }
            }
        };

        const auto solveAlongY = [&](auto read, auto low, auto high,
                                     auto write) {
            const auto in = std::span<const Scalar>{line.data(), _columns};
            const auto out = std::span<Scalar>{result.data(), _columns};
            for (std::size_t i = 0; i < _rows; ++i) {
                for (std::size_t j = 0; j < _columns; ++j) {
                    line[j] = read(i, j);
                }
                ySystem.Solve(in, low(i), high(i), out);
                for (std::size_t j = 0; j < _columns; ++j) {
                    write(i, j, out[j]);
                }
            }
        };

        solveAlongX(
            [&](auto i, auto j) { return V(i, j); },
            endData(alongX.low, slopes.xLow),
            endData(alongX.high, slopes.xHigh),
            [&](auto i, auto j, auto value) { _mx[i * _columns + j] = value; });
        solveAlongY(
            [&](auto i, auto j) { return V(i, j); },
            endData(alongY.low, slopes.yLow),
            endData(alongY.high, slopes.yHigh),
            [&](auto i, auto j, auto value) { _my[i * _columns + j] = value; });

        // The mixed term is the first-axis solve applied to the second-axis
        // curvatures, which is what makes this a true tensor product rather
        // than two independent one-dimensional fits. A clamped first-axis end
        // needs d/dx of those curvatures there, which is the second-axis
        // curvature of that edge's slopes. The second axis is not clamped
        // when the first is, so its solve needs no end data of its own.
        const auto edgeCurvature = [&](BoundaryCondition condition,
                                       std::span<const Scalar> edge) {
            std::vector<Scalar> curvature;
            if (condition == BoundaryCondition::Clamped) {
                curvature.resize(_columns);
                ySystem.Solve(edge, std::span<Scalar>{curvature});
            }
            return curvature;
        };
        const auto lowCurvature = edgeCurvature(alongX.low, slopes.xLow);
        const auto highCurvature = edgeCurvature(alongX.high, slopes.xHigh);

        solveAlongX(
            [&](auto i, auto j) { return My(i, j); },
            endData(alongX.low, std::span<const Scalar>{lowCurvature}),
            endData(alongX.high, std::span<const Scalar>{highCurvature}),
            [&](auto i, auto j, auto value) {
                _mxy[i * _columns + j] = value;
            });
    }
};

/// Borrow lvalue containers and own rvalue ones.
template <std::ranges::viewable_range X, std::ranges::viewable_range Y,
          std::ranges::viewable_range V, typename... Rest>
BicubicSpline(X &&, Y &&, V &&, Rest...)
    -> BicubicSpline<std::views::all_t<X>, std::views::all_t<Y>,
                     std::views::all_t<V>>;

} // namespace Interpolation

#endif // INTERPOLATION_BICUBIC_SPLINE_HPP
