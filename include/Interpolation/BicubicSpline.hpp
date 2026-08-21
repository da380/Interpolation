#ifndef INTERPOLATION_BICUBIC_SPLINE_HPP
#define INTERPOLATION_BICUBIC_SPLINE_HPP

#include <cstddef>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Grid.hpp>
#include <Interpolation/Samples.hpp>

namespace Interpolation {

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
 * segment formula the 1D class uses, and costs no allocation. The natural
 * condition is applied on all four edges.
 *
 * Construction takes ranges, borrowing lvalues and owning rvalues, with the
 * values flat and row-major: element `(i, j)` at `i * size(y) + j`. Queries
 * outside the grid continue the nearest edge cell.
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
     * @brief Construct from two axes and a flat, row-major value range.
     * @param x First-axis abscissae, strictly increasing, at least two.
     * @param y Second-axis abscissae, strictly increasing, at least two.
     * @param values `size(x) * size(y)` values in row-major order.
     * @throws std::invalid_argument if an axis is too short or not strictly
     *         increasing, or the value count does not match the axes.
     */
    BicubicSpline(XView x, YView y, VView values)
        : _x{std::move(x)}, _y{std::move(y)}, _v{std::move(values)} {
        Detail::ValidateGrid(_x, _y, _v, 2, "BicubicSpline");
        Solve();
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

    void Solve() {
        _rows = static_cast<std::size_t>(std::ranges::size(_x));
        _columns = static_cast<std::size_t>(std::ranges::size(_y));

        _mx.assign(_rows * _columns, Scalar{});
        _my.assign(_rows * _columns, Scalar{});
        _mxy.assign(_rows * _columns, Scalar{});

        // Copy the axes into contiguous buffers once: a view is not
        // guaranteed to be contiguous, and the solver takes spans.
        std::vector<Real> xNodes(_rows);
        std::vector<Real> yNodes(_columns);
        for (std::size_t i = 0; i < _rows; ++i) {
            xNodes[i] = _x[i];
        }
        for (std::size_t j = 0; j < _columns; ++j) {
            yNodes[j] = _y[j];
        }

        const auto longest = std::max(_rows, _columns);
        std::vector<Scalar> line(longest);
        std::vector<Scalar> result(longest);
        std::vector<Real> sub(longest), diag(longest), super(longest);

        const auto solveAlongX = [&](auto read, auto write) {
            for (std::size_t j = 0; j < _columns; ++j) {
                for (std::size_t i = 0; i < _rows; ++i) {
                    line[i] = read(i, j);
                }
                Detail::NaturalCurvatures<Real, Scalar>(
                    std::span<const Real>{xNodes},
                    std::span<const Scalar>{line.data(), _rows},
                    std::span<Scalar>{result.data(), _rows},
                    std::span<Real>{sub.data(), _rows},
                    std::span<Real>{diag.data(), _rows},
                    std::span<Real>{super.data(), _rows});
                for (std::size_t i = 0; i < _rows; ++i) {
                    write(i, j, result[i]);
                }
            }
        };

        const auto solveAlongY = [&](auto read, auto write) {
            for (std::size_t i = 0; i < _rows; ++i) {
                for (std::size_t j = 0; j < _columns; ++j) {
                    line[j] = read(i, j);
                }
                Detail::NaturalCurvatures<Real, Scalar>(
                    std::span<const Real>{yNodes},
                    std::span<const Scalar>{line.data(), _columns},
                    std::span<Scalar>{result.data(), _columns},
                    std::span<Real>{sub.data(), _columns},
                    std::span<Real>{diag.data(), _columns},
                    std::span<Real>{super.data(), _columns});
                for (std::size_t j = 0; j < _columns; ++j) {
                    write(i, j, result[j]);
                }
            }
        };

        solveAlongX(
            [&](auto i, auto j) { return V(i, j); },
            [&](auto i, auto j, auto value) { _mx[i * _columns + j] = value; });
        solveAlongY(
            [&](auto i, auto j) { return V(i, j); },
            [&](auto i, auto j, auto value) { _my[i * _columns + j] = value; });
        // The mixed term is the first-axis solve applied to the second-axis
        // curvatures, which is what makes this a true tensor product rather
        // than two independent one-dimensional fits.
        solveAlongX([&](auto i, auto j) { return My(i, j); },
                    [&](auto i, auto j, auto value) {
                        _mxy[i * _columns + j] = value;
                    });
    }
};

/// Borrow lvalue containers and own rvalue ones.
template <std::ranges::viewable_range X, std::ranges::viewable_range Y,
          std::ranges::viewable_range V>
BicubicSpline(X &&, Y &&,
              V &&) -> BicubicSpline<std::views::all_t<X>, std::views::all_t<Y>,
                                     std::views::all_t<V>>;

} // namespace Interpolation

#endif // INTERPOLATION_BICUBIC_SPLINE_HPP
