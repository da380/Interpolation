#ifndef INTERPOLATION_BILINEAR_HPP
#define INTERPOLATION_BILINEAR_HPP

#include <cstddef>
#include <ranges>
#include <utility>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Grid.hpp>
#include <Interpolation/Samples.hpp>

namespace Interpolation {

/**
 * @brief Bilinear interpolation on a rectilinear grid.
 *
 * The interpolant is the tensor product of the one-dimensional linear
 * interpolant along each axis, so it is built from the same segment formula
 * rather than a restatement of it. Being a tensor product, a mixed partial
 * derivative is just the two one-dimensional rules applied in turn.
 *
 * Construction takes ranges, borrowing lvalues and owning rvalues, and the
 * values are a flat row-major range: element `(i, j)` at `i * size(y) + j`.
 * Queries outside the grid continue the nearest edge cell.
 *
 * @tparam XView View over the first-axis abscissae.
 * @tparam YView View over the second-axis abscissae.
 * @tparam VView View over the flat, row-major values.
 */
template <typename XView, typename YView, typename VView>
    requires RealRange<XView> && RealRange<YView> &&
             RealOrComplexRange<VView> && std::ranges::view<XView> &&
             std::ranges::view<YView> && std::ranges::view<VView>
class Bilinear {
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
    Bilinear(XView x, YView y, VView values)
        : _x{std::move(x)}, _y{std::move(y)}, _v{std::move(values)} {
        Detail::ValidateGrid(_x, _y, _v, 2, "Bilinear");
    }

    /** @brief Number of nodes along the first axis. */
    std::size_t Rows() const {
        return static_cast<std::size_t>(std::ranges::size(_x));
    }

    /** @brief Number of nodes along the second axis. */
    std::size_t Columns() const {
        return static_cast<std::size_t>(std::ranges::size(_y));
    }

    /**
     * @brief Evaluate the interpolant or a mixed partial derivative.
     *
     * The interpolant is bilinear, so any order above one in either variable
     * is identically zero.
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

        const auto hy = _y[j + 1] - _y[j];
        const auto ty = y - _y[j];

        // Collapse the second axis first, then the first. Separability makes
        // this exact for the mixed derivative too.
        const auto lower = Detail::LinearPiece<Ny, Real, Scalar>(
            hy, ty, At(i, j), At(i, j + 1));
        const auto upper = Detail::LinearPiece<Ny, Real, Scalar>(
            hy, ty, At(i + 1, j), At(i + 1, j + 1));

        return Detail::LinearPiece<Nx, Real, Scalar>(_x[i + 1] - _x[i],
                                                     x - _x[i], lower, upper);
    }

    /** @brief Evaluate the interpolant; the same as `Evaluate<0, 0>`. */
    Scalar operator()(Real x, Real y) const { return Evaluate<0, 0>(x, y); }

  private:
    XView _x;
    YView _y;
    VView _v;

    Scalar At(std::size_t i, std::size_t j) const {
        return _v[i * Columns() + j];
    }
};

/// Borrow lvalue containers and own rvalue ones.
template <std::ranges::viewable_range X, std::ranges::viewable_range Y,
          std::ranges::viewable_range V>
Bilinear(X &&, Y &&,
         V &&) -> Bilinear<std::views::all_t<X>, std::views::all_t<Y>,
                           std::views::all_t<V>>;

} // namespace Interpolation

#endif // INTERPOLATION_BILINEAR_HPP
