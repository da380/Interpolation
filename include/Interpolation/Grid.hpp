#ifndef INTERPOLATION_GRID_HPP
#define INTERPOLATION_GRID_HPP

#include <array>
#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>

#include <Interpolation/Concepts.hpp>

namespace Interpolation {

/**
 * @brief Row-major extents of a rectilinear grid.
 *
 * This is the N-D seam. `std::mdspan` is the right long-term home for it, but
 * it needs GCC 15 or libc++ 18, neither of which is available here, so the
 * library carries this instead. It is deliberately small and deliberately
 * only does what the interpolators ask of it, so that swapping in `mdspan`
 * later is an internal change rather than a public one.
 *
 * @tparam Rank Number of dimensions.
 */
template <std::size_t Rank>
    requires(Rank > 0)
class Extents {
  public:
    /** @brief Construct from one extent per dimension. */
    constexpr explicit Extents(std::array<std::size_t, Rank> extents)
        : _extents{extents} {}

    /** @brief Number of dimensions. */
    static constexpr std::size_t Order() { return Rank; }

    /** @brief Extent along dimension `d`. */
    constexpr std::size_t Extent(std::size_t d) const { return _extents[d]; }

    /** @brief Total number of elements. */
    constexpr std::size_t Size() const {
        std::size_t total = 1;
        for (const auto extent : _extents) {
            total *= extent;
        }
        return total;
    }

    /** @brief Row-major offset of an index tuple. */
    constexpr std::size_t Index(std::array<std::size_t, Rank> indices) const {
        std::size_t offset = 0;
        for (std::size_t d = 0; d < Rank; ++d) {
            offset = offset * _extents[d] + indices[d];
        }
        return offset;
    }

  private:
    std::array<std::size_t, Rank> _extents;
};

/** @brief Extents of a two-dimensional grid. */
using Extents2D = Extents<2>;

/**
 * @brief A non-owning two-dimensional view over a flat range, row-major.
 *
 * Element `(i, j)` sits at offset `i * columns + j`, so the first index runs
 * over the abscissae of the first axis.
 *
 * @tparam VView View over the stored values.
 */
template <typename VView>
    requires RealOrComplexRange<VView> && std::ranges::view<VView>
class Grid2DView {
  public:
    /** @brief Value type, real or complex. */
    using Scalar = std::ranges::range_value_t<VView>;

    /** @brief Construct from a flat range and its extents. */
    constexpr Grid2DView(VView values, Extents2D extents)
        : _values{std::move(values)}, _extents{extents} {}

    /** @brief Number of rows, the extent along the first axis. */
    constexpr std::size_t Rows() const { return _extents.Extent(0); }

    /** @brief Number of columns, the extent along the second axis. */
    constexpr std::size_t Columns() const { return _extents.Extent(1); }

    /** @brief Element at row `i` and column `j`. */
    constexpr decltype(auto) operator()(std::size_t i, std::size_t j) const {
        return _values[_extents.Index({i, j})];
    }

  private:
    VView _values;
    Extents2D _extents;
};

namespace Detail {

/**
 * @brief Validate the axes and values of a rectilinear grid at construction.
 *
 * Same policy as the one-dimensional case: data errors are reported by
 * throwing, since construction is not the hot path, and both axes must be
 * strictly increasing.
 *
 * @param x First-axis abscissae.
 * @param y Second-axis abscissae.
 * @param values Flat, row-major values with `size(x) * size(y)` elements.
 * @param minimumNodes Fewest nodes each axis needs.
 * @param what Interpolator name, used in the exception message.
 * @throws std::invalid_argument if an axis is too short or not strictly
 *         increasing, or the value count does not match the axes.
 */
template <typename XRange, typename YRange, typename VRange>
void
ValidateGrid(const XRange &x, const YRange &y, const VRange &values,
             std::size_t minimumNodes, const char *what) {
    ValidateSamples(x, x, minimumNodes, what);
    ValidateSamples(y, y, minimumNodes, what);

    const auto nx = static_cast<std::size_t>(std::ranges::size(x));
    const auto ny = static_cast<std::size_t>(std::ranges::size(y));
    const auto nv = static_cast<std::size_t>(std::ranges::size(values));

    if (nv != nx * ny) {
        throw std::invalid_argument(
            std::string{what} + ": expected " + std::to_string(nx) + " x " +
            std::to_string(ny) + " = " + std::to_string(nx * ny) +
            " values, but was given " + std::to_string(nv));
    }
}

} // namespace Detail

} // namespace Interpolation

#endif // INTERPOLATION_GRID_HPP
