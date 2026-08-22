#ifndef INTERPOLATION_LAGRANGE_HPP
#define INTERPOLATION_LAGRANGE_HPP

#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Samples.hpp>

namespace Interpolation {

/**
 * @brief Lagrange cardinal basis on a fixed set of nodes.
 *
 * The basis is held in barycentric form. The weights
 * @f$ w_j = 1 / \prod_{k \neq j} (x_j - x_k) @f$ are computed once at
 * construction, in @f$ O(n^2) @f$, after which a basis value or its derivative
 * costs @f$ O(n) @f$. Evaluating directly from the product definition instead
 * costs @f$ O(n) @f$ per basis function and @f$ O(n^2) @f$ per derivative,
 * recomputing the same denominators on every call.
 *
 * Construction takes a range, borrowing an lvalue and owning an rvalue.
 *
 * @tparam XView View over real nodes.
 */
template <typename XView>
    requires RealRange<XView> && std::ranges::view<XView>
class LagrangeBasis {
  public:
    /** @brief Node precision. */
    using Real = std::ranges::range_value_t<XView>;

    /**
     * @brief Construct the cardinal basis on `x`.
     * @param x Strictly increasing nodes, at least one.
     * @throws std::invalid_argument if there are no nodes or they are not
     *         strictly increasing.
     */
    explicit LagrangeBasis(XView x) : _x{std::move(x)} {
        Detail::ValidateSamples(_x, _x, 1, "LagrangeBasis");
        ComputeWeights();
    }

    /** @brief Number of nodes. */
    std::size_t Size() const {
        return static_cast<std::size_t>(std::ranges::size(_x));
    }

    /** @brief Barycentric weight of node `i`. */
    Real Weight(std::size_t i) const { return _w[i]; }

    /** @brief Node `i`. */
    Real Node(std::size_t i) const { return _x[i]; }

    /**
     * @brief Evaluate the `i`th cardinal basis polynomial or its derivative.
     * @tparam N Derivative order, `0` or `1`.
     * @param i Zero-based basis index.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> Real Evaluate(std::size_t i, Real x) const {
        static_assert(N <= 1,
                      "LagrangeBasis evaluates only the value and the first "
                      "derivative of a cardinal basis function.");
        if constexpr (N == 1) {
            return Slope(i, x);
        } else {
            return Value(i, x);
        }
    }

    /** @brief Evaluate the `i`th basis function; the same as `Evaluate<0>`. */
    Real operator()(std::size_t i, Real x) const { return Value(i, x); }

  private:
    Real Value(std::size_t i, Real x) const {
        // Exactly on a node the barycentric quotient is 0/0, and the answer
        // is the Kronecker delta by definition of a cardinal basis.
        if (const auto node = NodeAt(x); node != NotANode()) {
            return node == i ? static_cast<Real>(1) : static_cast<Real>(0);
        }

        auto denominator = static_cast<Real>(0);
        for (std::size_t j = 0; j < Size(); ++j) {
            denominator += _w[j] / (x - _x[j]);
        }
        return (_w[i] / (x - _x[i])) / denominator;
    }

    Real Slope(std::size_t i, Real x) const {
        const auto n = Size();

        if (const auto k = NodeAt(x); k != NotANode()) {
            if (k != i) {
                return (_w[i] / _w[k]) / (_x[k] - _x[i]);
            }
            // Differentiating log L_i gives L_i'(x)/L_i(x) = sum 1/(x - x_j),
            // and L_i(x_i) = 1, so the diagonal is the sum itself.
            auto sum = static_cast<Real>(0);
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i) {
                    sum += static_cast<Real>(1) / (_x[i] - _x[j]);
                }
            }
            return sum;
        }

        auto denominator = static_cast<Real>(0);
        for (std::size_t j = 0; j < n; ++j) {
            denominator += _w[j] / (x - _x[j]);
        }
        const auto value = (_w[i] / (x - _x[i])) / denominator;

        // L_i'(x) = sum_j w_j (L_i(x) - delta_ij) / (x - x_j)^2 / D(x)
        auto sum = static_cast<Real>(0);
        for (std::size_t j = 0; j < n; ++j) {
            const auto gap = x - _x[j];
            const auto target = (j == i) ? value - static_cast<Real>(1) : value;
            sum += _w[j] * target / (gap * gap);
        }
        return sum / denominator;
    }

    XView _x;
    std::vector<Real> _w;

    static constexpr std::size_t NotANode() {
        return static_cast<std::size_t>(-1);
    }

    std::size_t NodeAt(Real x) const {
        for (std::size_t j = 0; j < Size(); ++j) {
            if (x == _x[j]) {
                return j;
            }
        }
        return NotANode();
    }

    void ComputeWeights() {
        const auto n = Size();
        _w.assign(n, static_cast<Real>(1));
        for (std::size_t j = 0; j < n; ++j) {
            auto product = static_cast<Real>(1);
            for (std::size_t k = 0; k < n; ++k) {
                if (k != j) {
                    product *= _x[j] - _x[k];
                }
            }
            _w[j] = static_cast<Real>(1) / product;
        }
    }
};

/// Borrow an lvalue container and own an rvalue one.
template <std::ranges::viewable_range X>
LagrangeBasis(X &&) -> LagrangeBasis<std::views::all_t<X>>;

/**
 * @brief Global polynomial interpolation through sampled values.
 *
 * Evaluation uses the second barycentric formula, which costs @f$ O(n) @f$
 * per query after an @f$ O(n^2) @f$ setup.
 *
 * Construction takes ranges, borrowing lvalues and owning rvalues. Evaluation
 * outside the node interval is polynomial extrapolation, which for a
 * high-degree interpolant diverges quickly.
 *
 * @tparam XView View over real abscissae.
 * @tparam YView View over real or complex ordinates.
 */
template <typename XView, typename YView>
    requires InterpolationRanges<XView, YView> && std::ranges::view<XView> &&
             std::ranges::view<YView>
class Lagrange {
  public:
    /** @brief Abscissa precision. */
    using Real = std::ranges::range_value_t<XView>;
    /** @brief Ordinate type, real or complex. */
    using Scalar = std::ranges::range_value_t<YView>;

    /**
     * @brief Construct from abscissa and ordinate ranges.
     * @param x Strictly increasing abscissae, at least one.
     * @param y Ordinates, the same length as `x`.
     * @throws std::invalid_argument if the ranges differ in length, are
     *         empty, or the abscissae are not strictly increasing.
     */
    Lagrange(XView x, YView y) : _y{std::move(y)}, _basis{std::move(x)} {
        // LagrangeBasis has already rejected empty and non-increasing nodes,
        // so only the pairing of the two ranges is left to check.
        const auto ny = static_cast<std::size_t>(std::ranges::size(_y));
        if (_basis.Size() != ny) {
            throw std::invalid_argument(
                "Lagrange: abscissa and ordinate ranges differ in length (" +
                std::to_string(_basis.Size()) + " and " + std::to_string(ny) +
                ")");
        }
    }

    /** @brief Number of interpolation nodes. */
    std::size_t Size() const { return _basis.Size(); }

    /**
     * @brief Evaluate the interpolating polynomial or its first derivative.
     *
     * @tparam N Derivative order, `0` or `1`.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> Scalar Evaluate(Real x) const {
        static_assert(N <= 1,
                      "Lagrange currently evaluates only the value and the "
                      "first derivative; higher barycentric derivatives need "
                      "a recurrence that is not implemented yet.");

        // Both branches accumulate over the nodes once. Summing the cardinal
        // basis functions instead would recompute the shared barycentric
        // denominator for every one of them, which is quadratic and throws
        // away the whole point of the barycentric form.
        const auto n = Size();
        const auto k = NodeIndex(x);

        if constexpr (N == 1) {
            if (k != NoNode()) {
                // p'(x_k) = sum_{j!=k} (w_j/w_k)(y_j - y_k)/(x_k - x_j)
                const auto wk = _basis.Weight(k);
                auto slope = Scalar{};
                for (std::size_t j = 0; j < n; ++j) {
                    if (j != k) {
                        slope += (_basis.Weight(j) / wk) * (_y[j] - _y[k]) /
                                 (_basis.Node(k) - _basis.Node(j));
                    }
                }
                return slope;
            }

            auto denominator = Real{};
            auto numerator = Scalar{};
            for (std::size_t j = 0; j < n; ++j) {
                const auto term = _basis.Weight(j) / (x - _basis.Node(j));
                denominator += term;
                numerator += term * _y[j];
            }
            const auto value = numerator / denominator;

            // p'(x) = sum_j w_j (p(x) - y_j) / (x - x_j)^2 / D(x)
            auto slope = Scalar{};
            for (std::size_t j = 0; j < n; ++j) {
                const auto gap = x - _basis.Node(j);
                slope += _basis.Weight(j) * (value - _y[j]) / (gap * gap);
            }
            return slope / denominator;
        } else {
            if (k != NoNode()) {
                return _y[k];
            }

            auto denominator = Real{};
            auto numerator = Scalar{};
            for (std::size_t j = 0; j < n; ++j) {
                const auto term = _basis.Weight(j) / (x - _basis.Node(j));
                denominator += term;
                numerator += term * _y[j];
            }
            return numerator / denominator;
        }
    }

    /** @brief Evaluate the interpolant; the same as `Evaluate<0>`. */
    Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    static constexpr std::size_t NoNode() {
        return static_cast<std::size_t>(-1);
    }

    std::size_t NodeIndex(Real x) const {
        for (std::size_t j = 0; j < Size(); ++j) {
            if (x == _basis.Node(j)) {
                return j;
            }
        }
        return NoNode();
    }

    YView _y;
    LagrangeBasis<XView> _basis;
};

/// Borrow lvalue containers and own rvalue ones.
template <std::ranges::viewable_range X, std::ranges::viewable_range Y>
Lagrange(X &&, Y &&) -> Lagrange<std::views::all_t<X>, std::views::all_t<Y>>;

} // namespace Interpolation

#endif // INTERPOLATION_LAGRANGE_HPP
