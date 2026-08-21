#ifndef INTERPOLATION_LAGRANGE_HPP
#define INTERPOLATION_LAGRANGE_HPP

#include <algorithm>
#include <cassert>
#include <concepts>
#include <iterator>
#include <numeric>
#include <ranges>
#include <vector>

#include <Interpolation/Concepts.hpp>

namespace Interpolation {

/**
 * @brief Lagrange cardinal basis on a fixed set of nodes.
 *
 * The object stores an iterator to the nodes rather than copying them. The
 * underlying container must remain alive and must not be reallocated while the
 * basis is in use.
 *
 * @tparam I Random-access iterator over real floating-point nodes.
 * @pre The range contains at least one strictly increasing node.
 */
template <RealIterator I> class LagrangePolynomial {
  public:
    /** @brief Scalar type of the interpolation nodes. */
    using value_t = std::iter_value_t<I>;

    /** @brief Default construction is not available without a node range. */
    LagrangePolynomial() = delete;

    /**
     * @brief Construct the cardinal basis for `[start, finish)`.
     * @param start Iterator to the first node.
     * @param finish Iterator one past the final node.
     */
    LagrangePolynomial(I start, I finish)
        : n{std::distance(start, finish)}, X{start} {
        assert(n > 0);
        assert(std::is_sorted(start, finish));
    }

    /**
     * @brief Evaluate one cardinal basis polynomial.
     * @param i Zero-based basis-function index.
     * @param x Query abscissa.
     * @return Value of the `i`th basis polynomial at `x`.
     * @pre `0 <= i <` the number of nodes.
     */
    value_t operator()(int i, value_t x) const {
        auto prod1 = static_cast<value_t>(1);
        auto prod2 = static_cast<value_t>(1);
        for (int j = 0; j < n; j++) {
            if (i != j) {
                prod1 *= x - X[j];
                prod2 *= X[i] - X[j];
            }
        }
        return prod1 / prod2;
    }

    /**
     * @brief Evaluate the derivative of one cardinal basis polynomial.
     * @param i Zero-based basis-function index.
     * @param x Query abscissa.
     * @return First derivative of the `i`th basis polynomial at `x`.
     * @pre `0 <= i <` the number of nodes.
     */
    value_t Derivative(int i, value_t x) const {
        auto hp = static_cast<value_t>(0);
        auto prod2 = static_cast<value_t>(1);
        for (int j = 0; j < n; j++) {
            if (i != j) {
                auto prod1 = static_cast<value_t>(1);
                for (int k = 0; k < n; k++) {
                    if (k != i && k != j)
                        prod1 *= x - X[k];
                }
                hp += prod1;
                prod2 *= X[i] - X[j];
            }
        }
        return hp / prod2;
    }

  private:
    std::ptrdiff_t n; // Number of nodes.
    I X;              // Iterator to the start of the nodes.
};

/**
 * @brief Global polynomial interpolation in the Lagrange cardinal basis.
 *
 * The object stores iterators rather than copying its samples. Both underlying
 * containers must remain alive and must not be reallocated while the
 * interpolator is in use. Evaluation outside the node interval is polynomial
 * extrapolation.
 *
 * @tparam xIter Random-access iterator over real abscissae.
 * @tparam yIter Random-access iterator over real or complex ordinates.
 * @pre The abscissa range contains at least one strictly increasing value.
 */
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
class Lagrange {
  public:
    /** @brief Scalar type used for abscissae. */
    using x_value_t = std::iter_value_t<xIter>;
    /** @brief Scalar type used for interpolated values. */
    using y_value_t = std::iter_value_t<yIter>;

    /**
     * @brief Construct an interpolator over non-owning sample ranges.
     * @param xS Iterator to the first abscissa.
     * @param xF Iterator one past the final abscissa.
     * @param yS Iterator to the ordinate corresponding to `xS`.
     */
    Lagrange(xIter xS, xIter xF, yIter yS)
        : xS{xS}, xF{xF}, yS{yS}, h{LagrangePolynomial(xS, xF)} {}

    /** @brief Return the number of interpolation nodes. */
    auto size() const { return std::distance(xS, xF); }

    /**
     * @brief Evaluate the interpolating polynomial.
     * @param x Query abscissa.
     * @return Interpolated or extrapolated ordinate.
     */
    y_value_t operator()(x_value_t x) const {
        auto y = static_cast<y_value_t>(0);
        for (int i = 0; i < size(); i++) {
            y += h(i, x) * yS[i];
        }
        return y;
    }

    /**
     * @brief Evaluate the first derivative of the interpolating polynomial.
     * @param x Query abscissa.
     * @return First derivative at `x`.
     */
    y_value_t Derivative(x_value_t x) const {
        auto yp = static_cast<y_value_t>(0);
        for (int i = 0; i < size(); i++) {
            yp += h.Derivative(i, x) * yS[i];
        }
        return yp;
    }

  private:
    xIter xS; // Iterator to the start of the x-values.
    xIter xF; // Iterator to the end of hte x-values.
    yIter yS; // Iterator to the start of the y-values.

    LagrangePolynomial<xIter> h; // Lagrange Polynomial for interpolation.
};

} // namespace Interpolation

#endif
