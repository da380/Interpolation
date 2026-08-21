#ifndef INTERPOLATION_TEST_CUBIC_SPLINE_GUARD_H
#define INTERPOLATION_TEST_CUBIC_SPLINE_GUARD_H

#include <Interpolation/Interpolation.hpp>
#include <algorithm>
#include <complex>
#include <cstddef>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

template <Interpolation::Real x_value_t, Interpolation::RealOrComplex y_value_t>
int
CubicSplineCheck() {
    using namespace Interpolation;

    // Make a random cubic polynomial.
    auto p = Polynomial1D<y_value_t>::Random(3);

    // set the number of sampling points randomly
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::uniform_int_distribution dint{5, 100};
    auto nSample = dint(gen);

    // Set the x values
    const x_value_t x1 = 0;
    const x_value_t x2 = 1;
    auto h = (x2 - x1) / static_cast<x_value_t>(nSample - 1);
    std::vector<x_value_t> x;
    std::generate_n(std::back_inserter(x), nSample,
                    [x1, h, m = 0]() mutable { return x1 + h * m++; });

    // add random shifts so that points are not equally spaced
    std::uniform_real_distribution<x_value_t> hDist(-0.1 * h, 0.1 * h);
    std::transform(std::next(x.begin()), std::prev(x.end()),
                   std::next(x.begin()),
                   [&](auto x) { return x + hDist(gen); });

    // Set the y-values
    std::vector<y_value_t> y;
    std::transform(x.begin(), x.end(), std::back_inserter(y),
                   [&](auto x) { return p(x); });

    // Form the interpolating function.
    auto f = CubicSpline(x.begin(), x.end(), y.begin(), CubicSplineBC::Clamped,
                         p.Derivative(x1), p.Derivative(x2));

    // Compare exact and interpolated values at randomly sampled points
    std::uniform_real_distribution<x_value_t> xDist{x1, x2};
    constexpr auto eps = 1000 * std::numeric_limits<x_value_t>::epsilon();
    const int nRandom = 100;
    int count = 0;
    while (count++ < nRandom) {
        auto xx = xDist(gen);
        x_value_t functionError = std::abs(f(xx) - p(xx));
        x_value_t derivativeError =
            std::abs(f.Derivative(xx) - p.Derivative(xx));
        if (functionError > eps)
            return 1;
        if (derivativeError > eps)
            return 1;
    }

    return 0;
}

namespace CubicSplineTest {

// Solve a general dense system by Gaussian elimination with partial pivoting.
//
// This is deliberately not the algorithm under test. CubicSpline exploits the
// fact that its system is tridiagonal and solves it by the Thomas algorithm
// without pivoting; this reference ignores the structure entirely and pivots,
// so agreement between the two is real evidence rather than the same code run
// twice. It previously called Eigen's fullPivLu, and was replaced in kind when
// the Eigen dependency was removed.
template <typename y_value_t>
std::vector<y_value_t>
SolveDense(std::vector<std::vector<y_value_t>> a, std::vector<y_value_t> b) {
    const auto n = a.size();
    for (std::size_t column = 0; column < n; ++column) {
        // Choose the pivot row by largest magnitude in this column.
        auto pivot = column;
        for (std::size_t row = column + 1; row < n; ++row) {
            if (std::abs(a[row][column]) > std::abs(a[pivot][column])) {
                pivot = row;
            }
        }
        if (std::abs(a[pivot][column]) == 0) {
            throw std::runtime_error("singular reference system");
        }
        std::swap(a[column], a[pivot]);
        std::swap(b[column], b[pivot]);

        for (std::size_t row = column + 1; row < n; ++row) {
            const auto factor = a[row][column] / a[column][column];
            if (std::abs(factor) == 0) {
                continue;
            }
            for (std::size_t k = column; k < n; ++k) {
                a[row][k] -= factor * a[column][k];
            }
            b[row] -= factor * b[column];
        }
    }

    // Back substitution.
    std::vector<y_value_t> x(n, y_value_t{});
    for (std::size_t i = n; i-- > 0;) {
        auto sum = b[i];
        for (std::size_t k = i + 1; k < n; ++k) {
            sum -= a[i][k] * x[k];
        }
        x[i] = sum / a[i][i];
    }
    return x;
}

template <typename y_value_t> class DenseReferenceSpline {
  public:
    using BC = Interpolation::CubicSplineBC;
    using Vector = std::vector<y_value_t>;
    using Matrix = std::vector<std::vector<y_value_t>>;

    DenseReferenceSpline(const std::vector<double> &x,
                         const std::vector<y_value_t> &y, BC left,
                         y_value_t ypl, BC right, y_value_t ypr)
        : x_{x}, y_{y} {
        const auto n = static_cast<std::ptrdiff_t>(x_.size());
        Matrix matrix(n, Vector(n, y_value_t{}));
        Vector rhs(n, y_value_t{});

        if (left == BC::Free) {
            matrix[0][0] = 1;
        } else {
            const auto h = x_[1] - x_[0];
            matrix[0][0] = h / 3.0;
            matrix[0][1] = h / 6.0;
            rhs[0] = (y_[1] - y_[0]) / h - ypl;
        }

        for (std::ptrdiff_t i = 1; i < n - 1; ++i) {
            const auto hPrev = x_[i] - x_[i - 1];
            const auto hNext = x_[i + 1] - x_[i];
            matrix[i][i - 1] = hPrev / 6.0;
            matrix[i][i] = (hPrev + hNext) / 3.0;
            matrix[i][i + 1] = hNext / 6.0;
            rhs[i] = (y_[i + 1] - y_[i]) / hNext - (y_[i] - y_[i - 1]) / hPrev;
        }

        if (right == BC::Free) {
            matrix[n - 1][n - 1] = 1;
        } else {
            const auto h = x_[n - 1] - x_[n - 2];
            matrix[n - 1][n - 2] = h / 6.0;
            matrix[n - 1][n - 1] = h / 3.0;
            rhs[n - 1] = ypr - (y_[n - 1] - y_[n - 2]) / h;
        }

        second_ = SolveDense(matrix, rhs);
    }

    y_value_t operator()(double value) const {
        const auto [lower, upper] = interval(value);
        const auto h = x_[upper] - x_[lower];
        const auto a = (x_[upper] - value) / h;
        const auto b = (value - x_[lower]) / h;
        return a * y_[lower] + b * y_[upper] +
               ((a * a * a - a) * second_[static_cast<std::size_t>(lower)] +
                (b * b * b - b) * second_[static_cast<std::size_t>(upper)]) *
                   h * h / 6.0;
    }

    y_value_t Derivative(double value) const {
        const auto [lower, upper] = interval(value);
        const auto h = x_[upper] - x_[lower];
        const auto a = (x_[upper] - value) / h;
        const auto b = (value - x_[lower]) / h;
        return (y_[upper] - y_[lower]) / h +
               h / 6.0 *
                   ((-3.0 * a * a + 1.0) *
                        second_[static_cast<std::size_t>(lower)] +
                    (3.0 * b * b - 1.0) *
                        second_[static_cast<std::size_t>(upper)]);
    }

  private:
    std::pair<std::ptrdiff_t, std::ptrdiff_t> interval(double value) const {
        auto iter = std::upper_bound(x_.begin(), x_.end(), value);
        if (iter == x_.begin())
            ++iter;
        if (iter == x_.end())
            --iter;
        const auto upper =
            static_cast<std::ptrdiff_t>(std::distance(x_.begin(), iter));
        return {upper - 1, upper};
    }

    std::vector<double> x_;
    std::vector<y_value_t> y_;
    Vector second_;
};

template <typename spline_t, typename y_value_t>
std::pair<y_value_t, y_value_t>
RecoverIntervalSecondDerivatives(const spline_t &spline, double xLower,
                                 double xUpper, const y_value_t &yLower,
                                 const y_value_t &yUpper) {
    const auto h = xUpper - xLower;
    const auto slope = (yUpper - yLower) / h;
    const auto a = 6.0 * (slope - spline.Derivative(xLower)) / h;
    const auto b = 6.0 * (spline.Derivative(xUpper) - slope) / h;
    return {(2.0 * a - b) / 3.0, (2.0 * b - a) / 3.0};
}

} // namespace CubicSplineTest

#endif // INTERPOLATION_TEST_CUBIC_SPLINE_GUARD_H
