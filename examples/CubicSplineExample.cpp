
#include <Interpolation/Interpolation.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numbers>
#include <numeric>
#include <random>
#include <vector>

int
main() {

    using namespace Interpolation;

    int n = 2;
    std::vector<double> x(n);
    std::vector<double> y(n);

    auto func = [](double x) { return std::sin(x); };

    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i) / static_cast<double>(n - 1);
        y[i] = func(x[i]);
    }

    CubicSpline f(x.begin(), x.end(), y.begin());

    int m = 50;
    std::ofstream file("Linear.out");
    for (int i = 0; i < m; i++) {
        auto x = static_cast<double>(i) / static_cast<double>(m - 1);
        file << x << " " << func(x) << " " << f(x) << std::endl;
    }

    /*
    // Set the type for x.
    using x_value_t = double;

    // Set the type for y (uncomment the desired option).
    using y_value_t = x_value_t;
    //  using y_value_t = std::complex<x_value_t>;

    // Make a random polynomial.
    auto p = Polynomial1D<y_value_t>::Random(3);

    // Set the function arrays.
    std::vector<x_value_t> x;
    std::vector<y_value_t> y;
    std::vector<y_value_t> yp;
    const int n = 2;
    const x_value_t x1 = 0;
    const x_value_t x2 = std::numbers::pi_v<x_value_t>;
    auto dx = (x2 - x1) / static_cast<x_value_t>(n - 1);
    std::generate_n(std::back_inserter(x), n,
                    [x1, dx, m = 0]() mutable { return x1 + dx * m++; });
    std::transform(x.begin(), x.end(), std::back_inserter(y),
                   [&](auto x) { return p(x); });
    std::transform(x.begin(), x.end(), std::back_inserter(yp),
                   [&](auto x) { return p.Derivative(x); });

    // Form the cubic spline.
    auto f =
        CubicSpline(x.begin(), x.end(), y.begin(), CubicSplineBC::Clamped,
                    p.Derivative(x1), CubicSplineBC::Clamped, p.Derivative(x2));

    // Compare exact and interpolated values at randomly sampled points
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::uniform_real_distribution<x_value_t> d{x1, x2};
    int count = 0;
    auto maxErr = static_cast<x_value_t>(0);
    auto maxDerivErr = static_cast<x_value_t>(0);
    while (count++ < 100) {
        auto xx = d(gen);
        auto err = std::abs(f(xx) - p(xx));
        if (err > maxErr)
            maxErr = err;
        auto derivErr = std::abs(f.Derivative(xx) - p.Derivative(xx));
        if (derivErr > maxDerivErr)
            maxDerivErr = derivErr;
    }
    std::cout << "Maximum interpolation error for function values = " << maxErr
              << std::endl;

    std::cout << "Maximum interpolation error for derivative values = "
              << maxDerivErr << std::endl;

    */
}
