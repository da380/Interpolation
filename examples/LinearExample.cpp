
#include <Interpolation/Interpolation.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

int
main() {
    using namespace Interpolation;

    int n = 2;
    std::vector<double> x(n);
    std::vector<std::complex<double>> y(n);

    auto func = [](double x) { return std::sin(x); };

    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i) / static_cast<double>(n - 1);
        y[i] = func(x[i]);
    }

    Linear f(x, y);

    int m = 50;
    std::ofstream file("Linear.out");
    for (int i = 0; i < m; i++) {
        auto x = static_cast<double>(i) / static_cast<double>(m - 1);
        file << x << " " << func(x) << " " << f(x) << std::endl;
    }
}
