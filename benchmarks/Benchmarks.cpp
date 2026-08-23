// Measures the two algorithmic changes made in phase 2 against the
// implementations they replaced, so the speedups are numbers rather than
// complexity arguments.
//
// The reference implementations below are deliberate reconstructions of the
// old code: Lagrange evaluated straight from the product definition, and the
// spline system solved densely with pivoting. They exist only here.
//
// This is a hand-run harness, not a statistical one. It reports the best of
// several repetitions, which is the least noisy simple summary on a busy
// machine, and it prints the work done so the optimiser cannot discard it.

#include <Interpolation/Interpolation.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <numeric>
#include <random>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

// Keeps a computed value observable so the optimiser cannot delete the work
// being timed.
volatile double g_sink = 0.0;

template <typename F>
double
BestOf(int repetitions, F &&f) {
    double best = 1.0e300;
    for (int r = 0; r < repetitions; ++r) {
        const auto start = Clock::now();
        f();
        const auto elapsed =
            std::chrono::duration<double, std::milli>(Clock::now() - start)
                .count();
        best = std::min(best, elapsed);
    }
    return best;
}

// The Lagrange evaluation this replaced: O(n) per basis function and O(n^2)
// per basis derivative, recomputing the denominators on every call.
double
NaiveLagrange(const std::vector<double> &x, const std::vector<double> &y,
              double query) {
    const auto n = x.size();
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double numerator = 1.0;
        double denominator = 1.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i) {
                numerator *= query - x[j];
                denominator *= x[i] - x[j];
            }
        }
        total += (numerator / denominator) * y[i];
    }
    return total;
}

double
NaiveLagrangeDerivative(const std::vector<double> &x,
                        const std::vector<double> &y, double query) {
    const auto n = x.size();
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double slope = 0.0;
        double denominator = 1.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i) {
                double product = 1.0;
                for (std::size_t k = 0; k < n; ++k) {
                    if (k != i && k != j) {
                        product *= query - x[k];
                    }
                }
                slope += product;
                denominator *= x[i] - x[j];
            }
        }
        total += (slope / denominator) * y[i];
    }
    return total;
}

// A dense solve with partial pivoting, standing in for the sparse
// factorisation the spline used before: O(n^3) against the Thomas O(n).
std::vector<double>
DenseSolveSplineSystem(const std::vector<double> &x,
                       const std::vector<double> &y) {
    const auto n = x.size();
    // Flat storage rather than a vector of vectors: contiguous, one
    // allocation, and row swaps move indices instead of whole containers.
    std::vector<double> a(n * n, 0.0);
    std::vector<double> b(n, 0.0);
    const auto at = [n](std::size_t r, std::size_t c) { return r * n + c; };

    a[at(0, 0)] = 1.0;
    a[at(n - 1, n - 1)] = 1.0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        const auto hPrev = x[i] - x[i - 1];
        const auto hNext = x[i + 1] - x[i];
        a[at(i, i - 1)] = hPrev / 6.0;
        a[at(i, i)] = (hPrev + hNext) / 3.0;
        a[at(i, i + 1)] = hNext / 6.0;
        b[i] = (y[i + 1] - y[i]) / hNext - (y[i] - y[i - 1]) / hPrev;
    }

    for (std::size_t col = 0; col < n; ++col) {
        auto pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            if (std::abs(a[at(row, col)]) > std::abs(a[at(pivot, col)])) {
                pivot = row;
            }
        }
        if (pivot != col) {
            for (std::size_t k = 0; k < n; ++k) {
                std::swap(a[at(col, k)], a[at(pivot, k)]);
            }
            std::swap(b[col], b[pivot]);
        }
        for (std::size_t row = col + 1; row < n; ++row) {
            const auto factor = a[at(row, col)] / a[at(col, col)];
            if (factor == 0.0) {
                continue;
            }
            for (std::size_t k = col; k < n; ++k) {
                a[at(row, k)] -= factor * a[at(col, k)];
            }
            b[row] -= factor * b[col];
        }
    }

    std::vector<double> solution(n, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        auto sum = b[i];
        for (std::size_t k = i + 1; k < n; ++k) {
            sum -= a[at(i, k)] * solution[k];
        }
        solution[i] = sum / a[at(i, i)];
    }
    return solution;
}

std::pair<std::vector<double>, std::vector<double>>
MakeSamples(std::size_t n, std::uint64_t seed) {
    std::mt19937_64 gen{seed};
    std::uniform_real_distribution<double> jitter{0.05, 1.0};
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(n);
    y.reserve(n);
    double position = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        position += jitter(gen);
        x.push_back(position);
        y.push_back(std::sin(position) + 0.1 * position);
    }
    return {std::move(x), std::move(y)};
}

// A negative agreement value means the comparison was skipped: see the note
// under the table.
void
Row(const std::string &what, std::size_t n, double before, double after,
    double agreement) {
    std::printf("  %-26s %6zu %12.3f %12.3f %9.1fx   ", what.c_str(), n, before,
                after, before / after);
    if (agreement < 0.0) {
        std::printf("%s\n", "not compared");
    } else {
        std::printf("%.3g\n", agreement);
    }
}

// Largest relative difference between the two implementations over the
// queries. A difference of accumulated sums hides both cancellation and the
// fact that the replaced implementation loses precision long before ours.
double
RelativeDifference(double a, double b) {
    const auto scale = std::max({std::abs(a), std::abs(b), 1.0});
    return std::abs(a - b) / scale;
}

} // namespace

int
main() {
    constexpr int repetitions = 5;

    std::printf("Interpolation benchmarks (best of %d, milliseconds)\n\n",
                repetitions);
    std::printf("  %-26s %6s %12s %12s %10s   %s\n", "case", "n", "replaced",
                "current", "speedup", "max rel diff");
    std::printf("  %s\n", std::string(88, '-').c_str());

    // Lagrange evaluation.
    for (const std::size_t n : {32u, 128u, 512u}) {
        const auto [x, y] = MakeSamples(n, 20260821ull);
        const Interpolation::Lagrange interpolant{x, y};

        const int queries = 2000;
        const double step = (x.back() - x.front()) / queries;

        const auto before = BestOf(repetitions, [&] {
            double total = 0.0;
            for (int q = 0; q < queries; ++q) {
                total += NaiveLagrange(x, y, x.front() + q * step);
            }
            g_sink = total;
        });

        const auto after = BestOf(repetitions, [&] {
            double total = 0.0;
            for (int q = 0; q < queries; ++q) {
                total += interpolant(x.front() + q * step);
            }
            g_sink = total;
        });

        const auto worst = n > 32 ? -1.0 : [&] {
            double w = 0.0;
            for (int q = 0; q < queries; ++q) {
                const auto at = x.front() + q * step;
                w = std::max(w, RelativeDifference(NaiveLagrange(x, y, at),
                                                   interpolant(at)));
            }
            return w;
        }();
        Row("Lagrange value", n, before, after, worst);
    }

    // Lagrange differentiation, where the old form was cubic.
    for (const std::size_t n : {16u, 32u, 64u}) {
        const auto [x, y] = MakeSamples(n, 20260821ull);
        const Interpolation::Lagrange interpolant{x, y};

        const int queries = 400;
        const double step = (x.back() - x.front()) / queries;

        const auto before = BestOf(repetitions, [&] {
            double total = 0.0;
            for (int q = 0; q < queries; ++q) {
                total += NaiveLagrangeDerivative(x, y, x.front() + q * step);
            }
            g_sink = total;
        });

        const auto after = BestOf(repetitions, [&] {
            double total = 0.0;
            for (int q = 0; q < queries; ++q) {
                total += interpolant.Evaluate<1>(x.front() + q * step);
            }
            g_sink = total;
        });

        const auto worst = n > 32 ? -1.0 : [&] {
            double w = 0.0;
            for (int q = 0; q < queries; ++q) {
                const auto at = x.front() + q * step;
                w = std::max(
                    w, RelativeDifference(NaiveLagrangeDerivative(x, y, at),
                                          interpolant.Evaluate<1>(at)));
            }
            return w;
        }();
        Row("Lagrange derivative", n, before, after, worst);
    }

    // Spline construction, which is where the solve happens.
    for (const std::size_t n : {64u, 256u, 1024u}) {
        const auto [x, y] = MakeSamples(n, 20260821ull);

        const auto before = BestOf(repetitions, [&] {
            const auto solution = DenseSolveSplineSystem(x, y);
            g_sink = std::accumulate(solution.begin(), solution.end(), 0.0);
        });

        const auto after = BestOf(repetitions, [&] {
            const Interpolation::CubicSpline spline{x, y};
            g_sink = spline.Evaluate<2>(x[n / 2]);
        });

        const auto dense = DenseSolveSplineSystem(x, y);
        const Interpolation::CubicSpline spline{x, y};
        double worst = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            worst = std::max(
                worst, RelativeDifference(dense[i], spline.Evaluate<2>(x[i])));
        }
        Row("CubicSpline construction", n, before, after, worst);
    }

    // One grid, many datasets: the shape a radial differentiation operator
    // has. The replaced form is a spline per line and a binary search per
    // node; the current one factorises the grid once and sweeps the nodes.
    for (const std::size_t n : {9u, 33u, 129u}) {
        constexpr std::size_t lineCount = 512;
        const auto [x, unused] = MakeSamples(n, 20260823ull);
        (void) unused;

        std::vector<std::vector<double>> lines;
        lines.reserve(lineCount);
        for (std::size_t line = 0; line < lineCount; ++line) {
            const auto [ignored, y] =
                MakeSamples(n, 20260823ull + 7919ull * line);
            (void) ignored;
            lines.push_back(y);
        }

        const auto before = BestOf(repetitions, [&] {
            double total = 0.0;
            for (const auto &line : lines) {
                const Interpolation::CubicSpline spline{x, line};
                for (std::size_t i = 0; i < n; ++i) {
                    total += spline.Evaluate<1>(x[i]);
                }
            }
            g_sink = total;
        });

        const auto after = BestOf(repetitions, [&] {
            // The factorisation is inside the timed region, so it is paid for
            // exactly once and the comparison stays honest.
            const Interpolation::CubicSplineSystem system{x};
            std::vector<double> curvature(n), derivative(n);
            double total = 0.0;
            for (const auto &line : lines) {
                system.Solve(line, std::span<double>{curvature});
                system.EvaluateAtNodes<1>(line, curvature,
                                          std::span<double>{derivative});
                for (std::size_t i = 0; i < n; ++i) {
                    total += derivative[i];
                }
            }
            g_sink = total;
        });

        const Interpolation::CubicSplineSystem system{x};
        std::vector<double> curvature(n), derivative(n);
        double worst = 0.0;
        for (const auto &line : lines) {
            const Interpolation::CubicSpline spline{x, line};
            system.Solve(line, std::span<double>{curvature});
            system.EvaluateAtNodes<1>(line, curvature,
                                      std::span<double>{derivative});
            for (std::size_t i = 0; i < n; ++i) {
                worst = std::max(worst,
                                 RelativeDifference(derivative[i],
                                                    spline.Evaluate<1>(x[i])));
            }
        }
        Row("d/dx at nodes, 512 lines", n, before, after, worst);
    }

    std::printf(
        "\n  The last column is the largest relative difference between the "
        "two\n  implementations, and should be at rounding level.\n\n"
        "  Node counts for Lagrange are kept small on purpose. High-degree\n"
        "  polynomial interpolation is ill-conditioned whatever the "
        "algorithm, and\n  the replaced form, which multiplies out every "
        "node difference, loses all\n  precision well before the "
        "barycentric one does; comparing them at large\n  n measures "
        "conditioning rather than speed.\n");
    return 0;
}
