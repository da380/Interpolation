// 11 - One grid, many datasets
//
// A cubic spline's matrix depends on the nodes alone; only the right-hand side
// carries the ordinates. So a caller with many datasets on one grid should
// factorise once and solve many times, and ask for the derivative at every
// node in one sweep rather than searching for segments it already knows.
//
// That is what CubicSplineSystem is for, and it is the shape a radial
// differentiation operator has: one grid, one line of data per column of a
// field, applied over and over.

#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/CubicSplineSystem.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <span>
#include <vector>

int
main() {
    using namespace Interpolation;
    using Complex = std::complex<double>;

    // A radial grid, and a field of complex coefficients sampled on it: one
    // line per coefficient, all sharing the same nodes.
    constexpr std::size_t nodeCount = 65;
    constexpr std::size_t lineCount = 2000;

    std::vector<double> radius;
    for (std::size_t i = 0; i < nodeCount; ++i) {
        // Deliberately uneven, thickening towards the surface.
        const double t = static_cast<double>(i) / (nodeCount - 1);
        radius.push_back(3480.0 + 2891.0 * t * t);
    }

    std::vector<std::vector<Complex>> lines(lineCount,
                                            std::vector<Complex>(nodeCount));
    for (std::size_t line = 0; line < lineCount; ++line) {
        const double k = 1.0 + 0.01 * static_cast<double>(line);
        for (std::size_t i = 0; i < nodeCount; ++i) {
            const double s = (radius[i] - 3480.0) / 2891.0;
            lines[line][i] = Complex{std::sin(k * s), std::cos(k * s)};
        }
    }

    std::printf("A field of %zu complex lines on one grid of %zu nodes.\n\n",
                lineCount, nodeCount);

    // The matrix is assembled and eliminated once, here. Everything after
    // this point is arithmetic on a right-hand side.
    const auto system = CubicSplineSystem{radius, BoundaryCondition::NotAKnot};

    // Two buffers for the whole sweep, not two per line. Solve and
    // EvaluateAtNodes both write into storage the caller owns, so the loop
    // below allocates nothing at all.
    std::vector<Complex> curvature(nodeCount);
    std::vector<Complex> derivative(nodeCount);

    const auto start = std::chrono::steady_clock::now();
    double checksum = 0.0;
    for (const auto &line : lines) {
        system.Solve(line, std::span{curvature});
        system.EvaluateAtNodes<1>(line, curvature, std::span{derivative});
        checksum += std::abs(derivative[nodeCount / 2]);
    }
    const auto elapsed = std::chrono::steady_clock::now() - start;

    std::printf("  d/dr at every node of every line: %.2f ms\n",
                std::chrono::duration<double, std::milli>(elapsed).count());
    std::printf("  checksum %.6f\n\n", checksum);

    // The same numbers, the long way round: a fresh spline per line, and a
    // query per node. Both re-do work the grid already determined.
    const auto slowStart = std::chrono::steady_clock::now();
    double slowChecksum = 0.0;
    for (const auto &line : lines) {
        const CubicSpline spline{radius, line, BoundaryCondition::NotAKnot,
                                 Complex{}, Complex{}};
        for (std::size_t i = 0; i < nodeCount; ++i) {
            const auto d = spline.Evaluate<1>(radius[i]);
            if (i == nodeCount / 2) {
                slowChecksum += std::abs(d);
            }
        }
    }
    const auto slowElapsed = std::chrono::steady_clock::now() - slowStart;

    std::printf("  a spline and a search per line:   %.2f ms\n",
                std::chrono::duration<double, std::milli>(slowElapsed).count());
    std::printf("  checksum %.6f\n\n", slowChecksum);

    // Identical answers: this is a reorganisation, not an approximation.
    std::printf("The two agree to %.1e.\n\n",
                std::abs(checksum - slowChecksum));

    // A spline that has already been built hands its factorisation back, for
    // the case where the second dataset turns up after the first fit.
    const CubicSpline first{radius, lines[0], BoundaryCondition::NotAKnot,
                            Complex{}, Complex{}};
    const auto reused = first.SplineSystem().Curvatures(lines[1]);
    std::printf("Reusing one spline's system on another line gives %zu "
                "curvatures without rebuilding the matrix.\n",
                reused.size());

    return 0;
}
