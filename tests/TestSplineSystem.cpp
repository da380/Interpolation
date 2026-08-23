#include <gtest/gtest.h>

#include <Interpolation/Interpolation.hpp>
#include <cmath>
#include <complex>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "AllocationCounter.hpp"
#include "TestUtilities.h"

// The reusable spline factorisation, and the nodal evaluation path that goes
// with it. Both exist because a caller with many ordinate sets on one grid —
// the components of a field, an ensemble, a line of spectral coefficients —
// should pay for the matrix once and never binary-search for a segment it
// already knows.

namespace {

using Interpolation::BoundaryCondition;
using Interpolation::CubicSpline;
using Interpolation::CubicSplineSystem;
using Interpolation::Side;
using Interpolation::SolveTridiagonal;
using Interpolation::TridiagonalFactorization;

// A diagonally dominant tridiagonal system, so the pivot-free elimination is
// legitimate and the reference solve is well conditioned.
struct System {
    std::vector<double> sub, diag, super;
};

System
DominantSystem(std::size_t n, std::uint64_t seed) {
    std::mt19937_64 engine{seed};
    std::uniform_real_distribution<double> off{-1.0, 1.0};
    std::uniform_real_distribution<double> extra{0.5, 1.5};

    System s{std::vector<double>(n), std::vector<double>(n),
             std::vector<double>(n)};
    for (std::size_t i = 0; i < n; ++i) {
        s.sub[i] = i == 0 ? 0.0 : off(engine);
        s.super[i] = i + 1 == n ? 0.0 : off(engine);
        s.diag[i] = std::abs(s.sub[i]) + std::abs(s.super[i]) + extra(engine);
    }
    return s;
}

template <typename Scalar>
std::vector<Scalar>
RandomVector(std::size_t n, std::uint64_t seed) {
    std::mt19937_64 engine{seed};
    std::uniform_real_distribution<double> value{-2.0, 2.0};
    std::vector<Scalar> v(n);
    for (auto &entry : v) {
        if constexpr (Interpolation::Complex<Scalar>) {
            entry = Scalar{value(engine), value(engine)};
        } else {
            entry = value(engine);
        }
    }
    return v;
}

// Nodes with genuinely uneven spacing: a uniform grid hides the not-a-knot
// pivot problem the slope formulation exists to avoid.
const std::vector<double> kNodes{0.0, 0.4, 1.5, 1.7, 2.9, 4.4, 4.6, 6.0};

std::vector<double>
SampledOn(const std::vector<double> &x) {
    std::vector<double> y;
    y.reserve(x.size());
    for (const auto value : x) {
        y.push_back(std::exp(-0.3 * value) * std::sin(1.7 * value) +
                    0.2 * value);
    }
    return y;
}

std::vector<std::complex<double>>
ComplexSampledOn(const std::vector<double> &x) {
    std::vector<std::complex<double>> y;
    y.reserve(x.size());
    for (const auto value : x) {
        y.emplace_back(std::cos(1.1 * value), std::exp(-0.2 * value));
    }
    return y;
}

} // namespace

// --- TridiagonalFactorization -----------------------------------------------

TEST(TridiagonalFactorization, AgreesWithTheOneShotSolve) {
    const auto seed = InterpolationTest::ReportedSeed(20260823);

    for (const std::size_t n : {1u, 2u, 3u, 9u, 64u}) {
        const auto system = DominantSystem(n, seed + n);
        const auto rhs = RandomVector<double>(n, seed + 1000 + n);

        // The one-shot solve consumes its diagonal, so it gets a copy.
        auto diagCopy = system.diag;
        auto reference = rhs;
        SolveTridiagonal<double, double>(std::span<const double>{system.sub},
                                         std::span<double>{diagCopy},
                                         std::span<const double>{system.super},
                                         std::span<double>{reference});

        const TridiagonalFactorization<double> factorization{
            std::span<const double>{system.sub},
            std::span<const double>{system.diag},
            std::span<const double>{system.super}};
        ASSERT_EQ(factorization.Size(), n);

        auto solution = rhs;
        factorization.Solve(std::span<double>{solution});

        for (std::size_t i = 0; i < n; ++i) {
            SCOPED_TRACE("n = " + std::to_string(n) + ", row " +
                         std::to_string(i));
            EXPECT_NEAR(solution[i], reference[i],
                        1.0e-12 * (1.0 + std::abs(reference[i])));
        }
    }
}

TEST(TridiagonalFactorization, OneFactorisationServesManyRightHandSides) {
    const auto seed = InterpolationTest::ReportedSeed(20260823);
    constexpr std::size_t n = 32;
    const auto system = DominantSystem(n, seed);

    const TridiagonalFactorization<double> factorization{
        std::span<const double>{system.sub},
        std::span<const double>{system.diag},
        std::span<const double>{system.super}};

    // The point of the class: solving repeatedly must not degrade, because
    // nothing the factorisation holds is consumed by a solve.
    for (int k = 0; k < 8; ++k) {
        const auto rhs = RandomVector<double>(n, seed + 7919 * k);

        auto diagCopy = system.diag;
        auto reference = rhs;
        SolveTridiagonal<double, double>(std::span<const double>{system.sub},
                                         std::span<double>{diagCopy},
                                         std::span<const double>{system.super},
                                         std::span<double>{reference});

        auto solution = rhs;
        factorization.Solve(std::span<double>{solution});
        for (std::size_t i = 0; i < n; ++i) {
            SCOPED_TRACE("solve " + std::to_string(k));
            EXPECT_NEAR(solution[i], reference[i],
                        1.0e-12 * (1.0 + std::abs(reference[i])));
        }
    }
}

TEST(TridiagonalFactorization, ARealMatrixActsOnAComplexRightHandSide) {
    const auto seed = InterpolationTest::ReportedSeed(20260823);
    constexpr std::size_t n = 24;
    const auto system = DominantSystem(n, seed);

    const TridiagonalFactorization<double> factorization{
        std::span<const double>{system.sub},
        std::span<const double>{system.diag},
        std::span<const double>{system.super}};

    const auto rhs = RandomVector<std::complex<double>>(n, seed + 11);
    auto solution = rhs;
    factorization.Solve(std::span<std::complex<double>>{solution});

    // The matrix is real, so the real and imaginary parts must solve
    // independently. That is exactly why the class is not templated on the
    // right-hand side type.
    std::vector<double> realPart(n), imagPart(n);
    for (std::size_t i = 0; i < n; ++i) {
        realPart[i] = rhs[i].real();
        imagPart[i] = rhs[i].imag();
    }
    factorization.Solve(std::span<double>{realPart});
    factorization.Solve(std::span<double>{imagPart});

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(solution[i].real(), realPart[i], 1.0e-13);
        EXPECT_NEAR(solution[i].imag(), imagPart[i], 1.0e-13);
    }
}

TEST(TridiagonalFactorization, RejectsMalformedSystems) {
    const std::vector<double> three{1.0, 1.0, 1.0};
    const std::vector<double> two{1.0, 1.0};
    const std::vector<double> none{};

    EXPECT_THROW(
        (TridiagonalFactorization<double>{std::span<const double>{none},
                                          std::span<const double>{none},
                                          std::span<const double>{none}}),
        std::invalid_argument);

    EXPECT_THROW(
        (TridiagonalFactorization<double>{std::span<const double>{two},
                                          std::span<const double>{three},
                                          std::span<const double>{three}}),
        std::invalid_argument);

    // A zero pivot is reported where it is met, rather than left to produce
    // infinities inside somebody's solve.
    const std::vector<double> zeroDiag{0.0, 1.0, 1.0};
    EXPECT_THROW(
        (TridiagonalFactorization<double>{std::span<const double>{three},
                                          std::span<const double>{zeroDiag},
                                          std::span<const double>{three}}),
        std::invalid_argument);

    // Dominant, so this one factorises: the check under test is the length of
    // the right-hand side, not the matrix.
    const std::vector<double> dominant{4.0, 4.0, 4.0};
    const TridiagonalFactorization<double> factorization{
        std::span<const double>{three}, std::span<const double>{dominant},
        std::span<const double>{three}};
    std::vector<double> wrongLength(2);
    EXPECT_THROW(factorization.Solve(std::span<double>{wrongLength}),
                 std::invalid_argument);
}

// --- CubicSplineSystem ------------------------------------------------------

template <typename Scalar>
void
CheckSystemReproducesTheSpline(const std::vector<double> &x,
                               const std::vector<Scalar> &y,
                               BoundaryCondition left, Scalar leftDerivative,
                               BoundaryCondition right,
                               Scalar rightDerivative) {
    const CubicSpline spline{
        x, y, left, leftDerivative, right, rightDerivative};
    const CubicSplineSystem system{x, left, right};

    std::vector<Scalar> curvature(x.size());
    system.Solve(y, leftDerivative, rightDerivative,
                 std::span<Scalar>{curvature});

    const auto reference = spline.Curvatures();
    ASSERT_EQ(curvature.size(), reference.size());
    for (std::size_t i = 0; i < curvature.size(); ++i) {
        SCOPED_TRACE("node " + std::to_string(i));
        // Bit-for-bit: the spline is solved by this very system, so anything
        // else would mean two assemblies had crept back in.
        EXPECT_EQ(curvature[i], reference[i]);
    }
}

TEST(CubicSplineSystem, ReproducesTheSplineForEveryBoundaryCondition) {
    const auto y = SampledOn(kNodes);

    CheckSystemReproducesTheSpline(kNodes, y, BoundaryCondition::Natural, 0.0,
                                   BoundaryCondition::Natural, 0.0);
    CheckSystemReproducesTheSpline(kNodes, y, BoundaryCondition::NotAKnot, 0.0,
                                   BoundaryCondition::NotAKnot, 0.0);
    CheckSystemReproducesTheSpline(kNodes, y, BoundaryCondition::Clamped, -0.7,
                                   BoundaryCondition::Clamped, 1.3);
    // Mixed ends: the derivative given for the Natural end is ignored.
    CheckSystemReproducesTheSpline(kNodes, y, BoundaryCondition::Clamped, -0.7,
                                   BoundaryCondition::Natural, 99.0);
    CheckSystemReproducesTheSpline(kNodes, y, BoundaryCondition::Natural, 99.0,
                                   BoundaryCondition::Clamped, 1.3);
}

TEST(CubicSplineSystem, ReproducesTheSplineForComplexOrdinates) {
    const auto y = ComplexSampledOn(kNodes);
    using C = std::complex<double>;

    CheckSystemReproducesTheSpline(kNodes, y, BoundaryCondition::Natural, C{},
                                   BoundaryCondition::Natural, C{});
    CheckSystemReproducesTheSpline(kNodes, y, BoundaryCondition::NotAKnot, C{},
                                   BoundaryCondition::NotAKnot, C{});
    CheckSystemReproducesTheSpline(kNodes, y, BoundaryCondition::Clamped,
                                   C{0.4, -1.1}, BoundaryCondition::Clamped,
                                   C{-0.2, 0.9});
}

TEST(CubicSplineSystem, OneFactorisationSolvesManyOrdinateSets) {
    const auto seed = InterpolationTest::ReportedSeed(20260823);
    const CubicSplineSystem system{kNodes, BoundaryCondition::NotAKnot};

    std::vector<double> curvature(kNodes.size());
    for (int line = 0; line < 16; ++line) {
        const auto y = RandomVector<double>(kNodes.size(), seed + 31 * line);

        system.Solve(y, std::span<double>{curvature});

        // Against a spline built from scratch on the same data, which is the
        // implementation GSHTrans is replacing.
        const CubicSpline spline{kNodes, y, BoundaryCondition::NotAKnot, 0.0,
                                 0.0};
        const auto reference = spline.Curvatures();
        for (std::size_t i = 0; i < curvature.size(); ++i) {
            SCOPED_TRACE("line " + std::to_string(line) + ", node " +
                         std::to_string(i));
            EXPECT_EQ(curvature[i], reference[i]);
        }
    }
}

TEST(CubicSplineSystem, SolvingAllocatesNothing) {
    if constexpr (!InterpolationTest::AllocationCountingEnabled()) {
        GTEST_SKIP() << "allocation counting stands down under sanitizers";
    } else {
        {
            const auto before = InterpolationTest::AllocationCount();
            void *probe = ::operator new(64);
            ::operator delete(probe);
            ASSERT_GT(InterpolationTest::AllocationCount(), before)
                << "the allocation counter is not observing allocations";
        }

        // This is the claim that makes the class worth having: the matrix is
        // built once, and a solve after that is arithmetic and nothing else.
        for (const auto condition :
             {BoundaryCondition::Natural, BoundaryCondition::NotAKnot}) {
            const CubicSplineSystem system{kNodes, condition};
            const auto y = SampledOn(kNodes);
            std::vector<double> curvature(kNodes.size());
            std::vector<double> nodal(kNodes.size());

            const auto before = InterpolationTest::AllocationCount();
            for (int line = 0; line < 64; ++line) {
                system.Solve(y, std::span<double>{curvature});
                system.EvaluateAtNodes<1>(y, curvature,
                                          std::span<double>{nodal});
            }
            const auto after = InterpolationTest::AllocationCount();

            EXPECT_EQ(after, before)
                << "solving allocated " << (after - before) << " times";
        }
    }
}

TEST(CubicSplineSystem, CurvaturesReturnsWhatSolveWrites) {
    const auto y = SampledOn(kNodes);
    const CubicSplineSystem system{kNodes, BoundaryCondition::Natural};

    std::vector<double> curvature(kNodes.size());
    system.Solve(y, std::span<double>{curvature});
    const auto returned = system.Curvatures(y);

    ASSERT_EQ(returned.size(), curvature.size());
    for (std::size_t i = 0; i < curvature.size(); ++i) {
        EXPECT_EQ(returned[i], curvature[i]);
    }
}

TEST(CubicSplineSystem, ReportsTheGridItWasBuiltOn) {
    const CubicSplineSystem system{kNodes, BoundaryCondition::Clamped,
                                   BoundaryCondition::Natural};
    ASSERT_EQ(system.Size(), kNodes.size());
    EXPECT_EQ(system.LeftCondition(), BoundaryCondition::Clamped);
    EXPECT_EQ(system.RightCondition(), BoundaryCondition::Natural);
    for (std::size_t i = 0; i < kNodes.size(); ++i) {
        EXPECT_EQ(system.Node(i), kNodes[i]);
    }
    for (std::size_t i = 0; i + 1 < kNodes.size(); ++i) {
        EXPECT_EQ(system.Spacing(i), kNodes[i + 1] - kNodes[i]);
    }
}

TEST(CubicSplineSystem, RejectsMalformedInput) {
    const std::vector<double> decreasing{0.0, 2.0, 1.0, 3.0};
    EXPECT_THROW((CubicSplineSystem{decreasing}), std::invalid_argument);

    const std::vector<double> tooShort{0.0};
    EXPECT_THROW((CubicSplineSystem{tooShort}), std::invalid_argument);

    // NotAKnot needs four nodes and both ends.
    const std::vector<double> three{0.0, 1.0, 2.0};
    EXPECT_THROW((CubicSplineSystem{three, BoundaryCondition::NotAKnot}),
                 std::invalid_argument);
    EXPECT_THROW((CubicSplineSystem{kNodes, BoundaryCondition::NotAKnot,
                                    BoundaryCondition::Natural}),
                 std::invalid_argument);

    const CubicSplineSystem system{kNodes};
    const auto y = SampledOn(kNodes);
    std::vector<double> curvature(kNodes.size());
    std::vector<double> wrongLength(kNodes.size() - 1);

    EXPECT_THROW(system.Solve(wrongLength, std::span<double>{curvature}),
                 std::invalid_argument);
    EXPECT_THROW(system.Solve(y, std::span<double>{wrongLength}),
                 std::invalid_argument);

    // A Clamped end has no derivative to use unless one is handed over.
    const CubicSplineSystem clamped{kNodes, BoundaryCondition::Clamped};
    EXPECT_THROW(clamped.Solve(y, std::span<double>{curvature}),
                 std::invalid_argument);
    EXPECT_NO_THROW(clamped.Solve(y, 0.5, -0.5, std::span<double>{curvature}));
}

TEST(CubicSplineSystem, ASplineExposesItsOwnFactorisation) {
    const auto y = SampledOn(kNodes);
    const CubicSpline spline{kNodes, y};

    // Having fitted one line, a caller who finds more data on the same grid
    // reuses the matrix rather than assembling it again.
    const auto other = ComplexSampledOn(kNodes);
    const auto curvature = spline.SplineSystem().Curvatures(other);

    const CubicSpline reference{kNodes, other};
    const auto expected = reference.Curvatures();
    ASSERT_EQ(curvature.size(), expected.size());
    for (std::size_t i = 0; i < curvature.size(); ++i) {
        EXPECT_EQ(curvature[i], expected[i]);
    }
}
