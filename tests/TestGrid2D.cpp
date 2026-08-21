#include <gtest/gtest.h>

#include <Interpolation/Interpolation.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TestUtilities.h"

namespace {

// Row-major sampling of f on the tensor grid x by y.
template <typename F>
std::vector<double>
SampleGrid(const std::vector<double> &x, const std::vector<double> &y, F f) {
    std::vector<double> values;
    values.reserve(x.size() * y.size());
    for (const auto xi : x) {
        for (const auto yj : y) {
            values.push_back(f(xi, yj));
        }
    }
    return values;
}

const std::vector<double> kX{0.0, 0.7, 1.9, 3.0};
const std::vector<double> kY{0.0, 1.1, 2.0, 2.6, 4.0};

} // namespace

TEST(Extents, IndexesRowMajorAndReportsSize) {
    const Interpolation::Extents2D extents{{3, 5}};
    EXPECT_EQ(extents.Order(), 2u);
    EXPECT_EQ(extents.Extent(0), 3u);
    EXPECT_EQ(extents.Extent(1), 5u);
    EXPECT_EQ(extents.Size(), 15u);
    EXPECT_EQ(extents.Index({0, 0}), 0u);
    EXPECT_EQ(extents.Index({0, 4}), 4u);
    EXPECT_EQ(extents.Index({1, 0}), 5u);
    EXPECT_EQ(extents.Index({2, 4}), 14u);
}

TEST(Grid2DView, ReadsThroughRowMajorOffsets) {
    const std::vector<double> flat{0, 1, 2, 3, 4, 5};
    const Interpolation::Grid2DView view{std::views::all(flat),
                                         Interpolation::Extents2D{{2, 3}}};
    EXPECT_EQ(view.Rows(), 2u);
    EXPECT_EQ(view.Columns(), 3u);
    EXPECT_DOUBLE_EQ(view(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(view(0, 2), 2.0);
    EXPECT_DOUBLE_EQ(view(1, 1), 4.0);
}

TEST(Bilinear, ReproducesABilinearSurfaceExactly) {
    const auto v = SampleGrid(kX, kY, [](double a, double b) {
        return 2.0 + 3.0 * a - 1.5 * b + 0.75 * a * b;
    });
    const Interpolation::Bilinear f{kX, kY, v};

    for (const auto &[px, py] :
         {std::pair{0.4, 0.5}, std::pair{1.5, 3.2}, std::pair{2.8, 1.0},
          std::pair{0.0, 0.0}, std::pair{3.0, 4.0}}) {
        SCOPED_TRACE(px * 100 + py);
        InterpolationTest::ExpectScaledNear(
            f(px, py), 2.0 + 3.0 * px - 1.5 * py + 0.75 * px * py);
    }
}

TEST(Bilinear, PartialDerivativesFollowFromSeparability) {
    const auto v = SampleGrid(kX, kY, [](double a, double b) {
        return 3.0 * a - 1.5 * b + 0.75 * a * b;
    });
    const Interpolation::Bilinear f{kX, kY, v};

    const double px = 1.5;
    const double py = 2.3;
    InterpolationTest::ExpectScaledNear(f.Evaluate<1, 0>(px, py),
                                        3.0 + 0.75 * py);
    InterpolationTest::ExpectScaledNear(f.Evaluate<0, 1>(px, py),
                                        -1.5 + 0.75 * px);
    InterpolationTest::ExpectScaledNear(f.Evaluate<1, 1>(px, py), 0.75);

    // The interpolant is bilinear, so anything above first order vanishes.
    InterpolationTest::ExpectScaledNear(f.Evaluate<2, 0>(px, py), 0.0);
    InterpolationTest::ExpectScaledNear(f.Evaluate<0, 2>(px, py), 0.0);
}

TEST(BicubicSpline, ReproducesABilinearSurfaceExactly) {
    const auto v = SampleGrid(kX, kY, [](double a, double b) {
        return 2.0 + 3.0 * a - 1.5 * b + 0.75 * a * b;
    });
    const Interpolation::BicubicSpline f{kX, kY, v};

    for (const auto &[px, py] :
         {std::pair{0.4, 0.5}, std::pair{1.5, 3.2}, std::pair{2.8, 1.0}}) {
        SCOPED_TRACE(px * 100 + py);
        InterpolationTest::ExpectScaledNear(
            f(px, py), 2.0 + 3.0 * px - 1.5 * py + 0.75 * px * py);
        InterpolationTest::ExpectScaledNear(f.Evaluate<1, 1>(px, py), 0.75);
    }
}

TEST(BicubicSpline, PassesThroughEveryGridPoint) {
    const auto v = SampleGrid(kX, kY, [](double a, double b) {
        return std::sin(a) * std::cos(b) + 0.3 * a * b;
    });
    const Interpolation::BicubicSpline f{kX, kY, v};

    for (std::size_t i = 0; i < kX.size(); ++i) {
        for (std::size_t j = 0; j < kY.size(); ++j) {
            SCOPED_TRACE(i * 10 + j);
            InterpolationTest::ExpectScaledNear(f(kX[i], kY[j]),
                                                v[i * kY.size() + j]);
        }
    }
}

TEST(BicubicSpline, ReducesToTheOneDimensionalSplineAlongAConstantAxis) {
    // With the second axis carrying no variation, the tensor product must
    // collapse onto the ordinary one-dimensional spline in the first axis.
    const auto v =
        SampleGrid(kX, kY, [](double a, double) { return std::sin(a); });
    const Interpolation::BicubicSpline surface{kX, kY, v};

    std::vector<double> line;
    for (const auto xi : kX) {
        line.push_back(std::sin(xi));
    }
    const Interpolation::CubicSpline curve{kX, line};

    for (const double px : {0.2, 0.7, 1.4, 2.5, 3.0}) {
        SCOPED_TRACE(px);
        InterpolationTest::ExpectScaledNear(surface(px, 2.0), curve(px));
        InterpolationTest::ExpectScaledNear(surface.Evaluate<1, 0>(px, 2.0),
                                            curve.Evaluate<1>(px));
        InterpolationTest::ExpectScaledNear(surface.Evaluate<2, 0>(px, 2.0),
                                            curve.Evaluate<2>(px));
    }
}

TEST(BicubicSpline, IsFarMoreAccurateThanBilinearInTheInterior) {
    const auto exact = [](double a, double b) {
        return std::sin(a) * std::cos(b);
    };

    const int n = 33;
    std::vector<double> x;
    std::vector<double> y;
    for (int i = 0; i < n; ++i) {
        x.push_back(3.0 * i / (n - 1));
        y.push_back(3.0 * i / (n - 1));
    }
    const auto v = SampleGrid(x, y, exact);

    const Interpolation::Bilinear bilinear{x, y, v};
    const Interpolation::BicubicSpline bicubic{x, y, v};

    double bilinearError = 0.0;
    double bicubicError = 0.0;
    const int m = 41;
    for (int a = 0; a < m; ++a) {
        for (int b = 0; b < m; ++b) {
            // Interior only. The natural condition on all four edges leaves
            // an O(h^2) error there, which would otherwise dominate and hide
            // the interior order of the scheme.
            const double px = 0.6 + 1.8 * a / (m - 1);
            const double py = 0.6 + 1.8 * b / (m - 1);
            const double truth = exact(px, py);
            bilinearError =
                std::max(bilinearError, std::abs(bilinear(px, py) - truth));
            bicubicError =
                std::max(bicubicError, std::abs(bicubic(px, py) - truth));
        }
    }

    EXPECT_LT(bicubicError, 1.0e-5);
    EXPECT_GT(bilinearError / bicubicError, 100.0)
        << "bilinear " << bilinearError << " bicubic " << bicubicError;
}

TEST(BicubicSpline, HandlesComplexValues) {
    using Complex = std::complex<double>;
    std::vector<Complex> v;
    for (const auto xi : kX) {
        for (const auto yj : kY) {
            v.emplace_back(xi * yj, xi - yj);
        }
    }
    const Interpolation::BicubicSpline f{kX, kY, v};

    for (const auto &[px, py] : {std::pair{0.4, 0.5}, std::pair{2.0, 3.0}}) {
        SCOPED_TRACE(px * 100 + py);
        InterpolationTest::ExpectScaledNear(f(px, py),
                                            Complex{px * py, px - py});
    }
}

TEST(Grid2D, RejectsBadInput) {
    const auto v = SampleGrid(kX, kY, [](double a, double b) { return a * b; });

    // Value count that does not match the axes.
    const std::vector<double> shortValues(v.begin(), v.end() - 1);
    EXPECT_THROW((Interpolation::Bilinear{kX, kY, shortValues}),
                 std::invalid_argument);
    EXPECT_THROW((Interpolation::BicubicSpline{kX, kY, shortValues}),
                 std::invalid_argument);

    // A repeated node on an axis.
    const std::vector<double> repeated{0.0, 1.0, 1.0, 2.0};
    const auto rv =
        SampleGrid(repeated, kY, [](double a, double b) { return a * b; });
    EXPECT_THROW((Interpolation::Bilinear{repeated, kY, rv}),
                 std::invalid_argument);
    EXPECT_THROW((Interpolation::BicubicSpline{repeated, kY, rv}),
                 std::invalid_argument);

    // An axis with too few nodes.
    const std::vector<double> single{1.0};
    const auto sv =
        SampleGrid(single, kY, [](double a, double b) { return a * b; });
    EXPECT_THROW((Interpolation::Bilinear{single, kY, sv}),
                 std::invalid_argument);
}

TEST(Grid2D, BorrowsLvaluesAndOwnsRvalues) {
    const auto make = [] {
        std::vector<double> x{0.0, 1.0, 2.0};
        std::vector<double> y{0.0, 1.0, 2.0};
        std::vector<double> v;
        for (const auto xi : x) {
            for (const auto yj : y) {
                v.push_back(xi * yj);
            }
        }
        return Interpolation::Bilinear{std::move(x), std::move(y),
                                       std::move(v)};
    };

    // Every source container is gone by the time this is evaluated.
    const auto owning = make();
    InterpolationTest::ExpectScaledNear(owning(0.5, 0.5), 0.25);
    EXPECT_EQ(owning.Rows(), 3u);
    EXPECT_EQ(owning.Columns(), 3u);
}

TEST(Grid2D, ExtrapolatesByContinuingTheEdgeCell) {
    const auto v = SampleGrid(kX, kY, [](double a, double b) { return a * b; });
    const Interpolation::Bilinear f{kX, kY, v};

    // Outside the grid the nearest cell's polynomial continues, which for a
    // bilinear surface remains exact.
    InterpolationTest::ExpectScaledNear(f(-1.0, 1.0), -1.0);
    InterpolationTest::ExpectScaledNear(f(4.0, 1.0), 4.0);
    InterpolationTest::ExpectScaledNear(f(1.0, -1.0), -1.0);
}
