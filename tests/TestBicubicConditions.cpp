#include <gtest/gtest.h>

#include <Interpolation/Interpolation.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "TestUtilities.h"

// End conditions chosen per axis on the bicubic spline: one axis clamped with
// its own slope at every node along the edge, the other periodic, or any
// other pairing short of clamping both. The reference throughout is the
// spline of splines — a one-dimensional spline along every line of the first
// axis, then one through what they give along the second — which shares no
// assembly with the tensor product and must agree with it to rounding when
// the tensor product is right.

namespace {

using Interpolation::AxisConditions;
using Interpolation::BicubicSpline;
using Interpolation::BoundaryCondition;
using Interpolation::CubicSpline;
using Interpolation::EdgeSlopes;

constexpr double kPi = std::numbers::pi;

// Colatitude-like and longitude-like axes, both deliberately uneven, each
// with its closing node so that either may be made periodic.
const std::vector<double> kX{0.0, 0.3, 0.9, 1.2, 1.8, 2.3, 2.7, 3.0, kPi};
const std::vector<double> kY{0.0, 0.4, 1.1, 1.5, 2.2, 2.9,      3.3,
                             4.0, 4.4, 5.1, 5.6, 6.0, 2.0 * kPi};

double
Field(double x, double y) {
    return std::cos(x) * (1.0 + 0.5 * std::sin(y)) +
           0.2 * std::cos(2.0 * y) * std::sin(x);
}

double
FieldDx(double x, double y) {
    return -std::sin(x) * (1.0 + 0.5 * std::sin(y)) +
           0.2 * std::cos(2.0 * y) * std::cos(x);
}

double
FieldDy(double x, double y) {
    return 0.5 * std::cos(x) * std::cos(y) -
           0.4 * std::sin(2.0 * y) * std::sin(x);
}

// A grid sampled from `f`, row-major, with the final row and the final
// column copied from the first wherever that axis is to be periodic.
template <typename F>
std::vector<double>
Sampled(const std::vector<double> &x, const std::vector<double> &y, F f,
        bool closeX, bool closeY) {
    const auto rows = x.size();
    const auto columns = y.size();
    std::vector<double> v(rows * columns);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < columns; ++j) {
            v[i * columns + j] = f(x[i], y[j]);
        }
    }
    if (closeX) {
        for (std::size_t j = 0; j < columns; ++j) {
            v[(rows - 1) * columns + j] = v[j];
        }
    }
    if (closeY) {
        for (std::size_t i = 0; i < rows; ++i) {
            v[i * columns + columns - 1] = v[i * columns];
        }
    }
    return v;
}

// One slope per node along an edge, closed if the edge runs along a periodic
// axis.
template <typename F>
std::vector<double>
EdgeOf(const std::vector<double> &along, F slope, bool close) {
    std::vector<double> s;
    s.reserve(along.size());
    for (const auto value : along) {
        s.push_back(slope(value));
    }
    if (close) {
        s.back() = s.front();
    }
    return s;
}

bool
Periodic(AxisConditions conditions) {
    return conditions.low == BoundaryCondition::Periodic;
}

bool
Clamped(BoundaryCondition condition) {
    return condition == BoundaryCondition::Clamped;
}

// The reference. Every line of the first axis is splined under the first
// axis's conditions and read at qx; the second axis is splined through those
// readings. For an Nx-th derivative in x the lines are differentiated first,
// which is right because d^Nx/dx^Nx S(qx, .) is itself a spline in y with
// those nodal values — and, on a clamped y-edge, with the x-spline of that
// edge's slopes, differentiated likewise, as its end derivatives.
template <std::size_t Nx, std::size_t Ny>
double
SplineOfSplines(const std::vector<double> &x, const std::vector<double> &y,
                const std::vector<double> &v, AxisConditions alongX,
                AxisConditions alongY, const EdgeSlopes<double> &slopes,
                double qx, double qy) {
    const auto rows = x.size();
    const auto columns = y.size();
    const auto across = [&](const std::vector<double> &line, double low,
                            double high) {
        const CubicSpline spline{x, line, alongX.low, low, alongX.high, high};
        return spline.template Evaluate<Nx>(qx);
    };

    std::vector<double> readings(columns);
    std::vector<double> line(rows);
    for (std::size_t j = 0; j < columns; ++j) {
        for (std::size_t i = 0; i < rows; ++i) {
            line[i] = v[i * columns + j];
        }
        readings[j] = across(line, Clamped(alongX.low) ? slopes.xLow[j] : 0.0,
                             Clamped(alongX.high) ? slopes.xHigh[j] : 0.0);
    }

    // Only one axis may be clamped, so the edge slopes of the second axis
    // are splined under the first axis's conditions with no end data.
    const auto edge = [&](BoundaryCondition condition,
                          std::span<const double> slope) {
        if (!Clamped(condition)) {
            return 0.0;
        }
        return across(std::vector<double>(slope.begin(), slope.end()), 0.0,
                      0.0);
    };
    const CubicSpline along{y,           readings,
                            alongY.low,  edge(alongY.low, slopes.yLow),
                            alongY.high, edge(alongY.high, slopes.yHigh)};
    return along.template Evaluate<Ny>(qy);
}

struct Case {
    std::string name;
    AxisConditions alongX;
    AxisConditions alongY;
};

std::vector<Case>
Pairings() {
    using enum BoundaryCondition;
    return {
        {"clamped x, periodic y", {Clamped, Clamped}, {Periodic, Periodic}},
        {"periodic x, clamped y", {Periodic, Periodic}, {Clamped, Clamped}},
        {"clamped-natural x, not-a-knot y",
         {Clamped, Natural},
         {NotAKnot, NotAKnot}},
        {"not-a-knot x, natural-clamped y",
         {NotAKnot, NotAKnot},
         {Natural, Clamped}},
        {"periodic x, periodic y", {Periodic, Periodic}, {Periodic, Periodic}},
    };
}

// Grid data and slopes for one pairing, closed wherever it must be.
struct Data {
    std::vector<double> values;
    std::vector<double> xLow, xHigh, yLow, yHigh;

    EdgeSlopes<double> Slopes() const {
        return {.xLow = xLow, .xHigh = xHigh, .yLow = yLow, .yHigh = yHigh};
    }
};

Data
DataFor(const Case &c, const std::vector<double> &x,
        const std::vector<double> &y) {
    const auto closeX = Periodic(c.alongX);
    const auto closeY = Periodic(c.alongY);
    return {
        Sampled(x, y, Field, closeX, closeY),
        EdgeOf(
            y, [&](double b) { return FieldDx(x.front(), b); }, closeY),
        EdgeOf(
            y, [&](double b) { return FieldDx(x.back(), b); }, closeY),
        EdgeOf(
            x, [&](double a) { return FieldDy(a, y.front()); }, closeX),
        EdgeOf(
            x, [&](double a) { return FieldDy(a, y.back()); }, closeX),
    };
}

template <std::size_t Nx, std::size_t Ny>
void
ExpectAgreement(const auto &spline, const Data &data, const Case &c, double a,
                double b) {
    SCOPED_TRACE("order " + std::to_string(Nx) + ", " + std::to_string(Ny));
    const auto actual = spline.template Evaluate<Nx, Ny>(a, b);
    const auto expected = SplineOfSplines<Nx, Ny>(
        kX, kY, data.values, c.alongX, c.alongY, data.Slopes(), a, b);
    if constexpr (Nx + Ny < 4) {
        InterpolationTest::ExpectScaledNear(actual, expected);
    } else {
        // The mixed fourth derivative has divided by four spacings, the
        // smallest here 0.14, and past the grid the edge cell's linear terms
        // multiply that rounding again. Worst seen: 4e-12, at the far
        // extrapolated corner. Evaluated in long double the two routes agree
        // to 1e-15 even there, so this is rounding and not a difference.
        EXPECT_NEAR(actual, expected, 2.0e-11);
    }
}

} // namespace

TEST(BicubicConditions, AgreesWithTheSplineOfSplinesUnderEveryPairing) {
    std::mt19937_64 engine{InterpolationTest::ReportedSeed(20261018)};
    std::uniform_real_distribution<double> qx{-0.3, kPi + 0.3};
    std::uniform_real_distribution<double> qy{-1.0, 2.0 * kPi + 1.0};

    for (const auto &c : Pairings()) {
        SCOPED_TRACE(c.name);
        const auto data = DataFor(c, kX, kY);
        const BicubicSpline spline{kX,       kY,       data.values,
                                   c.alongX, c.alongY, data.Slopes()};

        for (int k = 0; k < 60; ++k) {
            // Queries range past the grid on purpose: a periodic axis reduces
            // them and any other continues the edge cell, and the reference
            // does the same thing through the one-dimensional class.
            const auto a = qx(engine);
            const auto b = qy(engine);
            SCOPED_TRACE(std::to_string(a) + ", " + std::to_string(b));
            ExpectAgreement<0, 0>(spline, data, c, a, b);
            ExpectAgreement<1, 0>(spline, data, c, a, b);
            ExpectAgreement<0, 1>(spline, data, c, a, b);
            ExpectAgreement<1, 1>(spline, data, c, a, b);
            ExpectAgreement<2, 2>(spline, data, c, a, b);
        }
    }
}

TEST(BicubicConditions, PassesThroughEveryGridPoint) {
    for (const auto &c : Pairings()) {
        SCOPED_TRACE(c.name);
        const auto data = DataFor(c, kX, kY);
        const BicubicSpline spline{kX,       kY,       data.values,
                                   c.alongX, c.alongY, data.Slopes()};
        for (std::size_t i = 0; i < kX.size(); ++i) {
            for (std::size_t j = 0; j < kY.size(); ++j) {
                InterpolationTest::ExpectScaledNear(
                    spline(kX[i], kY[j]), data.values[i * kY.size() + j]);
            }
        }
    }
}

TEST(BicubicConditions, AClampedEdgeTakesItsSlopesAtAndBetweenTheNodes) {
    using enum BoundaryCondition;
    const Case c{"", {Clamped, Clamped}, {Periodic, Periodic}};
    const auto data = DataFor(c, kX, kY);
    const BicubicSpline spline{kX,       kY,       data.values,
                               c.alongX, c.alongY, data.Slopes()};

    // At every column the given slope, and between the columns the periodic
    // spline through them: the edge derivative is a spline in y like any
    // other, rather than whatever the assembly happened to leave there.
    const CubicSpline low{kY, data.xLow, Periodic, 0.0, Periodic, 0.0};
    const CubicSpline high{kY, data.xHigh, Periodic, 0.0, Periodic, 0.0};
    for (std::size_t j = 0; j < kY.size(); ++j) {
        SCOPED_TRACE(j);
        InterpolationTest::ExpectScaledNear(
            spline.Evaluate<1, 0>(kX.front(), kY[j]), data.xLow[j]);
        InterpolationTest::ExpectScaledNear(
            spline.Evaluate<1, 0>(kX.back(), kY[j]), data.xHigh[j]);
    }
    for (const double b : {0.1, 0.7, 2.5, 3.7, 5.9, 6.2}) {
        SCOPED_TRACE(b);
        InterpolationTest::ExpectScaledNear(
            spline.Evaluate<1, 0>(kX.front(), b), low(b));
        InterpolationTest::ExpectScaledNear(spline.Evaluate<1, 0>(kX.back(), b),
                                            high(b));
    }
}

TEST(BicubicConditions,
     IsTwiceContinuouslyDifferentiableAcrossThePeriodicSeam) {
    using enum BoundaryCondition;
    const Case c{"", {Clamped, Clamped}, {Periodic, Periodic}};
    const auto data = DataFor(c, kX, kY);
    const BicubicSpline spline{kX,       kY,       data.values,
                               c.alongX, c.alongY, data.Slopes()};

    const auto delta = 1.0e-10;
    const auto below = kY.back() - delta;
    for (const double a : {0.0, 0.5, 1.6, 2.8, kPi}) {
        SCOPED_TRACE(a);
        EXPECT_NEAR((spline.Evaluate<0, 0>(a, below)),
                    (spline.Evaluate<0, 0>(a, delta)), 1.0e-8);
        EXPECT_NEAR((spline.Evaluate<0, 1>(a, below)),
                    (spline.Evaluate<0, 1>(a, delta)), 1.0e-8);
        EXPECT_NEAR((spline.Evaluate<0, 2>(a, below)),
                    (spline.Evaluate<0, 2>(a, delta)), 1.0e-8);
        EXPECT_NEAR((spline.Evaluate<1, 1>(a, below)),
                    (spline.Evaluate<1, 1>(a, delta)), 1.0e-8);

        // And whole periods away it is the same surface.
        for (const int turns : {-2, 1, 3}) {
            const auto b = 1.3 + turns * (kY.back() - kY.front());
            InterpolationTest::ExpectScaledNear(spline(a, b), spline(a, 1.3));
        }
    }
}

TEST(BicubicConditions, ConvergesAtFourthOrder) {
    using enum BoundaryCondition;
    // The shape planetmodel asks for: clamped across the poles with the exact
    // slope at every meridian, periodic in longitude. Worst error over a
    // fixed sample of the whole domain; rates observed: 4.19, 3.96, 4.04.
    std::mt19937_64 engine{InterpolationTest::ReportedSeed(4)};
    std::uniform_real_distribution<double> qx{0.0, kPi};
    std::uniform_real_distribution<double> qy{0.0, 2.0 * kPi};
    std::vector<std::pair<double, double>> sample(400);
    for (auto &[a, b] : sample) {
        a = qx(engine);
        b = qy(engine);
    }

    const Case c{"", {Clamped, Clamped}, {Periodic, Periodic}};
    std::vector<double> worst;
    for (const std::size_t n : {8, 16, 32, 64}) {
        std::vector<double> x(n + 1), y(2 * n + 1);
        for (std::size_t i = 0; i <= n; ++i) {
            x[i] = kPi * static_cast<double>(i) / static_cast<double>(n);
        }
        for (std::size_t j = 0; j <= 2 * n; ++j) {
            y[j] = kPi * static_cast<double>(j) / static_cast<double>(n);
        }
        const auto data = DataFor(c, x, y);
        const BicubicSpline spline{x,        y,        data.values,
                                   c.alongX, c.alongY, data.Slopes()};
        double error = 0.0;
        for (const auto &[a, b] : sample) {
            error = std::max(error, std::abs(spline(a, b) - Field(a, b)));
        }
        worst.push_back(error);
    }
    for (std::size_t i = 0; i + 1 < worst.size(); ++i) {
        SCOPED_TRACE(i);
        const auto rate = std::log2(worst[i] / worst[i + 1]);
        EXPECT_GT(rate, 3.8);
        EXPECT_LT(rate, 4.4);
    }
}

TEST(BicubicConditions, HandlesComplexValues) {
    using enum BoundaryCondition;
    const Case c{"", {Clamped, Clamped}, {Periodic, Periodic}};
    const auto real = DataFor(c, kX, kY);
    // A second real data set for the imaginary part. The claim is algebraic,
    // so any data will do that still closes, and an affine map keeps it so.
    auto imaginary = DataFor(c, kX, kY);
    for (auto *array : {&imaginary.values, &imaginary.xLow, &imaginary.xHigh}) {
        for (auto &entry : *array) {
            entry = 0.5 * entry - 0.25;
        }
    }

    const auto combine = [](const std::vector<double> &re,
                            const std::vector<double> &im) {
        std::vector<std::complex<double>> z(re.size());
        for (std::size_t k = 0; k < z.size(); ++k) {
            z[k] = {re[k], im[k]};
        }
        return z;
    };
    const auto values = combine(real.values, imaginary.values);
    const auto xLow = combine(real.xLow, imaginary.xLow);
    const auto xHigh = combine(real.xHigh, imaginary.xHigh);

    const BicubicSpline spline{
        kX,
        kY,
        values,
        c.alongX,
        c.alongY,
        EdgeSlopes<std::complex<double>>{.xLow = xLow, .xHigh = xHigh}};
    const BicubicSpline re{kX,       kY,       real.values,
                           c.alongX, c.alongY, real.Slopes()};
    const BicubicSpline im{kX,       kY,       imaginary.values,
                           c.alongX, c.alongY, imaginary.Slopes()};

    for (const auto &[a, b] : {std::pair{0.2, 0.3}, std::pair{1.7, 4.1},
                               std::pair{3.0, 6.2}, std::pair{-0.1, 7.0}}) {
        SCOPED_TRACE(std::to_string(a) + ", " + std::to_string(b));
        InterpolationTest::ExpectScaledNear(
            spline.Evaluate<1, 1>(a, b),
            std::complex<double>{re.Evaluate<1, 1>(a, b),
                                 im.Evaluate<1, 1>(a, b)});
    }
}

TEST(BicubicConditions, OneConditionOnEveryEdgeMayBePeriodic) {
    // The short form with Periodic is a spline on the torus; it must be the
    // same object the long form builds.
    using enum BoundaryCondition;
    const Case c{"", {Periodic, Periodic}, {Periodic, Periodic}};
    const auto data = DataFor(c, kX, kY);
    const BicubicSpline shortForm{kX, kY, data.values, Periodic};
    const BicubicSpline longForm{kX, kY, data.values, c.alongX, c.alongY};
    for (const auto &[a, b] : {std::pair{0.2, 0.3}, std::pair{1.7, 4.1},
                               std::pair{-4.0, 9.0}, std::pair{7.5, -2.2}}) {
        SCOPED_TRACE(std::to_string(a) + ", " + std::to_string(b));
        EXPECT_EQ((shortForm.Evaluate<1, 1>(a, b)),
                  (longForm.Evaluate<1, 1>(a, b)));
    }
}

// Refused, and refused by this class in its own words, rather than by a
// solve it delegated to finding the same fault later.
template <typename Make>
void
ExpectRefusedByName(Make make) {
    try {
        make();
        ADD_FAILURE() << "expected std::invalid_argument";
    } catch (const std::invalid_argument &error) {
        EXPECT_TRUE(
            std::string_view{error.what()}.starts_with("BicubicSpline:"))
            << error.what();
    }
}

TEST(BicubicConditions, RefusesWhatItCannotUse) {
    using enum BoundaryCondition;
    const Case clampedX{"", {Clamped, Clamped}, {Periodic, Periodic}};
    const auto data = DataFor(clampedX, kX, kY);
    const auto slopes = data.Slopes();
    const auto build = [&](const std::vector<double> &values,
                           AxisConditions alongX, AxisConditions alongY,
                           const EdgeSlopes<double> &edges) {
        return [&, alongX, alongY, edges] {
            (void) BicubicSpline{kX, kY, values, alongX, alongY, edges};
        };
    };

    // Clamping both axes would need corner data this class does not take.
    const auto both =
        DataFor({"", {Clamped, Clamped}, {Clamped, Clamped}}, kX, kY);
    ExpectRefusedByName(build(both.values, {Clamped, Clamped},
                              {Clamped, Natural}, both.Slopes()));

    // A clamped edge needs a slope at every node along it.
    auto shortEdge = slopes;
    shortEdge.xHigh = shortEdge.xHigh.first(kY.size() - 1);
    ExpectRefusedByName(
        build(data.values, clampedX.alongX, clampedX.alongY, shortEdge));
    ExpectRefusedByName(
        build(data.values, clampedX.alongX, clampedX.alongY, {}));

    // Values, and the slopes along a periodic axis, must close exactly.
    auto open = data.values;
    open[2 * kY.size() + kY.size() - 1] += 1.0e-12;
    ExpectRefusedByName(build(open, clampedX.alongX, clampedX.alongY, slopes));
    auto openEdge = data.xLow;
    openEdge.back() = std::nextafter(openEdge.back(), 10.0);
    auto openSlopes = slopes;
    openSlopes.xLow = openEdge;
    ExpectRefusedByName(
        build(data.values, clampedX.alongX, clampedX.alongY, openSlopes));

    // On the torus, the final row closes as well as the final column.
    const auto torus =
        DataFor({"", {Periodic, Periodic}, {Periodic, Periodic}}, kX, kY);
    EXPECT_NO_THROW((BicubicSpline{kX, kY, torus.values, Periodic}));
    auto openRow = torus.values;
    openRow[(kX.size() - 1) * kY.size() + 3] += 1.0e-12;
    ExpectRefusedByName(
        [&] { (void) BicubicSpline{kX, kY, openRow, Periodic}; });

    // Periodic and NotAKnot are conditions on a whole axis.
    ExpectRefusedByName(
        build(data.values, clampedX.alongX, {Periodic, Natural}, slopes));
    ExpectRefusedByName(
        build(data.values, {NotAKnot, Clamped}, clampedX.alongY, slopes));

    // Three nodes, the closing one included, is too few for a period.
    const std::vector<double> few{0.0, 1.0, 2.0};
    const auto fewValues = Sampled(kX, few, Field, false, true);
    ExpectRefusedByName([&] {
        (void) BicubicSpline{kX, few, fewValues,
                             AxisConditions{Natural, Natural},
                             AxisConditions{Periodic, Periodic}};
    });
    // Natural needs only two, whatever the other axis does.
    const std::vector<double> two{0.0, 1.0};
    const auto twoValues = Sampled(kX, two, Field, true, false);
    EXPECT_NO_THROW(
        (BicubicSpline{kX, two, twoValues, AxisConditions{Periodic, Periodic},
                       AxisConditions{Natural, Natural}}));
}
