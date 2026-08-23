#include <gtest/gtest.h>

#include <Interpolation/Interpolation.hpp>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TestUtilities.h"

namespace {

// A layered profile with a genuine jump at 3.0, flagged in the data by the
// doubled abscissa, which is the convention the samples actually use.
const std::vector<double> kRadius{0.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0};
const std::vector<double> kDensity{8.0, 8.2, 8.6, 9.0, 12.0, 12.4, 12.9, 13.5};

auto
LayeredModel() {
    return Interpolation::SplitAtRepeats(kRadius, kDensity, [](auto r, auto d) {
        return Interpolation::CubicSpline{std::move(r), std::move(d)};
    });
}

// Stands in for a solver: it needs only something it can evaluate and the
// interval to work over, which is what a PieceView carries.
template <typename Layer>
double
MidpointIntegral(const Layer &layer, int steps) {
    const double h = layer.Width() / steps;
    double total = 0.0;
    for (int i = 0; i < steps; ++i) {
        total += layer(layer.Lower() + h * (i + 0.5)) * h;
    }
    return total;
}

} // namespace

TEST(Piecewise, SplitsWhereAnAbscissaRepeats) {
    const auto model = LayeredModel();

    EXPECT_EQ(model.PieceCount(), 2u);
    InterpolationTest::ExpectScaledNear(model.Lower(), 0.0);
    InterpolationTest::ExpectScaledNear(model.Upper(), 6.0);
    InterpolationTest::ExpectScaledNear(model.Breakpoint(1), 3.0);

    InterpolationTest::ExpectScaledNear(model.Piece(0).Lower(), 0.0);
    InterpolationTest::ExpectScaledNear(model.Piece(0).Upper(), 3.0);
    InterpolationTest::ExpectScaledNear(model.Piece(1).Lower(), 3.0);
    InterpolationTest::ExpectScaledNear(model.Piece(1).Upper(), 6.0);
    InterpolationTest::ExpectScaledNear(model.Piece(1).Width(), 3.0);
}

TEST(Piecewise, LimitsExposeTheJumpAtAnInterface) {
    const auto model = LayeredModel();

    const auto [below, above] = model.Limits(3.0);
    InterpolationTest::ExpectScaledNear(below, 9.0);
    InterpolationTest::ExpectScaledNear(above, 12.0);
    InterpolationTest::ExpectScaledNear(above - below, 3.0);

    // Away from a breakpoint the two sides agree exactly.
    for (const double r : {0.5, 1.5, 2.5, 4.5, 5.5}) {
        SCOPED_TRACE(r);
        const auto [left, right] = model.Limits(r);
        EXPECT_DOUBLE_EQ(left, right);
    }
}

TEST(Piecewise, EvaluationIsRightContinuousByDefault) {
    const auto model = LayeredModel();

    // The piece starting at the breakpoint owns it.
    InterpolationTest::ExpectScaledNear(model(3.0), 12.0);
    InterpolationTest::ExpectScaledNear(
        model.EvaluateFrom<0>(3.0, Interpolation::Side::Right), 12.0);
    InterpolationTest::ExpectScaledNear(
        model.EvaluateFrom<0>(3.0, Interpolation::Side::Left), 9.0);

    EXPECT_EQ(model.IndexOf(3.0, Interpolation::Side::Right), 1u);
    EXPECT_EQ(model.IndexOf(3.0, Interpolation::Side::Left), 0u);

    // Outside the domain the nearest piece answers, as the interpolators do.
    EXPECT_EQ(model.IndexOf(-1.0), 0u);
    EXPECT_EQ(model.IndexOf(99.0), 1u);
}

TEST(Piecewise, AnExtractedPieceIsUsableOnItsOwn) {
    const auto model = LayeredModel();
    const auto outer = model.Piece(1);

    // It evaluates exactly as the whole object does inside its interval.
    for (const double r : {3.2, 4.0, 5.5, 6.0}) {
        SCOPED_TRACE(r);
        InterpolationTest::ExpectScaledNear(outer(r), model(r));
        InterpolationTest::ExpectScaledNear(outer.Evaluate<1>(r),
                                            model.Evaluate<1>(r));
    }

    EXPECT_TRUE(outer.Contains(3.0));
    EXPECT_TRUE(outer.Contains(6.0));
    EXPECT_FALSE(outer.Contains(2.9));

    // And it carries the interval a solver needs. Integrating the outer layer
    // by midpoint must agree with the analytic integral of its spline.
    // The extracted layer is itself antidifferentiable, so its exact
    // integral needs no copy of the underlying spline.
    const auto exact =
        Interpolation::Primitive(outer).Integral(outer.Lower(), outer.Upper());
    EXPECT_NEAR(MidpointIntegral(outer, 20000), exact, 1.0e-5);
}

TEST(Piecewise, PieceAtSelectsBySideAndPiecesEnumerates) {
    const auto model = LayeredModel();

    InterpolationTest::ExpectScaledNear(
        model.PieceAt(3.0, Interpolation::Side::Left).Upper(), 3.0);
    InterpolationTest::ExpectScaledNear(
        model.PieceAt(3.0, Interpolation::Side::Right).Lower(), 3.0);

    const auto all = model.Pieces();
    ASSERT_EQ(all.size(), 2u);
    InterpolationTest::ExpectScaledNear(all[0].Lower(), 0.0);
    InterpolationTest::ExpectScaledNear(all[1].Upper(), 6.0);
}

TEST(Piecewise, ComposesWithTheAlgebraThroughRef) {
    using namespace Interpolation;
    const auto model = LayeredModel();

    static_assert(Function1D<decltype(model)>);
    static_assert(Function1D<decltype(model.Piece(0))>);

    // The operators take operands by value, and a model whose pieces own
    // their samples cannot be copied, so it enters the algebra by reference.
    const auto doubled = Ref(model) * 2.0;
    const auto slope = Derivative(Ref(model));

    for (const double r : {1.5, 4.5}) {
        SCOPED_TRACE(r);
        InterpolationTest::ExpectScaledNear(doubled(r), 2.0 * model(r));
        InterpolationTest::ExpectScaledNear(slope(r), model.Evaluate<1>(r));
    }

    // An extracted piece is cheap to copy, so it needs no wrapper.
    const auto scaledPiece = model.Piece(1) * 3.0;
    InterpolationTest::ExpectScaledNear(scaledPiece(4.0), 3.0 * model(4.0));
}

TEST(Piecewise, SamplesWithoutRepeatsGiveASinglePiece) {
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const std::vector<double> y{0.0, 1.0, 4.0, 9.0};
    const auto single = Interpolation::SplitAtRepeats(x, y, [](auto a, auto b) {
        return Interpolation::CubicSpline{std::move(a), std::move(b)};
    });

    EXPECT_EQ(single.PieceCount(), 1u);
    const Interpolation::CubicSpline plain{x, y};
    for (const double q : {0.5, 1.5, 2.5}) {
        SCOPED_TRACE(q);
        InterpolationTest::ExpectScaledNear(single(q), plain(q));
    }
}

TEST(Piecewise, ConstructsDirectlyFromBreakpointsAndPieces) {
    const std::vector<double> xa{0.0, 1.0};
    const std::vector<double> ya{0.0, 1.0};
    const std::vector<double> xb{1.0, 2.0};
    const std::vector<double> yb{5.0, 6.0};

    using Piece =
        Interpolation::Linear<std::ranges::ref_view<const std::vector<double>>,
                              std::ranges::ref_view<const std::vector<double>>>;
    std::vector<Piece> pieces{Piece{xa, ya}, Piece{xb, yb}};

    const Interpolation::Piecewise<Piece> f{{0.0, 1.0, 2.0}, std::move(pieces)};
    EXPECT_EQ(f.PieceCount(), 2u);

    const auto [below, above] = f.Limits(1.0);
    InterpolationTest::ExpectScaledNear(below, 1.0);
    InterpolationTest::ExpectScaledNear(above, 5.0);
}

TEST(Piecewise, RejectsMalformedInput) {
    using Piece =
        Interpolation::Linear<std::ranges::ref_view<const std::vector<double>>,
                              std::ranges::ref_view<const std::vector<double>>>;
    const std::vector<double> xa{0.0, 1.0};
    const std::vector<double> ya{0.0, 1.0};

    // No pieces at all.
    EXPECT_THROW((Interpolation::Piecewise<Piece>{{0.0}, {}}),
                 std::invalid_argument);

    // Breakpoints not one longer than the pieces.
    {
        std::vector<Piece> one{Piece{xa, ya}};
        EXPECT_THROW(
            (Interpolation::Piecewise<Piece>{{0.0, 1.0, 2.0}, std::move(one)}),
            std::invalid_argument);
    }

    // Breakpoints not strictly increasing.
    {
        std::vector<Piece> two{Piece{xa, ya}, Piece{xa, ya}};
        EXPECT_THROW(
            (Interpolation::Piecewise<Piece>{{0.0, 1.0, 1.0}, std::move(two)}),
            std::invalid_argument);
    }
}

TEST(Piecewise, SplitAtRepeatsRejectsMalformedSamples) {
    const auto build = [](auto a, auto b) {
        return Interpolation::CubicSpline{std::move(a), std::move(b)};
    };

    // Lengths differ.
    {
        const std::vector<double> x{0.0, 1.0, 2.0};
        const std::vector<double> y{0.0, 1.0};
        EXPECT_THROW(Interpolation::SplitAtRepeats(x, y, build),
                     std::invalid_argument);
    }

    // Decreasing abscissae are an error, not a breakpoint.
    {
        const std::vector<double> x{0.0, 2.0, 1.0, 3.0};
        const std::vector<double> y{0.0, 1.0, 2.0, 3.0};
        EXPECT_THROW(Interpolation::SplitAtRepeats(x, y, build),
                     std::invalid_argument);
    }

    // A tripled abscissa leaves the intended pieces ambiguous.
    {
        const std::vector<double> x{0.0, 1.0, 1.0, 1.0, 2.0, 3.0};
        const std::vector<double> y{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        EXPECT_THROW(Interpolation::SplitAtRepeats(x, y, build),
                     std::invalid_argument);
    }

    // A run too short for the piece type is rejected by the builder itself.
    {
        const std::vector<double> x{0.0, 0.0, 1.0, 2.0};
        const std::vector<double> y{1.0, 2.0, 3.0, 4.0};
        EXPECT_THROW(Interpolation::SplitAtRepeats(x, y, build),
                     std::invalid_argument);
    }
}

TEST(FunctionRef, BorrowsWithoutCopying) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const std::vector<double> y{0.0, 1.0, 8.0, 27.0};
    const CubicSpline s{x, y};

    const auto borrowed = Ref(s);
    static_assert(Function1D<decltype(borrowed)>);
    for (const double q : {0.5, 1.5, 2.5}) {
        SCOPED_TRACE(q);
        InterpolationTest::ExpectScaledNear(borrowed(q), s(q));
        InterpolationTest::ExpectScaledNear(borrowed.Evaluate<2>(q),
                                            s.Evaluate<2>(q));
    }
    EXPECT_EQ(&borrowed.Function(), &s);
}

TEST(CubicSpline, NotAKnotReproducesACubicExactly) {
    // Nonuniform, with a deliberately large spacing ratio.
    const std::vector<double> x{0.0, 0.31, 1.4, 1.55, 3.2, 5.0, 5.05, 7.0};
    const auto cubic = [](double v) {
        return 2.0 - 0.5 * v + 0.25 * v * v + 1.5 * v * v * v;
    };
    std::vector<double> y;
    for (const auto v : x) {
        y.push_back(cubic(v));
    }

    const Interpolation::CubicSpline notAKnot{
        x,
        y,
        Interpolation::BoundaryCondition::NotAKnot,
        0.0,
        Interpolation::BoundaryCondition::NotAKnot,
        0.0};
    const Interpolation::CubicSpline natural{x, y};

    double worstNotAKnot = 0.0;
    double worstNatural = 0.0;
    for (int i = 0; i <= 700; ++i) {
        const double q = 7.0 * i / 700.0;
        worstNotAKnot =
            std::max(worstNotAKnot, std::abs(notAKnot(q) - cubic(q)));
        worstNatural = std::max(worstNatural, std::abs(natural(q) - cubic(q)));
    }
    EXPECT_LT(worstNotAKnot, 1.0e-9);
    EXPECT_GT(worstNatural, 1.0);

    // Derivatives are exact too.
    InterpolationTest::ExpectScaledNear(notAKnot.Evaluate<1>(2.5),
                                        -0.5 + 0.5 * 2.5 + 4.5 * 2.5 * 2.5);
    InterpolationTest::ExpectScaledNear(notAKnot.Evaluate<2>(2.5),
                                        0.5 + 9.0 * 2.5);
}

TEST(CubicSpline, NotAKnotMustBeUsedAtBothEnds) {
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const std::vector<double> y{0.0, 1.0, 8.0, 27.0, 64.0};

    EXPECT_THROW((Interpolation::CubicSpline{
                     x, y, Interpolation::BoundaryCondition::NotAKnot, 0.0,
                     Interpolation::BoundaryCondition::Natural, 0.0}),
                 std::invalid_argument);

    // And it needs four nodes.
    const std::vector<double> few{0.0, 1.0, 2.0};
    const std::vector<double> fewY{0.0, 1.0, 8.0};
    EXPECT_THROW((Interpolation::CubicSpline{
                     few, fewY, Interpolation::BoundaryCondition::NotAKnot, 0.0,
                     Interpolation::BoundaryCondition::NotAKnot, 0.0}),
                 std::invalid_argument);
}

// --- The conventions a consumer depends on ----------------------------------
//
// GSHTrans has adopted this class's breakpoint semantics wholesale: its
// RadialGrid carries an element partition represented as breakpoints in
// exactly this sense, so that the two libraries cannot disagree about what
// happens at a core-mantle boundary. The tests above exercise the class; these
// state the contract, so that a change to it fails here and is noticed rather
// than discovered downstream.

namespace {

using Contract =
    Interpolation::Linear<std::ranges::ref_view<const std::vector<double>>,
                          std::ranges::ref_view<const std::vector<double>>>;

// Two pieces meeting at 1.0 with different values there: a genuine jump, of
// the kind a material discontinuity produces.
const std::vector<double> kLowerX{0.0, 1.0};
const std::vector<double> kLowerY{0.0, 2.0};
const std::vector<double> kUpperX{1.0, 2.0};
const std::vector<double> kUpperY{10.0, 14.0};

Interpolation::Piecewise<Contract>
Discontinuous() {
    std::vector<Contract> pieces{Contract{kLowerX, kLowerY},
                                 Contract{kUpperX, kUpperY}};
    return Interpolation::Piecewise<Contract>{{0.0, 1.0, 2.0},
                                              std::move(pieces)};
}

} // namespace

TEST(PiecewiseContract, ContinuityAcrossABreakpointIsNotChecked) {
    // Construction must succeed: a jump is data, not an error. This is the
    // clause a consumer relies on most, because enforcing continuity here
    // would make a layered earth model unrepresentable.
    const auto f = Discontinuous();

    ASSERT_EQ(f.PieceCount(), 2u);
    const auto [below, above] = f.Limits(1.0);
    EXPECT_DOUBLE_EQ(below, 2.0);
    EXPECT_DOUBLE_EQ(above, 10.0);
    EXPECT_NE(below, above);
}

TEST(PiecewiseContract, ThePiecesTileTheDomainWithNoGaps) {
    const auto f = Discontinuous();

    EXPECT_DOUBLE_EQ(f.Lower(), f.Breakpoint(0));
    EXPECT_DOUBLE_EQ(f.Upper(), f.Breakpoint(f.PieceCount()));

    for (std::size_t k = 0; k < f.PieceCount(); ++k) {
        SCOPED_TRACE("piece " + std::to_string(k));
        const auto piece = f.Piece(k);
        EXPECT_DOUBLE_EQ(piece.Lower(), f.Breakpoint(k));
        EXPECT_DOUBLE_EQ(piece.Upper(), f.Breakpoint(k + 1));
        // No gap: each piece begins exactly where the last one ended.
        if (k > 0) {
            EXPECT_DOUBLE_EQ(f.Piece(k - 1).Upper(), piece.Lower());
        }
    }
}

TEST(PiecewiseContract, PieceKOwnsTheHalfOpenIntervalFromTheLeft) {
    const auto f = Discontinuous();

    // Right-continuity: landing exactly on an interior breakpoint is answered
    // by the piece starting there.
    EXPECT_EQ(f.IndexOf(1.0), 1u);
    EXPECT_EQ(f.IndexOf(1.0, Interpolation::Side::Right), 1u);
    EXPECT_EQ(f.IndexOf(1.0, Interpolation::Side::Left), 0u);
    EXPECT_DOUBLE_EQ(f(1.0), 10.0);

    // The interior of a piece is unambiguous, whichever side is asked for.
    for (const auto x : {0.25, 0.75, 1.25, 1.75}) {
        SCOPED_TRACE("x = " + std::to_string(x));
        EXPECT_EQ(f.IndexOf(x, Interpolation::Side::Left),
                  f.IndexOf(x, Interpolation::Side::Right));
        const auto [left, right] = f.Limits(x);
        EXPECT_DOUBLE_EQ(left, right);
    }

    // The outer breakpoints have only one piece to answer from.
    EXPECT_EQ(f.IndexOf(0.0, Interpolation::Side::Left), 0u);
    EXPECT_EQ(f.IndexOf(2.0, Interpolation::Side::Right), 1u);
}

TEST(PiecewiseContract, LimitsReportsTheLeftValueFirst) {
    const auto f = Discontinuous();

    const auto limits = f.Limits(1.0);
    EXPECT_DOUBLE_EQ(limits.first,
                     f.EvaluateFrom<0>(1.0, Interpolation::Side::Left));
    EXPECT_DOUBLE_EQ(limits.second,
                     f.EvaluateFrom<0>(1.0, Interpolation::Side::Right));

    // The ordering is load-bearing for a caller reading a jump off a
    // boundary, so it is stated rather than left to be inferred.
    EXPECT_LT(limits.first, limits.second);
}

TEST(PiecewiseContract, DerivativesFollowTheSameConvention) {
    const auto f = Discontinuous();

    // Slopes 2 below and 4 above, so the derivative jumps too. Evaluate is
    // right-continuous at every order, not only at order zero.
    const auto [belowSlope, aboveSlope] = f.Limits<1>(1.0);
    EXPECT_DOUBLE_EQ(belowSlope, 2.0);
    EXPECT_DOUBLE_EQ(aboveSlope, 4.0);
    EXPECT_DOUBLE_EQ(f.Evaluate<1>(1.0), aboveSlope);
}

TEST(PiecewiseContract, OutsideTheDomainTheEndPieceContinues) {
    const auto f = Discontinuous();

    // Not clamped and not an error: the end piece is extrapolated, which is
    // what a solver stepping slightly past the last breakpoint needs.
    EXPECT_EQ(f.IndexOf(-1.0), 0u);
    EXPECT_EQ(f.IndexOf(3.0), f.PieceCount() - 1);
    EXPECT_DOUBLE_EQ(f(-1.0), -2.0);
    EXPECT_DOUBLE_EQ(f(3.0), 18.0);
}
