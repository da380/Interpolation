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
