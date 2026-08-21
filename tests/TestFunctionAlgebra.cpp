#include <gtest/gtest.h>

#include <Interpolation/Function.hpp>
#include <Interpolation/Interpolation.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <new>
#include <vector>

#include "AllocationCounter.hpp"
#include "TestUtilities.h"

namespace {

std::vector<double>
Cubes(const std::vector<double> &x) {
    std::vector<double> y;
    for (const auto value : x) {
        y.push_back(value * value * value);
    }
    return y;
}

} // namespace

TEST(FunctionAlgebra, InterpolatorsAndNodesModelFunction1D) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const auto y = Cubes(x);

    static_assert(
        Function1D<Linear<std::ranges::ref_view<const std::vector<double>>,
                          std::ranges::ref_view<const std::vector<double>>>>);
    static_assert(Function1D<Polynomial<double>>);
    static_assert(Function1D<Polynomial<std::complex<double>>>);
    static_assert(Function1D<Constant<double, double>>);
    static_assert(Function1D<Identity<double>>);

    const CubicSpline s{x, y};
    static_assert(Function1D<decltype(s)>);
    static_assert(Function1D<decltype(Derivative(s))>);
    static_assert(Function1D<decltype(s * s)>);
    static_assert(Function1D<decltype(s + 2.0)>);
    SUCCEED();
}

TEST(FunctionAlgebra, RoadmapExpressionMatchesTheAnalyticAnswer) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const auto y = Cubes(x);
    const CubicSpline s{x, y};

    const auto g = Derivative(s) * s + 2.0;

    for (const double query : {0.25, 1.0, 1.5, 2.75, 3.5}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(
            g(query), s.Evaluate<1>(query) * s(query) + 2.0);
    }
}

TEST(FunctionAlgebra, EvaluationDoesNotAllocate) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const auto y = Cubes(x);
    const CubicSpline s{x, y};

    if constexpr (!InterpolationTest::AllocationCountingEnabled()) {
        GTEST_SKIP() << "allocation counting is disabled under sanitizers";
    } else {

        // Guard against a vacuous pass: if the global operator new replacement
        // were not linked in, the counter would never move and the check below
        // would succeed without measuring anything.
        {
            // A direct call to operator new rather than a new-expression: the
            // standard lets a compiler elide the allocation for the latter, and
            // Clang does at -O2, which would make this guard fail spuriously.
            const auto before = InterpolationTest::AllocationCount();
            void *probe = ::operator new(64);
            ::operator delete(probe);
            ASSERT_GT(InterpolationTest::AllocationCount(), before)
                << "the allocation counter is not observing allocations";
        }

        // Build the expression first: construction may allocate, evaluation
        // must not. Borrowing lvalue ranges means the nodes hold views, not
        // data.
        const auto g = Derivative(s) * s + 2.0;
        volatile double sink = 0.0;

        const auto before = InterpolationTest::AllocationCount();
        for (int i = 0; i < 1000; ++i) {
            sink = sink + g(0.001 * i);
        }
        const auto after = InterpolationTest::AllocationCount();

        EXPECT_EQ(after, before) << "evaluating the expression allocated "
                                 << (after - before) << " times";
    }
}

TEST(FunctionAlgebra, NodesOwnTheirOperandsSoTemporariesDoNotDangle) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const auto y = Cubes(x);

    // Every operand here is a temporary. Nodes holding references would
    // dangle immediately; holding by value is what makes this valid.
    const auto g = [&] {
        return Interpolation::CubicSpline{x, y} *
                   Interpolation::CubicSpline{x, y} +
               Interpolation::Linear{x, y};
    }();

    const CubicSpline s{x, y};
    const Linear l{x, y};
    for (const double query : {0.5, 1.5, 2.5}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(g(query),
                                            s(query) * s(query) + l(query));
    }
}

TEST(FunctionAlgebra, ProductObeysLeibnizAtHigherOrders) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const auto y = Cubes(x);
    const CubicSpline s{x, y};

    const auto p = s * s;
    const double at = 1.5;

    // (f g)'' = f'' g + 2 f' g' + f g''
    InterpolationTest::ExpectScaledNear(
        p.Evaluate<2>(at), 2.0 * (s.Evaluate<2>(at) * s(at) +
                                  s.Evaluate<1>(at) * s.Evaluate<1>(at)));

    // (f g)''' = f''' g + 3 f'' g' + 3 f' g'' + f g'''
    InterpolationTest::ExpectScaledNear(
        p.Evaluate<3>(at), 2.0 * s.Evaluate<3>(at) * s(at) +
                               6.0 * s.Evaluate<2>(at) * s.Evaluate<1>(at));
}

TEST(FunctionAlgebra, DerivativeOfADerivativeIsTheSecondDerivative) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const auto y = Cubes(x);
    const CubicSpline s{x, y};

    const auto second = Derivative(Derivative(s));
    const auto alsoSecond = Derivative<2>(s);

    for (const double query : {0.5, 1.5, 2.5, 3.5}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(second(query),
                                            s.Evaluate<2>(query));
        InterpolationTest::ExpectScaledNear(alsoSecond(query),
                                            s.Evaluate<2>(query));
    }
}

TEST(FunctionAlgebra, CompositionUsesTheChainRule) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const auto y = Cubes(x);
    const CubicSpline s{x, y};

    // g(x) = x/2, so f(g(x)) = s(x/2) and its derivative is s'(x/2)/2.
    const auto halved = Compose(s, Identity<double>{} * 0.5);
    for (const double query : {0.5, 2.0, 3.5, 6.0}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(halved(query), s(0.5 * query));
        InterpolationTest::ExpectScaledNear(halved.Evaluate<1>(query),
                                            0.5 * s.Evaluate<1>(0.5 * query));
    }
}

TEST(FunctionAlgebra, PrimitiveIntegratesExactlyWhereItShould) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 0.7, 1.9, 2.4, 3.0};
    const auto y = Cubes(x);

    // A clamped spline reproduces cubic data exactly, so its integral over
    // the sampled interval must be the exact one.
    const CubicSpline s{x, y, BoundaryCondition::Clamped, 0.0, 27.0};
    const auto S = Primitive(s);
    InterpolationTest::ExpectScaledNear(S.Integral(0.0, 3.0), 81.0 / 4.0);
    InterpolationTest::ExpectScaledNear(S(0.0), 0.0);

    // Differentiating the primitive returns the original function.
    const auto back = Derivative(S);
    for (const double query : {0.3, 1.3, 2.2, 2.9}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(back(query), s(query));
    }

    // The piecewise-linear interpolant integrates to the trapezoid rule.
    const Linear l{x, y};
    const auto L = Primitive(l);
    double trapezoid = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        trapezoid += (x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2.0;
    }
    InterpolationTest::ExpectScaledNear(L.Integral(0.0, 3.0), trapezoid);

    // Akima is exact for this data too, since its pieces are cubic and the
    // endpoint slopes are recovered from the data.
    const AkimaSpline a{x, y};
    const auto A = Primitive(a);
    EXPECT_NEAR(A.Integral(0.0, 3.0), 81.0 / 4.0, 1.0);
}

TEST(FunctionAlgebra, PolynomialTakesTheClosedFormPrimitivePath) {
    using namespace Interpolation;
    const Polynomial<double> p{0.0, 0.0, 0.0, 1.0}; // x^3
    const auto P = Primitive(p);
    static_assert(!PiecewisePolynomial1D<Polynomial<double>>);
    static_assert(Antidifferentiable1D<Polynomial<double>>);

    InterpolationTest::ExpectScaledNear(P.Integral(0.0, 3.0), 81.0 / 4.0);
    InterpolationTest::ExpectScaledNear(Derivative(P)(2.0), p(2.0));
}

TEST(FunctionAlgebra, ComplexOrdinatesFlowThroughTheAlgebra) {
    using namespace Interpolation;
    using Complex = std::complex<double>;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const std::vector<Complex> y{
        {0.0, 0.0}, {1.0, 1.0}, {8.0, 2.0}, {27.0, 3.0}};

    const CubicSpline s{x, y};
    const auto g = Derivative(s) * s + Complex{2.0, 1.0};

    for (const double query : {0.5, 1.5, 2.5}) {
        SCOPED_TRACE(query);
        const auto expected =
            s.Evaluate<1>(query) * s(query) + Complex{2.0, 1.0};
        InterpolationTest::ExpectScaledNear(g(query), expected);
    }

    // A real scalar combined with a complex function stays complex.
    const auto scaled = 2.0 * s;
    static_assert(std::is_same_v<typename decltype(scaled)::Scalar, Complex>);
    InterpolationTest::ExpectScaledNear(scaled(1.5), 2.0 * s(1.5));
}

TEST(FunctionAlgebra, QuotientDerivativesFollowTheReciprocalRecurrence) {
    using namespace Interpolation;
    // q(x) = 1/(1+x), whose nth derivative is (-1)^n n! / (1+x)^(n+1).
    const Polynomial<double> numerator{1.0};
    const Polynomial<double> denominator{1.0, 1.0};
    const auto q = numerator / denominator;

    const double x = 1.0;
    const auto exact = [x](int n) {
        double factorial = 1.0;
        for (int i = 2; i <= n; ++i) {
            factorial *= i;
        }
        return (n % 2 ? -1.0 : 1.0) * factorial / std::pow(1.0 + x, n + 1);
    };

    InterpolationTest::ExpectScaledNear(q.Evaluate<0>(x), exact(0));
    InterpolationTest::ExpectScaledNear(q.Evaluate<1>(x), exact(1));
    InterpolationTest::ExpectScaledNear(q.Evaluate<2>(x), exact(2));
    InterpolationTest::ExpectScaledNear(q.Evaluate<3>(x), exact(3));
    InterpolationTest::ExpectScaledNear(q.Evaluate<4>(x), exact(4));
    InterpolationTest::ExpectScaledNear(q.Evaluate<5>(x), exact(5));
}

TEST(FunctionAlgebra, CompositionDerivativesFollowFaaDiBruno) {
    using namespace Interpolation;

    // An affine inner function: h(x) = (2x + 1)^4.
    {
        const Polynomial<double> outer{0.0, 0.0, 0.0, 0.0, 1.0}; // u^4
        const Polynomial<double> inner{1.0, 2.0};                // 2x + 1
        const auto h = Compose(outer, inner);
        const double x = 1.0;
        const double u = 3.0;

        InterpolationTest::ExpectScaledNear(h.Evaluate<0>(x), std::pow(u, 4));
        InterpolationTest::ExpectScaledNear(h.Evaluate<1>(x),
                                            8.0 * std::pow(u, 3));
        InterpolationTest::ExpectScaledNear(h.Evaluate<2>(x), 48.0 * u * u);
        InterpolationTest::ExpectScaledNear(h.Evaluate<3>(x), 192.0 * u);
        InterpolationTest::ExpectScaledNear(h.Evaluate<4>(x), 384.0);
        InterpolationTest::ExpectScaledNear(h.Evaluate<5>(x), 0.0);
    }

    // A nonlinear inner function, where the chain rule alone is not enough:
    // h(x) = (x^2 + 1)^3 = x^6 + 3x^4 + 3x^2 + 1.
    {
        const Polynomial<double> outer{0.0, 0.0, 0.0, 1.0}; // u^3
        const Polynomial<double> inner{1.0, 0.0, 1.0};      // x^2 + 1
        const auto h = Compose(outer, inner);
        const Polynomial<double> expanded{1.0, 0.0, 3.0, 0.0, 3.0, 0.0, 1.0};

        for (const double x : {0.4, 1.0, 1.7}) {
            SCOPED_TRACE(x);
            InterpolationTest::ExpectScaledNear(h.Evaluate<0>(x),
                                                expanded.Evaluate<0>(x));
            InterpolationTest::ExpectScaledNear(h.Evaluate<1>(x),
                                                expanded.Evaluate<1>(x));
            InterpolationTest::ExpectScaledNear(h.Evaluate<2>(x),
                                                expanded.Evaluate<2>(x));
            InterpolationTest::ExpectScaledNear(h.Evaluate<3>(x),
                                                expanded.Evaluate<3>(x));
            InterpolationTest::ExpectScaledNear(h.Evaluate<4>(x),
                                                expanded.Evaluate<4>(x));
            InterpolationTest::ExpectScaledNear(h.Evaluate<5>(x),
                                                expanded.Evaluate<5>(x));
            InterpolationTest::ExpectScaledNear(h.Evaluate<6>(x),
                                                expanded.Evaluate<6>(x));
        }
    }
}

TEST(FunctionAlgebra, HigherDerivativesComposeAndStayAllocationFree) {
    using namespace Interpolation;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const auto y = Cubes(x);
    const CubicSpline s{x, y};

    // Differentiating a quotient of interpolants twice now compiles, where
    // before it was a static_assert.
    const auto ratio = s / (s + 10.0);
    const auto curvature = Derivative<2>(ratio);

    // Check it against a central difference of the first derivative.
    const double at = 1.5;
    const double h = 1.0e-5;
    const double byDifference =
        (ratio.Evaluate<1>(at + h) - ratio.Evaluate<1>(at - h)) / (2 * h);
    EXPECT_NEAR(curvature(at), byDifference, 1.0e-4);

    if constexpr (InterpolationTest::AllocationCountingEnabled()) {
        volatile double sink = 0.0;
        const auto before = InterpolationTest::AllocationCount();
        for (int i = 0; i < 500; ++i) {
            sink = sink + ratio.Evaluate<4>(0.5 + 0.001 * i);
        }
        const auto after = InterpolationTest::AllocationCount();
        EXPECT_EQ(after, before) << "high-order evaluation allocated "
                                 << (after - before) << " times";
    }
}

TEST(AnyFunction1D, HoldsDifferentKindsBehindOneType) {
    using namespace Interpolation;
    using Any = AnyFunction1D<double, double>;
    static_assert(Function1D<Any>);

    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const auto y = Cubes(x);

    const CubicSpline spline{x, y};
    const Polynomial<double> polynomial{1.0, 2.0, 3.0};
    const Linear linear{x, y};

    std::vector<Any> mixed;
    mixed.emplace_back(spline);
    mixed.emplace_back(polynomial);
    mixed.emplace_back(linear);

    InterpolationTest::ExpectScaledNear(mixed[0](1.5), spline(1.5));
    InterpolationTest::ExpectScaledNear(mixed[1](1.5), polynomial(1.5));
    InterpolationTest::ExpectScaledNear(mixed[2](1.5), linear(1.5));

    InterpolationTest::ExpectScaledNear(mixed[0].Evaluate<1>(1.5),
                                        spline.Evaluate<1>(1.5));
    InterpolationTest::ExpectScaledNear(mixed[1].Evaluate<2>(1.5),
                                        polynomial.Evaluate<2>(1.5));
}

TEST(AnyFunction1D, IsCopyableEvenWhenWhatItHoldsIsNot) {
    using namespace Interpolation;
    using Any = AnyFunction1D<double, double>;

    // An interpolator built from rvalues owns its samples and is move-only,
    // so a clone-on-copy erasure could not hold it. Sharing immutable state
    // can, and the copy is independent as far as the interface allows.
    Any owning{CubicSpline{std::vector<double>{0.0, 1.0, 2.0, 3.0},
                           std::vector<double>{0.0, 1.0, 8.0, 27.0}}};
    Any copy = owning;
    Any assigned{Polynomial<double>{1.0}};
    assigned = copy;

    InterpolationTest::ExpectScaledNear(copy(1.5), owning(1.5));
    InterpolationTest::ExpectScaledNear(assigned(1.5), owning(1.5));

    // It therefore enters the algebra by value, with no Ref needed.
    const auto g = Derivative(copy) * copy + 2.0;
    InterpolationTest::ExpectScaledNear(
        g(1.5), owning.Evaluate<1>(1.5) * owning(1.5) + 2.0);
}

TEST(AnyFunction1D, GivesPiecewiseMixedPieceKinds) {
    using namespace Interpolation;
    using Any = AnyFunction1D<double, double>;

    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const auto y = Cubes(x);

    // Piecewise knows nothing about erasure; it simply holds one type, which
    // here happens to be one that can hold anything.
    std::vector<Any> pieces;
    pieces.emplace_back(Linear{x, y});
    pieces.emplace_back(Polynomial<double>{100.0, 1.0});

    const Piecewise<Any> layered{{0.0, 3.0, 6.0}, std::move(pieces)};
    EXPECT_EQ(layered.PieceCount(), 2u);

    const Linear reference{x, y};
    InterpolationTest::ExpectScaledNear(layered(1.5), reference(1.5));
    InterpolationTest::ExpectScaledNear(layered(4.0), 104.0);

    const auto [below, above] = layered.Limits(3.0);
    InterpolationTest::ExpectScaledNear(below, 27.0);
    InterpolationTest::ExpectScaledNear(above, 103.0);

    // An extracted layer still works on its own.
    const auto outer = layered.Piece(1);
    InterpolationTest::ExpectScaledNear(outer.Lower(), 3.0);
    InterpolationTest::ExpectScaledNear(outer(5.0), 105.0);
}

TEST(AnyFunction1D, HonoursItsMaxOrder) {
    using namespace Interpolation;
    // Two orders is enough for a value, a slope and a curvature.
    using Any = AnyFunction1D<double, double, 2>;
    static_assert(Any::HighestOrder == 2);

    const Any p{Polynomial<double>{1.0, 2.0, 3.0}};
    const Polynomial<double> reference{1.0, 2.0, 3.0};
    InterpolationTest::ExpectScaledNear(p.Evaluate<2>(1.5),
                                        reference.Evaluate<2>(1.5));
    // p.Evaluate<3>(1.5) is a compile error rather than a wrong answer.
}
