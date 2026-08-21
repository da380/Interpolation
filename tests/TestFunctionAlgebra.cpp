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
