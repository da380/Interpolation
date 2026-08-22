// 07 - The function algebra
//
// Everything here models one concept, Function1D: it names an abscissa type
// and a value type, and answers Evaluate<N>. Interpolators do, polynomials
// do, and so does every expression built from them. That is what lets them be
// combined arithmetically instead of only evaluated.

#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Function.hpp>
#include <Interpolation/Polynomial.hpp>

#include <cstdio>
#include <vector>

int
main() {
    using namespace Interpolation;

    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const std::vector<double> y{0.0, 1.0, 8.0, 27.0, 64.0};
    const CubicSpline s{x, y};

    static_assert(Function1D<decltype(s)>);

    // Differentiation produces a function, not a number, so it composes.
    const auto g = Derivative(s) * s + 2.0;
    const auto dg = Derivative(g);

    static_assert(Function1D<decltype(g)>);

    std::printf("g = f' * f + 2, built from a spline\n\n");
    std::printf("      x        g(x)      g'(x)\n");
    for (const double q : {0.5, 1.5, 2.5, 3.5}) {
        std::printf("  %6.3f  %10.5f  %10.5f\n", q, g(q), dg(q));
    }

    // The nodes hold their operands by value, so an expression built entirely
    // from temporaries is fine. Holding references would dangle here, since
    // every subexpression is a temporary.
    const auto fromTemporaries =
        CubicSpline{x, y} * CubicSpline{x, y} + Polynomial<double>{1.0};
    std::printf("\n  Built from temporaries: %.5f at x = 1.5\n",
                fromTemporaries(1.5));

    // Compose applies one function to the result of another.
    const auto halved = Compose(s, Identity<double>{} * 0.5);
    std::printf("\n  Compose(s, x/2) at 3.0 = %.6f, and s(1.5) = %.6f\n",
                halved(3.0), s(1.5));

    // Derivatives of a product follow the general Leibniz rule to any order,
    // expanded at compile time. Check the second against the hand-written
    // form f''g + 2f'g' + fg''.
    const auto product = s * s;
    const double at = 1.5;
    const double byHand = 2.0 * (s.Evaluate<2>(at) * s(at) +
                                 s.Evaluate<1>(at) * s.Evaluate<1>(at));
    std::printf("\n  (s*s)'' at 1.5: algebra %.6f, by hand %.6f\n",
                product.Evaluate<2>(at), byHand);

    // Quotients and compositions differentiate to any order too.
    const Polynomial<double> one{1.0};
    const Polynomial<double> linear{1.0, 1.0};
    const auto reciprocal = one / linear; // 1 / (1 + x)
    std::printf("\n  1/(1+x) at x = 1, whose nth derivative is "
                "(-1)^n n!/(1+x)^(n+1):\n");
    std::printf("    computed  %9.6f %9.6f %9.6f %9.6f\n",
                reciprocal.Evaluate<0>(1.0), reciprocal.Evaluate<1>(1.0),
                reciprocal.Evaluate<2>(1.0), reciprocal.Evaluate<3>(1.0));
    std::printf("    exact      0.500000 -0.250000  0.250000 -0.375000\n");

    return 0;
}
