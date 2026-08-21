// 09 - Piecewise-continuous functions
//
// The interpolators require strictly increasing abscissae, so data with a
// genuine jump has nowhere to go. Piecewise is where the discontinuity is
// represented instead: as structure, not as a coincidence in the samples.

#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Function.hpp>
#include <Interpolation/Linear.hpp>
#include <Interpolation/Piecewise.hpp>
#include <Interpolation/Polynomial.hpp>

#include <cstdio>
#include <utility>
#include <vector>

// Stands in for a solver: it needs something it can evaluate and the interval
// to work over, which is exactly what an extracted piece carries.
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

int
main() {
    using namespace Interpolation;

    // A layered profile. The abscissa 3.0 appears twice, which is a common
    // convention in tabulated data for marking an interface, and the value
    // differs across it.
    const std::vector<double> radius{0.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0};
    const std::vector<double> density{8.0,  8.2,  8.6,  9.0,
                                      12.0, 12.4, 12.9, 13.5};

    auto model = SplitAtRepeats(radius, density, [](auto r, auto d) {
        return CubicSpline{std::move(r), std::move(d)};
    });

    std::printf("%zu layers over [%.1f, %.1f]\n", model.PieceCount(),
                model.Lower(), model.Upper());
    for (std::size_t k = 0; k < model.PieceCount(); ++k) {
        const auto layer = model.Piece(k);
        std::printf("  layer %zu on [%.1f, %.1f]\n", k, layer.Lower(),
                    layer.Upper());
    }

    // At an interface both one-sided values are usually wanted; a single
    // number silently picks one of them.
    const auto [below, above] = model.Limits(3.0);
    std::printf("\nAt the interface r = 3:\n");
    std::printf("  from below %.4f, from above %.4f, jump %.4f\n", below, above,
                above - below);
    std::printf("  plain evaluation is right-continuous: %.4f\n", model(3.0));
    std::printf("  and from the left explicitly:         %.4f\n",
                model.EvaluateFrom<0>(3.0, Side::Left));

    // Extracting a layer is usually more useful than evaluating through the
    // whole object: it can be handed to something that works over exactly
    // that interval, and it carries the interval with it.
    const auto outer = model.Piece(1);
    std::printf("\nThe outer layer on its own:\n");
    std::printf("  midpoint integral   %.6f\n", MidpointIntegral(outer, 20000));
    std::printf("  analytic integral   %.6f\n",
                Primitive(outer).Integral(outer.Lower(), outer.Upper()));

    // The whole model is a function too, so it composes. Its pieces own their
    // samples and are therefore move-only, so it enters the algebra by
    // reference, spelled explicitly.
    const auto gradient = Derivative(Ref(model));
    std::printf("\n  d(density)/dr just above the interface: %.6f\n",
                gradient(3.0001));

    // Layers need not be the same kind of function. AnyFunction1D erases the
    // type, and Piecewise holds it without knowing anything about erasure.
    using Any = AnyFunction1D<double, double>;
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const std::vector<double> y{0.0, 1.0, 8.0, 27.0};
    std::vector<Any> mixed;
    mixed.emplace_back(Linear{x, y});
    mixed.emplace_back(Polynomial<double>{100.0, 1.0});
    const Piecewise<Any> hybrid{{0.0, 3.0, 6.0}, std::move(mixed)};

    std::printf("\nMixed layer kinds, linear then polynomial:\n");
    std::printf("  at 1.5 %.4f (linear), at 4.0 %.4f (polynomial)\n",
                hybrid(1.5), hybrid(4.0));

    return 0;
}
