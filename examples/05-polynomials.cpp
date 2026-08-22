// 05 - Polynomials
//
// Polynomial is a value type with arithmetic, calculus and the same
// Evaluate<N> interface as the interpolators, so it stands in wherever one of
// them can.

#include <Interpolation/Polynomial.hpp>

#include <complex>
#include <cstdio>
#include <iostream>
#include <vector>

int
main() {
    using Interpolation::Polynomial;

    // Coefficients are given in ascending powers: 1 + 2x + 3x^2.
    const Polynomial<double> p{1.0, 2.0, 3.0};
    const Polynomial<double> q{0.0, 1.0};

    std::cout << "p        = " << p << " (degree " << p.Degree() << ")\n";
    std::cout << "q        = " << q << " (degree " << q.Degree() << ")\n";
    std::cout << "p + q    = " << p + q << "\n";
    std::cout << "p * q    = " << p * q << "\n";
    std::cout << "p - p    = " << p - p << "\n";

    // Degree stays honest after cancellation: trailing zeros are trimmed, so
    // a subtraction that removes the leading term reports the lower degree
    // rather than the length of the coefficient array.
    const auto cancelled = p - Polynomial<double>{0.0, 0.0, 3.0};
    std::cout << "p - 3x^2 = " << cancelled << " (degree " << cancelled.Degree()
              << ")\n";

    // Mixing types promotes rather than truncating: a complex scalar added to
    // a real polynomial gives a complex polynomial, instead of quietly
    // dropping the imaginary part.
    const auto promoted = p + std::complex<double>{0.0, 1.0};
    std::cout << "p + i    = " << promoted << "\n";

    std::printf("\nCalculus at x = 2:\n");
    std::printf("  p            = %10.4f\n", p(2.0));
    std::printf("  dp/dx        = %10.4f\n", p.Evaluate<1>(2.0));
    std::printf("  d2p/dx2      = %10.4f\n", p.Evaluate<2>(2.0));
    std::printf("  d3p/dx3      = %10.4f  (a quadratic has none)\n",
                p.Evaluate<3>(2.0));
    std::printf("  antiderivative %10.4f  (vanishing at zero)\n",
                p.Antiderivative(2.0));

    // Random polynomials take an explicit seed. Anything randomised should
    // report the seed it used, so that a failure can be replayed exactly
    // rather than being irreproducible.
    const std::uint64_t seed = 20260821;
    const auto r = Polynomial<double>::Random(3, seed);
    std::cout << "\nRandom(3, " << seed << ") = " << r << "\n";
    std::cout << "same seed again      = "
              << Polynomial<double>::Random(3, seed) << "\n";

    return 0;
}
