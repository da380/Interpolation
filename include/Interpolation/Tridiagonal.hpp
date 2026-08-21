#ifndef INTERPOLATION_TRIDIAGONAL_HPP
#define INTERPOLATION_TRIDIAGONAL_HPP

#include <cassert>
#include <cstddef>
#include <span>

namespace Interpolation::Detail {

/**
 * @brief Solve a tridiagonal system in place by the Thomas algorithm.
 *
 * The coefficients and the right-hand side are typed separately so that a real
 * matrix can act on a complex right-hand side without instantiating a complex
 * matrix. The three diagonals are indexed against the full matrix:
 *
 * - `sub[i]` is `A(i, i - 1)` for `i` in `[1, n)`; `sub[0]` is not read.
 * - `diag[i]` is `A(i, i)`.
 * - `super[i]` is `A(i, i + 1)` for `i` in `[0, n - 1)`; `super[n - 1]` is not
 *   read.
 *
 * No pivoting is performed. That is safe for the spline systems this solves,
 * which are strictly diagonally dominant: an interior row has off-diagonal
 * magnitude `(h[i - 1] + h[i]) / 6` against a diagonal of
 * `(h[i - 1] + h[i]) / 3`, and every boundary row is dominant by the same
 * factor of two.
 *
 * @param sub Sub-diagonal, length `n`.
 * @param diag Diagonal, length `n`. Overwritten.
 * @param super Super-diagonal, length `n`.
 * @param rhs Right-hand side, length `n`. Overwritten with the solution.
 * @pre All four spans have the same non-zero length.
 */
template <typename Real, typename Scalar>
void
SolveTridiagonal(std::span<const Real> sub, std::span<Real> diag,
                 std::span<const Real> super, std::span<Scalar> rhs) {
    const auto n = diag.size();
    assert(n > 0);
    assert(sub.size() == n && super.size() == n && rhs.size() == n);

    // Forward elimination.
    for (std::size_t i = 1; i < n; ++i) {
        const auto factor = sub[i] / diag[i - 1];
        diag[i] -= factor * super[i - 1];
        rhs[i] -= factor * rhs[i - 1];
    }

    // Back substitution.
    rhs[n - 1] /= diag[n - 1];
    for (std::size_t i = n - 1; i-- > 0;) {
        rhs[i] = (rhs[i] - super[i] * rhs[i + 1]) / diag[i];
    }
}

} // namespace Interpolation::Detail

#endif // INTERPOLATION_TRIDIAGONAL_HPP
