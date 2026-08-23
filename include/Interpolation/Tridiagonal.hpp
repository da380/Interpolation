#ifndef INTERPOLATION_TRIDIAGONAL_HPP
#define INTERPOLATION_TRIDIAGONAL_HPP

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <Interpolation/Concepts.hpp>

namespace Interpolation {

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
 * This is the one-shot form: `diag` is consumed by the elimination, so nothing
 * survives that a second right-hand side could reuse. For many right-hand
 * sides on one matrix use TridiagonalFactorization, which keeps the pivots.
 *
 * @param sub Sub-diagonal, length `n`.
 * @param diag Diagonal, length `n`. Overwritten.
 * @param super Super-diagonal, length `n`.
 * @param rhs Right-hand side, length `n`. Overwritten with the solution.
 * @throws std::invalid_argument if the spans are empty or differ in length.
 */
template <typename Real, typename Scalar>
void
SolveTridiagonal(std::span<const Real> sub, std::span<Real> diag,
                 std::span<const Real> super, std::span<Scalar> rhs) {
    const auto n = diag.size();
    if (n == 0) {
        throw std::invalid_argument("SolveTridiagonal: the system is empty");
    }
    if (sub.size() != n || super.size() != n || rhs.size() != n) {
        throw std::invalid_argument(
            "SolveTridiagonal: the diagonals and the right-hand side must all "
            "have the same length");
    }

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

/**
 * @brief A tridiagonal matrix reduced once, ready to act on many right-hand
 * sides.
 *
 * The elimination in SolveTridiagonal destroys the diagonal it works on, so
 * solving `k` right-hand sides on one matrix costs `k` eliminations. That is
 * the right trade when the matrix is built for a single solve, and the wrong
 * one whenever a matrix outlives its right-hand side — a spline system on a
 * fixed grid applied to many datasets, a tensor-product grid applied line by
 * line, an operator applied every iteration of a matrix-free solve.
 *
 * This class does the elimination once at construction and keeps what it
 * produces. Each Solve is then a forward sweep and a back substitution, and
 * because the reciprocals of the pivots are stored rather than the pivots
 * themselves, it performs no division at all.
 *
 * The matrix is real and the right-hand side may be real or complex, so
 * `Scalar` is a parameter of Solve rather than of the class: one factorisation
 * serves both.
 *
 * @code
 * const auto factorization =
 *     TridiagonalFactorization<double>{sub, diag, super};
 * for (auto& line : lines) {
 *     factorization.Solve(std::span{line});  // in place
 * }
 * @endcode
 *
 * No pivoting is performed, for the reason given on SolveTridiagonal. A zero
 * pivot is reported at construction rather than left to produce infinities
 * during a solve.
 *
 * @tparam R Real precision of the matrix entries.
 */
template <typename R>
    requires Real<R>
class TridiagonalFactorization {
  public:
    /** @brief Precision of the matrix entries. */
    using Real = R;

    /**
     * @brief Reduce a tridiagonal matrix.
     *
     * The diagonals are indexed against the full matrix, exactly as in
     * SolveTridiagonal: `sub[0]` and `super[n - 1]` are not read.
     *
     * @param sub Sub-diagonal, length `n`.
     * @param diag Diagonal, length `n`.
     * @param super Super-diagonal, length `n`.
     * @throws std::invalid_argument if the spans are empty, differ in length,
     *         or the elimination meets a zero pivot.
     */
    TridiagonalFactorization(std::span<const R> sub, std::span<const R> diag,
                             std::span<const R> super) {
        const auto n = diag.size();
        if (n == 0) {
            throw std::invalid_argument(
                "TridiagonalFactorization: the system is empty");
        }
        if (sub.size() != n || super.size() != n) {
            throw std::invalid_argument("TridiagonalFactorization: the three "
                                        "diagonals must have the same length");
        }

        _multiplier.assign(n, R{});
        _inversePivot.assign(n, R{});
        _superOverPivot.assign(n, R{});

        auto pivot = diag[0];
        CheckPivot(pivot, 0);
        _inversePivot[0] = R{1} / pivot;

        for (std::size_t i = 1; i < n; ++i) {
            // The multiplier is the entry the elimination would have written
            // into L; keeping it is what lets a later right-hand side skip
            // the elimination entirely.
            _multiplier[i] = sub[i] * _inversePivot[i - 1];
            pivot = diag[i] - _multiplier[i] * super[i - 1];
            CheckPivot(pivot, i);
            _inversePivot[i] = R{1} / pivot;
            _superOverPivot[i - 1] = super[i - 1] * _inversePivot[i - 1];
        }
    }

    /** @brief Order of the system. */
    std::size_t Size() const { return _inversePivot.size(); }

    /**
     * @brief Solve for one right-hand side, in place.
     *
     * @tparam Scalar Right-hand side type, real or complex.
     * @param rhs Right-hand side, length `Size()`. Overwritten with the
     *        solution.
     * @throws std::invalid_argument if `rhs` has the wrong length.
     */
    template <typename Scalar> void Solve(std::span<Scalar> rhs) const {
        const auto n = Size();
        if (rhs.size() != n) {
            throw std::invalid_argument(
                "TridiagonalFactorization::Solve: the right-hand side has "
                "length " +
                std::to_string(rhs.size()) + " but the system has order " +
                std::to_string(n));
        }

        for (std::size_t i = 1; i < n; ++i) {
            rhs[i] -= _multiplier[i] * rhs[i - 1];
        }

        rhs[n - 1] *= _inversePivot[n - 1];
        for (std::size_t i = n - 1; i-- > 0;) {
            rhs[i] =
                rhs[i] * _inversePivot[i] - _superOverPivot[i] * rhs[i + 1];
        }
    }

  private:
    std::vector<R> _multiplier;     // sub[i] / pivot[i - 1]; [0] unused.
    std::vector<R> _inversePivot;   // 1 / pivot[i].
    std::vector<R> _superOverPivot; // super[i] / pivot[i]; [n - 1] unused.

    static void CheckPivot(R pivot, std::size_t row) {
        if (pivot == R{}) {
            throw std::invalid_argument(
                "TridiagonalFactorization: the elimination met a zero pivot "
                "at row " +
                std::to_string(row) +
                "; the matrix is singular or needs pivoting");
        }
    }
};

/// Deduce the precision from the diagonals.
template <typename R>
TridiagonalFactorization(std::span<const R>, std::span<const R>,
                         std::span<const R>) -> TridiagonalFactorization<R>;

} // namespace Interpolation

#endif // INTERPOLATION_TRIDIAGONAL_HPP
