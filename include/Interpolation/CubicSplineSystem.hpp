#ifndef INTERPOLATION_CUBIC_SPLINE_SYSTEM_HPP
#define INTERPOLATION_CUBIC_SPLINE_SYSTEM_HPP

#include <cstddef>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Samples.hpp>
#include <Interpolation/Side.hpp>
#include <Interpolation/Tridiagonal.hpp>

namespace Interpolation {

/** @brief Endpoint conditions for CubicSpline and CubicSplineSystem. */
enum class BoundaryCondition {
    /** The endpoint second derivative is zero. */
    Natural,
    /** The endpoint first derivative is supplied by the caller. */
    Clamped,
    /**
     * @brief The third derivative is continuous across the first and last
     * interior knots.
     *
     * Imposes nothing false at the boundary, so a cubic is reproduced exactly
     * and the scheme stays fourth order right up to the ends, where Natural
     * costs an order. It is a condition on the whole system rather than on one
     * endpoint, so it must be used at both ends and needs at least four nodes.
     */
    NotAKnot,
    /**
     * @brief The spline and its first two derivatives are continuous across
     * the seam where the final node meets the first.
     *
     * The final node closes the period: it lies one period beyond the first,
     * and its ordinate must equal the first's exactly, so the grid a caller
     * already holds needs no special case. Queries are reduced into
     * `[x_0, x_n)`. Like NotAKnot it constrains the whole system, so it must
     * be used at both ends, and it needs at least four nodes counting the
     * closing one.
     */
    Periodic
};

namespace Detail {

/**
 * @brief Check that two endpoint conditions may be paired, and return the
 * fewest nodes they need.
 *
 * Shared by every class that takes a pair, so the rule and its message are
 * written once.
 *
 * @param left Condition at the first node.
 * @param right Condition at the final node.
 * @param what Class name, used in the exception message.
 * @throws std::invalid_argument if NotAKnot or Periodic is used at only one
 *         end.
 */
inline std::size_t
SplineConditionMinimum(BoundaryCondition left, BoundaryCondition right,
                       const char *what) {
    const auto wholeSystem = [&](BoundaryCondition condition,
                                 const char *name) {
        if ((left == condition) != (right == condition)) {
            throw std::invalid_argument(
                std::string{what} + ": " + name +
                " constrains the whole system rather than one endpoint, so "
                "it must be used at both ends or neither");
        }
    };
    wholeSystem(BoundaryCondition::NotAKnot, "NotAKnot");
    wholeSystem(BoundaryCondition::Periodic, "Periodic");

    return left == BoundaryCondition::NotAKnot ||
                   left == BoundaryCondition::Periodic
               ? 4
               : 2;
}

/**
 * @brief Check that the final ordinate of a periodic set equals the first.
 *
 * Exact equality is asked for, not closeness: the closing ordinate is meant
 * to be a copy, and one that differs by rounding is a sign it was computed
 * separately, which is worth hearing about.
 *
 * @throws std::invalid_argument if the two differ.
 */
template <typename YRange>
void
CheckPeriodicClosure(const YRange &y, const char *what) {
    const auto n = static_cast<std::size_t>(std::ranges::size(y));
    if (y[n - 1] != y[0]) {
        throw std::invalid_argument(
            std::string{what} +
            ": a periodic spline's final ordinate closes the period, so it "
            "must equal the first exactly");
    }
}

} // namespace Detail

/**
 * @brief The cubic-spline linear system for a fixed set of nodes, factorised
 * once and reusable across any number of ordinate sets.
 *
 * A cubic spline is determined by a tridiagonal system whose **matrix depends
 * on the nodes alone**; only the right-hand side carries the ordinates. So a
 * caller holding many datasets on one grid — the components of a vector field,
 * an ensemble of profiles, a field of spherical-harmonic coefficients sampled
 * along a radial line — can factorise once and solve many times. That case is
 * common enough to deserve a name, and naming it is what this class is for.
 *
 * It is also what CubicSpline is built on, so there is one assembly of the
 * spline equations in the library rather than one per caller.
 *
 * @code
 * const auto system = CubicSplineSystem{radii};       // factorised here
 * std::vector<std::complex<double>> curvature(radii.size());
 * for (const auto& line : coefficients) {
 *     system.Solve(line, std::span{curvature});       // no allocation
 *     // ... use the nodal second derivatives
 * }
 * @endcode
 *
 * Solve allocates nothing and is `const`, so a fixed system may be shared
 * across threads. The matrix is real even when the ordinates are complex, so
 * the ordinate type is a parameter of Solve rather than of the class: one
 * system serves both.
 *
 * ### Which system
 *
 * Natural and Clamped are assembled in the **nodal-curvature** formulation and
 * solved directly for what the evaluator wants.
 *
 * NotAKnot is assembled in the **nodal-slope** formulation and converted
 * afterwards. Written directly in curvatures its boundary row reaches outside
 * the tridiagonal band, and eliminating that entry produces a leading
 * coefficient of `h0^2 - h1^2`, which vanishes on a uniform grid. In slopes the
 * same condition gives the boundary row
 *
 * @f[
 * h_1 d_0 + (h_0 + h_1) d_1 =
 *   \frac{h_1 (2h_1 + 3h_0)\delta_0 + h_0^2 \delta_1}{h_0 + h_1},
 * @f]
 *
 * whose diagonal is @f$h_1 > 0@f$ for any spacing. The conversion back to
 * curvatures happens in place in the output, so it costs no scratch buffer.
 *
 * No pivoting is used. The interior rows are diagonally dominant by a factor
 * of two; the two not-a-knot boundary rows are not, but the first elimination
 * step turns the second pivot into @f$h_0 + h_1@f$, and the pivots stay
 * positive.
 *
 * Periodic is assembled in curvatures with every row interior, indices taken
 * modulo the @f$m = n - 1@f$ distinct nodes. That system is cyclic: it is
 * tridiagonal but for one entry in each far corner, both
 * @f$c = h_{m-1}/6@f$. It is solved by Sherman-Morrison. With
 * @f$\gamma = -A_{00}@f$, folding the corners into the first and last
 * diagonal entries,
 *
 * @f[
 * T = A - u v^T, \quad u = (\gamma, 0, \ldots, 0, c)^T, \quad
 * v = (1, 0, \ldots, 0, c/\gamma)^T,
 * @f]
 *
 * leaves a tridiagonal @f$T@f$ that is still strictly diagonally dominant.
 * It is factorised once, and the one auxiliary solution
 * @f$z = T^{-1} u@f$ is kept with it, so each Solve is a single tridiagonal
 * solve and a rank-one correction, still `const` and still free of
 * allocation. The correction divides by @f$1 + v^T z = \det A / \det T@f$,
 * which is positive since both matrices are.
 *
 * Construction takes a range of nodes, borrowing an lvalue and owning an
 * rvalue, exactly as the interpolators do.
 *
 * @tparam XView View over the real abscissae.
 */
template <typename XView>
    requires RealRange<XView> && std::ranges::view<XView>
class CubicSplineSystem {
  public:
    /** @brief Abscissa precision. */
    using Real = std::ranges::range_value_t<XView>;

    /**
     * @brief Assemble and factorise the system.
     * @param x Strictly increasing abscissae; at least two, or four when
     *        NotAKnot or Periodic is used. For Periodic the final node closes
     *        the period.
     * @param left Condition at the first node.
     * @param right Condition at the final node.
     * @throws std::invalid_argument if the abscissae are not strictly
     *         increasing, there are too few of them, or NotAKnot or Periodic
     *         is asked for at only one end.
     */
    CubicSplineSystem(XView x, BoundaryCondition left, BoundaryCondition right)
        : _x{std::move(x)}, _left{left}, _right{right},
          _notAKnot{left == BoundaryCondition::NotAKnot},
          _periodic{left == BoundaryCondition::Periodic},
          _factorization{Assemble(_x, left, right, _spacing)} {
        if (_periodic) {
            PrepareCorrection();
        }
    }

    /** @brief Assemble with the natural condition at both ends. */
    explicit CubicSplineSystem(XView x)
        : CubicSplineSystem{std::move(x), BoundaryCondition::Natural,
                            BoundaryCondition::Natural} {}

    /** @brief Assemble with the same condition type at both ends. */
    CubicSplineSystem(XView x, BoundaryCondition both)
        : CubicSplineSystem{std::move(x), both, both} {}

    /** @brief Number of nodes. */
    std::size_t Size() const {
        return static_cast<std::size_t>(std::ranges::size(_x));
    }

    /** @brief Abscissa of node `i`. */
    Real Node(std::size_t i) const { return _x[i]; }

    /** @brief The nodes themselves. */
    const XView &Nodes() const { return _x; }

    /** @brief Width of segment `i`, from node `i` to node `i + 1`. */
    Real Spacing(std::size_t i) const { return _spacing[i]; }

    /** @brief Condition imposed at the first node. */
    BoundaryCondition LeftCondition() const { return _left; }

    /** @brief Condition imposed at the final node. */
    BoundaryCondition RightCondition() const { return _right; }

    /**
     * @brief The period, `x_n - x_0`, when the condition is Periodic; empty
     * otherwise.
     */
    std::optional<Real> Period() const {
        if (!_periodic) {
            return std::nullopt;
        }
        return _x[Size() - 1] - _x[0];
    }

    /**
     * @brief Solve for the nodal second derivatives of one ordinate set.
     *
     * Neither endpoint may be Clamped, since no endpoint derivative is given;
     * use the four-argument overload for that.
     *
     * @param y Ordinates, one per node.
     * @param curvature Output, one per node. Overwritten.
     * @throws std::invalid_argument if either length is wrong, or an endpoint
     *         is Clamped.
     */
    template <typename YRange, typename Scalar>
        requires InterpolationRanges<XView, YRange>
    void Solve(const YRange &y, std::span<Scalar> curvature) const {
        if (_left == BoundaryCondition::Clamped ||
            _right == BoundaryCondition::Clamped) {
            throw std::invalid_argument(
                "CubicSplineSystem::Solve: a Clamped endpoint needs its "
                "prescribed derivative, so use the overload that takes one");
        }
        Solve(y, Scalar{}, Scalar{}, curvature);
    }

    /**
     * @brief Solve for the nodal second derivatives, giving the endpoint
     * derivatives a Clamped condition needs.
     *
     * The derivatives are ignored at an endpoint whose condition is not
     * Clamped, so a caller with a mixed pair need not care which is which.
     *
     * `y` and `curvature` must not overlap: the right-hand side is read as
     * the solution is written.
     *
     * @param y Ordinates, one per node.
     * @param leftDerivative First derivative at the first node.
     * @param rightDerivative First derivative at the final node.
     * @param curvature Output, one per node. Overwritten.
     * @throws std::invalid_argument if either length is wrong, or the system
     *         is Periodic and the final ordinate differs from the first.
     */
    template <typename YRange, typename Scalar>
        requires InterpolationRanges<XView, YRange>
    void Solve(const YRange &y, Scalar leftDerivative, Scalar rightDerivative,
               std::span<Scalar> curvature) const {
        const auto n = Size();
        CheckLength(static_cast<std::size_t>(std::ranges::size(y)), n,
                    "ordinate range");
        CheckLength(curvature.size(), n, "curvature output");

        if (_notAKnot) {
            SolveNotAKnot(y, curvature);
        } else if (_periodic) {
            Detail::CheckPeriodicClosure(y, "CubicSplineSystem::Solve");
            SolvePeriodic(y, curvature);
        } else {
            SolveCurvatureForm(y, leftDerivative, rightDerivative, curvature);
        }
    }

    /**
     * @brief Index of the segment that owns node `k`, from a chosen side.
     *
     * `Side::Right` agrees with the segment an ordinary query at that node
     * would select. On a Periodic system the seam is not an end, so the
     * closing node is answered from the first segment on its right, and the
     * first node from the last segment on its left.
     */
    std::size_t NodeSegment(std::size_t k, Side side = Side::Right) const {
        const auto n = Size();
        if (_periodic) {
            if (side == Side::Right && k + 1 == n) {
                return 0;
            }
            if (side == Side::Left && k == 0) {
                return n - 2;
            }
        }
        return Detail::NodeSegment(n, k, side);
    }

    /**
     * @brief Evaluate a solved spline, or its `N`th derivative, at every node.
     *
     * The other half of what a differentiation operator wants. Given the
     * ordinates and the curvatures Solve produced for them, this walks the
     * segments once and writes the nodal values — where the general query path
     * would pay a binary search per point to find segments it already knows.
     *
     * It takes the ordinates and curvatures rather than holding them, so a
     * caller sweeping many datasets over one grid keeps its own buffers and
     * this class stays a description of the grid.
     *
     * Values at a node agree from both sides, as do the first and second
     * derivatives; the third jumps, and `side` chooses which limit is
     * reported. The default is the right-hand one, matching the rest of the
     * library. Orders above three are identically zero.
     *
     * @tparam N Derivative order; `0` is the value itself.
     * @param y Ordinates, one per node.
     * @param curvature Nodal second derivatives, as returned by Solve.
     * @param out Output, one per node. Overwritten.
     * @param side Which adjoining segment to answer from at a node.
     * @throws std::invalid_argument if any length is wrong.
     */
    template <std::size_t N = 0, typename YRange, typename CRange,
              typename Scalar>
        requires InterpolationRanges<XView, YRange>
    void EvaluateAtNodes(const YRange &y, const CRange &curvature,
                         std::span<Scalar> out, Side side = Side::Right) const {
        const auto n = Size();
        CheckLength(static_cast<std::size_t>(std::ranges::size(y)), n,
                    "ordinate range");
        CheckLength(static_cast<std::size_t>(std::ranges::size(curvature)), n,
                    "curvature range");
        CheckLength(out.size(), n, "output range");

        for (std::size_t k = 0; k < n; ++k) {
            const auto i = NodeSegment(k, side);
            // Node k is the left end of its segment, or the right end. The
            // one case where neither index says which is the closing node of
            // a periodic system, answered from the start of the first
            // segment.
            const auto atLeftEnd =
                i == k || (_periodic && k + 1 == n && i == 0);
            out[k] = Detail::SplinePiece<N, Real, Scalar>(
                _spacing[i], atLeftEnd ? Real{} : _spacing[i], y[i], y[i + 1],
                curvature[i], curvature[i + 1]);
        }
    }

    /**
     * @brief Nodal second derivatives, returned in a fresh vector.
     *
     * The convenient form of Solve. It allocates, so the per-line inner loop
     * of an operator should use Solve with its own buffer instead.
     */
    template <typename YRange>
        requires InterpolationRanges<XView, YRange>
    std::vector<std::ranges::range_value_t<YRange>>
    Curvatures(const YRange &y) const {
        std::vector<std::ranges::range_value_t<YRange>> curvature(Size());
        Solve(y, std::span{curvature});
        return curvature;
    }

  private:
    XView _x;
    BoundaryCondition _left;
    BoundaryCondition _right;
    bool _notAKnot;
    bool _periodic;
    std::vector<Real> _spacing; // Segment widths; filled by Assemble.
    TridiagonalFactorization<Real> _factorization;

    // The Sherman-Morrison correction a Periodic system needs; empty and
    // unread otherwise.
    std::vector<Real> _correction; // z = T^{-1} u
    Real _cornerRatio{};           // c / gamma, the last entry of v
    Real _inverseDenominator{};    // 1 / (1 + v^T z)

    static void CheckLength(std::size_t given, std::size_t expected,
                            const char *what) {
        if (given != expected) {
            throw std::invalid_argument(
                std::string{"CubicSplineSystem::Solve: the "} + what +
                " has length " + std::to_string(given) +
                " but the system has " + std::to_string(expected) + " nodes");
        }
    }

    // Build the three diagonals, factorise them, and leave the segment widths
    // behind in `spacing`. Called from the member initialiser list, so it takes
    // the pieces it fills rather than reading half-built members.
    static TridiagonalFactorization<Real> Assemble(const XView &x,
                                                   BoundaryCondition left,
                                                   BoundaryCondition right,
                                                   std::vector<Real> &spacing) {
        Detail::ValidateSamples(
            x, x,
            Detail::SplineConditionMinimum(left, right, "CubicSplineSystem"),
            "CubicSplineSystem");

        const auto n = static_cast<std::size_t>(std::ranges::size(x));
        spacing.resize(n - 1);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            spacing[i] = x[i + 1] - x[i];
        }

        // A periodic system has one unknown fewer than there are nodes: the
        // closing node's curvature is the first's.
        const auto order = left == BoundaryCondition::Periodic ? n - 1 : n;
        std::vector<Real> sub(order, Real{}), diag(order, Real{}),
            super(order, Real{});
        if (left == BoundaryCondition::NotAKnot) {
            AssembleNotAKnot(spacing, sub, diag, super);
        } else if (left == BoundaryCondition::Periodic) {
            AssemblePeriodic(spacing, sub, diag, super);
        } else {
            AssembleCurvatureForm(spacing, left, right, sub, diag, super);
        }

        return TridiagonalFactorization<Real>{std::span<const Real>{sub},
                                              std::span<const Real>{diag},
                                              std::span<const Real>{super}};
    }

    // Interior rows express continuity of the first derivative, written in
    // nodal curvatures.
    static void
    AssembleCurvatureForm(const std::vector<Real> &h, BoundaryCondition left,
                          BoundaryCondition right, std::vector<Real> &sub,
                          std::vector<Real> &diag, std::vector<Real> &super) {
        const auto n = diag.size();
        constexpr auto oneThird = static_cast<Real>(1) / static_cast<Real>(3);
        constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);

        for (std::size_t i = 1; i + 1 < n; ++i) {
            sub[i] = oneSixth * h[i - 1];
            diag[i] = oneThird * (h[i - 1] + h[i]);
            super[i] = oneSixth * h[i];
        }

        // A Natural condition states directly that the end curvature is zero.
        if (left == BoundaryCondition::Natural) {
            diag[0] = 1;
        } else {
            diag[0] = oneThird * h[0];
            super[0] = oneSixth * h[0];
        }

        if (right == BoundaryCondition::Natural) {
            diag[n - 1] = 1;
        } else {
            sub[n - 1] = oneSixth * h[n - 2];
            diag[n - 1] = oneThird * h[n - 2];
        }
    }

    // Interior rows express continuity of the second derivative, written in
    // nodal slopes; the boundary rows are the not-a-knot condition.
    static void AssembleNotAKnot(const std::vector<Real> &h,
                                 std::vector<Real> &sub,
                                 std::vector<Real> &diag,
                                 std::vector<Real> &super) {
        const auto n = diag.size();

        for (std::size_t i = 1; i + 1 < n; ++i) {
            sub[i] = h[i];
            diag[i] = 2 * (h[i - 1] + h[i]);
            super[i] = h[i - 1];
        }

        diag[0] = h[1];
        super[0] = h[0] + h[1];

        sub[n - 1] = h[n - 2] + h[n - 3];
        diag[n - 1] = h[n - 3];
    }

    // Every row is interior, with indices taken modulo the m distinct nodes.
    // The two corner entries are left out here and folded into the first and
    // last diagonal entries, which leaves the tridiagonal T that is
    // factorised; PrepareCorrection supplies the rank-one remainder.
    static void AssemblePeriodic(const std::vector<Real> &h,
                                 std::vector<Real> &sub,
                                 std::vector<Real> &diag,
                                 std::vector<Real> &super) {
        const auto m = diag.size();
        constexpr auto oneThird = static_cast<Real>(1) / static_cast<Real>(3);
        constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);

        for (std::size_t i = 0; i < m; ++i) {
            const auto previous = h[i == 0 ? m - 1 : i - 1];
            sub[i] = oneSixth * previous;
            diag[i] = oneThird * (previous + h[i]);
            super[i] = oneSixth * h[i];
        }

        const auto corner = oneSixth * h[m - 1];
        const auto gamma = -diag[0];
        diag[0] -= gamma;
        diag[m - 1] -= corner * corner / gamma;
    }

    // Solve T z = u once, for u = (gamma, 0, ..., 0, c), and keep what every
    // later correction needs. Run from the constructor body, once the
    // factorisation it uses exists.
    void PrepareCorrection() {
        const auto m = Size() - 1;
        constexpr auto oneThird = static_cast<Real>(1) / static_cast<Real>(3);
        constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);
        const auto corner = oneSixth * _spacing[m - 1];
        const auto gamma = -oneThird * (_spacing[m - 1] + _spacing[0]);

        _correction.assign(m, Real{});
        _correction[0] = gamma;
        _correction[m - 1] = corner;
        _factorization.Solve(std::span<Real>{_correction});

        _cornerRatio = corner / gamma;
        _inverseDenominator =
            static_cast<Real>(1) / (static_cast<Real>(1) + _correction[0] +
                                    _cornerRatio * _correction[m - 1]);
    }

    template <typename YRange, typename Scalar>
    void SolveCurvatureForm(const YRange &y, Scalar leftDerivative,
                            Scalar rightDerivative,
                            std::span<Scalar> curvature) const {
        const auto n = Size();
        const auto &h = _spacing;

        for (std::size_t i = 1; i + 1 < n; ++i) {
            curvature[i] =
                (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1];
        }

        curvature[0] = _left == BoundaryCondition::Natural
                           ? Scalar{}
                           : (y[1] - y[0]) / h[0] - leftDerivative;
        curvature[n - 1] =
            _right == BoundaryCondition::Natural
                ? Scalar{}
                : rightDerivative - (y[n - 1] - y[n - 2]) / h[n - 2];

        _factorization.Solve(curvature);
    }

    template <typename YRange, typename Scalar>
    void SolvePeriodic(const YRange &y, std::span<Scalar> curvature) const {
        const auto m = Size() - 1;
        const auto &h = _spacing;
        const auto secant = [&](std::size_t i) {
            return (y[i + 1] - y[i]) / h[i];
        };

        curvature[0] = secant(0) - secant(m - 1);
        for (std::size_t i = 1; i < m; ++i) {
            curvature[i] = secant(i) - secant(i - 1);
        }

        const auto unknowns = curvature.first(m);
        _factorization.Solve(unknowns);

        const auto factor = (unknowns[0] + _cornerRatio * unknowns[m - 1]) *
                            _inverseDenominator;
        for (std::size_t i = 0; i < m; ++i) {
            unknowns[i] -= factor * _correction[i];
        }
        curvature[m] = curvature[0];
    }

    template <typename YRange, typename Scalar>
    void SolveNotAKnot(const YRange &y, std::span<Scalar> curvature) const {
        const auto n = Size();
        const auto &h = _spacing;
        const auto secant = [&](std::size_t i) {
            return (y[i + 1] - y[i]) / h[i];
        };

        // The right-hand side is built in the output and solved there, so the
        // nodal slopes arrive in `curvature` before being converted.
        for (std::size_t i = 1; i + 1 < n; ++i) {
            curvature[i] = static_cast<Real>(3) *
                           (h[i] * secant(i - 1) + h[i - 1] * secant(i));
        }
        curvature[0] = (h[1] * (2 * h[1] + 3 * h[0]) * secant(0) +
                        h[0] * h[0] * secant(1)) /
                       (h[0] + h[1]);
        curvature[n - 1] =
            (h[n - 3] * (2 * h[n - 3] + 3 * h[n - 2]) * secant(n - 2) +
             h[n - 2] * h[n - 2] * secant(n - 3)) /
            (h[n - 3] + h[n - 2]);

        _factorization.Solve(curvature);

        // Convert nodal slopes to nodal curvatures: on segment i the piece is
        // y_i + d_i t + c t^2 + e t^3, so S''(x_i) = 2c. The conversion is done
        // in place, which needs the final node computed first — it is the only
        // one that reads a slope the forward sweep has already overwritten.
        const auto last = (static_cast<Real>(2) * curvature[n - 2] +
                           static_cast<Real>(4) * curvature[n - 1] -
                           static_cast<Real>(6) * secant(n - 2)) /
                          h[n - 2];
        for (std::size_t i = 0; i + 1 < n; ++i) {
            curvature[i] =
                static_cast<Real>(2) *
                (static_cast<Real>(3) * secant(i) -
                 static_cast<Real>(2) * curvature[i] - curvature[i + 1]) /
                h[i];
        }
        curvature[n - 1] = last;
    }
};

/// Borrow an lvalue container of nodes and own an rvalue one.
template <std::ranges::viewable_range X, typename... Rest>
CubicSplineSystem(X &&, Rest...) -> CubicSplineSystem<std::views::all_t<X>>;

} // namespace Interpolation

#endif // INTERPOLATION_CUBIC_SPLINE_SYSTEM_HPP
