#ifndef INTERPOLATION_CUBIC_SPLINE_SYSTEM_HPP
#define INTERPOLATION_CUBIC_SPLINE_SYSTEM_HPP

#include <cstddef>
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
    NotAKnot
};

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
     *        NotAKnot is used.
     * @param left Condition at the first node.
     * @param right Condition at the final node.
     * @throws std::invalid_argument if the abscissae are not strictly
     *         increasing, there are too few of them, or NotAKnot is asked for
     *         at only one end.
     */
    CubicSplineSystem(XView x, BoundaryCondition left, BoundaryCondition right)
        : _x{std::move(x)}, _left{left}, _right{right},
          _notAKnot{left == BoundaryCondition::NotAKnot},
          _factorization{Assemble(_x, left, right, _spacing)} {}

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
     * @throws std::invalid_argument if either length is wrong.
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
        } else {
            SolveCurvatureForm(y, leftDerivative, rightDerivative, curvature);
        }
    }

    /**
     * @brief Index of the segment that owns node `k`, from a chosen side.
     *
     * `Side::Right` agrees with the segment an ordinary query at that node
     * would select.
     */
    std::size_t NodeSegment(std::size_t k, Side side = Side::Right) const {
        return Detail::NodeSegment(Size(), k, side);
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
            const auto i = Detail::NodeSegment(n, k, side);
            out[k] = Detail::SplinePiece<N, Real, Scalar>(
                _spacing[i], _x[k] - _x[i], y[i], y[i + 1], curvature[i],
                curvature[i + 1]);
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
    std::vector<Real> _spacing; // Segment widths; filled by Assemble.
    TridiagonalFactorization<Real> _factorization;

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
        const auto notAKnot = left == BoundaryCondition::NotAKnot;
        if (notAKnot != (right == BoundaryCondition::NotAKnot)) {
            throw std::invalid_argument(
                "CubicSplineSystem: NotAKnot constrains the whole system "
                "rather than one endpoint, so it must be used at both ends or "
                "neither");
        }
        Detail::ValidateSamples(x, x, notAKnot ? 4 : 2, "CubicSplineSystem");

        const auto n = static_cast<std::size_t>(std::ranges::size(x));
        spacing.resize(n - 1);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            spacing[i] = x[i + 1] - x[i];
        }

        std::vector<Real> sub(n, Real{}), diag(n, Real{}), super(n, Real{});
        if (notAKnot) {
            AssembleNotAKnot(spacing, sub, diag, super);
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
