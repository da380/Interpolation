#ifndef INTERPOLATION_SIDE_HPP
#define INTERPOLATION_SIDE_HPP

namespace Interpolation {

/**
 * @brief Which side of a breakpoint or node a query should be answered from.
 *
 * The library is right-continuous by default: a piece owns the half-open
 * interval `[lower, upper)`, so a query landing exactly on a breakpoint is
 * answered from the piece starting there. `Side` is how a caller asks for the
 * other limit.
 *
 * It lives in its own header because the convention is shared. `Piecewise`
 * uses it across breakpoints, where a genuine discontinuity may sit; the
 * one-dimensional interpolators use it at their own nodes, where the value
 * agrees from both sides but a high enough derivative need not — the third
 * derivative of a cubic spline and the first of a piecewise-linear
 * interpolant both jump there.
 */
enum class Side {
    /** The piece ending at the breakpoint. */
    Left,
    /** The piece starting at the breakpoint. */
    Right
};

} // namespace Interpolation

#endif // INTERPOLATION_SIDE_HPP
