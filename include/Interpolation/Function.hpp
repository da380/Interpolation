#ifndef INTERPOLATION_FUNCTION_HPP
#define INTERPOLATION_FUNCTION_HPP

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>

namespace Interpolation {

/**
 * @brief A real-argument function that can report derivatives of any order.
 *
 * Every interpolator in the library models this, and so does every node built
 * by the operators below, which is what makes them composable. The concept
 * asks for two type names and one operation:
 *
 * - `Real`, the abscissa precision;
 * - `Scalar`, the value type, real or complex;
 * - `Evaluate<N>(x)`, the `N`th derivative at `x`, with `N == 0` the value.
 *
 * `operator()` is `Evaluate<0>` and is provided by every model, but it is not
 * part of the requirement: a node only has to know how to differentiate.
 */
template <typename F>
concept Function1D = requires(const F &f, typename F::Real x) {
    typename F::Real;
    typename F::Scalar;
    requires ::Interpolation::Real<typename F::Real>;
    requires RealOrComplex<typename F::Scalar>;
    { f.template Evaluate<0>(x) } -> std::convertible_to<typename F::Scalar>;
    { f.template Evaluate<1>(x) } -> std::convertible_to<typename F::Scalar>;
};

/**
 * @brief A function made of polynomial pieces between known nodes.
 *
 * Antidifferentiation needs more than Evaluate<N>: it needs to know where the
 * pieces meet and what each piece integrates to. Interpolators expose that
 * through three small accessors rather than carrying an antiderivative that
 * most callers never ask for.
 */
template <typename F>
concept PiecewisePolynomial1D =
    requires(const F &f, std::size_t i, typename F::Real x) {
        requires Function1D<F>;
        { f.Size() } -> std::convertible_to<std::size_t>;
        { f.Node(i) } -> std::convertible_to<typename F::Real>;
        { f.Segment(x) } -> std::convertible_to<std::size_t>;
        { f.SegmentIntegral(i, x) } -> std::convertible_to<typename F::Scalar>;
    };

/** @brief A function that knows its own antiderivative in closed form. */
template <typename F>
concept Antidifferentiable1D = requires(const F &f, typename F::Real x) {
    requires Function1D<F>;
    { f.Antiderivative(x) } -> std::convertible_to<typename F::Scalar>;
};

/** @brief Two functions that share an abscissa type and can be combined. */
template <typename F, typename G>
concept CompatibleFunctions1D = requires() {
    requires Function1D<F>;
    requires Function1D<G>;
    requires std::same_as<typename F::Real, typename G::Real>;
};

namespace Detail {

/** @brief Binomial coefficient, evaluated at compile time. */
consteval std::size_t
Binomial(std::size_t n, std::size_t k) {
    if (k > n) {
        return 0;
    }
    std::size_t result = 1;
    for (std::size_t i = 0; i < k; ++i) {
        // Exact in this order: the partial product is always divisible.
        result = result * (n - i) / (i + 1);
    }
    return result;
}

/** @brief Result value type when two functions are combined. */
template <typename F, typename G>
using CombinedScalar =
    std::common_type_t<typename F::Scalar, typename G::Scalar>;

} // namespace Detail

// ---------------------------------------------------------------------------
// Leaves
// ---------------------------------------------------------------------------

/**
 * @brief A constant function.
 *
 * Scalars enter the algebra through this, so `f + 2.0` needs no separate
 * node type and no separate derivative rule.
 */
template <typename R, typename S>
    requires ::Interpolation::Real<R> && RealOrComplex<S>
class Constant {
  public:
    /** @brief Abscissa precision. */
    using Real = R;
    /** @brief Value type, real or complex. */
    using Scalar = S;

    /** @brief Construct the node, taking ownership of its operands. */
    constexpr explicit Constant(Scalar value) : _value{value} {}

    /**
     * @brief Evaluate the node or its `N`th derivative.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        (void) x;
        if constexpr (N == 0) {
            return _value;
        } else {
            return Scalar{};
        }
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    Scalar _value;
};

/** @brief The identity function, `x`. */
template <typename R>
    requires ::Interpolation::Real<R>
class Identity {
  public:
    /** @brief Abscissa precision. */
    using Real = R;
    /** @brief Value type, real or complex. */
    using Scalar = R;

    /**
     * @brief Evaluate the node or its `N`th derivative.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        if constexpr (N == 0) {
            return x;
        } else if constexpr (N == 1) {
            return static_cast<Scalar>(1);
        } else {
            return Scalar{};
        }
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }
};

// ---------------------------------------------------------------------------
// Nodes
//
// Every node stores its operands BY VALUE. That is the whole point: a node
// holding references would dangle the moment it was built from a temporary,
// which is exactly what `f + g` produces when either side is an expression.
// Interpolators are cheap to copy, since they hold views rather than data
// unless they were given rvalues.
// ---------------------------------------------------------------------------

/** @brief The pointwise negation of a function. */
template <Function1D F> class Negate {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = typename F::Scalar;

    /** @brief Construct the node, taking ownership of its operands. */
    constexpr explicit Negate(F f) : _f{std::move(f)} {}

    /**
     * @brief Evaluate the node or its `N`th derivative.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        return -_f.template Evaluate<N>(x);
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    F _f;
};

/** @brief The pointwise sum of two functions. */
template <typename F, typename G>
    requires CompatibleFunctions1D<F, G>
class Sum {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = Detail::CombinedScalar<F, G>;

    /** @brief Construct the node, taking ownership of its operands. */
    constexpr Sum(F f, G g) : _f{std::move(f)}, _g{std::move(g)} {}

    /** @brief Differentiation is linear, so this holds for every order. */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        return static_cast<Scalar>(_f.template Evaluate<N>(x)) +
               static_cast<Scalar>(_g.template Evaluate<N>(x));
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    F _f;
    G _g;
};

/** @brief The pointwise difference of two functions. */
template <typename F, typename G>
    requires CompatibleFunctions1D<F, G>
class Difference {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = Detail::CombinedScalar<F, G>;

    /** @brief Construct the node, taking ownership of its operands. */
    constexpr Difference(F f, G g) : _f{std::move(f)}, _g{std::move(g)} {}

    /**
     * @brief Evaluate the node or its `N`th derivative.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        return static_cast<Scalar>(_f.template Evaluate<N>(x)) -
               static_cast<Scalar>(_g.template Evaluate<N>(x));
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    F _f;
    G _g;
};

/** @brief The pointwise product of two functions. */
template <typename F, typename G>
    requires CompatibleFunctions1D<F, G>
class Product {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = Detail::CombinedScalar<F, G>;

    /** @brief Construct the node, taking ownership of its operands. */
    constexpr Product(F f, G g) : _f{std::move(f)}, _g{std::move(g)} {}

    /**
     * @brief Evaluate by the general Leibniz rule.
     *
     * @f$ (fg)^{(N)} = \sum_k \binom{N}{k} f^{(k)} g^{(N-k)} @f$, expanded at
     * compile time, so an arbitrary derivative order costs no branching and
     * no storage.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        return [&]<std::size_t... K>(std::index_sequence<K...>) {
            return ((static_cast<Real>(Detail::Binomial(N, K)) *
                     static_cast<Scalar>(_f.template Evaluate<K>(x)) *
                     static_cast<Scalar>(_g.template Evaluate<N - K>(x))) +
                    ...);
        }(std::make_index_sequence<N + 1>{});
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    F _f;
    G _g;
};

/** @brief The pointwise quotient of two functions. */
template <typename F, typename G>
    requires CompatibleFunctions1D<F, G>
class Quotient {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = Detail::CombinedScalar<F, G>;

    /** @brief Construct the node, taking ownership of its operands. */
    constexpr Quotient(F f, G g) : _f{std::move(f)}, _g{std::move(g)} {}

    /**
     * @brief Evaluate the quotient or its first derivative.
     *
     * Higher orders would need the general recurrence for derivatives of a
     * reciprocal. That is deliberately not guessed at here: the alternative
     * is silently returning a wrong number.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        static_assert(N <= 1,
                      "Quotient supports the value and the first derivative; "
                      "higher orders need a recurrence that is not "
                      "implemented.");
        const auto g = static_cast<Scalar>(_g.template Evaluate<0>(x));
        if constexpr (N == 0) {
            return static_cast<Scalar>(_f.template Evaluate<0>(x)) / g;
        } else {
            const auto f = static_cast<Scalar>(_f.template Evaluate<0>(x));
            const auto df = static_cast<Scalar>(_f.template Evaluate<1>(x));
            const auto dg = static_cast<Scalar>(_g.template Evaluate<1>(x));
            return (df * g - f * dg) / (g * g);
        }
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    F _f;
    G _g;
};

/**
 * @brief A function differentiated `K` times.
 *
 * This is the change that makes differentiation composable: it is a function
 * in its own right, so it can be added, multiplied and differentiated again,
 * rather than being a member that returns a number and ends the story.
 */
template <std::size_t K, Function1D F> class DerivativeNode {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = typename F::Scalar;

    /** @brief Construct the node, taking ownership of its operands. */
    constexpr explicit DerivativeNode(F f) : _f{std::move(f)} {}

    /**
     * @brief Evaluate the node or its `N`th derivative.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        return _f.template Evaluate<N + K>(x);
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    F _f;
};

/** @brief The composition `f(g(x))`. */
template <typename F, typename G>
    requires Function1D<F> && Function1D<G> &&
             std::convertible_to<typename G::Scalar, typename F::Real>
class Composition {
  public:
    /** @brief Abscissa precision. */
    using Real = typename G::Real;
    /** @brief Value type, real or complex. */
    using Scalar = typename F::Scalar;

    /** @brief Construct the node, taking ownership of its operands. */
    constexpr Composition(F f, G g) : _f{std::move(f)}, _g{std::move(g)} {}

    /**
     * @brief Evaluate the composition or its first derivative.
     *
     * Higher orders need Faa di Bruno's formula, which is not implemented.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        static_assert(N <= 1,
                      "Composition supports the value and the first "
                      "derivative; higher orders need Faa di Bruno's formula, "
                      "which is not implemented.");
        const auto inner =
            static_cast<typename F::Real>(_g.template Evaluate<0>(x));
        if constexpr (N == 0) {
            return _f.template Evaluate<0>(inner);
        } else {
            return _f.template Evaluate<1>(inner) *
                   static_cast<Scalar>(_g.template Evaluate<1>(x));
        }
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

  private:
    F _f;
    G _g;
};

/**
 * @brief The antiderivative of a piecewise-polynomial function.
 *
 * Vanishes at the first node. The cumulative integral at each node is
 * computed once, here, so an interpolator carries no antiderivative state
 * unless one is actually asked for. Evaluation then costs a segment lookup
 * and one polynomial.
 *
 * Its own derivatives are just the operand's, shifted down by one, so a
 * primitive composes like anything else in the algebra.
 */
template <PiecewisePolynomial1D F> class PiecewisePrimitive {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = typename F::Scalar;

    /** @brief Construct the node, taking ownership of its operands. */
    explicit PiecewisePrimitive(F f) : _f{std::move(f)} {
        const auto n = _f.Size();
        _cumulative.assign(n, Scalar{});
        for (std::size_t i = 0; i + 1 < n; ++i) {
            _cumulative[i + 1] =
                _cumulative[i] +
                _f.SegmentIntegral(i, _f.Node(i + 1) - _f.Node(i));
        }
    }

    /**
     * @brief Evaluate the node or its `N`th derivative.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> Scalar Evaluate(Real x) const {
        if constexpr (N == 0) {
            const auto i = _f.Segment(x);
            return _cumulative[i] + _f.SegmentIntegral(i, x - _f.Node(i));
        } else {
            return _f.template Evaluate<N - 1>(x);
        }
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    Scalar operator()(Real x) const { return Evaluate<0>(x); }

    /** @brief Definite integral over `[a, b]`. */
    Scalar Integral(Real a, Real b) const {
        return Evaluate<0>(b) - Evaluate<0>(a);
    }

  private:
    F _f;
    std::vector<Scalar> _cumulative;
};

/**
 * @brief The antiderivative of a function that knows its own.
 *
 * Used for closed-form functions such as Polynomial, where no cumulative
 * table is needed.
 */
template <Antidifferentiable1D F> class ClosedFormPrimitive {
  public:
    /** @brief Abscissa precision. */
    using Real = typename F::Real;
    /** @brief Value type, real or complex. */
    using Scalar = typename F::Scalar;

    /** @brief Construct the node, taking ownership of its operands. */
    explicit constexpr ClosedFormPrimitive(F f) : _f{std::move(f)} {}

    /**
     * @brief Evaluate the node or its `N`th derivative.
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        if constexpr (N == 0) {
            return _f.Antiderivative(x);
        } else {
            return _f.template Evaluate<N - 1>(x);
        }
    }

    /** @brief Evaluate the node; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

    /** @brief Definite integral over `[a, b]`. */
    constexpr Scalar Integral(Real a, Real b) const {
        return Evaluate<0>(b) - Evaluate<0>(a);
    }

  private:
    F _f;
};

// ---------------------------------------------------------------------------
// Free functions and operators
//
// These are constrained on Function1D, so they take part in overload
// resolution only for actual functions and cannot capture unrelated types.
// ---------------------------------------------------------------------------

/** @brief Differentiate a function `K` times, giving another function. */
template <std::size_t K = 1, Function1D F>
constexpr auto
Derivative(F f) {
    return DerivativeNode<K, F>{std::move(f)};
}

/**
 * @brief The antiderivative of `f`, as a function.
 *
 * For a piecewise-polynomial function this vanishes at the first node; for a
 * closed-form one it is whatever that type defines. Either way the result is
 * a Function1D, so it can be differentiated back, added, or composed.
 */
template <PiecewisePolynomial1D F>
auto
Primitive(F f) {
    return PiecewisePrimitive<F>{std::move(f)};
}

/** @copydoc Primitive */
template <Antidifferentiable1D F>
    requires(!PiecewisePolynomial1D<F>)
constexpr auto
Primitive(F f) {
    return ClosedFormPrimitive<F>{std::move(f)};
}

/** @brief Compose two functions, giving `f(g(x))`. */
template <typename F, typename G>
constexpr auto
Compose(F f, G g) {
    return Composition<F, G>{std::move(f), std::move(g)};
}

/** @brief Wrap a scalar as a constant function over the abscissae of `F`. */
template <Function1D F, typename S>
    requires RealOrComplex<S>
constexpr auto
AsConstant(S value) {
    return Constant<typename F::Real, S>{value};
}

/** @brief Negate a function pointwise. */
template <Function1D F>
constexpr auto
operator-(F f) {
    return Negate<F>{std::move(f)};
}

/** @brief Add two functions pointwise. */
template <typename F, typename G>
    requires CompatibleFunctions1D<F, G>
constexpr auto
operator+(F f, G g) {
    return Sum<F, G>{std::move(f), std::move(g)};
}

/** @brief Subtract two functions pointwise. */
template <typename F, typename G>
    requires CompatibleFunctions1D<F, G>
constexpr auto
operator-(F f, G g) {
    return Difference<F, G>{std::move(f), std::move(g)};
}

/** @brief Multiply two functions pointwise. */
template <typename F, typename G>
    requires CompatibleFunctions1D<F, G>
constexpr auto
operator*(F f, G g) {
    return Product<F, G>{std::move(f), std::move(g)};
}

/** @brief Divide two functions pointwise. */
template <typename F, typename G>
    requires CompatibleFunctions1D<F, G>
constexpr auto
operator/(F f, G g) {
    return Quotient<F, G>{std::move(f), std::move(g)};
}

/** @brief Add a scalar to a function. */
template <Function1D F, typename S>
    requires RealOrComplex<S>
constexpr auto
operator+(F f, S value) {
    return std::move(f) + AsConstant<F>(value);
}

/** @brief Add a function to a scalar. */
template <Function1D F, typename S>
    requires RealOrComplex<S>
constexpr auto
operator+(S value, F f) {
    return AsConstant<F>(value) + std::move(f);
}

/** @brief Subtract a scalar from a function. */
template <Function1D F, typename S>
    requires RealOrComplex<S>
constexpr auto
operator-(F f, S value) {
    return std::move(f) - AsConstant<F>(value);
}

/** @brief Subtract a function from a scalar. */
template <Function1D F, typename S>
    requires RealOrComplex<S>
constexpr auto
operator-(S value, F f) {
    return AsConstant<F>(value) - std::move(f);
}

/** @brief Scale a function by a scalar. */
template <Function1D F, typename S>
    requires RealOrComplex<S>
constexpr auto
operator*(F f, S value) {
    return std::move(f) * AsConstant<F>(value);
}

/** @brief Scale a function by a scalar (left). */
template <Function1D F, typename S>
    requires RealOrComplex<S>
constexpr auto
operator*(S value, F f) {
    return AsConstant<F>(value) * std::move(f);
}

/** @brief Divide a function by a scalar. */
template <Function1D F, typename S>
    requires RealOrComplex<S>
constexpr auto
operator/(F f, S value) {
    return std::move(f) / AsConstant<F>(value);
}

} // namespace Interpolation

#endif // INTERPOLATION_FUNCTION_HPP
