#ifndef INTERPOLATION_POLYNOMIAL_HPP
#define INTERPOLATION_POLYNOMIAL_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>

namespace Interpolation {

/**
 * @brief One-dimensional polynomial with coefficients in ascending power order.
 *
 * A coefficient sequence `{a0, a1, ..., an}` represents
 * `a0 + a1*x + ... + an*x^n`.
 *
 * @tparam T Real or complex floating-point coefficient type.
 */
template <typename T>
    requires RealOrComplex<T>
class Polynomial {
  public:
    /** @brief Scalar type used for coefficients and evaluation. */
    using value_type = T;
    /**
     * @brief Abscissa precision, so that Polynomial models Function1D.
     *
     * A complex polynomial is still a function of a real abscissa as far as
     * the algebra is concerned; its coefficients and values are what carry
     * the imaginary part.
     */
    using Real = RemoveComplex<T>;
    /** @brief Value type, real or complex. */
    using Scalar = T;
    /** @brief Mutable coefficient iterator. */
    using iterator = std::vector<T>::iterator;
    /** @brief Constant coefficient iterator. */
    using const_iterator = std::vector<T>::const_iterator;

    /** @brief Construct the zero polynomial with one coefficient equal to zero.
     */
    Polynomial() : _a{T{}} {}

    /** @brief Copy coefficients from a vector in ascending power order. */
    Polynomial(const std::vector<T> &a) : _a{a} {}
    /** @brief Move coefficients from a vector in ascending power order. */
    Polynomial(std::vector<T> &&a) : _a{std::move(a)} {}

    /** @brief Construct from coefficients in ascending power order. */
    Polynomial(std::initializer_list<T> list) : _a{std::vector<T>{list}} {}

    /** @brief Copy a polynomial. */
    Polynomial(const Polynomial &) = default;
    /** @brief Move a polynomial. */
    Polynomial(Polynomial &&) = default;

    /** @brief Copy-assign a polynomial with the same coefficient type. */
    Polynomial &operator=(const Polynomial &) = default;
    /** @brief Move-assign a polynomial with the same coefficient type. */
    Polynomial &operator=(Polynomial &&) = default;

    /**
     * @brief Copy a polynomial while converting its coefficient type.
     * @tparam FLOAT Source coefficient type.
     * @param rhs Source polynomial.
     */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial(const Polynomial<FLOAT> &rhs) {
        std::transform(rhs.cbegin(), rhs.cend(), std::back_inserter(_a),
                       [](auto x) { return static_cast<T>(x); });
    }

    /**
     * @brief Assign coefficients from a polynomial with a convertible type.
     * @tparam FLOAT Source coefficient type.
     * @param polinit Source polynomial.
     * @return This polynomial.
     */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial<T> &operator=(const Polynomial<FLOAT> &polinit) {
        std::vector<T> converted;
        converted.reserve(std::distance(polinit.cbegin(), polinit.cend()));
        std::transform(polinit.cbegin(), polinit.cend(),
                       std::back_inserter(converted),
                       [](auto value) { return static_cast<T>(value); });
        _a = std::move(converted);
        return *this;
    }

    /**
     * @brief Generate a real polynomial with normally distributed coefficients.
     * @param n Requested nonnegative degree.
     * @param seed Seed for the generator. Explicit so that a randomised test
     *        can report it and a failure can be replayed exactly.
     * @return Random polynomial of degree `n`.
     */
    static Polynomial Random(int n, std::uint64_t seed)
        requires ::Interpolation::Real<T>
    {
        std::mt19937_64 gen{seed};
        std::normal_distribution<T> d{};
        std::vector<T> a;
        std::generate_n(std::back_inserter(a), n + 1, [&] { return d(gen); });
        return Polynomial(std::move(a));
    }

    /**
     * @brief Generate a complex polynomial with normally distributed
     * components.
     * @param n Requested nonnegative degree.
     * @param seed Seed for the generator. Explicit so that a randomised test
     *        can report it and a failure can be replayed exactly.
     * @return Random polynomial of degree `n`.
     */
    static Polynomial Random(int n, std::uint64_t seed)
        requires ::Interpolation::Complex<T>
    {
        using S = typename T::value_type;
        std::mt19937_64 gen{seed};
        std::normal_distribution<S> d{};
        std::vector<T> a;
        std::generate_n(std::back_inserter(a), n + 1,
                        [&] { return T{d(gen), d(gen)}; });
        return Polynomial(std::move(a));
    }

    /**
     * @brief Draw a seed from the system entropy source.
     *
     * Randomised tests should seed deterministically and report the seed they
     * used, so that a failure can be replayed. Call this once, print it, and
     * pass it to Random().
     */
    static std::uint64_t RandomSeed() {
        std::random_device rd{};
        return (static_cast<std::uint64_t>(rd()) << 32) ^ rd();
    }

    /**
     * @brief Drop trailing zero coefficients.
     *
     * Degree() is derived from the coefficient count, so without this a
     * subtraction that cancels the leading terms leaves the polynomial
     * claiming a degree it no longer has.
     */
    Polynomial &Trim() {
        while (_a.size() > 1 && _a.back() == T{}) {
            _a.pop_back();
        }
        return *this;
    }

    /**
     * @brief Return the polynomial degree.
     *
     * Signed deliberately: an unsigned `size() - 1` yields SIZE_MAX for an
     * empty coefficient vector, and makes every loop over the coefficients a
     * signed/unsigned comparison.
     */
    int Degree() const { return static_cast<int>(_a.size()) - 1; }

    /** @brief Evaluate the polynomial at `x` using Horner's method. */
    T operator()(T x) const {
        return std::accumulate(_a.rbegin(), _a.rend(), static_cast<T>(0),
                               [x](auto p, auto a) { return p * x + a; });
    }

    /** @brief Evaluate the first derivative at `x`. */
    T Derivative(T x) const {
        return std::accumulate(_a.rbegin(), std::prev(_a.rend()),
                               static_cast<T>(0),
                               [x, m = Degree()](auto p, auto a) mutable {
                                   return p * x + static_cast<T>(m--) * a;
                               });
    }

    /**
     * @brief Evaluate the polynomial or its `N`th derivative.
     *
     * This is the spelling every interpolator in the library shares, so a
     * polynomial can stand in wherever one of them can.
     *
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> T Evaluate(T x) const {
        if constexpr (N == 0) {
            return (*this)(x);
        } else {
            // Differentiate N times by repeated coefficient shifting, which
            // costs N passes over the coefficients and no allocation beyond
            // the working vector.
            const auto degree = Degree();
            if (degree < static_cast<int>(N)) {
                return static_cast<T>(0);
            }
            std::vector<T> c{_a};
            for (std::size_t pass = 0; pass < N; ++pass) {
                for (std::size_t k = 1; k < c.size(); ++k) {
                    c[k - 1] = static_cast<T>(k) * c[k];
                }
                c.pop_back();
            }
            return std::accumulate(
                c.rbegin(), c.rend(), static_cast<T>(0),
                [x](auto acc, auto a) { return acc * x + a; });
        }
    }

    /**
     * @brief Evaluate the zero-constant antiderivative at `x`.
     * @return `sum(a[i] * x^(i+1) / (i+1))`.
     */
    T Primitive(T x) const {
        return std::accumulate(_a.rbegin(), _a.rend(), static_cast<T>(0),
                               [x, m = Degree() + 1](auto p, auto a) mutable {
                                   return p * x + a * x / static_cast<T>(m--);
                               });
    }

    /**
     * @brief The antiderivative vanishing at zero, spelled for Function1D.
     *
     * Named so that the free function Primitive() can find it through the
     * Antidifferentiable1D concept.
     * @param x Query abscissa.
     */
    T Antiderivative(T x) const { return Primitive(x); }

    /**
     * @brief Compatibility spelling for Primitive().
     * @deprecated Use Primitive().
     */
    [[deprecated("Use Primitive()")]] T Primative(T x) const {
        return Primitive(x);
    }

    /** @brief Return the definite integral over `[a, b]`. */
    T Integrate(T a, T b) const { return Primitive(b) - Primitive(a); }

    /** @brief Return a mutable iterator to the first coefficient. */
    auto begin() { return _a.begin(); }
    /** @brief Return a mutable iterator one past the final coefficient. */
    auto end() { return _a.end(); }

    /** @brief Return a constant iterator to the first coefficient. */
    auto cbegin() const { return _a.cbegin(); }
    /** @brief Return a constant iterator one past the final coefficient. */
    auto cend() const { return _a.cend(); }

    /**
     * @brief Return a coefficient value.
     * @param idx Zero-based power index.
     * @pre `idx` is within the stored coefficient range.
     */
    T operator[](int idx) { return _a[idx]; }

    /** @brief Return a copy of the coefficient vector. */
    std::vector<T> polycoeff() const { return _a; };

    /**
     * @brief Return one coefficient, padding powers outside the range with
     * zero.
     * @param i Power index.
     */
    T polycoeff(int i) const {
        if (i > this->Degree()) {
            return 0.0;
        } else if (i < 0) {
            return 0.0;
        } else {
            return this->_a[i];
        }
    };

    /** @brief Return the coefficient-wise additive inverse. */
    Polynomial<T> operator-() const {

        std::vector<T> myvec(this->Degree() + 1);
        for (int idx = 0; idx < this->Degree() + 1; ++idx) {
            myvec[idx] = -this->_a[idx];
        };
        Polynomial<T> myret{myvec};
        return myret;
    };

    /** @brief Add a scalar to the constant coefficient. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial<T> &operator+=(FLOAT b) {
        _a[0] += b;
        return *this;
    }

    /** @brief Subtract a scalar from the constant coefficient. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial<T> &operator-=(FLOAT b) {
        _a[0] -= b;
        return *this;
    }

    /** @brief Multiply every coefficient by a scalar. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial<T> &operator*=(FLOAT b) {
        for (int idx = 0; idx < this->Degree() + 1; ++idx) {
            _a[idx] *= b;
        };
        return *this;
    }

    /** @brief Divide every coefficient by a nonzero scalar. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial<T> &operator/=(FLOAT b) {
        for (int idx = 0; idx < this->Degree() + 1; ++idx) {
            _a[idx] /= b;
        };
        return *this;
    }

    /** @brief Add another polynomial coefficient-wise. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial<T> &operator+=(const Polynomial<FLOAT> &b) {

        bool sdeg = (this->Degree() < b.Degree());
        for (int idx = 0; idx < this->Degree() + 1; ++idx) {
            _a[idx] += b.polycoeff(idx);
        };
        if (sdeg) {
            for (int idx = this->Degree() + 1; idx < b.Degree() + 1; ++idx) {
                _a.push_back(b.polycoeff(idx));
            };
        };
        return *this;
    }

    /** @brief Subtract another polynomial coefficient-wise. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial<T> &operator-=(const Polynomial<FLOAT> &b) {

        bool sdeg = (this->Degree() < b.Degree());
        for (int idx = 0; idx < this->Degree() + 1; ++idx) {
            _a[idx] -= b.polycoeff(idx);
        };
        if (sdeg) {
            for (int idx = this->Degree() + 1; idx < b.Degree() + 1; ++idx) {
                _a.push_back(-b.polycoeff(idx));
            };
        };
        return *this;
    }

    /** @brief Multiply by another polynomial using coefficient convolution. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial<T> &operator*=(const Polynomial<FLOAT> &b) {
        int maxc = this->Degree() + b.Degree();

        std::vector<T> newcoeff(maxc + 1, 0.0);
        for (int idx = 0; idx < maxc + 1; ++idx) {
            for (int idx2 = 0; idx2 < idx + 1; ++idx2) {
                newcoeff[idx] +=
                    this->polycoeff(idx2) * b.polycoeff(idx - idx2);
            }
        };
        this->_a = std::move(newcoeff);
        return *this;
    }

    /**
     * @brief Write coefficients in ascending power order, separated by spaces.
     * @param os Destination stream.
     * @param obj Polynomial to write.
     * @return `os`.
     */
    friend std::ostream &operator<<(std::ostream &os,
                                    const Polynomial<T> &obj) {
        // Write obj to stream
        for (int idx = 0; idx < obj.Degree() + 1; ++idx) {
            os << obj.polycoeff(idx) << " ";
        }
        return os;
    };

    // ---------------------------------------------------------------------
    // Arithmetic.
    //
    // These are hidden friends: found by argument-dependent lookup on a
    // Polynomial operand and invisible otherwise, rather than templates at
    // global namespace scope that take part in every unrelated overload
    // resolution in the program.
    //
    // The result type is common_type_t of the two operands rather than the
    // left operand's type, so adding a complex scalar to a real polynomial
    // gives a complex polynomial instead of silently discarding the
    // imaginary part.
    // ---------------------------------------------------------------------

    /** @brief Add a scalar to the constant coefficient. */
    template <typename U>
        requires(RealOrComplex<U> || std::integral<U>)
    friend Polynomial<std::common_type_t<T, U>> operator+(const Polynomial &a,
                                                          const U &b) {
        Polynomial<std::common_type_t<T, U>> result{a};
        result += b;
        return result;
    }

    /** @brief Add a polynomial to a scalar. */
    template <typename U>
        requires(RealOrComplex<U> || std::integral<U>)
    friend Polynomial<std::common_type_t<T, U>> operator+(const U &b,
                                                          const Polynomial &a) {
        return a + b;
    }

    /** @brief Subtract a scalar from the constant coefficient. */
    template <typename U>
        requires(RealOrComplex<U> || std::integral<U>)
    friend Polynomial<std::common_type_t<T, U>> operator-(const Polynomial &a,
                                                          const U &b) {
        Polynomial<std::common_type_t<T, U>> result{a};
        result -= b;
        return result;
    }

    /** @brief Subtract a polynomial from a scalar. */
    template <typename U>
        requires(RealOrComplex<U> || std::integral<U>)
    friend Polynomial<std::common_type_t<T, U>> operator-(const U &b,
                                                          const Polynomial &a) {
        return -a + b;
    }

    /** @brief Multiply every coefficient by a scalar. */
    template <typename U>
        requires(RealOrComplex<U> || std::integral<U>)
    friend Polynomial<std::common_type_t<T, U>> operator*(const Polynomial &a,
                                                          const U &b) {
        Polynomial<std::common_type_t<T, U>> result{a};
        result *= b;
        return result;
    }

    /** @brief Multiply a scalar by a polynomial. */
    template <typename U>
        requires(RealOrComplex<U> || std::integral<U>)
    friend Polynomial<std::common_type_t<T, U>> operator*(const U &b,
                                                          const Polynomial &a) {
        return a * b;
    }

    /** @brief Divide every coefficient by a nonzero scalar. */
    template <typename U>
        requires(RealOrComplex<U> || std::integral<U>)
    friend Polynomial<std::common_type_t<T, U>> operator/(const Polynomial &a,
                                                          const U &b) {
        Polynomial<std::common_type_t<T, U>> result{a};
        result /= b;
        return result;
    }

    /** @brief Add two polynomials coefficient-wise, trimming the result. */
    template <typename U>
    friend Polynomial<std::common_type_t<T, U>>
    operator+(const Polynomial &a, const Polynomial<U> &b) {
        using R = std::common_type_t<T, U>;
        const auto degree = std::max(a.Degree(), b.Degree());
        std::vector<R> c(static_cast<std::size_t>(degree + 1), R{});
        for (int i = 0; i <= degree; ++i) {
            c[static_cast<std::size_t>(i)] =
                static_cast<R>(a.polycoeff(i)) + static_cast<R>(b.polycoeff(i));
        }
        return Polynomial<R>{std::move(c)}.Trim();
    }

    template <typename U>
    /** @brief Subtract two polynomials coefficient-wise, trimming the result.
     */
    friend Polynomial<std::common_type_t<T, U>>
    operator-(const Polynomial &a, const Polynomial<U> &b) {
        using R = std::common_type_t<T, U>;
        const auto degree = std::max(a.Degree(), b.Degree());
        std::vector<R> c(static_cast<std::size_t>(degree + 1), R{});
        for (int i = 0; i <= degree; ++i) {
            c[static_cast<std::size_t>(i)] =
                static_cast<R>(a.polycoeff(i)) - static_cast<R>(b.polycoeff(i));
        }
        return Polynomial<R>{std::move(c)}.Trim();
    }

    /** @brief Multiply two polynomials by coefficient convolution. */
    template <typename U>
    friend Polynomial<std::common_type_t<T, U>>
    operator*(const Polynomial &a, const Polynomial<U> &b) {
        using R = std::common_type_t<T, U>;
        if (a.Degree() < 0 || b.Degree() < 0) {
            return Polynomial<R>{};
        }
        const auto degree = a.Degree() + b.Degree();
        std::vector<R> c(static_cast<std::size_t>(degree + 1), R{});
        for (int i = 0; i <= degree; ++i) {
            for (int j = 0; j <= i; ++j) {
                c[static_cast<std::size_t>(i)] +=
                    static_cast<R>(a.polycoeff(j)) *
                    static_cast<R>(b.polycoeff(i - j));
            }
        }
        return Polynomial<R>{std::move(c)};
    }

  private:
    std::vector<T> _a; // Vector of polynomial coefficients.
};

} // namespace Interpolation

#endif // INTERPOLATION_POLYNOMIAL_HPP
