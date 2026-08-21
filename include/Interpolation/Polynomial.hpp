#ifndef INTERPOLATION_POLYNOMIAL_HPP
#define INTERPOLATION_POLYNOMIAL_HPP

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
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
requires RealOrComplexFloatingPoint<T>
class Polynomial1D {
  public:
    /** @brief Scalar type used for coefficients and evaluation. */
    using value_type = T;
    /** @brief Mutable coefficient iterator. */
    using iterator = std::vector<T>::iterator;
    /** @brief Constant coefficient iterator. */
    using const_iterator = std::vector<T>::const_iterator;

    /** @brief Construct the zero polynomial with one coefficient equal to zero. */
    Polynomial1D() : _a{T{}} {}

    /** @brief Copy coefficients from a vector in ascending power order. */
    Polynomial1D(const std::vector<T> &a) : _a{a} {}
    /** @brief Move coefficients from a vector in ascending power order. */
    Polynomial1D(std::vector<T> &&a) : _a{std::move(a)} {}

    /** @brief Construct from coefficients in ascending power order. */
    Polynomial1D(std::initializer_list<T> list) : _a{std::vector<T>{list}} {}

    /** @brief Copy a polynomial. */
    Polynomial1D(const Polynomial1D &) = default;
    /** @brief Move a polynomial. */
    Polynomial1D(Polynomial1D &&) = default;

    /** @brief Copy-assign a polynomial with the same coefficient type. */
    Polynomial1D &operator=(const Polynomial1D &) = default;
    /** @brief Move-assign a polynomial with the same coefficient type. */
    Polynomial1D &operator=(Polynomial1D &&) = default;

    /**
     * @brief Copy a polynomial while converting its coefficient type.
     * @tparam FLOAT Source coefficient type.
     * @param rhs Source polynomial.
     */
    template <typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
    Polynomial1D(const Polynomial1D<FLOAT> &rhs) {
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
        requires std::is_convertible_v<FLOAT, T> Polynomial1D<T>
    &operator=(const Polynomial1D<FLOAT> &polinit) {
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
     * @return Random polynomial of degree `n`.
     */
    static Polynomial1D Random(int n)
        requires RealFloatingPoint<T>
    {
        std::random_device rd{};
        std::mt19937_64 gen{rd()};
        std::normal_distribution<T> d{};
        std::vector<T> a;
        std::generate_n(std::back_inserter(a), n + 1, [&]() {
            if constexpr (RealFloatingPoint<T>) {
                return d(gen);
            }
        });
        return Polynomial1D(a);
    }

    /**
     * @brief Generate a complex polynomial with normally distributed components.
     * @param n Requested nonnegative degree.
     * @return Random polynomial of degree `n`.
     */
    static Polynomial1D Random(int n)
        requires ComplexFloatingPoint<T>
    {
        using S = typename T::value_type;
        std::random_device rd{};
        std::mt19937_64 gen{rd()};
        std::normal_distribution<S> d{};
        std::vector<T> a;
        std::generate_n(std::back_inserter(a), n + 1, [&]() {
            return T{d(gen), d(gen)};
        });
        return Polynomial1D(a);
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
     * @brief Return one coefficient, padding powers outside the range with zero.
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
    Polynomial1D<T> operator-() const {

        std::vector<T> myvec(this->Degree() + 1);
        for (int idx = 0; idx < this->Degree() + 1; ++idx) {
            myvec[idx] = -this->_a[idx];
        };
        Polynomial1D<T> myret{myvec};
        return myret;
    };

    /** @brief Add a scalar to the constant coefficient. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial1D<T> &operator+=(FLOAT b) {
        _a[0] += b;
        return *this;
    }

    /** @brief Subtract a scalar from the constant coefficient. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial1D<T> &operator-=(FLOAT b) {
        _a[0] -= b;
        return *this;
    }

    /** @brief Multiply every coefficient by a scalar. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial1D<T> &operator*=(FLOAT b) {
        for (int idx = 0; idx < this->Degree() + 1; ++idx) {
            _a[idx] *= b;
        };
        return *this;
    }

    /** @brief Divide every coefficient by a nonzero scalar. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial1D<T> &operator/=(FLOAT b) {
        for (int idx = 0; idx < this->Degree() + 1; ++idx) {
            _a[idx] /= b;
        };
        return *this;
    }

    /** @brief Add another polynomial coefficient-wise. */
    template <typename FLOAT>
        requires std::is_convertible_v<FLOAT, T>
    Polynomial1D<T> &operator+=(const Polynomial1D<FLOAT> &b) {

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
    Polynomial1D<T> &operator-=(const Polynomial1D<FLOAT> &b) {

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
    Polynomial1D<T> &operator*=(const Polynomial1D<FLOAT> &b) {
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
                                    const Polynomial1D<T> &obj) {
        // Write obj to stream
        for (int idx = 0; idx < obj.Degree() + 1; ++idx) {
            os << obj.polycoeff(idx) << " ";
        }
        return os;
    };

  private:
    std::vector<T> _a;   // Vector of polynomial coefficients.
};

}   // namespace Interpolation

/**
 * @brief Add a scalar to a polynomial's constant coefficient.
 * @param a Polynomial operand.
 * @param b Scalar operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator+(Interpolation::Polynomial1D<T> a, FLOAT b) {
    Interpolation::Polynomial1D<T> myval = a;
    myval += b;
    return myval;
};

/**
 * @brief Add a scalar to a polynomial's constant coefficient.
 * @param b Scalar operand.
 * @param a Polynomial operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator+(FLOAT b, Interpolation::Polynomial1D<T> a) {
    Interpolation::Polynomial1D<T> myval = a;
    myval += b;
    return myval;
};

/**
 * @brief Subtract a scalar from a polynomial's constant coefficient.
 * @param a Polynomial operand.
 * @param b Scalar operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator-(Interpolation::Polynomial1D<T> a, FLOAT b) {
    Interpolation::Polynomial1D<T> myval = a;
    myval -= b;
    return myval;
};

/**
 * @brief Subtract a polynomial from a scalar constant polynomial.
 * @param b Scalar operand.
 * @param a Polynomial operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator-(FLOAT b, Interpolation::Polynomial1D<T> a) {
    Interpolation::Polynomial1D<T> myval = -a;
    myval += b;
    return myval;
};

/**
 * @brief Multiply every polynomial coefficient by a scalar.
 * @param a Polynomial operand.
 * @param b Scalar operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator*(Interpolation::Polynomial1D<T> a, FLOAT b) {
    Interpolation::Polynomial1D<T> myval = a;
    myval *= b;
    return myval;
};

/**
 * @brief Multiply every polynomial coefficient by a scalar.
 * @param b Scalar operand.
 * @param a Polynomial operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator*(FLOAT b, Interpolation::Polynomial1D<T> a) {
    Interpolation::Polynomial1D<T> myval = a;
    myval *= b;
    return myval;
};

/**
 * @brief Divide every polynomial coefficient by a nonzero scalar.
 * @param a Polynomial operand.
 * @param b Scalar operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator/(Interpolation::Polynomial1D<T> a, FLOAT b) {
    Interpolation::Polynomial1D<T> myval = a;
    myval /= b;
    return myval;
};

/**
 * @brief Add two polynomials coefficient-wise.
 * @param a Left polynomial operand.
 * @param b Right polynomial operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator+(const Interpolation::Polynomial1D<T> &a,
          const Interpolation::Polynomial1D<FLOAT> &b) {
    int maxval = std::max(a.Degree(), b.Degree());
    std::vector<T> myvec(maxval + 1, 0.0);
    for (int idx = 0; idx < maxval + 1; ++idx) {
        myvec[idx] = a.polycoeff(idx) + b.polycoeff(idx);
    };
    Interpolation::Polynomial1D<T> myval{myvec};
    return myval;
};

/**
 * @brief Subtract two polynomials coefficient-wise.
 * @param a Left polynomial operand.
 * @param b Right polynomial operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator-(const Interpolation::Polynomial1D<T> &a,
          const Interpolation::Polynomial1D<FLOAT> &b) {
    int maxval = std::max(a.Degree(), b.Degree());
    std::vector<T> myvec(maxval + 1, 0.0);
    for (int idx = 0; idx < maxval + 1; ++idx) {
        myvec[idx] = a.polycoeff(idx) - b.polycoeff(idx);
    };
    Interpolation::Polynomial1D<T> myval{myvec};
    return myval;
};

/**
 * @brief Multiply two polynomials by coefficient convolution.
 * @param a Left polynomial operand.
 * @param b Right polynomial operand.
 */
template <typename T, typename FLOAT>
    requires std::is_convertible_v<FLOAT, T>
Interpolation::Polynomial1D<T>
operator*(const Interpolation::Polynomial1D<T> &a,
          const Interpolation::Polynomial1D<FLOAT> &b) {
    int maxc = a.Degree() + b.Degree();
    std::vector<T> newcoeff(maxc + 1, 0.0);
    for (int idx = 0; idx < maxc + 1; ++idx) {
        for (int idx2 = 0; idx2 < idx + 1; ++idx2) {
            newcoeff[idx] += a.polycoeff(idx2) * b.polycoeff(idx - idx2);
        }
    };
    Interpolation::Polynomial1D<T> myval{newcoeff};
    return myval;
};

#endif   // INTERPOLATION_POLYNOMIAL_HPP
