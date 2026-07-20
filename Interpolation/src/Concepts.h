#ifndef INTERPOLATION_CONCEPTS_GUARD_H
#define INTERPOLATION_CONCEPTS_GUARD_H

#include <complex>
#include <concepts>
#include <iterator>
#include <type_traits>

/**
 * @brief Interpolation algorithms and supporting mathematical utilities.
 */
namespace Interpolation {

/**
 * @brief Identifies complex numbers whose component type is floating point.
 * @tparam T Type to inspect.
 */
template <typename T> struct IsComplexFloatingPoint : public std::false_type {};

/**
 * @brief Specialization for `std::complex` values.
 * @tparam T Component type of the complex value.
 */
template <typename T>
struct IsComplexFloatingPoint<std::complex<T>>
    : public std::bool_constant<std::is_floating_point_v<T>> {};

/** @brief A real floating-point scalar type. */
template <typename T>
concept RealFloatingPoint = std::floating_point<T>;

/** @brief A `std::complex` type with floating-point components. */
template <typename T>
concept ComplexFloatingPoint =
    IsComplexFloatingPoint<std::remove_const_t<T>>::value;

/** @brief A supported real or complex interpolation ordinate type. */
template <typename T>
concept RealOrComplexFloatingPoint =
    RealFloatingPoint<T> or ComplexFloatingPoint<T>;

/** @brief A random-access iterator over real floating-point values. */
template <typename T>
concept RealFloatingPointIterator = requires() {
    requires std::random_access_iterator<T>;
    requires RealFloatingPoint<std::iter_value_t<T>>;
};

/** @brief A random-access iterator over complex floating-point values. */
template <typename T>
concept ComplexFloatingPointIterator = requires() {
    requires std::random_access_iterator<T>;
    requires ComplexFloatingPoint<std::iter_value_t<T>>;
};

/** @brief A random-access iterator over supported real or complex values. */
template <typename T>
concept RealOrComplexFloatingPointIterator = requires() {
    requires std::random_access_iterator<T>;
    requires RealOrComplexFloatingPoint<std::iter_value_t<T>>;
};

/**
 * @brief Compatible abscissa and ordinate iterators for interpolation.
 *
 * The abscissa iterator must contain real floating-point values. The ordinate
 * iterator may contain real or complex floating-point values, and the two
 * scalar types must support the arithmetic required by the interpolators.
 *
 * @tparam xIter Random-access iterator type for abscissae.
 * @tparam yIter Random-access iterator type for ordinates.
 */
template <typename xIter, typename yIter>
concept InterpolationIteratorPair = requires(xIter x, yIter y) {
    requires RealFloatingPointIterator<xIter>;
    requires RealOrComplexFloatingPointIterator<yIter>;
    requires std::convertible_to<std::iter_value_t<xIter>,
                                 std::iter_value_t<yIter>>;
    { (*x) + (*y) } -> std::convertible_to<std::iter_value_t<yIter>>;
    { (*x) * (*y) } -> std::convertible_to<std::iter_value_t<yIter>>;
    { (*y) / (*x) } -> std::convertible_to<std::iter_value_t<yIter>>;
};

}   // namespace Interpolation

#endif   //  INTERPOLATION_CONCEPTS_GUARD_H
