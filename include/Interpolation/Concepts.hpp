#ifndef INTERPOLATION_CONCEPTS_HPP
#define INTERPOLATION_CONCEPTS_HPP

#include <iterator>
#include <ranges>

#include <NumericConcepts/Iterators.hpp>
#include <NumericConcepts/Numeric.hpp>
#include <NumericConcepts/Ranges.hpp>

/**
 * @brief Interpolation algorithms and supporting mathematical utilities.
 */
namespace Interpolation {

/**
 * @brief A real floating-point scalar, used for abscissae.
 */
using NumericConcepts::Real;

/**
 * @brief A `std::complex` with floating-point components.
 */
using NumericConcepts::Complex;

/**
 * @brief A supported ordinate type: real or complex.
 */
using NumericConcepts::RealOrComplex;

/**
 * @brief The underlying real precision of a real or complex type.
 */
using NumericConcepts::RemoveComplex;

/**
 * @brief A range of abscissae.
 *
 * This refines `NumericConcepts::RealRange`, which requires only an
 * `input_range`. Interpolation locates a query by binary search and needs the
 * node count up front, so random access and a known size are both required.
 *
 * @tparam T Range type to check.
 */
template <typename T>
concept RealRange = requires() {
    requires NumericConcepts::RealRange<T>;
    requires std::ranges::random_access_range<T>;
    requires std::ranges::sized_range<T>;
};

/**
 * @brief A range of ordinates, real or complex, with the same refinement.
 * @tparam T Range type to check.
 */
template <typename T>
concept RealOrComplexRange = requires() {
    requires NumericConcepts::RealOrComplexRange<T>;
    requires std::ranges::random_access_range<T>;
    requires std::ranges::sized_range<T>;
};

/**
 * @brief Compatible abscissa and ordinate ranges for interpolation.
 *
 * The abscissae must be real. The ordinates may be real or complex, and the
 * two scalar types must support the arithmetic the interpolators perform.
 *
 * @tparam X Abscissa range type.
 * @tparam Y Ordinate range type.
 */
template <typename X, typename Y>
concept InterpolationRanges =
    requires(std::ranges::range_value_t<X> x, std::ranges::range_value_t<Y> y) {
        requires RealRange<X>;
        requires RealOrComplexRange<Y>;
        requires std::convertible_to<std::ranges::range_value_t<X>,
                                     std::ranges::range_value_t<Y>>;
        { x + y } -> std::convertible_to<std::ranges::range_value_t<Y>>;
        { x *y } -> std::convertible_to<std::ranges::range_value_t<Y>>;
        { y / x } -> std::convertible_to<std::ranges::range_value_t<Y>>;
    };

/** @brief A random-access iterator over real abscissae. */
template <typename T>
concept RealIterator = requires() {
    requires std::random_access_iterator<T>;
    requires Real<std::iter_value_t<T>>;
};

/** @brief A random-access iterator over complex ordinates. */
template <typename T>
concept ComplexIterator = requires() {
    requires std::random_access_iterator<T>;
    requires Complex<std::iter_value_t<T>>;
};

/** @brief A random-access iterator over real or complex ordinates. */
template <typename T>
concept RealOrComplexIterator = requires() {
    requires std::random_access_iterator<T>;
    requires RealOrComplex<std::iter_value_t<T>>;
};

/**
 * @brief Compatible abscissa and ordinate iterators for interpolation.
 *
 * @tparam xIter Random-access iterator type for abscissae.
 * @tparam yIter Random-access iterator type for ordinates.
 */
template <typename xIter, typename yIter>
concept InterpolationIteratorPair = requires(xIter x, yIter y) {
    requires RealIterator<xIter>;
    requires RealOrComplexIterator<yIter>;
    requires std::convertible_to<std::iter_value_t<xIter>,
                                 std::iter_value_t<yIter>>;
    { (*x) + (*y) } -> std::convertible_to<std::iter_value_t<yIter>>;
    { (*x) * (*y) } -> std::convertible_to<std::iter_value_t<yIter>>;
    { (*y) / (*x) } -> std::convertible_to<std::iter_value_t<yIter>>;
};

} // namespace Interpolation

#endif // INTERPOLATION_CONCEPTS_HPP
