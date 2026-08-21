#include <gtest/gtest.h>

#include <Interpolation/AkimaSpline.hpp>
#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Interpolation.hpp>
#include <Interpolation/Lagrange.hpp>
#include <Interpolation/Linear.hpp>
#include <Interpolation/Polynomial.hpp>
#include <NumericConcepts/Ranges.hpp>
#include <complex>
#include <list>
#include <vector>

using RealIterator = std::vector<double>::iterator;
using ConstRealIterator = std::vector<double>::const_iterator;
using ComplexIterator = std::vector<std::complex<double>>::iterator;
using IntegerIterator = std::vector<int>::iterator;
using ListIterator = std::list<double>::iterator;

static_assert(Interpolation::Real<float>);
static_assert(Interpolation::Real<double>);
static_assert(!Interpolation::Real<int>);
static_assert(Interpolation::Complex<std::complex<double>>);
static_assert(!Interpolation::Complex<double>);
static_assert(Interpolation::RealOrComplex<long double>);
static_assert(Interpolation::RealOrComplex<std::complex<long double>>);

static_assert(Interpolation::RealIterator<RealIterator>);
static_assert(Interpolation::RealIterator<ConstRealIterator>);
static_assert(!Interpolation::RealIterator<ComplexIterator>);
static_assert(!Interpolation::RealIterator<IntegerIterator>);
static_assert(!Interpolation::RealIterator<ListIterator>);
static_assert(Interpolation::ComplexIterator<ComplexIterator>);
static_assert(
    Interpolation::InterpolationIteratorPair<RealIterator, RealIterator>);
static_assert(
    Interpolation::InterpolationIteratorPair<RealIterator, ComplexIterator>);
static_assert(!Interpolation::InterpolationIteratorPair<ComplexIterator,
                                                        ComplexIterator>);
static_assert(
    !Interpolation::InterpolationIteratorPair<IntegerIterator, RealIterator>);

// The range concepts are refinements, not adoptions: NumericConcepts requires
// only an input_range, while interpolation needs random access and a size.
static_assert(Interpolation::RealRange<std::vector<double>>);
static_assert(
    Interpolation::RealOrComplexRange<std::vector<std::complex<double>>>);
static_assert(!Interpolation::RealRange<std::vector<int>>);
static_assert(!Interpolation::RealRange<std::list<double>>);
static_assert(NumericConcepts::RealRange<std::list<double>>,
              "the underlying concept should accept a list; ours should not");

static_assert(Interpolation::InterpolationRanges<std::vector<double>,
                                                 std::vector<double>>);
static_assert(Interpolation::InterpolationRanges<
              std::vector<double>, std::vector<std::complex<double>>>);
static_assert(
    !Interpolation::InterpolationRanges<std::vector<std::complex<double>>,
                                        std::vector<std::complex<double>>>);
static_assert(
    !Interpolation::InterpolationRanges<std::vector<int>, std::vector<double>>);

TEST(ConceptsAndHeaders, PublicForwardingHeadersCompileTogether) { SUCCEED(); }
