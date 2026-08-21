#include <gtest/gtest.h>

#include <Interpolation/AkimaSpline.hpp>
#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Interpolation.hpp>
#include <Interpolation/Lagrange.hpp>
#include <Interpolation/Linear.hpp>
#include <Interpolation/Polynomial.hpp>
#include <complex>
#include <list>
#include <vector>

using RealIterator = std::vector<double>::iterator;
using ConstRealIterator = std::vector<double>::const_iterator;
using ComplexIterator = std::vector<std::complex<double>>::iterator;
using IntegerIterator = std::vector<int>::iterator;
using ListIterator = std::list<double>::iterator;

static_assert(Interpolation::RealFloatingPoint<float>);
static_assert(Interpolation::RealFloatingPoint<double>);
static_assert(!Interpolation::RealFloatingPoint<int>);
static_assert(Interpolation::ComplexFloatingPoint<std::complex<double>>);
static_assert(!Interpolation::ComplexFloatingPoint<double>);
static_assert(Interpolation::RealOrComplexFloatingPoint<long double>);
static_assert(
    Interpolation::RealOrComplexFloatingPoint<std::complex<long double>>);

static_assert(Interpolation::RealFloatingPointIterator<RealIterator>);
static_assert(Interpolation::RealFloatingPointIterator<ConstRealIterator>);
static_assert(!Interpolation::RealFloatingPointIterator<ComplexIterator>);
static_assert(!Interpolation::RealFloatingPointIterator<IntegerIterator>);
static_assert(!Interpolation::RealFloatingPointIterator<ListIterator>);
static_assert(Interpolation::ComplexFloatingPointIterator<ComplexIterator>);
static_assert(
    Interpolation::InterpolationIteratorPair<RealIterator, RealIterator>);
static_assert(
    Interpolation::InterpolationIteratorPair<RealIterator, ComplexIterator>);
static_assert(!Interpolation::InterpolationIteratorPair<ComplexIterator,
                                                        ComplexIterator>);
static_assert(
    !Interpolation::InterpolationIteratorPair<IntegerIterator, RealIterator>);

TEST(ConceptsAndHeaders, PublicForwardingHeadersCompileTogether) { SUCCEED(); }
