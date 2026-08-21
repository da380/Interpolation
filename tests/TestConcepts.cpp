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

TEST(ConceptsAndHeaders, PublicForwardingHeadersCompileTogether) { SUCCEED(); }
