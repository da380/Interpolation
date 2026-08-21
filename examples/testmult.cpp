// #include <PlanetaryModel/All>
#include <Interpolation/Interpolation.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
// #include "PREM.h"

// class DummyClass {};

// template <typename T>
// requires PlanetaryModel::SphericalGeometryModel<T>
// void Dummy(T t) {
//   t.UpperRadius(2);
// }

int
main() {
    //   auto myprem = PREM();
    //   // //     for (int i = 0; i < 13; ++i){
    //   // // std::cout << myprem.Density(i)(0) << " " << myprem.VPH(i)(0) << "
    //   "
    //   <<
    //   // myprem.A(i)(0) << std::endl;
    //   // //     }
    //   //     std::cout << myprem.DensityNorm() << " "  << std::endl;
    //   static_assert(
    //       PlanetaryModel::HasNormalisationInformation<EarthConstants<double>
    //       >);
    //   static_assert(PlanetaryModel::SphericalGeometryModel<PREM<double, int>
    //   >); static_assert(PlanetaryModel::SphericalDensityModel<PREM<double,
    //   int>
    //   >); static_assert(PlanetaryModel::SphericalElasticModel<PREM<double,
    //   int>
    //   >);

    // polynomials:
    // Interpolation::Polynomial<double> vecpoly{1.0, 2.0};
    // Interpolation::Polynomial<double> vecpoly2 = vecpoly * 2.0;
    // vecpoly2 += vecpoly;
    // vecpoly2 +=1.0;
    // std::cout << vecpoly2.polycoeff(0) << vecpoly2.polycoeff(1) << std::endl;
    // std::cout << vecpoly2 << std::endl;

    Interpolation::Polynomial<double> vecpoly{1.0, 2.0};
    Interpolation::Polynomial<double> vecpoly3{1.0, 2.0};
    Interpolation::Polynomial<double> vecpolyadd = 2 + vecpoly;
    Interpolation::Polynomial<double> vecpolysub = 2.0 - vecpoly;
    Interpolation::Polynomial<double> vecpolymult = 2.0 * vecpoly;
    Interpolation::Polynomial<double> vecpolydiv = vecpoly / 2.0;
    Interpolation::Polynomial<double> vecpoly2{1.0, 2.0, 3.0};

    std::cout << "Initial: " << vecpoly << std::endl;
    std::cout << "Addition of 2: " << vecpolyadd << std::endl;
    std::cout << "Subtraction by 2: " << vecpolysub << std::endl;
    std::cout << "Multiplication by 2: " << vecpolymult << std::endl;
    std::cout << "Division by 2: " << vecpolydiv << std::endl << std::endl;

    std::cout << "Initial: " << vecpoly3 << std::endl;
    vecpoly3 += vecpoly2;
    std::cout << "Added: " << vecpoly2 << std::endl;
    std::cout << "Final: " << vecpoly3 << std::endl;
    vecpoly3 += vecpoly;
    std::cout << "Final: " << vecpoly3 << std::endl;
    vecpoly3 += 1.0;
    std::cout << "Final: " << vecpoly3 << std::endl;
    vecpoly3 -= 1.0;
    std::cout << "Final: " << vecpoly3 << std::endl;
    vecpoly3 *= 2.0;
    std::cout << "Final: " << vecpoly3 << std::endl;
    vecpoly3 /= 2.0;
    std::cout << "Final: " << vecpoly3 << std::endl;
    vecpoly3 += vecpoly;
    std::cout << "Final: " << vecpoly3 << std::endl;
    vecpoly3 -= vecpoly;
    std::cout << "Final: " << vecpoly3 << std::endl;

    std::cout << std::endl;
    std::cout << "Initial: " << vecpoly3 << std::endl;
    std::cout << "Multiplier: " << vecpoly << std::endl;
    vecpoly3 *= vecpoly;
    std::cout << "Final: " << vecpoly3 << std::endl;

    std::cout << std::endl;
    std::cout << "Initial: " << vecpoly << std::endl;
    std::cout << "Multiplier: " << vecpoly3 << std::endl;
    vecpoly *= vecpoly3;
    std::cout << "Final: " << vecpoly << std::endl;
    std::cout << std::endl;
    std::cout << "Initial: " << vecpoly << std::endl;
    std::cout << "Multiplier: " << vecpoly3 << std::endl;
    Interpolation::Polynomial<float> vecpoly4 = vecpoly3 * vecpoly;
    std::cout << "Final: " << vecpoly4 << std::endl;

    std::cout << "Full output: " << vecpoly << std::endl;
    std::cout << "Output []: " << vecpoly[0] << std::endl;
    //   vecpoly2 = vecpoly * 2.0;
}