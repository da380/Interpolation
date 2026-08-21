
#include <Interpolation/Interpolation.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

int main() {
  using namespace Interpolation;

  // Set the nodes
  int N = 3;
  double x1 = 0.0;
  double x2 = 1.0;
  double dx = (x2 - x1) / static_cast<double>(N - 1);
  std::vector<double> X(N);
  for (int i = 0; i < N; i++) {
    X[i] = x1 + i * dx;
  }

  // Set the Lagrange polynomial
  auto p = LagrangePolynomial(X.begin(), X.end());

  // Set the values for plotting
  int n = 100;
  dx = (x2 - x1) / static_cast<double>(n - 1);

  auto file = std::ofstream("Lagrange.out");
  for (int i = 0; i < n; i++) {
    double x = x1 + i * dx;
    int j = 1;
    file << x << " " << p(j, x) << " " << p.Derivative(j, x) << std::endl;
  }
}
