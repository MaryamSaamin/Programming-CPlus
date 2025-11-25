// main.cpp -- simple tests and small timing for InverseCumulativeNormal
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include "InverseCumulativeNormal.h"

int main() {
    quant::InverseCumulativeNormal icn; // mean=0, sigma=1

    std::cout << std::setprecision(10);

    // Scalar test points (including extreme tails)
    double xs[] = {1e-12, 1e-6, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0 - 1e-6, 1.0};
    std::cout << "Scalar tests:\n";
    for (double x : xs) {
        double z = icn(x);
        std::cout << "  x = " << std::setw(10) << std::defaultfloat << x
                  << "  (hex " << std::hexfloat << x << std::defaultfloat << ")"
                  << "  -> z = " << std::setw(14) << z << '\n';
    }

    // Vector test
    const double xin[] = {0.0001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.9999};
    double zout[std::size(xin)];

    icn(xin, zout, std::size(xin)); // run once to warm caches
    std::cout << "\nVector results:\n";
    for (std::size_t i = 0; i < std::size(xin); ++i) {
        std::cout << "  x = " << xin[i] << "  -> z = " << zout[i] << '\n';
    }

    // Symmetry check for a few values
    std::cout << "\nSymmetry checks (g(1-x) + g(x)):\n";
    double sym_xs[] = {1e-6, 0.01, 0.1, 0.3, 0.45};
    for (double x : sym_xs) {
        double a = icn(x);
        double b = icn(1.0 - x);
        std::cout << "  x=" << x << "  g(x)=" << a << "  g(1-x)=" << b
                  << "  sum=" << (a + b) << '\n';
    }

    // Quick timing for vector path (repeat to get measurable time)
    const std::size_t N = 1000000; // 1e6 elements
    std::vector<double> large_in(N, 0.123456); // arbitrary non-edge value
    std::vector<double> large_out(N, 0.0);

    const int repeats = 5;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < repeats; ++r) {
        icn(large_in.data(), large_out.data(), N);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "\nVector path timing: processed " << (N * repeats)
              << " elements in " << ms << " ms  (avg " << (ms / (N * repeats)) << " ms/elem)\n";

    return 0;
}
