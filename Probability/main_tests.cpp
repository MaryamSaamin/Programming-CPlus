// main_tests.cpp
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "InverseCumulativeNormal.h"

// Baseline bisection inversion (copy from your header if not public)
static double baseline_bisect(double x) {
    double lo = -12.0, hi = 12.0;
    if (x < 0.5) hi = 0.0; else lo = 0.0;
    for (int iter = 0; iter < 80; ++iter) {
        double mid = 0.5 * (lo + hi);
        double cdf = 0.5 * std::erfc(-mid * 0.7071067811865475244);
        if (cdf < x) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

int main() {
    using clk = std::chrono::high_resolution_clock;
    quant::InverseCumulativeNormal icn; // default mean=0,sigma=1

    std::cout << std::setprecision(12);

    // ----------------- Build grid -----------------
    std::vector<double> xs;
    // central dense linear grid
    const double xmin = 1e-12, xmax = 1.0 - 1e-12;
    const int Ncenter = 200000; // adjust for resolution
    for (int i = 0; i < Ncenter; ++i) {
        double t = double(i) / double(Ncenter-1);
        xs.push_back(1e-12 + t*(1 - 2e-12));
    }
    // add log-spaced tails for extra sampling
    for (int k = 1; k <= 16; ++k) {
        double v = std::pow(10.0, -k);
        xs.push_back(v);
        xs.push_back(1.0 - v);
    }
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());

    // ----------------- Correctness: round-trip -----------------
    std::vector<double> errs; errs.reserve(xs.size());
    for (double x : xs) {
        double z = icn(x);
        double fx = 0.5 * std::erfc(-z * 0.7071067811865475244); // Phi(z)
        errs.push_back(std::fabs(fx - x));
    }
    double maxe = *std::max_element(errs.begin(), errs.end());
    double meane = std::accumulate(errs.begin(), errs.end(), 0.0) / errs.size();
    std::sort(errs.begin(), errs.end());
    double p999 = errs[std::min<size_t>(errs.size()-1, size_t(errs.size()*0.999))];

    std::cout << "Round-trip (Phi(icn(x)) - x):\n";
    std::cout << "  max abs error = " << maxe << "\n";
    std::cout << "  mean abs error = " << meane << "\n";
    std::cout << "  99.9%ile abs error = " << p999 << "\n\n";

    // ----------------- Symmetry & monotonicity -----------------
    double max_sym = 0.0;
    for (double x : {1e-12,1e-9,1e-6,1e-4,1e-3,0.01,0.1,0.3,0.45}) {
        double a = icn(x);
        double b = icn(1.0 - x);
        max_sym = std::max(max_sym, std::fabs(a + b));
    }
    std::cout << "Max symmetry violation (g(1-x)+g(x)) sample = " << max_sym << "\n";

    // monotonicity check (strictly increasing)
    double min_diff = std::numeric_limits<double>::infinity();
    for (size_t i = 1; i < xs.size(); ++i) {
        double d = icn(xs[i]) - icn(xs[i-1]);
        if (d <= 0.0) { std::cout << "Monotonicity FAILED at x="<<xs[i-1]<<'\n'; break; }
        min_diff = std::min(min_diff, d);
    }
    std::cout << "Min positive step in g(x) (sampled grid) = " << min_diff << "\n\n";

    // ----------------- Derivative sanity -----------------
    // numeric derivative: central difference on modest grid
    double max_derr = 0.0;
    for (double x : {1e-12,1e-9,1e-6,1e-4,1e-3,0.01,0.1,0.3,0.5,0.7,0.9,0.99,0.9999}) {
        double h = std::max(1e-12, x*1e-6);
        double z1 = icn(std::max(x-h, 1e-300));
        double z2 = icn(std::min(x+h, 1.0-1e-300));
        double dnum = (z2 - z1) / ( (std::min(x+h,1.0-1e-300)) - (std::max(x-h,1e-300)) );
        double z = icn(x);
        double analytic = 1.0 / (0.398942280401432677939946 * std::exp(-0.5*z*z)); // 1/phi(z)
        double rel_err = std::fabs(dnum - analytic) / std::max(analytic, 1e-300);
        max_derr = std::max(max_derr, rel_err);
    }
    std::cout << "Max relative derivative error (sampled) = " << max_derr << "\n\n";

    // ----------------- Performance: scalar throughput -----------------
    // Warmup
    const size_t N = 2000000; // 2e6
    std::vector<double> xs_perf; xs_perf.reserve(N);
    for (size_t i=0;i<N;++i) xs_perf.push_back(0.123456 + 1e-6 * (i%100)); // non-edge
    double sink=0.0;
    // warm
    for (int r=0;r<3;++r) for (size_t i=0;i<N;++i) sink += icn(xs_perf[i]);
    // timed runs
    const int runs = 5;
    auto t0 = clk::now();
    for (int r=0;r<runs;++r) for (size_t i=0;i<N;++i) sink += icn(xs_perf[i]);
    auto t1 = clk::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double ns_per = (ms * 1e6) / (double(N) * runs);
    std::cout << "Scalar throughput: " << ns_per << " ns/eval (icn)\n";

    // Baseline bisection speed (same input)
    // warm
    for (int r=0;r<1;++r) for (size_t i=0;i<200000;++i) sink += baseline_bisect(xs_perf[i]);
    auto tb0 = clk::now();
    for (int r=0;r<1;++r) for (size_t i=0;i<200000;++i) sink += baseline_bisect(xs_perf[i]);
    auto tb1 = clk::now();
    double msb = std::chrono::duration<double, std::milli>(tb1 - tb0).count();
    double ns_b = (msb * 1e6) / 200000.0;
    std::cout << "Baseline bisection: " << ns_b << " ns/eval (small sample 2e5)\n";
    std::cout << "Speedup over baseline (rough) = " << (ns_b / ns_per) << "x\n\n";

    // ----------------- Vector path timing -----------------
    const size_t M = 1000000;
    std::vector<double> vin(M, 0.123456);
    std::vector<double> vout(M, 0.0);
    // warm
    icn(vin.data(), vout.data(), M);
    // timed
    auto tv0 = clk::now();
    icn(vin.data(), vout.data(), M);
    auto tv1 = clk::now();
    double msv = std::chrono::duration<double, std::milli>(tv1 - tv0).count();
    std::cout << "Vector path: " << msv << " ms for " << M << " elems  => " << (msv*1e6/M) << " ns/elem\n";

    // naive single-thread loop timing (no vectorization)
    std::vector<double> vout_naive(M, 0.0);
    auto tna0 = clk::now();
    for (size_t i=0;i<M;++i) vout_naive[i] = icn(vin[i]);
    auto tna1 = clk::now();
    double msna = std::chrono::duration<double, std::milli>(tna1 - tna0).count();
    std::cout << "Naive loop: " << msna << " ms => " << (msna*1e6/M) << " ns/elem\n";
    std::cout << "Vector speedup = " << (msna / msv) << "x\n\n";

    // Avoid optimizing out sink
    if (sink == 0.123456789) std::cout << "sink\n";

    return 0;
}
