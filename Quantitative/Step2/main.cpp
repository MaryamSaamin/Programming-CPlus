#include <iostream>
#include <iomanip>
#include "bs_call_price.h"

int main() {
    // --- Parameters ---
    double S     = 100.0;
    double K     = 100.0;
    double r     = 0.05;
    double q     = 0.02;
    double sigma = 0.2;
    double T     = 1.0;

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "=== Analytic ===\n";
    double delta_analytic = bs_delta_call(S, K, r, q, sigma, T);
    double gamma_analytic = bs_gamma(S, K, r, q, sigma, T);
    std::cout << "Delta: " << delta_analytic << "\n";
    std::cout << "Gamma: " << gamma_analytic << "\n\n";

    std::cout << "=== Forward Difference ===\n";
    double delta_fwd_val = delta_fwd(S, K, r, q, sigma, T, 1e-3);
    double gamma_fwd_val = gamma_fwd(S, K, r, q, sigma, T, 1e-2);
    std::cout << "Delta: " << delta_fwd_val << "\n";
    std::cout << "Gamma: " << gamma_fwd_val << "\n\n";

    std::cout << "=== Complex-Step (Pure Imaginary) ===\n";
    double delta_cs_val = delta_cs(S, K, r, q, sigma, T, 1e-8);
    double gamma_cs_val = gamma_cs(S, K, r, q, sigma, T, 1e-8);
    std::cout << "Delta: " << delta_cs_val << "\n";
    std::cout << "Gamma: " << gamma_cs_val << "\n\n";

    std::cout << "=== Complex-Step (45-degree ω = (1+i)/√2) ===\n";
    double gamma_cs45_val = gamma_cs45(S, K, r, q, sigma, T, 1e-6);
    std::cout << "Gamma (O(h^4)): " << gamma_cs45_val << "\n\n";

    // --- Compare errors ---
    std::cout << "=== Relative Errors (vs Analytic) ===\n";
    std::cout << "Delta (fwd): " << std::abs((delta_fwd_val - delta_analytic) / delta_analytic) << "\n";
    std::cout << "Gamma (fwd): " << std::abs((gamma_fwd_val - gamma_analytic) / gamma_analytic) << "\n";
    std::cout << "Delta (cs):  " << std::abs((delta_cs_val  - delta_analytic) / delta_analytic) << "\n";
    std::cout << "Gamma (cs):  " << std::abs((gamma_cs_val  - gamma_analytic) / gamma_analytic) << "\n";
    std::cout << "Gamma (cs45):" << std::abs((gamma_cs45_val - gamma_analytic) / gamma_analytic) << "\n";

     // --- Test Phi_t with small imaginary step ---
    std::cout << "=== Test Phi_t (complex branch) ===\n";
    cplx z_test(S, 1e-8);  // S + i h
    cplx Phi_z = Phi_t(z_test);
    std::cout << "Phi_t(S + i h) = " << Phi_z << " (imag part = h * phi(S))\n";

    return 0;
}
