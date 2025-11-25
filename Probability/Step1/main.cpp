#include <iostream>
#include "bs_call_price.h"

int main() {
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;

    std::cout << "=== Analytic ===\n";
    std::cout << "Delta: " << bs_delta_call(S,K,r,q,sigma,T) << "\n";
    std::cout << "Gamma: " << bs_gamma(S,K,r,q,sigma,T) << "\n";

    std::cout << "\n=== Forward ===\n";
    std::cout << "Delta: " << delta_fwd(S,K,r,q,sigma,T) << "\n";
    std::cout << "Gamma: " << gamma_fwd(S,K,r,q,sigma,T) << "\n";

    std::cout << "\n=== Complex-Step ===\n";
    std::cout << "Delta: " << delta_cs(S,K,r,q,sigma,T) << "\n";
    std::cout << "Gamma: " << gamma_cs(S,K,r,q,sigma,T) << "\n";

    return 0;
}