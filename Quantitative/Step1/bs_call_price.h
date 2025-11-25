#pragma once
#include <cmath>
#include <complex>
#include <algorithm>

using cplx = std::complex<double>;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// High-precision CDF
inline double Phi_real(double z) {
    return 0.5 * std::erfc(-z * 0.70710678118654752440);
}

// === CORE BS PRICE (S can be double or complex) ===
template<typename S_type>
inline auto bs_price_call_impl(S_type S, double K, double r, double q, double sigma, double tau)
{
    using std::exp; using std::log; using std::sqrt; using std::max; using std::abs;
    using result_t = decltype(S * 1.0);

    const result_t DF = exp(-r * tau);
    const result_t F  = S * exp((r - q) * tau);
    const double sigmaT = sigma * sqrt(max(tau, 0.0));

    if (sigmaT == 0.0) {
        double payoff = (std::real(F) > K) ? (std::real(F) - K) : 0.0;
        return result_t(std::real(DF) * payoff);
    }

    result_t ln_F_over_K;
    if (K > 0.0) {
        result_t ratio = F / K;
        if (std::abs(ratio - result_t(1.0)) < 1e-8) {
            ln_F_over_K = result_t(std::log1p(std::real(ratio) - 1.0));
        } else {
            ln_F_over_K = log(ratio);
        }
    } else {
        ln_F_over_K = log(F / K);
    }

    result_t d1 = (ln_F_over_K + 0.5 * sigma * sigma * tau) / sigmaT;
    result_t d2 = d1 - sigmaT;

    double Phi_d1 = Phi_real(std::real(d1));
    double Phi_d2 = Phi_real(std::real(d2));

    return DF * (F * Phi_d1 - K * Phi_d2);
}

// === OVERLOADS ===
inline double bs_price_call(double S, double K, double r, double q, double sigma, double T) {
    return std::real(bs_price_call_impl(S, K, r, q, sigma, T));
}

inline cplx bs_price_call(cplx S, double K, double r, double q, double sigma, double T) {
    return bs_price_call_impl(S, K, r, q, sigma, T);
}

// === ANALYTIC GREEKS ===
inline double bs_delta_call(double S, double K, double r, double q, double sigma, double T) {
    const double F = S * std::exp((r - q) * T);
    const double sigmaT = sigma * std::sqrt(std::max(T, 0.0));
    if (sigmaT == 0.0) return (F > K) ? std::exp(-q * T) : 0.0;
    const double d1 = (std::log(F / K) + 0.5 * sigma * sigma * T) / sigmaT;
    return std::exp(-q * T) * Phi_real(d1);
}

inline double bs_gamma(double S, double K, double r, double q, double sigma, double T) {
    const double F = S * std::exp((r - q) * T);
    const double sigmaT = sigma * std::sqrt(std::max(T, 0.0));
    if (sigmaT == 0.0) return 0.0;
    const double d1 = (std::log(F / K) + 0.5 * sigma * sigma * T) / sigmaT;
    const double phi_d1 = std::exp(-0.5 * d1 * d1) / std::sqrt(2.0 * M_PI);
    return std::exp(-q * T) * phi_d1 / (S * sigmaT);
}

// === FORWARD DIFF ===
inline double delta_fwd(double S, double K, double r, double q, double sigma, double T, double h = 1e-3) {
    return (bs_price_call(S + h, K, r, q, sigma, T) - bs_price_call(S, K, r, q, sigma, T)) / h;
}

inline double gamma_fwd(double S, double K, double r, double q, double sigma, double T, double h = 1e-2) {
    double C0 = bs_price_call(S, K, r, q, sigma, T);
    double C1 = bs_price_call(S + h, K, r, q, sigma, T);
    double C2 = bs_price_call(S + 2 * h, K, r, q, sigma, T);
    return (C2 - 2 * C1 + C0) / (h * h);
}

// === COMPLEX-STEP (ROBUST & ACCURATE) ===
inline double delta_cs(double S, double K, double r, double q, double sigma, double T, double h = 1e-20) {
    cplx S_c(S, h);
    return bs_price_call(S_c, K, r, q, sigma, T).imag() / h;
}

inline double gamma_cs(double S, double K, double r, double q, double sigma, double T, double h = 1e-8) {
    double C0 = bs_price_call(S, K, r, q, sigma, T);
    cplx S1(S,  h);
    cplx S2(S, -h);
    cplx C1 = bs_price_call(S1, K, r, q, sigma, T);
    cplx C2 = bs_price_call(S2, K, r, q, sigma, T);
    return (C1.real() + C2.real() - 2.0 * C0) / (h * h);
}