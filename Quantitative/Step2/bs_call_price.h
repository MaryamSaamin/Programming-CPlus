#pragma once
#include <cmath>
#include <complex>
#include <algorithm>

using cplx = std::complex<double>;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// High-precision CDF (real argument)
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
        // robust log for ratio near 1
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

// === ANALYTIC GREEKS (exact) ===
// Delta: exp(-qT) * Phi(d1)
// Gamma: exp(-qT) * phi(d1) / (S * sigma * sqrt(T))
// Compute phi(d1) via log phi to avoid underflow.
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
    // log-safe phi
    const double log_phi = -0.5 * d1 * d1 - 0.5 * std::log(2.0 * M_PI);
    const double phi_d1 = std::exp(log_phi);
    return std::exp(-q * T) * phi_d1 / (S * sigmaT);
}

// === FORWARD DIFF (classical) ===
inline double delta_fwd(double S, double K, double r, double q, double sigma, double T, double h = 1e-3) {
    return (bs_price_call(S + h, K, r, q, sigma, T) - bs_price_call(S, K, r, q, sigma, T)) / h;
}

// central-forward-style for gamma (as requested by instruction)
inline double gamma_fwd(double S, double K, double r, double q, double sigma, double T, double h = 1e-2) {
    double C0 = bs_price_call(S, K, r, q, sigma, T);
    double C1 = bs_price_call(S + h, K, r, q, sigma, T);
    double C2 = bs_price_call(S + 2 * h, K, r, q, sigma, T);
    return (C2 - 2.0 * C1 + C0) / (h * h);
}

// === COMPLEX-STEP (robust & accurate) ===
// First derivative (delta): Im[f(S + i h)] / h
inline double delta_cs(double S, double K, double r, double q, double sigma, double T, double h = 1e-8) {
    cplx S_c(S, h);
    return bs_price_call(S_c, K, r, q, sigma, T).imag() / h;
}

// Second derivative (gamma) - using the estimator derived from Re[f(S + i h)]
// f''(S) ≈ -2 * (Re f(S + i h) - f(S)) / h^2 (truncation error O(h^2))
// algebraically this equals - (Re f(S+ih) + Re f(S-ih) - 2 f(S)) / h^2
inline double gamma_cs(double S, double K, double r, double q, double sigma, double T, double h = 1e-8) {
    double C0 = bs_price_call(S, K, r, q, sigma, T);            // f(S)
    cplx S_plus(S, h);                                          // S + i h
    cplx S_minus(S, -h);                                        // S - i h
    double C_plus  = std::real(bs_price_call(S_plus, K, r, q, sigma, T)); 
    double C_minus = std::real(bs_price_call(S_minus, K, r, q, sigma, T));
    return (C_plus + C_minus - 2.0 * C0) / (h * h);  // ✅ no negative sign
}


// Optional: higher-order (45-degree) complex-step estimator for second derivative
// ω = e^{iπ/4} = (1+i)/sqrt(2)
// f''(x) ≈ Imag( f(x + h ω) + f(x - h ω) ) / h^2   (truncation O(h^4))
inline double gamma_cs45(double S, double K, double r, double q, double sigma, double T, double h = 1e-6) {
    const cplx omega = cplx(1.0,1.0) / std::sqrt(2.0);
    cplx z1 = cplx(S,0.0) + h * omega;
    cplx z2 = cplx(S,0.0) - h * omega;
    cplx f1 = bs_price_call(z1, K, r, q, sigma, T);
    cplx f2 = bs_price_call(z2, K, r, q, sigma, T);
    return (f1.imag() + f2.imag()) / (h * h);
}

// High-quality normal PDF (real)
inline double phi_real(double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI);
}

// Type-dispatched CDF Φt (works for real and complex)
template<typename T>
inline auto Phi_t(T z) {
    if constexpr (std::is_same_v<T, double>) {
        return Phi_real(z);  // real branch
    } else if constexpr (std::is_same_v<T, cplx>) {
        double zr = std::real(z);
        double zi = std::imag(z);
        return cplx(Phi_real(zr), zi * phi_real(zr));  // complex branch (first-order Taylor)
    }
}


