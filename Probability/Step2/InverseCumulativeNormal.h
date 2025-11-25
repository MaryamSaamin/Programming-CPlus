#pragma once
#include <cmath>
#include <cstddef>
#include <limits>
#include <algorithm>

namespace quant {

class InverseCumulativeNormal {
  public:
    explicit InverseCumulativeNormal(double average = 0.0, double sigma = 1.0)
    : average_(average), sigma_(sigma) {}

    // Scalar call: return average + sigma * Φ^{-1}(x)
    inline double operator()(double x) const {
        return average_ + sigma_ * standard_value(x);
    }

    // Vector overload: out[i] = average + sigma * Φ^{-1}(in[i]) for i in [0, n)
    inline void operator()(const double* in, double* out, std::size_t n) const {
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = average_ + sigma_ * standard_value(in[i]);
        }
    }

    // Standardized value: inverse CDF with average=0, sigma=1.
    static inline double standard_value(double x) {
        // Defensive edge cases
        if (x <= 0.0) return -std::numeric_limits<double>::infinity();
        if (x >= 1.0) return  std::numeric_limits<double>::infinity();

        // Piecewise: central fast rational, tails fallback to bisection (robust)
        if (x >= x_low_ && x <= x_high_) {
            double z = central_rational(x);
        #ifdef ICN_ENABLE_HALLEY_REFINEMENT
            z = halley_refine(z, x);
        #endif
            return z;
        } else {
            // Use the safe bisection baseline for tails to guarantee correct sign & magnitude.
            // (This is slower for extreme tails; replace later with a rational tail for speed.)
            return invert_bisect(x);
        }
    }

  private:
    // ---- Baseline numerics (stable bisection retained for tails) ----------

    // Standard normal pdf
    static inline double phi(double z) {
        constexpr double INV_SQRT_2PI =
            0.398942280401432677939946059934381868475858631164934657; // 1/sqrt(2π)
        return INV_SQRT_2PI * std::exp(-0.5 * z * z);
    }

    // Standard normal cdf using erfc: Φ(z) = 0.5 * erfc(-z/√2)
    static inline double Phi(double z) {
        constexpr double INV_SQRT_2 =
            0.707106781186547524400844362104849039284835937688474036588; // 1/√2
        return 0.5 * std::erfc(-z * INV_SQRT_2);
    }

    // Crude but reliable invert via bisection; used here for tails (robust).
    static inline double invert_bisect(double x) {
        double lo = -12.0;
        double hi =  12.0;
        if (x < 0.5) {
            hi = 0.0;
        } else {
            lo = 0.0;
        }
        for (int iter = 0; iter < 80; ++iter) {
            double mid = 0.5 * (lo + hi);
            double cdf = Phi(mid);
            if (cdf < x) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    // ---- Central-region rational (Acklam) --------------------------------

    static inline double central_rational(double x) {
        // Coefficients (Acklam), ordered for Horner with increasing power:
        constexpr double a[] = {
            -3.969683028665376e+01,
             2.209460984245205e+02,
            -2.759285104469687e+02,
             1.383577518672690e+02,
            -3.066479806614716e+01,
             2.506628277459239e+00
        };
        constexpr double b[] = {
            -5.447609879822406e+01,
             1.615858368580409e+02,
            -1.556989798598866e+02,
             6.680131188771972e+01,
            -1.328068155288572e+01
        };

        double q = x - 0.5;
        double r = q * q;

        // Horner: numerator = a0 + a1*r + a2*r^2 + ...
        double num = a[0];
        num = num * r + a[1];
        num = num * r + a[2];
        num = num * r + a[3];
        num = num * r + a[4];
        num = num * r + a[5];

        // Horner: denominator = b0 + b1*r + ... + 1*r^5 (we add the trailing +1)
        double den = b[0];
        den = den * r + b[1];
        den = den * r + b[2];
        den = den * r + b[3];
        den = den * r + b[4];
        den = den * r + 1.0;

        return q * (num / den);
    }

#ifdef ICN_ENABLE_HALLEY_REFINEMENT
    // One-step Halley refinement (3rd order). Usually brings result to full double precision.
    static inline double halley_refine(double z, double x) {
        const double f = Phi(z);
        const double p = phi(z);
        const double r = (f - x) / std::max(p, std::numeric_limits<double>::min());
        const double denom = 1.0 - 0.5 * z * r;
        return z - r / (denom != 0.0 ? denom
                                     : std::copysign(std::numeric_limits<double>::infinity(), denom));
    }
#endif

    // ---- State & constants ------------------------------------------------

    double average_, sigma_;

    // Region split: keep Acklam's recommended split for central vs tail
    static constexpr double x_low_  = 0.02425;         // ~ Φ(-2.0)
    static constexpr double x_high_ = 1.0 - x_low_;
};

} // namespace quant
