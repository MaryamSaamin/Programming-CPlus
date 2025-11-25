#pragma once
#include <cmath>
#include <cstddef>
#include <limits>

namespace quant {

class InverseCumulativeNormal {
public:
    explicit InverseCumulativeNormal(double average = 0.0, double sigma = 1.0)
      : average_(average), sigma_(sigma) {}

    inline double operator()(double x) const {
        return average_ + sigma_ * standard_value(x);
    }

    inline void operator()(const double* in, double* out, std::size_t n) const {
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = average_ + sigma_ * standard_value(in[i]);
        }
    }

    static inline double standard_value(double x) {
        if (x <= 0.0) return -std::numeric_limits<double>::infinity();
        if (x >= 1.0) return  std::numeric_limits<double>::infinity();

        if (x >= x_low_ && x <= x_high_) {
            double z = central_rational(x);
#ifdef ICN_ENABLE_HALLEY_REFINEMENT
            z = halley_refine(z, x);
#endif
            return z;
        } else {
            const bool left = (x < 0.5);
            double m = left ? x : (1.0 - x);
            if (m <= 0.0) return left ? -std::numeric_limits<double>::infinity()
                                      :  std::numeric_limits<double>::infinity();
            double zpos = std::abs(tail_rational_pos(m));
            double z = left ? -zpos : zpos;
#ifdef ICN_ENABLE_HALLEY_REFINEMENT
            z = halley_refine(z, x);
#endif
            return z;
        }
    }

private:
    static inline double phi(double z) {
        constexpr double INV_SQRT_2PI = 0.398942280401432677939946059934381868475858631164934657;
        return INV_SQRT_2PI * std::exp(-0.5 * z * z);
    }
    static inline double Phi(double z) {
        constexpr double INV_SQRT_2 = 0.707106781186547524400844362104849039284835937688474036588;
        return 0.5 * std::erfc(-z * INV_SQRT_2);
    }

    // Central region (Acklam) coefficients (m=6,n=6)
    static inline double central_rational(double x) {
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

        double num = a[0];
        num = num * r + a[1];
        num = num * r + a[2];
        num = num * r + a[3];
        num = num * r + a[4];
        num = num * r + a[5];

        double den = b[0];
        den = den * r + b[1];
        den = den * r + b[2];
        den = den * r + b[3];
        den = den * r + b[4];
        den = den * r + 1.0;

        return q * (num / den);
    }

    // Tail rational (Acklam) coefficients (p=6,q=4)
    static inline double tail_rational_pos(double m) {
        constexpr double c[] = {
            -7.784894002430293e-03,
            -3.223964580411365e-01,
            -2.400758277161838e+00,
            -2.549732539343734e+00,
             4.374664141464968e+00,
             2.938163982698783e+00
        };
        constexpr double d[] = {
             7.784695709041462e-03,
             3.224671290700398e-01,
             2.445134137142996e+00,
             3.754408661907416e+00
        };

        const double q = std::sqrt(-2.0 * std::log(m));

        double num = c[0];
        num = num * q + c[1];
        num = num * q + c[2];
        num = num * q + c[3];
        num = num * q + c[4];
        num = num * q + c[5];

        double den = d[0];
        den = den * q + d[1];
        den = den * q + d[2];
        den = den * q + d[3];
        den = den * q + 1.0;

        return num / den;
    }

#ifdef ICN_ENABLE_HALLEY_REFINEMENT
    static inline double stable_residual(double z, double x) {
        const double p = phi(z);
        const double f = Phi(z);
        if (f > 1e-12 && (1.0 - f) > 1e-12) {
            return (f - x) / std::max(p, std::numeric_limits<double>::min());
        }
        const double absz = std::fabs(z);
        double phi_ratio;
        if (absz > 1e2) {
            phi_ratio = 1.0 / absz;
        } else {
            const double arg = -z * 0.7071067811865475244;
            double erfc_val = std::erfc(arg);
            double logPhi;
            if (erfc_val > 0.0) logPhi = std::log(0.5) + std::log(erfc_val);
            else logPhi = std::log(p) - std::log(absz + 1.0);
            double logphi = std::log(std::max(p, std::numeric_limits<double>::min()));
            phi_ratio = std::exp(logPhi - logphi);
        }
        const double inv_phi = 1.0 / std::max(p, std::numeric_limits<double>::min());
        return phi_ratio - x * inv_phi;
    }
    static inline double halley_refine(double z, double x) {
        const double r = stable_residual(z, x);
        const double denom = 1.0 - 0.5 * z * r;
        return z - r / (denom != 0.0 ? denom
                                     : std::copysign(std::numeric_limits<double>::infinity(), denom));
    }
#endif

    double average_, sigma_;
    static constexpr double x_low_  = 0.02425;
    static constexpr double x_high_ = 1.0 - x_low_;
};

} // namespace quant
