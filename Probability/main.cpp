#include <iostream> 
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>
#include "bs_call_price.h" // your BS price header

using namespace std;
using cplx = complex<double>;

int main() {
    // === Scenario parameters ===
    struct Scenario {
        string name;
        double S, K, r, q, sigma, T;
    };
    vector<Scenario> scenarios = {
        {"scenario1", 100.0, 100.0, 0.0, 0.0, 0.20, 1.0},      // ATM
        {"scenario2", 100.0, 100.0, 0.0, 0.0, 0.01, 1.0/365.0} // Near-expiry low-vol
    };

    // Relative step sizes to sweep
    vector<double> h_rel_list = {1e-2, 1e-4, 1e-6, 1e-8};

    for (auto &sc : scenarios) {
        ofstream fout("bs_fd_vs_complex_" + sc.name + ".csv");
        fout << fixed << setprecision(8);
        fout << "h_rel,h,"
             << "Delta_analytic,Delta_fd,Delta_cs,err_D_fd,err_D_cs,"
             << "Gamma_analytic,Gamma_fd,Gamma_cs_real,Gamma_cs_45,"
             << "err_G_fd,err_G_cs_real,err_G_cs_45\n";

        // Analytic Delta & Gamma using correct names
        double delta_analytic = bs_delta_call(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T);
        double gamma_analytic = bs_gamma(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T);

        for (double h_rel : h_rel_list) {
            double h = h_rel * sc.S;

            // Forward difference
            double delta_fd = (bs_price_call_t(sc.S + h, sc.K, sc.r, sc.q, sc.sigma, sc.T) -
                               bs_price_call_t(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T)) / h;

            double gamma_fd = (bs_price_call_t(sc.S + h, sc.K, sc.r, sc.q, sc.sigma, sc.T) -
                               2.0 * bs_price_call_t(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T) +
                               bs_price_call_t(sc.S - h, sc.K, sc.r, sc.q, sc.sigma, sc.T)) / (h*h);

            // Complex-step Delta
            cplx S_c(sc.S, h);
            double delta_cs = imag(bs_price_call_t(S_c, sc.K, sc.r, sc.q, sc.sigma, sc.T)) / h;

            // Complex-step Gamma (real-part)
            double gamma_cs_real = -2.0 * (real(bs_price_call_t(S_c, sc.K, sc.r, sc.q, sc.sigma, sc.T)) -
                                           bs_price_call_t(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T)) / (h*h);

            // 45-degree complex-step
            cplx omega = cplx(1.0/sqrt(2.0), 1.0/sqrt(2.0));
            double gamma_cs_45 = imag(
                bs_price_call_t(sc.S + h*omega, sc.K, sc.r, sc.q, sc.sigma, sc.T) +
                bs_price_call_t(sc.S - h*omega, sc.K, sc.r, sc.q, sc.sigma, sc.T)
            ) / (h*h);

            // Absolute errors
            double err_D_fd = abs(delta_fd - delta_analytic);
            double err_D_cs = abs(delta_cs - delta_analytic);
            double err_G_fd = abs(gamma_fd - gamma_analytic);
            double err_G_cs_real = abs(gamma_cs_real - gamma_analytic);
            double err_G_cs_45 = abs(gamma_cs_45 - gamma_analytic);

            // CSV line
            fout << h_rel << "," << h << ","
                 << delta_analytic << "," << delta_fd << "," << delta_cs << "," << err_D_fd << "," << err_D_cs << ","
                 << gamma_analytic << "," << gamma_fd << "," << gamma_cs_real << "," << gamma_cs_45 << ","
                 << err_G_fd << "," << err_G_cs_real << "," << err_G_cs_45 << "\n";
        }

        fout.close();
        cout << "CSV written for " << sc.name << endl;
    }

    return 0;
}
